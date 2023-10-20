mzCalc <- function(mim, name, ionization=c('pos','neg'), adduct) {
  if(missing(adduct)) adduct <- ifelse(ionization=='pos','M+H','M-H')
  if(missing(name)) name <- mim
  if(ionization=='pos') {
    if(!exists('addPos')) data('adductLists', package = 'CMSITools')
    adductList <- addPos
    adductList <- rbind(adductList, adductList[1,])
    nAdd <- nrow(adductList)
    adductList$Ion.name[nAdd] <- 'M-H2O+H'
    adductList$Ion.mass[nAdd] <- NA
    adductList$Charge[nAdd] <- 1
    adductList$Mass[nAdd] <- 17.0038
  } else {
    if(!exists('addNeg')) data('adductLists', package = 'CMSITools')
    adductList <- addNeg
  }
  if(adduct[1]=='all') adduct <- adductList$Ion.name
  if(adduct[1]=='std') {
    if(ionization=='pos') {
      adduct <- c('M+H','M+Na','M+K','M+NH4','M+2H','M+CH3OH+H','M+ACN+H','M-H2O+H')
    } else {
      adduct <- c('M-H','M-H2O-H','M-2H','M+Na-2H','M+K-2H','M+Cl','M+FA-H','M+Hac-H')
    }
  }
  adductList <- adductList[adductList$Ion.name%in%adduct,]
  reportList <- list()
  for (m in 1:length(mim)) {
    reportList[[m]] <- data.frame(name=name[m], mim=mim[m], polarity=ionization, adduct=adductList$Ion.name, mz=adductList$Mass+mim[m]/abs(adductList$Charge))
  }
  return(do.call(rbind, reportList))
}

featureInfo <- function(names, start = 1, separator = '_', mzcol = 1) {
  names <- substring(names, start)
  info <- strsplit(names, separator)
  info <- do.call(rbind, info)
  featInfo <- data.frame(mz = info[,mzcol] %>% as.numeric(), 
                         rt = info[,-mzcol] %>% as.numeric())
}


extractFeats <- function(mz, name, PT, ppm = 10) {
  if(missing(name)) name <- mz
  nLookup <- length(mz)
  if (nLookup!=length(name)) stop('different lengths of mz and names in extractFeats()')
  featInfo <- featureInfo(colnames(PT))
  hits <- list()
  r <- 0
  for (i in 1:nLookup) {
    lookup <- mz[i]
    dmz <- ppm * lookup * 1e-6
    whichHit <- which(abs(featInfo$mz-lookup) < dmz)
    if(length(whichHit > 0)) {
      r <- r + 1
      diff <- featInfo$mz[whichHit] - lookup
      hits[[r]] <- data.frame(name = name[i],
                              lookup = lookup, 
                              feat_nr = whichHit, 
                              feat_name = colnames(PT)[whichHit], 
                              feat_mz = featInfo$mz[whichHit], 
                              feat_rt = featInfo$rt[whichHit],
                              dmz = diff,
                              ppm = 1e6 * diff / lookup)
    }
  }
  hits <- do.call(rbind, hits)
}

pairedTest <- function(extracts, PT, meta = C_meta) {
  nFeat <- nrow(extracts)
  pFix <- pRand <- numeric(nFeat)
  meanInt <- minInt <- maxInt <- numeric(nFeat)
  FC <- FCGeom <- numeric(nFeat)
  for (i in 1:nFeat) {
    cat('\nHit:', extracts$name[i], '- Feature:', extracts$feat_name[i], '- Nr:', extracts$feat_nr[i])
    featDF <- data.frame(feature = PT[,extracts$feat_nr[i]], 
                         pair = as.factor(meta$caseset2), 
                         casecontrol = as.factor(meta$casecontrol))
    featDF <- featDF[order(featDF$casecontrol, featDF$pair),]
    # t1 <- proc.time()[3]
    # mod <- glm(feature ~ casecontrol + pair, 
    #            data = featDF)
    # aov <- anova(mod, test='F')
    # t1 <- proc.time()[3] - t1
    # pFix[i] <- aov$`Pr(>F)`[2]
    # t2 <- proc.time()[3]
    mod <- lme4::lmer(feature ~ casecontrol + (1 | pair),
                      data = featDF)
    aov <- car::Anova(mod)
    # t2 <- proc.time()[3] - t2
    pRand[i] <- aov$`Pr(>Chisq)`
    meanInt[i] <- mean(featDF$feature)
    minInt[i] <- min(featDF$feature)
    maxInt[i] <- max(featDF$feature)
    FC[i] <- mean(featDF$feature[featDF$casecontrol==1]) / mean(featDF$feature[featDF$casecontrol==0])
    FCGeom[i] <- log(featDF$feature[featDF$casecontrol==1] / featDF$feature[featDF$casecontrol==0]) %>% mean %>% exp
  }
  # extracts$pFix <- pFix
  extracts$meanInt <- meanInt
  extracts$minInt <- minInt
  extracts$maxInt <- maxInt
  extracts$pRand <- pRand
  extracts$FC <- FC
  extracts$FCGeom <- FCGeom
  return(extracts)
}

getCoelutes <- function(name, featName, dRT = 2, rLim = 0.7) {
  if(missing(name)) name <- featName
  nFeat <- length(featName)
  coElutes <- list()
  for (i in 1:nFeat) {
    polarity <- substring(name[i],nchar(name[i])-2) %>% tolower
    if (polarity == 'neg') PT <- RN else PT <- RP
    featInfo <- featureInfo(colnames(PT))
    featNr <- which(colnames(PT) == featName[i])
    lookup_mz <- featInfo$mz[featNr]
    lookup_rt <- featInfo$rt[featNr]
    whichCoelutes <- which(abs(featInfo$rt-lookup_rt) < dRT)
    whichCoelutes <- c(featNr, whichCoelutes[!whichCoelutes%in%featNr])
    coElutes[[i]] <- data.frame(featNr = whichCoelutes, 
                                featName = colnames(PT)[whichCoelutes],
                                mz = featInfo$mz[whichCoelutes],
                                rt = featInfo$rt[whichCoelutes],
                                r = cor(PT[,whichCoelutes])[,1],
                                fc = colMeans(PT[,whichCoelutes]) / mean(PT[,featNr]),
                                meanInt = colMeans(PT[,whichCoelutes]))
    coElutes[[i]] <- coElutes[[i]][coElutes[[i]]$r>rLim,]
  }
  return(coElutes)
}
