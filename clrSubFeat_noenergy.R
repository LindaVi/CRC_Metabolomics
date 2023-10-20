# Functions for performing CLR on continuous data (assuming linearity)


# Convenience function to do stratified CLR by site
clrSubFeat_noenergy <- function(scores, ...) {
  clrSite <- function(name, scores, site, ...) {
    siteclr <- clrFeat_noenergy(scores, site = site, ...)
    if (nrow(siteclr)==0) return(NULL) else return(data.frame(site=name, siteclr))
  }
  CLR <- rbind(
    clrSite(name='1+2: Colon', scores, site = 1:2, ...),
    clrSite(name='1: Proximal colon', scores, site = 1, ...),
    clrSite(name='2: Distal colon', scores, site = 2, ...),
    clrSite(name='3: Rectum', scores, site = 3, ...)
  ) 
  return(CLR)
}

# main CLR function
clrFeat_noenergy <- function(scores, scoreName, verbose = F, meta = dd2, meas = 1, site = c(1:3, 9999), sex = 1:2, p = 0.1, noAlco = FALSE, alcoBinary=FALSE, scale=TRUE, scaleBySex=FALSE, round = F) {
  library(survival)
  
  if(missing(scoreName)) scoreName <- colnames(scores)
  
  # selection of individuals in analysis (based on measurement occation and site) -> "include" variable
  sitePair <- meta$caseset2[meta$lokal3_new%in%site]
  include <- meta$measurement==meas & meta$caseset2%in%sitePair & meta$gender%in%sex
  
  # Limit observations to those included
  scores <- scores[include,, drop=F] 
  gender <- meta$gender[include] # 1 for male and 2 for female
  if (scale) {
    if(scaleBySex) {
      if (sum(gender==1)>0) scores[gender==1,] <- scale(scores[gender==1,])
      if (sum(gender==2)>0) scores[gender==2,] <- scale(scores[gender==2,])
    } else {
      scores <- scale(scores)
    }
  }
  pairs <- meta$caseset2[include]
  casestat <- meta$casecontrol[include]
  cleanVarFact <- function(var) { # Some different namings for different types of missingness -> harmonized to '9999' for all missingness
    var[var==6666] <- 9999
    var[var==8888] <- 9999
    var[is.na(var)] <- 9999
    return(var)
  }
  age <- meta$age[include]
  bmi <- meta$bmi[include]
  year <- meta$pyear[include]
  # Physical activity needs some TLC
  pa <- meta$g6[include] # Main PA variable
  pa <- StatTools::switchNames(pa, c(1,2,3,4,5,8888,9999),c(1,2,3,3,4,8888,9999))
  pa_mon <- meta$monica_motion_fritid_86_09[include] # Auxiliary PA assessment
  pa_mon <- StatTools::switchNames(pa_mon, c(1,2,3,4,5,6,8888,9999),c(1,2,3,3,4,4,8888,9999)) # Convert 6-level scale to 4-level scale
  pa[pa==8888] <- pa_mon[pa==8888] # Use auxiliary PA, when main variable is missing
  pa <- pa %>% cleanVarFact %>% as.factor()
  energy <- meta$ensum1[include]
  smoke <- meta$sm_status[include]
  smoke <- StatTools::switchNames(smoke, c(1,2,3,4,5,6666,9999), c(3,2,1,3,2,9999,9999)) # Convert 5-level scale to 3-level scale (current, former, never)
  smoke <- smoke %>% cleanVarFact  %>% as.factor()
  edu <- meta$utbild[include] 
  edu <- StatTools::switchNames(edu, c(1,2,3,4,9999), c(1,2,2,3,9999)) # Convert 4-level scale to 3-level scale
  edu <- edu %>% cleanVarFact  %>% as.factor()
  marital <- meta$civil[include] %>% cleanVarFact %>% as.factor()
  alcoSum <- meta$alkosum1[include]
  # Rearrange alcohol to by-gender 
  alco <- character(sum(include))
  alcoMedMale <- median(alcoSum[alcoSum!=0 & gender==1], na.rm = T)
  alcoMedFem <- median(alcoSum[alcoSum!=0 & gender==2], na.rm = T)
  alco[gender==1 & alcoSum <= alcoMedMale] <- 'low'
  alco[gender==1 & alcoSum > alcoMedMale] <- 'high'
  alco[gender==2 & alcoSum <= alcoMedFem] <- 'low'
  alco[gender==2 & alcoSum > alcoMedFem] <- 'high'
  alco[alcoSum==0] <- 'abstain' # Abstainers are separate category
  if (alcoBinary) alco <- ifelse(alco=='abstain', 'abstain', 'consume')
  if(verbose) print(table(alco))
  
  # Allocate model output objects
  models <- list()
  pvect <- numeric(ncol(scores))
  or <- character(ncol(scores))
  est <- se <- n <- numeric(ncol(scores))
  
  # Perform CLR for all exposures (columns in score DF)
  for (i in 1:ncol(scores)) {
    
    # Prepare dataframe for analysis
    myDF <- data.frame(score=scores[,i], casestat, pairs, bmi, energy, pa, edu, smoke, alco)
    # Keep complete cases only
    myDF <- myDF[complete.cases(myDF),]
    # Sort out complete case-control pairs (incomplete data often only affects either case or control in the pairs)
    doubles <- myDF$pairs[duplicated(myDF$pairs)]
    myDF <- myDF[myDF$pairs%in%doubles,]
    
    # Separate models for whether to adjust for alcohol intake (e.g. alcohol exposure shouldn't adjust for alcohol)
    if (noAlco) {
      clrMod <- clogit(casestat ~ score + bmi + pa + edu + smoke + strata(pairs), data=myDF)
    } else {
      clrMod <- clogit(casestat ~ score + bmi + pa + edu + smoke + alco + strata(pairs), data=myDF)
    }
    
    # Extract model outputs
    models[[i]] <- summary(clrMod)
    if(verbose) print(models[[i]])
    coef <- models[[i]]$coefficients
    pvect[i] <- coef[1,5]
    est[i] <- coef[1,1]
    se[i] <- coef[1,3]
    or[i] <- paste0(exp(coef[1,1]) %>% sprintf('%.3f',.), ' (',exp(coef[1,1]-1.96*coef[1,3]) %>% sprintf('%.3f',.),'-',exp(coef[1,1]+1.96*coef[1,3]) %>% sprintf('%.3f',.),')')
    n[i] <- clrMod$n
  }
  
  # round off p-values
  if (round) pvect <- pvect %>% sprintf('%.3f',.)
  # Rearrange to dataframe and limit to those with p < limit
  clrDF <- data.frame(variable=scoreName, n=n, est=est, se=se, OR=or, p=pvect)
  clrDF <- clrDF[clrDF$p < p, ]
  
  # return output
  return(clrDF)
}

