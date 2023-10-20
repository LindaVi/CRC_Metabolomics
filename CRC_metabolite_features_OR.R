rm(list=ls())

sinkfile <- file('CRC_metabolite_features_OR_sink.txt', open = "wt")
sink(sinkfile)
sink(sinkfile, type = 'message')

# Load data objects
load(file='combined_data_4_MB_at_BL.rda')
# Not available in Repo due to GDPR
# C_meta (containing all outcome variables and covariates)

# Load original metabolite data00
load(file='combAllFeatures.rda') # Data wrangling in "TargetedHypothesesMBDataWrangl.r"
# Not available in Repo due to GDPR
# RN: peak table in Reversed phase - Negative (4461 analytical features - before RAMCLUST) - row-wise matched to C_meta
# > head(colnames(RN))
# [1] "170.833153089859_31.3572528845699" "168.836022973106_31.3587694975369" "172.830741241854_31.3803517340893" "90.9334456616951_31.6270713280501"
# [5] "158.978557797786_32.3358836381859" "365.885545476401_34.4323109227088"
# RP: peak table in Reversed phase - Positive (4808 analytical features - before RAMCLUST) - row-wise matched to C_meta
# > head(colnames(RP))
# [1] "194.11744504438_6.10549748003088"  "118.086413576838_10.3610150064102" "222.112377003487_27.1242139955774" "118.08642234546_27.6089883203931" 
# [5] "194.117459424836_27.6338992725433" "224.128109788562_27.9615758286661"

# Load candidates from literature data - for which to perform replication
# Supplied in Repo
load(file = 'targeted.rda')
targ$hMIM[targ$hMIM == 0] <- NA
colnames(targ)
# [1] "Metabolite.name"                                    "HMDB.number"                                       
# [3] "Reported.mass"                                      "MIM"                                               
# [5] "RT..Sec."                                           "mode"                                              
# [7] "Adduct"                                             "Prominent.MS2.fragments.with.possible.fragment.IDs"
# [9] "Study"                                              "Prospective"                                       
# [11] "Biospecimen"                                        "Chemical.formula"                                  
# [13] "Kommentar"                                          "hMIM"                                              
# [15] "hWarn"                                              "hError"


# Extract candidates with masses 
# (neg_mz are for candidates with reported m/z in negative mode)
# (pos_mz are for candidates with reported m/z in positive mode)
# (calc are Monoisotopic masses (MIM) for candidates reported without m/z, but matched by metabolite ID)
neg_mz <- tolower(targ$mode)=='neg' & !is.na(targ$Reported.mass)
pos_mz <- tolower(targ$mode)=='pos' & !is.na(targ$Reported.mass)
calc <- !is.na(targ$MIM)
sum(neg_mz) # 127 candidates with negative m/z
sum(pos_mz) # 59 candidates with positive m/z
sum(calc) # 287 candidates with MIM -> search masses with potential adducts

# Extract corresponding candidate datasets
neg_mz <- targ[neg_mz,]
pos_mz <- targ[pos_mz,]
calc <- targ[calc,]

# Load functions to find hits in data potentially corresponding to candidates
source(file = 'getHits.R') # Supplied in Repo
ppm <- 10 # Tolerance for finding candidate hits in actual peak tables
# Potential hits (m/z) in Negative 
negHits <- extractFeats(mz = neg_mz$Reported.mass, name = paste0(neg_mz$Metabolite.name, '_mzNeg'), PT = RN, ppm = ppm)
negHits <- pairedTest(extracts = negHits, PT = RN) # Simple t-test between groups for analytical feature - potential hits
# Potential hits (m/z) in Positive
posHits <- extractFeats(mz = pos_mz$Reported.mass, name = paste0(pos_mz$Metabolite.name, '_mzPos'), PT = RP, ppm = ppm)
posHits <- pairedTest(extracts = posHits, PT = RP)
# Potential hits (Calculated mass) in Positive 
candCalcPos <- mzCalc(name = calc$Metabolite.name, mim = calc$MIM, ionization = 'pos', adduct = 'std')
CalcPosHits <- extractFeats(mz = candCalcPos$mz, name = paste(candCalcPos$name, candCalcPos$adduct, 'calc_pos', sep = '_'), PT = RP, ppm = ppm)
CalcPosHits <- pairedTest(extracts = CalcPosHits, PT = RP)
# Potential hits (Calculated mass) in Negative 
candCalcNeg <- mzCalc(name = calc$Metabolite.name, mim = calc$MIM, ionization = 'neg', adduct = 'std')
CalcNegHits <- extractFeats(mz = candCalcNeg$mz, name = paste(candCalcNeg$name, candCalcNeg$adduct, 'calc_neg', sep = '_'), PT = RN, ppm = ppm)
CalcNegHits <- pairedTest(extracts = CalcNegHits, PT = RN)

# Combine all potential hits
allHits <- rbind(negHits,
                 CalcNegHits,
                 posHits,
                 CalcPosHits)

head(allHits)
# name   lookup feat_nr                         feat_name  feat_mz  feat_rt           dmz       ppm   meanInt    minInt   maxInt      pRand       FC    FCGeom
# 1 D-Ribulose 5-phosphate_mzNeg 210.9931    1145   210.994364807223_155.4644989408 210.9944 155.4645  0.0012648072  5.994543 62277.055  667.2734 222907.3 0.03345278 1.059206 1.0603475
# 2   D-Ribose 5-phosphate_mzNeg 210.9934    1145   210.994364807223_155.4644989408 210.9944 155.4645  0.0009648072  4.572689 62277.055  667.2734 222907.3 0.03345278 1.059206 1.0603475
# 3   Xylulose 5-phosphate_mzNeg 210.9929    1145   210.994364807223_155.4644989408 210.9944 155.4645  0.0014648072  6.942448 62277.055  667.2734 222907.3 0.03345278 1.059206 1.0603475
# 4     Ribose 1-phosphate_mzNeg 210.9933    1145   210.994364807223_155.4644989408 210.9944 155.4645  0.0010648072  5.046640 62277.055  667.2734 222907.3 0.03345278 1.059206 1.0603475
# 5   Glycyl-Phenylalanine_mzNeg 221.0840    1300 221.081934345434_255.683367831992 221.0819 255.6834 -0.0020656546 -9.343302  9089.697  142.7036 117293.9 0.05893916 0.919889 0.9875029
# 6   Glycyl-Phenylalanine_mzNeg 221.0840    1490  221.08260176329_311.458860482799 221.0826 311.4589 -0.0013982367 -6.324459 23413.125 1688.6178 313729.0 0.49818899 1.019562 1.0034467

dim(allHits) # 1051 x 14
# Find and remove duplicate hits
duplicated(allHits$feat_name) %>% sum # 332
allHits <- allHits[!duplicated(allHits$feat_name),]
allHits$fdrRand <- p.adjust(allHits$pRand, method='fdr')

####################################
# Conditional logistic regression for potential hits

# sub-analyses by site
source(file = 'clrSubFeat_noenergy.R') # Submitted to Repo

####################################
# Variables from metabolomics

p <- 1
scaleBySex <- TRUE
polarity <- substring(allHits$name,nchar(allHits$name)-2) %>% tolower

# Extract metabolites-of-interest from positive and negative modes respectively
M_sub <- list()
for (i in 1:nrow(allHits)) {
  if (polarity[i] == 'neg') {
    M_sub[[i]] <- RN[,allHits$feat_nr[i]]
  } else {
    M_sub[[i]] <- RP[,allHits$feat_nr[i]]
  }
}
M_sub <- do.call(cbind, M_sub) # 1051 variables in matrix
colnames(M_sub) <- allHits$name

# Perform actual CLR
featCLR <- clrFeat_noenergy(scores=M_sub, meta = C_meta, p = p, scaleBySex = scaleBySex)
colnames(featCLR)
featCLR$fdr <- p.adjust(featCLR$p, method = 'fdr')
allHits$pCLR <- featCLR$p
allHits$fdrCLR <- featCLR$fdr
allHits$OR <- featCLR$OR

# Selected hits from CLR
selHits <- allHits[allHits$pCLR<0.05,] # 61 metabolite features x 17 columns
selHits

sink(type = 'message')
sink()