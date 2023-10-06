# DESCRIPTION:
# Classification of cases vs controls using the MUVR algorithm (random forest )


rm(list=ls())
#Load required librarys
library(MUVR)
library(doParallel)


# REQUIRED INPUT
# MB_data - a data.frame were each row corresponds to one sample and each column to a metabolic feature
# MB_Conf_data - a data.frame were the possible confounders bmi, alkohol, smoking, physical activity, education and diabetes is added to the metabolic features 
# OutC - data frame where one column is named casecontrol and corresponds to disease status (0 or 1) and another is named caseset and indicates which samples that are matched (integers) 

MB_data <- readRDS(file='Path_to_Data/MB_data.rds')
MB_conf_data <- readRDS(file='Path_to_Data/MB_conf_data.rds')
OutC <- readRDS(file='Path_to_Data/Outcome.rds')



#Set parameters for random forest analysis

md <- 'RF'
vR <- 0.85
nrep <- 30
nO <- 6


# Set up for parallel processing
nCore <- detectCores()-1

cl=makeCluster(nCore)
registerDoParallel(cl)

#Random forest analysis with built in variable selection based on 1) metabolomic featutes  2) metabolomic features together with bmi, smoking, alkohol, education,physical activity and diabetes variables
Res_MB_data <- qMUVR(X=MB_data, Y=OutC$casecontrol,ID=OutC$caseset, method=md, varRatio = vR, nrep = nrep, nOuter = nO, nCore=nCore)
Res_MB_Conf_data <- qMUVR(X=MB_conf_data, Y=OutC$casecontrol,ID=OutC$caseset, method=md, varRatio = vR, nrep = nrep, nOuter = nO, nCore=nCore)

stopCluster(cl)



  
#Balanced error rates for classification based on minimal optimal model
  BER_MB <- BER(Res_MB_data$inData$Y,Res_MB_data$yClass[,1]) 
  BER_MB_Cobf <- BER(Res_MB_Conf_data$inData$Y,Res_MB_Conf_data$yClass[,1])
  
  






