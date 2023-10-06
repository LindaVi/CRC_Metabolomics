
#Load required packages
library(survival)
library(dplyr)

#################################

#load datasets
#MB_data is a data frame where each coulmn is a metabolic feature and each row corresponds to one sample
#Conf is a data frame with confounder variables containing the columns bmi, smoking, diabetes, Alkohol, phys_act
#OutC is a data fram with two columns; casecontrol (0 or 1) which describes CRC status and caseset which is an integer describing which samples that are matched

MB_data <- readRDS(file='Path_to_Data/MBdata.rds')
Conf <- readRDS(file='Path_to_Data/Confounders.rds')
OutC <- readRDS(file='Path_to_Data/Outcome.rds')


#Initialize vectors
clogit_v <- vector()
clogit_v2 <- vector()


Nr_mb <- dim(MB_data)[2]
for(i in 1:Nr_mb)
{
  
  #Conditional logistic regression with metabolomic features
  clog_temp <- clogit(OutC$casecontrol ~ MB_data[,i] +strata (OutC$caseset))
  clogit_v[i] <-summary(clog_temp)$coefficients[5]
  
  
  #Conditional logistic regression with metabolomic features and possible confounders
  clog_temp2 <- clogit(OutC$casecontrol ~ MB_data[,i]+Conf$bmi + Conf$smoking +Conf$diabetes + Conf$alkohol + Conf$education + Conf$phys_act   +strata (OutC$caseset))
  clogit_v2[i] <- summary(clog_temp2)$coefficients[1,5]
  

}

#Cut offs for significance (nominal p <0.05 and false discovery rate (FDR) controlled at 0.25)
SL <- 0.05
BH_SL <- 0.25


#Number of significant features
sum(clogit_v<SL)
sum(clogit_v2<SL)


P_clogit <- p.adjust(clogit_v, method="BH")
sum(P_clogit < BH_SL)

P_clogit2 <- p.adjust(clogit_v2, method="BH")
sum(P_clogit2 < BH_SL)


#look at Odds ratios and 95% confidence intervals for features with FDR <0.25

ind_sign <- adjP_clogit<BH_SL
MB_sign <- MB_data[,ind_sign]

mod1_OR <- vector()
mod2_OR <- vector()
mod1_CI <- vector()
mod2_CI <- vector()


for(n in 1:dim(MB_sign)[2])
{
  
  clog_model <- clogit(OutC$casecontrol ~ MB_sign[,n] +strata (OutC$caseset))
  mod1_OR[n]<-exp(coef(clog_model))
  S <- summary(clog_model)
  mod1_CI[n] <- paste0('(',exp(S$coefficients[1,1]-1.96*S$coefficients[1,3]) %>% round(.,6),'-',exp(S$coefficients[1,1]+1.96*S$coefficients[1,3]) %>% round(.,6),')')
  
  clog_model2 <- clogit(OutC$casecontrol ~ MB_sign[,n]+Conf$bmi + Conf$smoking +Conf$diabetes + Conf$alkohol + Conf$education + Conf$phys_act  +strata (OutC$caseset))
  mod2_OR[n]<-exp(coef(clog_model2))[1]
  S2 <- summary(clog_model2)
  mod2_CI[n] <- paste0('(',exp(S2$coefficients[1,1]-1.96*S2$coefficients[1,3]) %>% round(.,6),'-',exp(S2$coefficients[1,1]+1.96*S2$coefficients[1,3]) %>% round(.,6),')')
  
}



#Make data frame with name of significant features, p-values, OR and CI for CLR models without and without adjustment for possible confounders
Signif_feat_names <- colnames(MB_sign)
P_nominal <- clogit_v[ind_sign]
P_FDR <- adjP_clogit[ind_sign]
P_model2_nominal <- clogit_v2[ind_sign]
P_model2_FDR <- adjP_clogit2[ind_sign]


Signif_feat <- cbind(Signif_feat_names, P_nominal,P_FDR,P_model2_nominal, P_model2_FDR, paste(mod1_OR,mod1_CI, sep=""), paste(mod2_OR,mod2_CI, sep=""))
colnames(Signif_feat) <- c("Feature", "Nominal P univariable","FDR P univariable","Nominal P multivariable","FDR P multivariable","OR (95 % CI) univariable", "OR(95% CI ) multivariable")
Signif_feat

