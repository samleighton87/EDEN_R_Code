#load required libraries
library(readr)
library(caret)
library(doParallel)
library(foreign)
library(haven)
library(psych)
library(RANN)
library(pROC)
library(clinfun)
library(combinat)
library(gtools)
library(DescTools)
library(plyr)
library(matrixStats)

#enable multicore which roughly halves time for analysis runs 
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

options(max.print=1000000)

#don't use scientific notation (revert back with options(scipen=0)
options(scipen=999)
options(digits = 5)

#REMISSION
#
#load study data from your current working directory - setwd() if required
#Set up study data as per my example_remission.csv template file
data_Y1_Rem_OPUS_all = read_csv("OPUS/example_remission.csv")
#Outcome is factor
data_Y1_Rem_OPUS_all$M12_PANSS_Period_Rem = as.factor(data_Y1_Rem_OPUS_all$M12_PANSS_Period_Rem)
#Remove any rows with missing outcome columns (keeping missing data in predictor columns)
data_Y1_Rem_OPUS_all = data_Y1_Rem_OPUS_all[which(!is.na(data_Y1_Rem_OPUS_all$M12_PANSS_Period_Rem)),]
#Temporarily remove outcome column to standardise the predictor columns
data_Y1_Rem_OPUS = data_Y1_Rem_OPUS_all[ ,!(colnames(data_Y1_Rem_OPUS_all) %in% c("M12_PANSS_Period_Rem"))]
#standardise the columns before building model so that it doesn't matter if our variables are scaled differently
preProcValues = preProcess(data_Y1_Rem_OPUS, method = c("center", "scale"))
data_Y1_Rem_OPUS_Stand = predict(preProcValues, data_Y1_Rem_OPUS)
#Add factor outcome back in
data_Y1_Rem_OPUS_Stand$M12_PANSS_Period_Rem = data_Y1_Rem_OPUS_all$M12_PANSS_Period_Rem

#load my glm model which predicts 1 year PANSS period remission (mod_Y1_Rem_EDEN_all) from current wd
load("OPUS/mod_Y1_Rem_OPUS.rda")
#Have a look at the model coefficients if you want
summary(mod_Y1_Rem_EDEN_all)

#Test our EDEN model on your OPUS data with the same preprocessing as I did (knnImpute, k=5) - information held in model object
Y1_Rem_OPUS_glm_result <- predict(mod_Y1_Rem_EDEN_all, data_Y1_Rem_OPUS_Stand, type = "prob", na.action = na.pass)

#Make a ROC curve object
Y1_Rem_OPUS_glm_result_roc = roc(predictor = Y1_Rem_OPUS_glm_result$Yes, response = data_Y1_Rem_OPUS_Stand$M12_PANSS_Period_Rem=="Yes")

plot(Y1_Rem_OPUS_glm_result_roc)

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(Y1_Rem_OPUS_glm_result$Yes, data_Y1_Rem_OPUS_Stand$M12_PANSS_Period_Rem=="Yes")
#95CI = 1.96*SE
1.96*sqrt(roc.area.test(Y1_Rem_OPUS_glm_result$Yes, data_Y1_Rem_OPUS_Stand$M12_PANSS_Period_Rem=="Yes")$var)

#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
Y1_Rem_OPUS_glm_result_auc_null = NULL
#takes about a minute
for(i in seq (1:10001))
{
  Y1_Rem_OPUS_perm = permute(data_Y1_Rem_OPUS_Stand$M12_PANSS_Period_Rem)
  Y1_Rem_OPUS_glm_result_auc_null = c(Y1_Rem_OPUS_glm_result_auc_null, roc(predictor = Y1_Rem_OPUS_glm_result$Yes, response = Y1_Rem_OPUS_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_Rem_OPUS_glm_result_auc_null >= Y1_Rem_OPUS_glm_result_roc$auc))/10001

#Get other metrics for the point on the roc curve equivalent to Youden's index
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
set.seed(987)
#Doesn't work with test data as not enough data points
ci.coords(Y1_Rem_OPUS_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

#GAF
#
#load study data from your current working directory - setwd() if required
#Set up study data as per my example_GAF.csv template file
data_Y1_GAF_OPUS_all = read_csv("OPUS/example_GAF.csv")
#Outcome is factor
data_Y1_GAF_OPUS_all$M12_GAF_Binary = as.factor(data_Y1_GAF_OPUS_all$M12_GAF_Binary)
#Remove any rows with missing outcome columns (keeping missing data in predictor columns)
data_Y1_GAF_OPUS_all = data_Y1_GAF_OPUS_all[which(!is.na(data_Y1_GAF_OPUS_all$M12_GAF_Binary)),]
#Temporarily remove outcome column to standardise the predictor columns
data_Y1_GAF_OPUS = data_Y1_GAF_OPUS_all[ ,!(colnames(data_Y1_GAF_OPUS_all) %in% c("M12_GAF_Binary"))]
#standardise the columns before building model so that it doesn't matter if our variables are scaled differently
preProcValues = preProcess(data_Y1_GAF_OPUS, method = c("center", "scale"))
data_Y1_GAF_OPUS_Stand = predict(preProcValues, data_Y1_GAF_OPUS)
#Add factor outcome back in
data_Y1_GAF_OPUS_Stand$M12_GAF_Binary = data_Y1_GAF_OPUS_all$M12_GAF_Binary

#load my glm model which predicts 1 year PANSS period remission (mod_Y1_GAF_EDEN_all) from current wd
load("OPUS/mod_Y1_GAF_OPUS.rda")
#Have a look at the model coefficients if you want
summary(mod_Y1_GAF_EDEN_all)

#Test our EDEN model on your OPUS data with the same preprocessing as I did (knnImpute, k=5) - information held in model object
Y1_GAF_OPUS_glm_result <- predict(mod_Y1_GAF_EDEN_all, data_Y1_GAF_OPUS_Stand, type = "prob", na.action = na.pass)

#Make a ROC curve object
Y1_GAF_OPUS_glm_result_roc = roc(predictor = Y1_GAF_OPUS_glm_result$Yes, response = data_Y1_GAF_OPUS_Stand$M12_GAF_Binary=="Yes")

plot(Y1_GAF_OPUS_glm_result_roc)

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(Y1_GAF_OPUS_glm_result$Yes, data_Y1_GAF_OPUS_Stand$M12_GAF_Binary=="Yes")
#95CI = 1.96*SE
1.96*sqrt(roc.area.test(Y1_GAF_OPUS_glm_result$Yes, data_Y1_GAF_OPUS_Stand$M12_GAF_Binary=="Yes")$var)

#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
Y1_GAF_OPUS_glm_result_auc_null = NULL
#takes about a minute
for(i in seq (1:10001))
{
  Y1_GAF_OPUS_perm = permute(data_Y1_GAF_OPUS_Stand$M12_GAF_Binary)
  Y1_GAF_OPUS_glm_result_auc_null = c(Y1_GAF_OPUS_glm_result_auc_null, roc(predictor = Y1_GAF_OPUS_glm_result$Yes, response = Y1_GAF_OPUS_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_GAF_OPUS_glm_result_auc_null >= Y1_GAF_OPUS_glm_result_roc$auc))/10001

#Get other metrics for the point on the roc curve equivalent to Youden's index
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
set.seed(987)
#Doesn't work with test data as not enough data points
ci.coords(Y1_GAF_OPUS_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

#EET
#
#load study data from your current working directory - setwd() if required
#Set up study data as per my example_EET.csv template file
data_Y1_EET_OPUS_all = read_csv("OPUS/example_EET.csv")
#Outcome is factor
data_Y1_EET_OPUS_all$M12_EET = as.factor(data_Y1_EET_OPUS_all$M12_EET)
#Remove any rows with missing outcome columns (keeping missing data in predictor columns)
data_Y1_EET_OPUS_all = data_Y1_EET_OPUS_all[which(!is.na(data_Y1_EET_OPUS_all$M12_EET)),]
#Temporarily remove outcome column to standardise the predictor columns
data_Y1_EET_OPUS = data_Y1_EET_OPUS_all[ ,!(colnames(data_Y1_EET_OPUS_all) %in% c("M12_EET"))]
#standardise the columns before building model so that it doesn't matter if our variables are scaled differently
preProcValues = preProcess(data_Y1_EET_OPUS, method = c("center", "scale"))
data_Y1_EET_OPUS_Stand = predict(preProcValues, data_Y1_EET_OPUS)
#Add factor outcome back in
data_Y1_EET_OPUS_Stand$M12_EET = data_Y1_EET_OPUS_all$M12_EET

#load my glm model which predicts 1 year PANSS period remission (mod_Y1_EET_EDEN_all) from current wd
load("OPUS/mod_Y1_EET_OPUS.rda")
#Have a look at the model coefficients if you want
summary(mod_Y1_EET_EDEN_all)

#Test our EDEN model on your OPUS data with the same preprocessing as I did (knnImpute, k=5) - information held in model object
Y1_EET_OPUS_glm_result <- predict(mod_Y1_EET_EDEN_all, data_Y1_EET_OPUS_Stand, type = "prob", na.action = na.pass)

#Make a ROC curve object
Y1_EET_OPUS_glm_result_roc = roc(predictor = Y1_EET_OPUS_glm_result$Yes, response = data_Y1_EET_OPUS_Stand$M12_EET=="Yes")

plot(Y1_EET_OPUS_glm_result_roc)

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(Y1_EET_OPUS_glm_result$Yes, data_Y1_EET_OPUS_Stand$M12_EET=="Yes")
#95CI = 1.96*SE
1.96*sqrt(roc.area.test(Y1_EET_OPUS_glm_result$Yes, data_Y1_EET_OPUS_Stand$M12_EET=="Yes")$var)

#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
Y1_EET_OPUS_glm_result_auc_null = NULL
#takes about a minute
for(i in seq (1:10001))
{
  Y1_EET_OPUS_perm = permute(data_Y1_EET_OPUS_Stand$M12_EET)
  Y1_EET_OPUS_glm_result_auc_null = c(Y1_EET_OPUS_glm_result_auc_null, roc(predictor = Y1_EET_OPUS_glm_result$Yes, response = Y1_EET_OPUS_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_EET_OPUS_glm_result_auc_null >= Y1_EET_OPUS_glm_result_roc$auc))/10001

#Get other metrics for the point on the roc curve equivalent to Youden's index
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
set.seed(987)
#Doesn't work with test data as not enough data points
ci.coords(Y1_EET_OPUS_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

#QoL
#
#load study data from your current working directory - setwd() if required
#Set up study data as per my example_QoL.csv template file
data_Y1_QoL_OPUS_all = read_csv("OPUS/example_QoL.csv")
#Outcome is factor
data_Y1_QoL_OPUS_all$M12_QoL = as.factor(data_Y1_QoL_OPUS_all$M12_QoL)
#Remove any rows with missing outcome columns (keeping missing data in predictor columns)
data_Y1_QoL_OPUS_all = data_Y1_QoL_OPUS_all[which(!is.na(data_Y1_QoL_OPUS_all$M12_QoL)),]
#Temporarily remove outcome column to standardise the predictor columns
data_Y1_QoL_OPUS = data_Y1_QoL_OPUS_all[ ,!(colnames(data_Y1_QoL_OPUS_all) %in% c("M12_QoL"))]
#standardise the columns before building model so that it doesn't matter if our variables are scaled differently
preProcValues = preProcess(data_Y1_QoL_OPUS, method = c("center", "scale"))
data_Y1_QoL_OPUS_Stand = predict(preProcValues, data_Y1_QoL_OPUS)
#Add factor outcome back in
data_Y1_QoL_OPUS_Stand$M12_QoL = data_Y1_QoL_OPUS_all$M12_QoL

#load my glm model which predicts 1 year PANSS period remission (mod_Y1_QoL_EDEN_all) from current wd
load("OPUS/mod_Y1_QoL_OPUS.rda")
#Have a look at the model coefficients if you want
summary(mod_Y1_QoL_EDEN_all)

#Test our EDEN model on your OPUS data with the same preprocessing as I did (knnImpute, k=5) - information held in model object
Y1_QoL_OPUS_glm_result <- predict(mod_Y1_QoL_EDEN_all, data_Y1_QoL_OPUS_Stand, type = "prob", na.action = na.pass)

#Make a ROC curve object
Y1_QoL_OPUS_glm_result_roc = roc(predictor = Y1_QoL_OPUS_glm_result$Yes, response = data_Y1_QoL_OPUS_Stand$M12_QoL=="Yes")

plot(Y1_QoL_OPUS_glm_result_roc)

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(Y1_QoL_OPUS_glm_result$Yes, data_Y1_QoL_OPUS_Stand$M12_QoL=="Yes")
#95CI = 1.96*SE
1.96*sqrt(roc.area.test(Y1_QoL_OPUS_glm_result$Yes, data_Y1_QoL_OPUS_Stand$M12_QoL=="Yes")$var)

#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
Y1_QoL_OPUS_glm_result_auc_null = NULL
#takes about a minute
for(i in seq (1:10001))
{
  Y1_QoL_OPUS_perm = permute(data_Y1_QoL_OPUS_Stand$M12_QoL)
  Y1_QoL_OPUS_glm_result_auc_null = c(Y1_QoL_OPUS_glm_result_auc_null, roc(predictor = Y1_QoL_OPUS_glm_result$Yes, response = Y1_QoL_OPUS_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_QoL_OPUS_glm_result_auc_null >= Y1_QoL_OPUS_glm_result_roc$auc))/10001

#Get other metrics for the point on the roc curve equivalent to Youden's index
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
set.seed(987)
#Doesn't work with test data as not enough data points
ci.coords(Y1_QoL_OPUS_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))
