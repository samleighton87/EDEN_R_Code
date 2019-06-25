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

#load study data from your current working directory - setwd() if required
#Set up study data as per my exampleData.csv template file
data_Y1_Rem_CRISPFEP_all = read_csv("C:/Users/sam_l/Desktop/EDEN_Pred/Manuscript/NestedCV/CRISPFEP/crispfep.csv")
#Define factors
data_Y1_Rem_CRISPFEP_all$M12_PANSS_Period_Rem = as.factor(data_Y1_Rem_CRISPFEP_all$M12_PANSS_Period_Rem)
data_Y1_Rem_CRISPFEP_all$Site = as.factor(data_Y1_Rem_CRISPFEP_all$Site)

#Remove any rows with missing outcome columns (keeping missing data in predictor columns)
data_Y1_Rem_CRISPFEP_all = data_Y1_Rem_CRISPFEP_all[which(!is.na(data_Y1_Rem_CRISPFEP_all$M12_PANSS_Period_Rem)),]

#take out outcome and site
data_Y1_Rem_CRISPFEP = data_Y1_Rem_CRISPFEP_all[ ,!(colnames(data_Y1_Rem_CRISPFEP_all) %in% c("M12_PANSS_Period_Rem", "Site"))]

#standardise the columns before trying model so that it doesn't matter if our variables are scaled differently
preProcValues = preProcess(data_Y1_Rem_CRISPFEP, method = c("center", "scale"))
data_Y1_Rem_CRISPFEP_Stand = predict(preProcValues, data_Y1_Rem_CRISPFEP)
#Add factor outcome back in
data_Y1_Rem_CRISPFEP_Stand$M12_PANSS_Period_Rem = data_Y1_Rem_CRISPFEP_all$M12_PANSS_Period_Rem

#Test our EDEN model on Scottish data with the same preprocessing as I did (knnImpute, k=5) - information held in model object
Y1_Rem_CRISPFEP_glm_result <- predict(mod_Y1_Rem_EDEN_CF_all, data_Y1_Rem_CRISPFEP_Stand, type = "prob", na.action = na.pass)

#Make a ROC curve object
Y1_Rem_CRISPFEP_glm_result_roc = roc(predictor = Y1_Rem_CRISPFEP_glm_result$Yes, response = data_Y1_Rem_CRISPFEP_Stand$M12_PANSS_Period_Rem=="Yes", ci=T)

plot(Y1_Rem_CRISPFEP_glm_result_roc)

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(Y1_Rem_CRISPFEP_glm_result$Yes, data_Y1_Rem_CRISPFEP_Stand$M12_PANSS_Period_Rem=="Yes")
#95CI = 1.96*SE
(roc.area.test(Y1_Rem_CRISPFEP_glm_result$Yes, data_Y1_Rem_CRISPFEP_Stand$M12_PANSS_Period_Rem=="Yes")$area) - (1.96*sqrt(roc.area.test(Y1_Rem_CRISPFEP_glm_result$Yes, data_Y1_Rem_CRISPFEP_Stand$M12_PANSS_Period_Rem=="Yes")$var))
(roc.area.test(Y1_Rem_CRISPFEP_glm_result$Yes, data_Y1_Rem_CRISPFEP_Stand$M12_PANSS_Period_Rem=="Yes")$area) + (1.96*sqrt(roc.area.test(Y1_Rem_CRISPFEP_glm_result$Yes, data_Y1_Rem_CRISPFEP_Stand$M12_PANSS_Period_Rem=="Yes")$var))

#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
Y1_Rem_CRISPFEP_glm_result_auc_null = NULL
#takes about a minute
for(i in seq (1:10001))
{
  Y1_Rem_CRISPFEP_perm = permute(data_Y1_Rem_CRISPFEP_Stand$M12_PANSS_Period_Rem)
  Y1_Rem_CRISPFEP_glm_result_auc_null = c(Y1_Rem_CRISPFEP_glm_result_auc_null, roc(predictor = Y1_Rem_CRISPFEP_glm_result$Yes, response = Y1_Rem_CRISPFEP_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_Rem_CRISPFEP_glm_result_auc_null >= Y1_Rem_CRISPFEP_glm_result_roc$auc))/10001

#Get other metrics for the point on the roc curve equivalent to Youden's index
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
set.seed(987)
ci.coords(Y1_Rem_CRISPFEP_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

#make ROC diagram
png("figure2A.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(Y1_Rem_CRISPFEP_glm_result_roc, auc.polygon=TRUE, grid=c(0.1, 0.1))
text(0.3,0.3,"AUC = 0·680\n(0·587, 0·773)\np=0·0004")
dev.off()

#load study data from your current working directory - setwd() if required
data_Y1_EET_CRISPFEP_all = read_csv("C:/Users/sam_l/Desktop/EDEN_Pred/Manuscript/NestedCV/CRISPFEP/crispfep2.csv")
#Define factors
data_Y1_EET_CRISPFEP_all$M12_EET = as.factor(data_Y1_EET_CRISPFEP_all$M12_EET)
data_Y1_EET_CRISPFEP_all$Site = as.factor(data_Y1_EET_CRISPFEP_all$Site)

#Remove any rows with missing outcome columns (keeping missing data in predictor columns)
data_Y1_EET_CRISPFEP_all = data_Y1_EET_CRISPFEP_all[which(!is.na(data_Y1_EET_CRISPFEP_all$M12_EET)),]

#take out outcome and site
data_Y1_EET_CRISPFEP = data_Y1_EET_CRISPFEP_all[ ,!(colnames(data_Y1_EET_CRISPFEP_all) %in% c("M12_EET", "Site"))]

#standardise the columns before trying model so that it doesn't matter if our variables are scaled differently
preProcValues = preProcess(data_Y1_EET_CRISPFEP, method = c("center", "scale"))
data_Y1_EET_CRISPFEP_Stand = predict(preProcValues, data_Y1_EET_CRISPFEP)
#Add factor outcome back in
data_Y1_EET_CRISPFEP_Stand$M12_EET = data_Y1_EET_CRISPFEP_all$M12_EET

#Test our EDEN model on Scottish data with the same preprocessing as I did (knnImpute, k=5) - information held in model object
Y1_EET_CRISPFEP_glm_result <- predict(mod_Y1_EET_EDEN_CF_all, data_Y1_EET_CRISPFEP_Stand, type = "prob", na.action = na.pass)

#Make a ROC curve object
Y1_EET_CRISPFEP_glm_result_roc = roc(predictor = Y1_EET_CRISPFEP_glm_result$Yes, response = data_Y1_EET_CRISPFEP_Stand$M12_EET=="Yes")

plot(Y1_EET_CRISPFEP_glm_result_roc)

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(Y1_EET_CRISPFEP_glm_result$Yes, data_Y1_EET_CRISPFEP_Stand$M12_EET=="Yes")
#95CI = 1.96*SE
(roc.area.test(Y1_EET_CRISPFEP_glm_result$Yes, data_Y1_EET_CRISPFEP_Stand$M12_EET=="Yes")$area) - (1.96*sqrt(roc.area.test(Y1_EET_CRISPFEP_glm_result$Yes, data_Y1_EET_CRISPFEP_Stand$M12_EET=="Yes")$var))
(roc.area.test(Y1_EET_CRISPFEP_glm_result$Yes, data_Y1_EET_CRISPFEP_Stand$M12_EET=="Yes")$area) + (1.96*sqrt(roc.area.test(Y1_EET_CRISPFEP_glm_result$Yes, data_Y1_EET_CRISPFEP_Stand$M12_EET=="Yes")$var))

#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
Y1_EET_CRISPFEP_glm_result_auc_null = NULL
#takes about a minute
for(i in seq (1:10001))
{
  Y1_EET_CRISPFEP_perm = permute(data_Y1_EET_CRISPFEP_Stand$M12_EET)
  Y1_EET_CRISPFEP_glm_result_auc_null = c(Y1_EET_CRISPFEP_glm_result_auc_null, roc(predictor = Y1_EET_CRISPFEP_glm_result$Yes, response = Y1_EET_CRISPFEP_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_EET_CRISPFEP_glm_result_auc_null >= Y1_EET_CRISPFEP_glm_result_roc$auc))/10001

#Get other metrics for the point on the roc curve equivalent to Youden's index
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
set.seed(987)
ci.coords(Y1_EET_CRISPFEP_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

#make ROC diagram
png("figure2B.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(Y1_EET_CRISPFEP_glm_result_roc, auc.polygon=TRUE, grid=c(0.1, 0.1))
text(0.3,0.3,"AUC = 0·867\n(0·805, 0·930)\np<0·0001")
dev.off()

#load study data from your current working directory - setwd() if required
data_Y1_QoL_CRISPFEP_all = read_csv("C:/Users/sam_l/Desktop/EDEN_Pred/Manuscript/NestedCV/CRISPFEP/crispfep3.csv")
#Define factors
data_Y1_QoL_CRISPFEP_all$M12_QoL = as.factor(data_Y1_QoL_CRISPFEP_all$M12_QoL)
data_Y1_QoL_CRISPFEP_all$Site = as.factor(data_Y1_QoL_CRISPFEP_all$Site)

#Remove any rows with missing outcome columns (keeping missing data in predictor columns)
data_Y1_QoL_CRISPFEP_all = data_Y1_QoL_CRISPFEP_all[which(!is.na(data_Y1_QoL_CRISPFEP_all$M12_QoL)),]

#take out outcome and site
data_Y1_QoL_CRISPFEP = data_Y1_QoL_CRISPFEP_all[ ,!(colnames(data_Y1_QoL_CRISPFEP_all) %in% c("M12_QoL", "Site"))]

#standardise the columns before trying model so that it doesn't matter if our variables are scaled differently
preProcValues = preProcess(data_Y1_QoL_CRISPFEP, method = c("center", "scale"))
data_Y1_QoL_CRISPFEP_Stand = predict(preProcValues, data_Y1_QoL_CRISPFEP)
#Add factor outcome back in
data_Y1_QoL_CRISPFEP_Stand$M12_QoL = data_Y1_QoL_CRISPFEP_all$M12_QoL

#Test our EDEN model on your Scottish data with the same preprocessing as I did (knnImpute, k=5) - information held in model object
Y1_QoL_CRISPFEP_glm_result <- predict(mod_Y1_QoL_EDEN_CF_all, data_Y1_QoL_CRISPFEP_Stand, type = "prob", na.action = na.pass)

#Make a ROC curve object
Y1_QoL_CRISPFEP_glm_result_roc = roc(predictor = Y1_QoL_CRISPFEP_glm_result$Yes, response = data_Y1_QoL_CRISPFEP_Stand$M12_QoL=="Yes")

plot(Y1_QoL_CRISPFEP_glm_result_roc)

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(Y1_QoL_CRISPFEP_glm_result$Yes, data_Y1_QoL_CRISPFEP_Stand$M12_QoL=="Yes")
#95CI = 1.96*SE
(roc.area.test(Y1_QoL_CRISPFEP_glm_result$Yes, data_Y1_QoL_CRISPFEP_Stand$M12_QoL=="Yes")$area) - (1.96*sqrt(roc.area.test(Y1_QoL_CRISPFEP_glm_result$Yes, data_Y1_QoL_CRISPFEP_Stand$M12_QoL=="Yes")$var))
(roc.area.test(Y1_QoL_CRISPFEP_glm_result$Yes, data_Y1_QoL_CRISPFEP_Stand$M12_QoL=="Yes")$area) + (1.96*sqrt(roc.area.test(Y1_QoL_CRISPFEP_glm_result$Yes, data_Y1_QoL_CRISPFEP_Stand$M12_QoL=="Yes")$var))

#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
Y1_QoL_CRISPFEP_glm_result_auc_null = NULL
#takes about a minute
for(i in seq (1:10001))
{
  Y1_QoL_CRISPFEP_perm = permute(data_Y1_QoL_CRISPFEP_Stand$M12_QoL)
  Y1_QoL_CRISPFEP_glm_result_auc_null = c(Y1_QoL_CRISPFEP_glm_result_auc_null, roc(predictor = Y1_QoL_CRISPFEP_glm_result$Yes, response = Y1_QoL_CRISPFEP_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_QoL_CRISPFEP_glm_result_auc_null >= Y1_QoL_CRISPFEP_glm_result_roc$auc))/10001

#Get other metrics for the point on the roc curve equivalent to Youden's index
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
set.seed(987)
ci.coords(Y1_QoL_CRISPFEP_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

#make ROC diagram
png("figure2C.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(Y1_QoL_CRISPFEP_glm_result_roc, auc.polygon=TRUE, grid=c(0.1, 0.1))
text(0.3,0.3,"AUC = 0·679\n(0·522, 0·836)\np=0·0368")
dev.off()
