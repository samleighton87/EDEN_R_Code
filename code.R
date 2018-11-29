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

#enable multicore (windows) which roughly halves time for analysis runs !!!!WINDOWS ONLY!!!!!!!
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

options(max.print=1000000)

#don't use scientific notation (revert back with options(scipen=0)
options(scipen=999)
options(digits = 5)

#load study data
preproc_all = read_csv("eden_preproc_lancet.csv")

preproc_all$Site = as.factor(preproc_all$Site)
preproc_all$Sex = as.factor(preproc_all$Sex)
preproc_all$Ethnicity = as.factor(preproc_all$Ethnicity)
preproc_all$AthiestAgnostic = as.factor(preproc_all$AthiestAgnostic)
preproc_all$Country_Birth = as.factor(preproc_all$Country_Birth)
preproc_all$Fluency = as.factor(preproc_all$Fluency)
preproc_all$Marital_Status = as.factor(preproc_all$Marital_Status)
preproc_all$Living_Status = as.factor(preproc_all$Living_Status)
preproc_all$Housing_Type = as.factor(preproc_all$Housing_Type)
preproc_all$BL_Paid_Employ = as.factor(preproc_all$BL_Paid_Employ)
preproc_all$BL_Vol_Employ = as.factor(preproc_all$BL_Vol_Employ)
preproc_all$BL_Education = as.factor(preproc_all$BL_Education)
preproc_all$BL_EET = as.factor(preproc_all$BL_EET)
preproc_all$Benefits_Any = as.factor(preproc_all$Benefits_Any)
preproc_all$Main_Income_Source = as.factor(preproc_all$Main_Income_Source)
preproc_all$Nature_Baseline_Admission = as.factor(preproc_all$Nature_Baseline_Admission)
preproc_all$AdmittedLast12Month_Acute_Psych = as.factor(preproc_all$AdmittedLast12Month_Acute_Psych)
preproc_all$AdmittedLast12Month_Psych_Rehab = as.factor(preproc_all$AdmittedLast12Month_Psych_Rehab)
preproc_all$AdmittedLast12Month_Long_Stay = as.factor(preproc_all$AdmittedLast12Month_Long_Stay)
preproc_all$AdmittedLast12Month_Emerg_Crisis = as.factor(preproc_all$AdmittedLast12Month_Emerg_Crisis)
preproc_all$AdmittedLast12Month_Gen_Med = as.factor(preproc_all$AdmittedLast12Month_Gen_Med)
preproc_all$CPN_Contact_Last3Mon = as.factor(preproc_all$CPN_Contact_Last3Mon)
preproc_all$Help_FriendRelOutsideHome_Last3Mon = as.factor(preproc_all$Help_FriendRelOutsideHome_Last3Mon)
preproc_all$Help_FriendRelAroundHouse_Last3Mon = as.factor(preproc_all$Help_FriendRelAroundHouse_Last3Mon)
preproc_all$Help_FriendRelPersonal_Last3Mon = as.factor(preproc_all$Help_FriendRelPersonal_Last3Mon)
preproc_all$Help_FriendRelChild_Last3Mon = as.factor(preproc_all$Help_FriendRelChild_Last3Mon)
preproc_all$Help_FriendRelAny_Last3Mon = as.factor(preproc_all$Help_FriendRelAny_Last3Mon)
preproc_all$Previous_Antipsychotic = as.factor(preproc_all$Previous_Antipsychotic)
preproc_all$Past_Drug_Use = as.factor(preproc_all$Past_Drug_Use)
preproc_all$Known_Cannabis = as.factor(preproc_all$Known_Cannabis)
preproc_all$Known_Amphet = as.factor(preproc_all$Known_Amphet)
preproc_all$Known_Coke = as.factor(preproc_all$Known_Coke)
preproc_all$Known_Ecstasy = as.factor(preproc_all$Known_Ecstasy)
preproc_all$Known_Mushrooms = as.factor(preproc_all$Known_Mushrooms)
preproc_all$Known_Opiates = as.factor(preproc_all$Known_Opiates)
preproc_all$Known_LSD = as.factor(preproc_all$Known_LSD)
preproc_all$Known_Khat = as.factor(preproc_all$Known_Khat)
preproc_all$Known_Ketamine = as.factor(preproc_all$Known_Ketamine)
preproc_all$Known_Solvents = as.factor(preproc_all$Known_Solvents)
preproc_all$Possible_Dev_Disorder = as.factor(preproc_all$Possible_Dev_Disorder)
preproc_all$Head_Injury = as.factor(preproc_all$Head_Injury)
preproc_all$Previous_Seconday_Care = as.factor(preproc_all$Previous_Seconday_Care)
preproc_all$PW1_First_Contact = as.factor(preproc_all$PW1_First_Contact)
preproc_all$PW1_Who_Suggested_Care = as.factor(preproc_all$PW1_Who_Suggested_Care)
preproc_all$PW1_Appointment_Attended = as.factor(preproc_all$PW1_Appointment_Attended)
preproc_all$Help_Sought = as.factor(preproc_all$Help_Sought)
preproc_all$Mother_Tongue = as.factor(preproc_all$Mother_Tongue)
preproc_all$Any_Dep_Bipolar = as.factor(preproc_all$Any_Dep_Bipolar)
preproc_all$Any_Schizo_Psych = as.factor(preproc_all$Any_Schizo_Psych)
preproc_all$Contact_Just_Services_Last3Mon = as.factor(preproc_all$Contact_Just_Services_Last3Mon)
preproc_all$BL_Leisure = as.factor(preproc_all$BL_Leisure)
preproc_all$BL_Sport = as.factor(preproc_all$BL_Sport)
preproc_all$BL_Socialising = as.factor(preproc_all$BL_Socialising)
preproc_all$BL_Resting = as.factor(preproc_all$BL_Resting)
preproc_all$BL_Hobbies = as.factor(preproc_all$BL_Hobbies)
preproc_all$BL_Childcare = as.factor(preproc_all$BL_Childcare)
preproc_all$BL_8Hours_Sleep = as.factor(preproc_all$BL_8Hours_Sleep)
preproc_all$BL_Housework = as.factor(preproc_all$BL_Housework)
preproc_all$BL_MostSeriousHarm_How = as.factor(preproc_all$BL_MostSeriousHarm_How)
preproc_all$BL_MostSeriousHarm_Premeditation = as.factor(preproc_all$BL_MostSeriousHarm_Premeditation)
preproc_all$BL_MostSeriousViolence_Victim = as.factor(preproc_all$BL_MostSeriousViolence_Victim)
preproc_all$BL_MostSeriousViolence_VictimGender = as.factor(preproc_all$BL_MostSeriousViolence_VictimGender)
preproc_all$BL_MostSeriousViolence_VictimHarm = as.factor(preproc_all$BL_MostSeriousViolence_VictimHarm)
preproc_all$BL_MostSeriousViolence_VictimInjury = as.factor(preproc_all$BL_MostSeriousViolence_VictimInjury)

#outcomes are all binary
preproc_all$M12_PANSS_Period_Rem = as.factor(preproc_all$M12_PANSS_Period_Rem)
preproc_all$M12_GAF_Binary = as.factor(preproc_all$M12_GAF_Binary)
preproc_all$M12_EET = as.factor(preproc_all$M12_EET)
preproc_all$M12_EQ5D_UK_TTO_Index_Binary = as.factor(preproc_all$M12_EQ5D_UK_TTO_Index_Binary)

preproc_all$FUP_PANSS_Peroid_Rem = as.factor(preproc_all$FUP_PANSS_Peroid_Rem)
preproc_all$FUP_GAF_Binary = as.factor(preproc_all$FUP_GAF_Binary)
preproc_all$FUP_EET = as.factor(preproc_all$FUP_EET)
preproc_all$FUP_EQ5D_UK_TTO_Index_Binary = as.factor(preproc_all$FUP_EQ5D_UK_TTO_Index_Binary)

#get birm
birm_prep = preproc_all[ which(preproc_all$Site=="Birmingham EIS"| preproc_all$Site=="East Anglia Norfolk EIS"), ]
#get others
others_prep = preproc_all[ -which(preproc_all$Site=="Birmingham EIS" | preproc_all$Site=="East Anglia Norfolk EIS"), ]

#get columns I want birm
birm_prep = birm_prep[ ,!(colnames(birm_prep) %in% c("Site"))]
#get columns I want others
others_prep = others_prep[ ,!(colnames(others_prep) %in% c("Site"))]

# set tolerance level in custom method (tol = ...) and reset control2 below
toleranceSam <- function (x, metric, tol = 7.5, maximize) 
{
  index <- 1:nrow(x)
  if (!maximize) {
    best <- min(x[, metric])
    perf <- (x[, metric] - best)/best * 100
    candidates <- index[perf < tol]
    bestIter <- min(candidates)
  }
  else {
    best <- max(x[, metric])
    perf <- (x[, metric] - best)/best * -100
    candidates <- index[perf < tol]
    bestIter <- min(candidates)
  }
  bestIter
}

#change to use one SE Rule - simplest model within one se of best performance Breiman et al. Replace with toleranceSam for custom 3% tolerance
control <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary, selectionFunction ="oneSE")
control2 <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary, selectionFunction ="toleranceSam")
control3 <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary)


#optimise over alpha 0.5 but range of lambda choosing onese option
#if results in model with <=10 predictors, good, if not optimise min tolerance which results in <=10 predictors.
#https://stats.stackexchange.com/questions/138569/why-is-lambda-within-one-standard-error-from-the-minimum-is-a-recommended-valu?rq=1
#https://stats.stackexchange.com/questions/123205/relation-between-the-tuning-parameter-lambda-parameter-estimates-beta-i-a
tuneGrid = expand.grid(alpha=0.5, lambda = 10 ^ seq(-0.3, -5, length = 100))

##
##
##Y1 Rem (for 6 months) built on others and tested on birm
##
#data.frame to model predictors vs Y1 remission
birm_prep_Y1_Rem = birm_prep[ ,!(colnames(birm_prep) %in% c("M12_GAF_Binary","M12_EET","M12_EQ5D_UK_TTO_Index_Binary","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]
others_prep_Y1_Rem = others_prep[ ,!(colnames(others_prep) %in% c("M12_GAF_Binary","M12_EET","M12_EQ5D_UK_TTO_Index_Binary","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]

#omit nas from outcome only
birm_prep_Y1_Rem = birm_prep_Y1_Rem[which(!is.na(birm_prep_Y1_Rem$M12_PANSS_Period_Rem)),]
others_prep_Y1_Rem = others_prep_Y1_Rem[which(!is.na(others_prep_Y1_Rem$M12_PANSS_Period_Rem)),]

#make replicable
set.seed(987)
#build glmnet model with knn imputation on birm data using oneSE rule (optimal)
others_prep_Y1_Rem_glm_mod <- train(M12_PANSS_Period_Rem ~ ., data=others_prep_Y1_Rem, method="glmnet", metric="ROC", tuneGrid = tuneGrid, preProc = c("nzv","zv","center", "scale","knnImpute"), trControl=control, na.action = na.pass)

#Test model on others data applying the same preprocessing as the train (model) object
others_prep_Y1_Rem_glm_result <- predict(others_prep_Y1_Rem_glm_mod, birm_prep_Y1_Rem, type = "prob", na.action = na.pass)

others_prep_Y1_Rem_glm_result_roc = roc(predictor = others_prep_Y1_Rem_glm_result$Yes, response = birm_prep_Y1_Rem$M12_PANSS_Period_Rem=="Yes")

others_prep_Y1_Rem_glm_result_roc$auc
plot(others_prep_Y1_Rem_glm_result_roc)
coords(others_prep_Y1_Rem_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

#make ROC diagram
png("figure1.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(others_prep_Y1_Rem_glm_result_roc, auc.polygon=TRUE, grid=c(0.1, 0.1))
segments(x0=0.7006, x1=0.7006, y0=0.7218, y1=(1-0.7006), col = "#00BFC4", lwd=2)
points(x=0.7006,y=0.7218, pch=4, col="#00BFC4", lwd=2)
abline(a=1,b=-1, lty=2)
text(0.8,0.9,"Sens = 0.722\nSpec = 0.701")
text(0.3,0.2,"AUC = 0.750")
dev.off()

coef(others_prep_Y1_Rem_glm_mod$finalModel, others_prep_Y1_Rem_glm_mod$bestTune$lambda)

#Make odds ratio diagram
#final model coefs
others_coef_Y1_Rem = coef(others_prep_Y1_Rem_glm_mod$finalModel, others_prep_Y1_Rem_glm_mod$bestTune$lambda)
names = dimnames(others_coef_Y1_Rem)[[1]]
others_coef_Y1_Rem = others_coef_Y1_Rem[-1]
names = names[-1]
others_coef_Y1_Rem = others_coef_Y1_Rem[c(9,58,70,83,85,86,101,104,136)]
names = names[c(9,58,70,83,85,86,101,104,136)]
# rename
names = c("Living in Own Home/Parents Home", "Duration of untreated psychosis", "PAS Late Adolescence - Sociability & Withdrawal","PAS General - Highest Level of Functioning","PAS General - Degree of Interest in Life","PAS General - Energy Level", "PANSS P3 - Hallucinatory Behaviour","PANSS P6 - Suspiciousness","GAF Total")
group = ifelse(others_coef_Y1_Rem>0,"+","-")

others_odds_ratio_Y1_Rem = exp(others_coef_Y1_Rem)
df = data.frame(names, others_coef_Y1_Rem, group, others_odds_ratio_Y1_Rem)
df$names <- factor(df$names, levels = df$names[order(df$others_odds_ratio_Y1_Rem)])
png("figure2.png", res = 300, width = 40, height = ((9/10)*15), units = 'cm')
ggplot(df,aes(x=names,y=others_odds_ratio_Y1_Rem,fill=group))+scale_y_continuous(limits = c(0.75, 1.25),breaks = c(0.75, 1, 1.25), trans = scales::log_trans())+geom_bar(stat="identity")+coord_flip()+theme(legend.position="none")+xlab(NULL)+ylab("Odds Ratio (log scale)")+theme(axis.text.y = element_text(face="bold",size=12))+geom_text(aes(label=sprintf("%0.3f", round(others_odds_ratio_Y1_Rem, digits = 3))), y=ifelse(others_coef_Y1_Rem< 0.03 & others_coef_Y1_Rem > 0, others_coef_Y1_Rem+0.017, ifelse(others_coef_Y1_Rem > -0.03 & others_coef_Y1_Rem < 0 ,others_coef_Y1_Rem-0.017,others_coef_Y1_Rem/2)), colour = ifelse(others_coef_Y1_Rem< 0.03 & others_coef_Y1_Rem > 0 ,"#00BFC4", ifelse(others_coef_Y1_Rem> -0.03 & others_coef_Y1_Rem < 0 ,"#F8766D","white")), fontface = "bold") 
dev.off()


#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
others_prep_Y1_Rem_glm_result_auc_null = NULL
for(i in seq (1:10000))
{
  birm_perm = permute(birm_prep_Y1_Rem$M12_PANSS_Period_Rem)
  others_prep_Y1_Rem_glm_result_auc_null = c(others_prep_Y1_Rem_glm_result_auc_null, roc(predictor = others_prep_Y1_Rem_glm_result$Yes, response = birm_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(others_prep_Y1_Rem_glm_result_auc_null >= others_prep_Y1_Rem_glm_result_roc$auc))/10000

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(others_prep_Y1_Rem_glm_result$Yes, birm_prep_Y1_Rem$M12_PANSS_Period_Rem=="Yes")
#95CI = 1.96*SE=1.96*S.D./sqrt(n)
(1.96*sqrt(roc.area.test(others_prep_Y1_Rem_glm_result$Yes, birm_prep_Y1_Rem$M12_PANSS_Period_Rem=="Yes")$var)/sqrt(length(others_prep_Y1_Rem_glm_result$Yes)))

#
#data.frame to test 1 year remission model on birm FUP remission data
#
birm_prep_FUP_Rem = birm_prep[ ,!(colnames(birm_prep) %in% c("M12_PANSS_Period_Rem","M12_GAF_Binary","M12_EET","M12_EQ5D_UK_TTO_Index_Binary","FUP_GAF_Binary","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]

#omit nas from outcome only
birm_prep_FUP_Rem = birm_prep_FUP_Rem[which(!is.na(birm_prep_FUP_Rem$FUP_PANSS_Peroid_Rem)),]

#Test model on others data applying the same preprocessing as the train (model) object
others_prep_Y1_FUP_Rem_glm_result <- predict(others_prep_Y1_Rem_glm_mod, birm_prep_FUP_Rem, type = "prob", na.action = na.pass)

others_prep_Y1_FUP_Rem_glm_result_roc = roc(predictor = others_prep_Y1_FUP_Rem_glm_result$Yes, response = birm_prep_FUP_Rem$FUP_PANSS_Peroid_Rem=="Yes")

others_prep_Y1_FUP_Rem_glm_result_roc$auc

plot(others_prep_Y1_FUP_Rem_glm_result_roc)

coords(others_prep_Y1_FUP_Rem_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

#make ROC diagram
png("figure3.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(others_prep_Y1_FUP_Rem_glm_result_roc, auc.polygon=TRUE, grid=c(0.1, 0.1))
segments(x0=0.6344, x1=0.6344, y0=0.7576, y1=(1-0.6344), col = "#00BFC4", lwd=2)
points(x=0.6344,y=0.7576, pch=4, col="#00BFC4", lwd=2)
abline(a=1,b=-1, lty=2)
text(0.8,0.9,"Sens = 0.758\nSpec = 0.634")
text(0.3,0.2,"AUC = 0.699")
dev.off()

#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
others_prep_Y1_FUP_Rem_glm_result_auc_null = NULL
for(i in seq (1:10000))
{
  birm_perm = permute(birm_prep_FUP_Rem$FUP_PANSS_Peroid_Rem)
  others_prep_Y1_FUP_Rem_glm_result_auc_null = c(others_prep_Y1_FUP_Rem_glm_result_auc_null, roc(predictor = others_prep_Y1_FUP_Rem_glm_result$Yes, response = birm_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(others_prep_Y1_FUP_Rem_glm_result_auc_null >= others_prep_Y1_FUP_Rem_glm_result$auc))/10000

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(others_prep_Y1_FUP_Rem_glm_result$Yes, birm_prep_FUP_Rem$FUP_PANSS_Peroid_Rem=="Yes")
#95CI = 1.96*SE=1.96*S.D./sqrt(n)
(1.96*sqrt(roc.area.test(others_prep_Y1_FUP_Rem_glm_result$Yes, birm_prep_FUP_Rem$FUP_PANSS_Peroid_Rem=="Yes")$var)/sqrt(length(others_prep_Y1_FUP_Rem_glm_result$Yes)))

##
##
##Y1 GAF built on others and tested on birm
##
#data.frame to model predictors vs Y1 GAF
birm_prep_Y1_GAF = birm_prep[ ,!(colnames(birm_prep) %in% c("M12_PANSS_Period_Rem","M12_EET","M12_EQ5D_UK_TTO_Index_Binary","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]
others_prep_Y1_GAF = others_prep[ ,!(colnames(others_prep) %in% c("M12_PANSS_Period_Rem","M12_EET","M12_EQ5D_UK_TTO_Index_Binary","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]

#omit nas from outcome only
birm_prep_Y1_GAF = birm_prep_Y1_GAF[which(!is.na(birm_prep_Y1_GAF$M12_GAF_Binary)),]
others_prep_Y1_GAF = others_prep_Y1_GAF[which(!is.na(others_prep_Y1_GAF$M12_GAF_Binary)),]

set.seed(987)
#optimised tolerance to 1.75 above
others_prep_Y1_GAF_glm_mod <- train(M12_GAF_Binary ~ ., data=others_prep_Y1_GAF, method="glmnet", metric="ROC", tuneGrid = tuneGrid, preProc = c("nzv","zv","center", "scale","knnImpute"), trControl=control2, na.action = na.pass)

#Test model on others data applying the same preprocessing as the train (model) object
others_prep_Y1_GAF_glm_result <- predict(others_prep_Y1_GAF_glm_mod, birm_prep_Y1_GAF, type = "prob", na.action = na.pass)

others_prep_Y1_GAF_glm_result_roc = roc(predictor = others_prep_Y1_GAF_glm_result$Yes, response = birm_prep_Y1_GAF$M12_GAF_Binary=="Yes")
others_prep_Y1_GAF_glm_result_roc$auc
plot(others_prep_Y1_GAF_glm_result_roc)
coords(others_prep_Y1_GAF_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

coef(others_prep_Y1_GAF_glm_mod$finalModel, others_prep_Y1_GAF_glm_mod$bestTune$lambda)

#roc plot
png("figure4.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(others_prep_Y1_GAF_glm_result_roc, auc.polygon=TRUE, grid=c(0.1, 0.1))
segments(x0=0.8479, x1=0.8479, y0=0.5175, y1=(1-0.8479), col = "#00BFC4", lwd=2)
points(x=0.8479,y=0.5175, pch=4, col="#00BFC4", lwd=2)
abline(a=1,b=-1, lty=2)
text(0.8,0.9,"Sens = 0.518\nSpec = 0.848")
text(0.3,0.2,"AUC = 0.744")
dev.off()

#coef diagram
#final model coefs
others_coef_Y1_GAF = coef(others_prep_Y1_GAF_glm_mod$finalModel, others_prep_Y1_GAF_glm_mod$bestTune$lambda)
names = dimnames(others_coef_Y1_GAF)[[1]]
others_coef_Y1_GAF = others_coef_Y1_GAF[-1]
names = names[-1]
#non-zero coefs
others_coef_Y1_GAF = others_coef_Y1_GAF[c(13,21,22,44,69,71,75,84,104,138)]
names = names[c(13,21,22,44,69,71,75,84,104,138)]
# rename
names = c("Qualification Level", "Main Source of Income - Salary/Wage","Main Source of Income - State Benefits","Previous Amphetamine Use","PAS Early Adolescence - Adaption to School","PAS Late Adolescence - Sociability & Withdrawal", "PAS Late Adolescence - Social Sexual Aspects","PAS General - Highest Level of Functioning","PANSS P4 - Excitement","GAF Total")
group = ifelse(others_coef_Y1_GAF>0,"+","-")

others_odds_ratio_Y1_GAF = exp(others_coef_Y1_GAF)
df = data.frame(names, others_coef_Y1_GAF, group, others_odds_ratio_Y1_GAF)
df$names <- factor(df$names, levels = df$names[order(df$others_odds_ratio_Y1_GAF)])
png("figure5.png", res = 300, width = 40, height = 15, units = 'cm')
ggplot(df,aes(x=names,y=others_odds_ratio_Y1_GAF,fill=group))+scale_y_continuous(limits = c(0.75, 1.25),breaks = c(0.75, 1, 1.25), trans = scales::log_trans())+geom_bar(stat="identity")+coord_flip()+theme(legend.position="none")+xlab(NULL)+ylab("Odds Ratio (log scale)")+theme(axis.text.y = element_text(face="bold",size=12))+geom_text(aes(label=sprintf("%0.3f", round(others_odds_ratio_Y1_GAF, digits = 3))), y=ifelse(others_coef_Y1_GAF< 0.03 & others_coef_Y1_GAF > 0, others_coef_Y1_GAF+0.017, ifelse(others_coef_Y1_GAF > -0.03 & others_coef_Y1_GAF < 0 ,others_coef_Y1_GAF-0.017,others_coef_Y1_GAF/2)), colour = ifelse(others_coef_Y1_GAF< 0.03 & others_coef_Y1_GAF > 0 ,"#00BFC4", ifelse(others_coef_Y1_GAF> -0.03 & others_coef_Y1_GAF < 0 ,"#F8766D","white")), fontface = "bold") 
dev.off()


#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
others_prep_Y1_GAF_glm_result_auc_null = NULL
for(i in seq (1:10000))
{
  birm_perm = permute(birm_prep_Y1_GAF$M12_GAF_Binary)
  others_prep_Y1_GAF_glm_result_auc_null = c(others_prep_Y1_GAF_glm_result_auc_null, roc(predictor = others_prep_Y1_GAF_glm_result$Yes, response = birm_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(others_prep_Y1_GAF_glm_result_auc_null >= others_prep_Y1_GAF_glm_result_roc$auc))/10000

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(others_prep_Y1_GAF_glm_result$Yes, birm_prep_Y1_GAF$M12_GAF_Binary=="Yes")
#95CI = 1.96*SE=1.96*S.D./sqrt(n)
(1.96*sqrt(roc.area.test(others_prep_Y1_GAF_glm_result$Yes, birm_prep_Y1_GAF$M12_GAF_Binary=="Yes")$var)/sqrt(length(others_prep_Y1_GAF_glm_result$Yes)))


#
#data.frame to test 1 year GAF model on birm FUP GAF data
#
birm_prep_FUP_GAF = birm_prep[ ,!(colnames(birm_prep) %in% c("M12_PANSS_Period_Rem","M12_GAF_Binary","M12_EET","M12_EQ5D_UK_TTO_Index_Binary","FUP_PANSS_Peroid_Rem","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]

#omit nas from outcome only
birm_prep_FUP_GAF = birm_prep_FUP_GAF[which(!is.na(birm_prep_FUP_GAF$FUP_GAF_Binary)),]


#Test model on others data applying the same preprocessing as the train (model) object
others_prep_Y1_FUP_GAF_glm_result = predict(others_prep_Y1_GAF_glm_mod, birm_prep_FUP_GAF, type = "prob", na.action = na.pass)

others_prep_Y1_FUP_GAF_glm_result_roc = roc(predictor = others_prep_Y1_FUP_GAF_glm_result$Yes, response = birm_prep_FUP_GAF$FUP_GAF_Binary=="Yes")

others_prep_Y1_FUP_GAF_glm_result_roc$auc

plot(others_prep_Y1_FUP_GAF_glm_result_roc)

coords(others_prep_Y1_FUP_GAF_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

#roc plot
png("figure6.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(others_prep_Y1_FUP_GAF_glm_result_roc, auc.polygon=TRUE, grid=c(0.1, 0.1))
segments(x0=0.5781, x1=0.5781, y0=0.6935, y1=(1-0.5781), col = "#00BFC4", lwd=2)
points(x=0.5781,y=0.6935, pch=4, col="#00BFC4", lwd=2)
abline(a=1,b=-1, lty=2)
text(0.8,0.9,"Sens = 0.694\nSpec = 0.578")
text(0.3,0.2,"AUC = 0.615")
dev.off()

#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
others_prep_FUP_Y1_GAF_glm_result_auc_null = NULL
for(i in seq (1:10000))
{
  birm_perm = permute(birm_prep_FUP_GAF$FUP_GAF_Binary)
  others_prep_FUP_Y1_GAF_glm_result_auc_null = c(others_prep_FUP_Y1_GAF_glm_result_auc_null, roc(predictor = others_prep_Y1_FUP_GAF_glm_result$Yes, response = birm_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(others_prep_FUP_Y1_GAF_glm_result_auc_null >= others_prep_Y1_FUP_GAF_glm_result_roc$auc))/10000

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(others_prep_Y1_FUP_GAF_glm_result$Yes, birm_prep_FUP_GAF$FUP_GAF_Binary=="Yes")
#95CI = 1.96*SE=1.96*S.D./sqrt(n)
(1.96*sqrt(roc.area.test(others_prep_Y1_FUP_GAF_glm_result$Yes, birm_prep_FUP_GAF$FUP_GAF_Binary=="Yes")$var)/sqrt(length(others_prep_Y1_FUP_GAF_glm_result$Yes)))


#
##Y1 EET built on others and tested on birm
#
#
birm_prep_Y1_EET = birm_prep[ ,!(colnames(birm_prep) %in% c("M12_PANSS_Period_Rem","M12_GAF_Binary","M12_EQ5D_UK_TTO_Index_Binary","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]
others_prep_Y1_EET = others_prep[ ,!(colnames(others_prep) %in% c("M12_PANSS_Period_Rem","M12_GAF_Binary","M12_EQ5D_UK_TTO_Index_Binary","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]

#omit nas from outcome only
birm_prep_Y1_EET = birm_prep_Y1_EET[which(!is.na(birm_prep_Y1_EET$M12_EET)),]
others_prep_Y1_EET = others_prep_Y1_EET[which(!is.na(others_prep_Y1_EET$M12_EET)),]

set.seed(987)
#build glmnet model with knn imputation on birm data using oneSE rule (optimal)
others_prep_Y1_EET_glm_mod <- train(M12_EET ~ ., data=others_prep_Y1_EET, method="glmnet", metric="ROC", tuneGrid = tuneGrid, preProc = c("nzv","zv","center", "scale","knnImpute"), trControl=control, na.action = na.pass)

others_prep_Y1_EET_glm_result <- predict(others_prep_Y1_EET_glm_mod, birm_prep_Y1_EET, type = "prob", na.action = na.pass)

others_prep_Y1_EET_glm_result_roc = roc(predictor = others_prep_Y1_EET_glm_result$Yes, response = birm_prep_Y1_EET$M12_EET=="Yes")

others_prep_Y1_EET_glm_result_roc$auc

plot(others_prep_Y1_EET_glm_result_roc)

coords(others_prep_Y1_EET_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

#final model coefs
coef(others_prep_Y1_EET_glm_mod$finalModel, others_prep_Y1_EET_glm_mod$bestTune$lambda)

#roc
png("figure7.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(others_prep_Y1_EET_glm_result_roc, auc.polygon=TRUE, grid=c(0.1, 0.1))
segments(x0=0.8432, x1=0.8432, y0=0.5081, y1=(1-0.8432), col = "#00BFC4", lwd=2)
points(x=0.8432,y=0.5081, pch=4, col="#00BFC4", lwd=2)
abline(a=1,b=-1, lty=2)
text(0.8,0.9,"Sens = 0.508\nSpec = 0.843")
text(0.3,0.2,"AUC = 0.725")
dev.off()

#final model coefs figure
others_coef_Y1_EET = coef(others_prep_Y1_EET_glm_mod$finalModel, others_prep_Y1_EET_glm_mod$bestTune$lambda)
names = dimnames(others_coef_Y1_EET)[[1]]
others_coef_Y1_EET = others_coef_Y1_EET[-1]
names = names[-1]
others_coef_Y1_EET = others_coef_Y1_EET[c(13,18,22,45,80,81,139,141)]
names = names[c(13,18,22,45,80,81,139,141)]
# rename
names = c("Qualification Level","In Employment, Education or Training","Main Source of Income - State Benefits","Previous Amphetamine Use", "PAS General - Education", "PAS General - Employed/At School", "GAF Total","GAF Disability Total")
group = ifelse(others_coef_Y1_EET>0,"+","-")

others_odds_ratio_Y1_EET = exp(others_coef_Y1_EET)
df = data.frame(names, others_coef_Y1_EET, group, others_odds_ratio_Y1_EET)
df$names <- factor(df$names, levels = df$names[order(df$others_odds_ratio_Y1_EET)])
png("figure8.png", res = 300, width = 39, height = ((8/10)*15), units = 'cm')
ggplot(df,aes(x=names,y=others_odds_ratio_Y1_EET,fill=group))+scale_y_continuous(limits = c(0.75, 1.25),breaks = c(0.75, 1, 1.25), trans = scales::log_trans())+geom_bar(stat="identity")+coord_flip()+theme(legend.position="none")+xlab(NULL)+ylab("Odds Ratio (odds scale)")+theme(axis.text.y = element_text(face="bold",size=12))+geom_text(aes(label=sprintf("%0.3f", round(others_odds_ratio_Y1_EET, digits = 3))), y=ifelse(others_coef_Y1_EET< 0.03 & others_coef_Y1_EET > 0, others_coef_Y1_EET+0.017, ifelse(others_coef_Y1_EET > -0.03 & others_coef_Y1_EET < 0 ,others_coef_Y1_EET-0.017,others_coef_Y1_EET/2)), colour = ifelse(others_coef_Y1_EET< 0.03 & others_coef_Y1_EET > 0 ,"#00BFC4", ifelse(others_coef_Y1_EET> -0.03 & others_coef_Y1_EET < 0 ,"#F8766D","white")), fontface = "bold") 
dev.off()

#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
others_prep_Y1_EET_glm_result_auc_null = NULL
for(i in seq (1:10000))
{
  birm_perm = permute(birm_prep_Y1_EET$M12_EET)
  others_prep_Y1_EET_glm_result_auc_null = c(others_prep_Y1_EET_glm_result_auc_null, roc(predictor = others_prep_Y1_EET_glm_result$Yes, response = birm_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(others_prep_Y1_EET_glm_result_auc_null >= others_prep_Y1_EET_glm_result_roc$auc))/10000

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(others_prep_Y1_EET_glm_result$Yes, birm_prep_Y1_EET$M12_EET=="Yes")
#95CI = 1.96*SE=1.96*S.D./sqrt(n)
(1.96*sqrt(roc.area.test(others_prep_Y1_EET_glm_result$Yes, birm_prep_Y1_EET$M12_EET=="Yes")$var)/sqrt(length(others_prep_Y1_EET_glm_result$Yes)))

#
#data.frame to test 1 year EET model on birm FUP EET data
#
birm_prep_FUP_EET = birm_prep[ ,!(colnames(birm_prep) %in% c("M12_PANSS_Period_Rem","M12_EET","M12_GAF_Binary","M12_EQ5D_UK_TTO_Index_Binary","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EQ5D_UK_TTO_Index_Binary"))]

#omit nas from outcome only
birm_prep_FUP_EET = birm_prep_FUP_EET[which(!is.na(birm_prep_FUP_EET$FUP_EET)),]

#Test model on others data applying the same preprocessing as the train (model) object
others_prep_Y1_FUP_EET_glm_result <- predict(others_prep_Y1_EET_glm_mod, birm_prep_FUP_EET, type = "prob", na.action = na.pass)

others_prep_Y1_FUP_EET_glm_result_roc = roc(predictor = others_prep_Y1_FUP_EET_glm_result$Yes, response = birm_prep_FUP_EET$FUP_EET=="Yes")

others_prep_Y1_FUP_EET_glm_result_roc$auc

plot(others_prep_Y1_FUP_EET_glm_result_roc)

#roc
png("figure9.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(others_prep_Y1_FUP_EET_glm_result_roc, auc.polygon=TRUE, grid=c(0.1, 0.1))
segments(x0=0.5872, x1=0.5872, y0=0.8481, y1=(1-0.5872), col = "#00BFC4", lwd=2)
points(x=0.5872,y=0.8481, pch=4, col="#00BFC4", lwd=2)
abline(a=1,b=-1, lty=2)
text(0.8,0.9,"Sens = 0.8481\nSpec = 0.587")
text(0.3,0.2,"AUC = 0.739")
dev.off()

#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
others_prep_Y1_FUP_EET_glm_result_auc_null = NULL
for(i in seq (1:10000))
{
  birm_perm = permute(birm_prep_FUP_EET$FUP_EET)
  others_prep_Y1_FUP_EET_glm_result_auc_null = c(others_prep_Y1_FUP_EET_glm_result_auc_null, roc(predictor = others_prep_Y1_FUP_EET_glm_result$Yes, response = birm_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(others_prep_Y1_FUP_EET_glm_result_auc_null >= others_prep_Y1_FUP_EET_glm_result_roc$auc))/10000

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(others_prep_Y1_FUP_EET_glm_result$Yes, birm_prep_FUP_EET$FUP_EET=="Yes")
#95CI = 1.96*SE=1.96*S.D./sqrt(n)
(1.96*sqrt(roc.area.test(others_prep_Y1_FUP_EET_glm_result$Yes, birm_prep_FUP_EET$FUP_EET=="Yes")$var)/sqrt(length(others_prep_Y1_FUP_EET_glm_result$Yes)))

#
##Y1 EQ5D UK TTO Index binary (>median which is 0.848) built on others and tested on birm
#
#
birm_prep_Y1_EQ5D_TTO = birm_prep[ ,!(colnames(birm_prep) %in% c("M12_PANSS_Period_Rem","M12_GAF_Binary","M12_EET","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]
others_prep_Y1_EQ5D_TTO = others_prep[ ,!(colnames(others_prep) %in% c("M12_PANSS_Period_Rem","M12_GAF_Binary","M12_EET","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]

#omit nas from outcome only
birm_prep_Y1_EQ5D_TTO = birm_prep_Y1_EQ5D_TTO[which(!is.na(birm_prep_Y1_EQ5D_TTO$M12_EQ5D_UK_TTO_Index_Binary)),]
others_prep_Y1_EQ5D_TTO = others_prep_Y1_EQ5D_TTO[which(!is.na(others_prep_Y1_EQ5D_TTO$M12_EQ5D_UK_TTO_Index_Binary)),]

set.seed(987)
#build glmnet model with knn imputation on birm data
#oneSE does not result in a sparse model, tolerance 7.5 is optimal
others_prep_Y1_EQ5D_TTO_glm_mod <- train(M12_EQ5D_UK_TTO_Index_Binary ~ ., data=others_prep_Y1_EQ5D_TTO, method="glmnet", metric="ROC", tuneGrid = tuneGrid, preProc = c("nzv","zv","center", "scale","knnImpute"), trControl=control2, na.action = na.pass)

others_prep_Y1_EQ5D_TTO_glm_result <- predict(others_prep_Y1_EQ5D_TTO_glm_mod, birm_prep_Y1_EQ5D_TTO, type = "prob", na.action = na.pass)

others_prep_Y1_EQ5D_TTO_glm_result_roc = roc(predictor = others_prep_Y1_EQ5D_TTO_glm_result$Yes, response = birm_prep_Y1_EQ5D_TTO$M12_EQ5D_UK_TTO_Index_Binary=="Yes")

others_prep_Y1_EQ5D_TTO_glm_result_roc$auc

plot(others_prep_Y1_EQ5D_TTO_glm_result_roc)

coords(others_prep_Y1_EQ5D_TTO_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

#final model coefs
coef(others_prep_Y1_EQ5D_TTO_glm_mod$finalModel, others_prep_Y1_EQ5D_TTO_glm_mod$bestTune$lambda)

#roc
png("figure10.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(others_prep_Y1_EQ5D_TTO_glm_result_roc, auc.polygon=TRUE, grid=c(0.1, 0.1))
segments(x0=0.5249, x1=0.5249, y0=0.8069, y1=(1-0.5249), col = "#00BFC4", lwd=2)
points(x=0.5249,y=0.8069, pch=4, col="#00BFC4", lwd=2)
abline(a=1,b=-1, lty=2)
text(0.8,0.9,"Sens = 0.807\nSpec = 0.525")
text(0.3,0.2,"AUC = 0.702")
dev.off()

#final model coefs figure
others_coef_Y1_EQ5D_TTO = coef(others_prep_Y1_EQ5D_TTO_glm_mod$finalModel, others_prep_Y1_EQ5D_TTO_glm_mod$bestTune$lambda)
names = dimnames(others_coef_Y1_EQ5D_TTO)[[1]]
others_coef_Y1_EQ5D_TTO = others_coef_Y1_EQ5D_TTO[-1]
names = names[-1]
others_coef_Y1_EQ5D_TTO = others_coef_Y1_EQ5D_TTO[c(9,82,84,87,102,136,144,146,162)]
names = names[c(9,82,84,87,102,136,144,146,162)]
# rename
names = c("Living in Own Home/Parents Home","PAS General - Job Change/Interrupted School\n Attendance", "PAS General - Highest Level of Functioning","PAS General - Energy Level", "PANSS P3 - Hallucinatory Behaviour", "Calgary Depression Scale Total","EQ5D - Anxiety/Depression","EQ5D - UK TTO Index", "Never Self-Harmed/Attempted Suicide")
group = ifelse(others_coef_Y1_EQ5D_TTO>0,"+","-")

others_odds_ratio_Y1_EQ5D_TTO = exp(others_coef_Y1_EQ5D_TTO)
df = data.frame(names, others_coef_Y1_EQ5D_TTO, group, others_odds_ratio_Y1_EQ5D_TTO)
df$names <- factor(df$names, levels = df$names[order(df$others_odds_ratio_Y1_EQ5D_TTO)])
png("figure11.png", res = 300, width = 39.5, height =((9/10)* 15), units = 'cm')
ggplot(df,aes(x=names,y=others_odds_ratio_Y1_EQ5D_TTO,fill=group))+scale_y_continuous(limits = c(0.75, 1.25),breaks = c(0.75, 1, 1.25), trans = scales::log_trans())+geom_bar(stat="identity")+coord_flip()+theme(legend.position="none")+xlab(NULL)+ylab("Odds Ratio (log scale)")+theme(axis.text.y = element_text(face="bold",size=12))+geom_text(aes(label=sprintf("%0.3f", round(others_odds_ratio_Y1_EQ5D_TTO, digits = 3))), y=ifelse(others_coef_Y1_EQ5D_TTO< 0.03 & others_coef_Y1_EQ5D_TTO > 0, others_coef_Y1_EQ5D_TTO+0.017, ifelse(others_coef_Y1_EQ5D_TTO > -0.03 & others_coef_Y1_EQ5D_TTO < 0 ,others_coef_Y1_EQ5D_TTO-0.017,others_coef_Y1_EQ5D_TTO/2)), colour = ifelse(others_coef_Y1_EQ5D_TTO< 0.03 & others_coef_Y1_EQ5D_TTO > 0 ,"#00BFC4", ifelse(others_coef_Y1_EQ5D_TTO> -0.03 & others_coef_Y1_EQ5D_TTO < 0 ,"#F8766D","white")), fontface = "bold") 
dev.off()

#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
others_prep_Y1_EQ5D_TTO_glm_result_auc_null = NULL
for(i in seq (1:10000))
{
  birm_perm = permute(birm_prep_Y1_EQ5D_TTO$M12_EQ5D_UK_TTO_Index_Binary)
  others_prep_Y1_EQ5D_TTO_glm_result_auc_null = c(others_prep_Y1_EQ5D_TTO_glm_result_auc_null, roc(predictor = others_prep_Y1_EQ5D_TTO_glm_result$Yes, response = birm_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(others_prep_Y1_EQ5D_TTO_glm_result_auc_null >= others_prep_Y1_EQ5D_TTO_glm_result_roc$auc))/10000

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(others_prep_Y1_EQ5D_TTO_glm_result$Yes, birm_prep_Y1_EQ5D_TTO$M12_EQ5D_UK_TTO_Index_Binary=="Yes")
#95CI = 1.96*SE=1.96*S.D./sqrt(n)
(1.96*sqrt(roc.area.test(others_prep_Y1_EQ5D_TTO_glm_result$Yes, birm_prep_Y1_EQ5D_TTO$M12_EQ5D_UK_TTO_Index_Binary=="Yes")$var)/sqrt(length(others_prep_Y1_EQ5D_TTO_glm_result$Yes)))

#
#data.frame to test 1 year EQ5D TTO model on birm FUP EQ5D TTO data
#
birm_prep_FUP_EQ5D_TTO = birm_prep[ ,!(colnames(birm_prep) %in% c("M12_PANSS_Period_Rem","M12_EET","M12_GAF_Binary","M12_EQ5D_UK_TTO_Index_Binary","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EET"))]

#omit nas from outcome only
birm_prep_FUP_EQ5D_TTO = birm_prep_FUP_EQ5D_TTO[which(!is.na(birm_prep_FUP_EQ5D_TTO$FUP_EQ5D_UK_TTO_Index_Binary)),]

#Test model on others data applying the same preprocessing as the train (model) object
others_prep_Y1_FUP_EQ5D_TTO_glm_result <- predict(others_prep_Y1_EQ5D_TTO_glm_mod, birm_prep_FUP_EQ5D_TTO, type = "prob", na.action = na.pass)

others_prep_Y1_FUP_EQ5D_TTO_glm_result_roc = roc(predictor = others_prep_Y1_FUP_EQ5D_TTO_glm_result$Yes, response = birm_prep_FUP_EQ5D_TTO$FUP_EQ5D_UK_TTO_Index_Binary=="Yes")

others_prep_Y1_FUP_EQ5D_TTO_glm_result_roc$auc

plot(others_prep_Y1_FUP_EQ5D_TTO_glm_result_roc)

coords(others_prep_Y1_FUP_EQ5D_TTO_glm_result_roc,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

#roc
png("figure12.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(others_prep_Y1_FUP_EQ5D_TTO_glm_result_roc, auc.polygon=TRUE, grid=c(0.1, 0.1))
segments(x0=0.2857, x1=0.2857, y0=0.9101, y1=(1-0.2857), col = "#00BFC4", lwd=2)
points(x=0.2857,y=0.9101, pch=4, col="#00BFC4", lwd=2)
abline(a=1,b=-1, lty=2)
text(0.8,0.9,"Sens = 0.9101\nSpec = 0.286")
text(0.3,0.2,"AUC = 0.572")
dev.off()

#permutation test of outcomes https://www.quora.com/How-is-statistical-significance-determined-for-ROC-curves-and-AUC-values
set.seed(987)
others_prep_Y1_FUP_EQ5D_TTO_glm_result_auc_null = NULL
for(i in seq (1:10000))
{
  birm_perm = permute(birm_prep_FUP_EQ5D_TTO$FUP_EQ5D_UK_TTO_Index_Binary)
  others_prep_Y1_FUP_EQ5D_TTO_glm_result_auc_null = c(others_prep_Y1_FUP_EQ5D_TTO_glm_result_auc_null, roc(predictor = others_prep_Y1_FUP_EQ5D_TTO_glm_result$Yes, response = birm_perm=="Yes")$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(others_prep_Y1_FUP_EQ5D_TTO_glm_result_auc_null >= others_prep_Y1_FUP_EQ5D_TTO_glm_result_roc$auc))/10000

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(others_prep_Y1_FUP_EQ5D_TTO_glm_result$Yes, birm_prep_FUP_EQ5D_TTO$FUP_EQ5D_UK_TTO_Index_Binary=="Yes")
#95CI = 1.96*SE=1.96*S.D./sqrt(n)
(1.96*sqrt(roc.area.test(others_prep_Y1_FUP_EQ5D_TTO_glm_result$Yes, birm_prep_FUP_EQ5D_TTO$FUP_EQ5D_UK_TTO_Index_Binary=="Yes")$var)/sqrt(length(others_prep_Y1_FUP_EQ5D_TTO_glm_result$Yes)))

##END
