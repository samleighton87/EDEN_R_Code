# METHODS
#
# -Dummy code categorical variables.
# -Participants with missing outcome excluded.
# -Predictor variables with >20% missing data excluded (addressing previous criticism)
# -Predictor variables with zero or near zero variance excluded.
#
# -Model training 10 fold cross validation over a grid of alpha and lambda 
# choosing sparsest model within 1 standard error of best training performance 
#
# -After each cross validation split, predictor variables are standardised and 
# missing data is imputed via single imputation using k nearest neighbours (k=5).
#
# -Internal-External validation via a leave one-site out cross-validation. 14 sites 
# so each loop, leave one site out, train model on remaining 13 sites as described above, 
# trained model tested on left out one site. Repeat 14 times, combining predictions to 
# get an average ROC AUC performance for all 14 models tested on each left out site. 
# Confirm significance by permutation testing. Look at stability of feature selection.
#
# -Get shared predictor variables from among the top predictor variables in LOSOCV
# (those shared across all 14 LOSOCV models), rebuild models with only shared variables
# using GLM (no improvement with more complex models)
#
# -Look at correlation between outcomes, look at correlation between model probability
# outputs, look at performance of models on each outcome.


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

#https://github.com/nogueirs/JMLR2018
getStability <- function(X,alpha=0.05) {
  ## the input X is a binary matrix of size M*d where:
  ## M is the number of bootstrap replicates
  ## d is the total number of features
  ## alpha is the level of significance (e.g. if alpha=0.05, we will get 95% confidence intervals)
  ## it's an optional argument and is set to 5% by default
  ### first we compute the stability
  
  M<-nrow(X)
  d<-ncol(X)
  hatPF<-colMeans(X)
  kbar<-sum(hatPF)
  v_rand=(kbar/d)*(1-kbar/d)
  stability<-1-(M/(M-1))*mean(hatPF*(1-hatPF))/v_rand ## this is the stability estimate
  
  ## then we compute the variance of the estimate
  ki<-rowSums(X)
  phi_i<-rep(0,M)
  for(i in 1:M){ 
    phi_i[i]<-(1/v_rand)*((1/d)*sum(X[i,]*hatPF)-(ki[i]*kbar)/d^2-(stability/2)*((2*kbar*ki[i])/d^2-ki[i]/d-kbar/d+1))
  }
  phi_bar=mean(phi_i)
  var_stab=(4/M^2)*sum((phi_i-phi_bar)^2) ## this is the variance of the stability estimate
  
  ## then we calculate lower and upper limits of the confidence intervals
  z<-qnorm(1-alpha/2) # this is the standard normal cumulative inverse at a level 1-alpha/2
  upper<-stability+z*sqrt(var_stab) ## the upper bound of the (1-alpha) confidence interval
  lower<-stability-z*sqrt(var_stab) ## the lower bound of the (1-alpha) confidence interval
  
  return(list("stability"=stability,"variance"=var_stab,"lower"=lower,"upper"=upper))
  
}

#enable multicore which roughly halves time for analysis runs 
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

options(max.print=1000000)

#don't use scientific notation (revert back with options(scipen=0)
options(scipen=999)
options(digits = 5)

#load study data
preproc_all = read_csv("eden_preproc.csv")

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

#Not using four year outcomes as too high percentage missing data
preproc_all$FUP_PANSS_Peroid_Rem = as.factor(preproc_all$FUP_PANSS_Peroid_Rem)
preproc_all$FUP_GAF_Binary = as.factor(preproc_all$FUP_GAF_Binary)
preproc_all$FUP_EET = as.factor(preproc_all$FUP_EET)
preproc_all$FUP_EQ5D_UK_TTO_Index_Binary = as.factor(preproc_all$FUP_EQ5D_UK_TTO_Index_Binary)

#Calculate Correlation Phi between each of the 4 binary outcomes - this is what happens when you compute Pearsons for binary variables
EDEN_Corr = preproc_all[ ,(colnames(preproc_all) %in% c("M12_PANSS_Period_Rem","M12_GAF_Binary","M12_EET","M12_EQ5D_UK_TTO_Index_Binary"))]
EDEN_Corr = data.matrix(EDEN_Corr)
print(corr.test(EDEN_Corr), short = F)

#Remission Y1 data
prep_Y1_Rem = preproc_all[ ,!(colnames(preproc_all) %in% c("Site","M12_PANSS_Period_Rem","M12_GAF_Binary","M12_EET","M12_EQ5D_UK_TTO_Index_Binary","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]
#Dummy code (not outcome or site)
dummies = dummyVars(~ ., data = prep_Y1_Rem, fullRank = T)
prep_Y1_Rem <- data.frame(predict(dummies, newdata = prep_Y1_Rem))
#Add factor site and outcome back in
prep_Y1_Rem$Site = preproc_all$Site
prep_Y1_Rem$M12_PANSS_Period_Rem = preproc_all$M12_PANSS_Period_Rem
#remove na from outcome
prep_Y1_Rem = prep_Y1_Rem[which(!is.na(prep_Y1_Rem$M12_PANSS_Period_Rem)),]
#remove columns with more than 20% missing data
prep_Y1_Rem = prep_Y1_Rem[, colMeans(is.na(prep_Y1_Rem)) <= .2]
#remove zero and near zero variance columns
nzv_cols <- nearZeroVar(prep_Y1_Rem)
if(length(nzv_cols) > 0) prep_Y1_Rem <- prep_Y1_Rem[, -nzv_cols]

#Get shared predictor columns data for OPUS external validation
prep_Y1_Rem_OPUS = prep_Y1_Rem[,c(92,125,54,66,14,91,100,
                                  10,149,45,136,
                                  68,112,122,21,
                                  33,17,152)]
#Rename the columns to something easier
colnames(prep_Y1_Rem_OPUS) = c("P3","GAF","DUP","PASClientLASocWith","QualLevel","P2","N4",
                               "HousingOwnParents","SelfHarmKnifeRazor","LSD","AnyLeisure",
                               "PASClientLASocSex","G9","InsightSymptoms","IncomeSalary",
                               "HelpHouseFriendRel","InEducation","NeverSelfHarmed")
#standardise the columns before building model
preProcValues = preProcess(prep_Y1_Rem_OPUS, method = c("center", "scale"))
prep_Y1_Rem_OPUS_Stand = predict(preProcValues, prep_Y1_Rem_OPUS)
#Add factor site and outcome back in
prep_Y1_Rem_OPUS_Stand$Site = prep_Y1_Rem$Site
prep_Y1_Rem_OPUS_Stand$M12_PANSS_Period_Rem = prep_Y1_Rem$M12_PANSS_Period_Rem

#Get shared data for CRISPFEP (Scottish datasets) external validation
prep_Y1_Rem_CRISPFEP = prep_Y1_Rem[,c(92,54,26,66,14,
                                      91,100,10,31,
                                      68,112,122,56,
                                      21,33,52,17)]
#Rename the columns to something easier and the same as CRISPFEP columns
colnames(prep_Y1_Rem_CRISPFEP) = c("P3","DUP","AdmissionVol","PASClientLASocWith","QualLevel",
                                   "P2","N4","HousingOwnParents","CPN",
                                   "PASClientLASocSex","G9","InsightSymptoms","HelpSeek",
                                   "IncomeSalary","HelpFriendRel","FamilySuggestCare", "InEducation")
#standardise the columns before building model
preProcValues = preProcess(prep_Y1_Rem_CRISPFEP, method = c("center", "scale"))
prep_Y1_Rem_CRISPFEP_Stand = predict(preProcValues, prep_Y1_Rem_CRISPFEP)
#Add factor site and outcome back in
prep_Y1_Rem_CRISPFEP_Stand$Site = prep_Y1_Rem$Site
prep_Y1_Rem_CRISPFEP_Stand$M12_PANSS_Period_Rem = prep_Y1_Rem$M12_PANSS_Period_Rem

#GAF Y1 data
prep_Y1_GAF = preproc_all[ ,!(colnames(preproc_all) %in% c("Site","M12_PANSS_Period_Rem","M12_GAF_Binary","M12_EET","M12_EQ5D_UK_TTO_Index_Binary","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]
#Dummy code (not outcome or site)
dummies = dummyVars(~ ., data = prep_Y1_GAF, fullRank = T)
prep_Y1_GAF <- data.frame(predict(dummies, newdata = prep_Y1_GAF))
#Add factor site and outcome back in
prep_Y1_GAF$Site = preproc_all$Site
prep_Y1_GAF$M12_GAF_Binary = preproc_all$M12_GAF_Binary
#remove na from outcome
prep_Y1_GAF = prep_Y1_GAF[which(!is.na(prep_Y1_GAF$M12_GAF_Binary)),]
#remove columns with more than 20% missing data
prep_Y1_GAF = prep_Y1_GAF[, colMeans(is.na(prep_Y1_GAF)) <= .2]
#remove zero and near zero variance columns
nzv_cols <- nearZeroVar(prep_Y1_GAF)
if(length(nzv_cols) > 0) prep_Y1_GAF <- prep_Y1_GAF[, -nzv_cols]

#Get shared data for OPUS
prep_Y1_GAF_OPUS = prep_Y1_GAF[,c(126,21,69,128,14,
                                  22,92,40,4,67,93,
                                  97,101,33,115,127,90,
                                  35)]
#Rename the columns to something easier
colnames(prep_Y1_GAF_OPUS) = c("GAF","IncomeSalary","PASClientLASocSex","GAFDisability","QualLevel",
                               "IncomeBenefits","P2","Amphetamine","AthiestAgnostic","PASClientLASocWith","P3",
                               "P7","N4","HelpHouseFriendRel","G11","GAFSymptoms","CriminalJustice",
                               "HelpAnyFriendRel")
#standardise the columns before building model
preProcValues = preProcess(prep_Y1_GAF_OPUS, method = c("center", "scale"))
prep_Y1_GAF_OPUS_Stand = predict(preProcValues, prep_Y1_GAF_OPUS)
#Add factor site and outcome back in
prep_Y1_GAF_OPUS_Stand$Site = prep_Y1_GAF$Site
prep_Y1_GAF_OPUS_Stand$M12_GAF_Binary = prep_Y1_GAF$M12_GAF_Binary

#EET Y1 data
prep_Y1_EET = preproc_all[ ,!(colnames(preproc_all) %in% c("Site","M12_PANSS_Period_Rem","M12_GAF_Binary","M12_EET","M12_EQ5D_UK_TTO_Index_Binary","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]
#Dummy code (not outcome or site)
dummies = dummyVars(~ ., data = prep_Y1_EET, fullRank = T)
prep_Y1_EET <- data.frame(predict(dummies, newdata = prep_Y1_EET))
#Add factor site and outcome back in
prep_Y1_EET$Site = preproc_all$Site
prep_Y1_EET$M12_EET = preproc_all$M12_EET
#remove na from outcome
prep_Y1_EET = prep_Y1_EET[which(!is.na(prep_Y1_EET$M12_EET)),]
#remove columns with more than 20% missing data
prep_Y1_EET = prep_Y1_EET[, colMeans(is.na(prep_Y1_EET)) <= .2]
#remove zero and near zero variance columns
nzv_cols <- nearZeroVar(prep_Y1_EET)
if(length(nzv_cols) > 0) prep_Y1_EET <- prep_Y1_EET[, -nzv_cols]

#Get shared data for OPUS
prep_Y1_EET_OPUS = prep_Y1_EET[,c(18,14,129,21,17,
                                  22,140,3,126,127,93,
                                  118,151,56,59,
                                  128,10,104,
                                  33,36,91,116)]
#Rename the columns to something easier
colnames(prep_Y1_EET_OPUS) = c("EET","QualLevel","GAFDisability","IncomeSalary","InEducation",
                               "IncomeBenefits","AnySport","White","Depression","GAF","P2",
                               "G13","SelfHarmKnifeRazor","DUP","PASClientChSocWith",
                               "GAFSymptoms","HousingOwnParents","N6",
                               "HelpOutsideHomeFriendRel","HelpAnyFriendRel","CriminalJustice","G11")
#standardise the columns before building model
preProcValues = preProcess(prep_Y1_EET_OPUS, method = c("center", "scale"))
prep_Y1_EET_OPUS_Stand = predict(preProcValues, prep_Y1_EET_OPUS)
#Add factor site and outcome back in
prep_Y1_EET_OPUS_Stand$Site = prep_Y1_EET$Site
prep_Y1_EET_OPUS_Stand$M12_EET = prep_Y1_EET$M12_EET

#Get shared data for CRISPFEP
prep_Y1_EET_CRISPFEP = prep_Y1_EET[,c(18,14,129,21,17,
                                      22,146,3,52,
                                      126,93,118,56,59,
                                      2,10,54,104,
                                      33,11,36,116)]
#Rename the columns to something easier and the same as CRISPFEP columns
colnames(prep_Y1_EET_CRISPFEP) = c("BL_EET","QualLevel","Disability","IncomeSalary","InEducation",
                                   "IncomeBenefits","Childcare","WhiteBrit","FirstContactPolice",
                                   "Depression","P2","G13","DUP","PASClientChSocWith",
                                   "Pakistani","HousingOwnParents","FamilySuggestCare","N6",
                                   "HelpOutside","Renting","HelpFriendRel","G11")
#standardise the columns before building model
preProcValues = preProcess(prep_Y1_EET_CRISPFEP, method = c("center", "scale"))
prep_Y1_EET_CRISPFEP_Stand = predict(preProcValues, prep_Y1_EET_CRISPFEP)
#Add factor site and outcome back in
prep_Y1_EET_CRISPFEP_Stand$Site = prep_Y1_EET$Site
prep_Y1_EET_CRISPFEP_Stand$M12_EET = prep_Y1_EET$M12_EET

#EQ5D TTO Y1 data
prep_Y1_EQ5D_TTO = preproc_all[ ,!(colnames(preproc_all) %in% c("Site","M12_PANSS_Period_Rem","M12_GAF_Binary","M12_EET","M12_EQ5D_UK_TTO_Index_Binary","FUP_PANSS_Peroid_Rem","FUP_GAF_Binary","FUP_EET","FUP_EQ5D_UK_TTO_Index_Binary"))]
#Dummy code (not outcome or site)
dummies = dummyVars(~ ., data = prep_Y1_EQ5D_TTO, fullRank = T)
prep_Y1_EQ5D_TTO <- data.frame(predict(dummies, newdata = prep_Y1_EQ5D_TTO))
#Add factor site and outcome back in
prep_Y1_EQ5D_TTO$Site = preproc_all$Site
prep_Y1_EQ5D_TTO$M12_EQ5D_UK_TTO_Index_Binary = preproc_all$M12_EQ5D_UK_TTO_Index_Binary
#remove na from outcome
prep_Y1_EQ5D_TTO = prep_Y1_EQ5D_TTO[which(!is.na(prep_Y1_EQ5D_TTO$M12_EQ5D_UK_TTO_Index_Binary)),]
#remove columns with more than 20% missing data
prep_Y1_EQ5D_TTO = prep_Y1_EQ5D_TTO[, colMeans(is.na(prep_Y1_EQ5D_TTO)) <= .2]
#remove zero and near zero variance columns
nzv_cols <- nearZeroVar(prep_Y1_EQ5D_TTO)
if(length(nzv_cols) > 0) prep_Y1_EQ5D_TTO <- prep_Y1_EQ5D_TTO[, -nzv_cols]

#Get shared data for OPUS
prep_Y1_QoL_OPUS = prep_Y1_EQ5D_TTO[,c(40,93,62,10,
                                       66,89,14,21,115,
                                       92,95,110,126,41,67,
                                       22,8,150,97,
                                       151,103,33,55,1,
                                       45,39)]
#Rename the columns to something easier
colnames(prep_Y1_QoL_OPUS) = c("Amphetamine","P3","PASClientEASocWith","HousingOwnParents",
                               "PASClientEASocSex","EduLevel","QualLevel","IncomeSalary","G11",
                               "P2","P5","G6","GAF","Cocaine","PASClientLASocWith",
                               "IncomeBenefits","LivingParents","SelfHarmKnifeRazor","P7",
                               "NeverSelfHarmed","N6","HelpHouseFriendRel","DUP","Male",
                               "LSD","Cannabis")
#standardise the columns before building model
preProcValues = preProcess(prep_Y1_QoL_OPUS, method = c("center", "scale"))
prep_Y1_QoL_OPUS_Stand = predict(preProcValues, prep_Y1_QoL_OPUS)
#Add factor site and outcome back in
prep_Y1_QoL_OPUS_Stand$Site = prep_Y1_EQ5D_TTO$Site
prep_Y1_QoL_OPUS_Stand$M12_QoL = prep_Y1_EQ5D_TTO$M12_EQ5D_UK_TTO_Index_Binary

#Get shared data for CRISPFEP
prep_Y1_QoL_CRISPFEP = prep_Y1_EQ5D_TTO[,c(106,93,62,10,66,
                                           89,51,21,115,92,95,
                                           110,41,67,22,96,
                                           8,97,103,11,55,102,1,39)]
#Rename the columns to something easier and the same as CRISPFEP columns
colnames(prep_Y1_QoL_CRISPFEP) = c("G2","P3","PASClientEASocWith","HousingOwnParents","PASClientEASocSex",
                                   "EduLevel","FirstContactPolice","IncomeSalary","G11","P2","P5",
                                   "G6","Cocaine","PASClientLASocWith","IncomeBenefits","P6",
                                   "LivingParents","P7","N6","Renting","DUP","N5","Male","Cannabis")
#standardise the columns before building model
preProcValues = preProcess(prep_Y1_QoL_CRISPFEP, method = c("center", "scale"))
prep_Y1_QoL_CRISPFEP_Stand = predict(preProcValues, prep_Y1_QoL_CRISPFEP)
#Add factor site and outcome back in
prep_Y1_QoL_CRISPFEP_Stand$Site = prep_Y1_EQ5D_TTO$Site
prep_Y1_QoL_CRISPFEP_Stand$M12_QoL = prep_Y1_EQ5D_TTO$M12_EQ5D_UK_TTO_Index_Binary

#Takes 10 mins
#do a 10 fold cv for training as we are nesting this in a 14 site cv (as per BJPsych https://doi.org/10.1192/bjp.2018.122)
control <- trainControl(method="cv", number=10, classProbs=TRUE, summaryFunction=twoClassSummary, selectionFunction ="oneSE")

##############################################
#Y1 Remission Outcome and Model
#
results_Y1_Rem = list()
mods_Y1_Rem = list()

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  
  #set up leave one site out cv
  test_set = prep_Y1_Rem[ which(prep_Y1_Rem$Site == levels(prep_Y1_Rem$Site)[i]), ]
  train_set = prep_Y1_Rem[ -which(prep_Y1_Rem$Site == levels(prep_Y1_Rem$Site)[i]), ]
  
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  train_set = train_set[ ,!(colnames(train_set) %in% c("Site"))]
 
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_PANSS_Period_Rem) 
  
  #train model over grid of 10 by 10 lambda alpha, knn impute, standardise
  mod <- train(M12_PANSS_Period_Rem ~ ., data=train_set, method="glmnet", metric="ROC", tuneLength = 10, preProc = c("center", "scale","knnImpute"), trControl=control,na.action = na.pass)
  result$pred = predict(mod, test_set, type = "prob", na.action = na.pass)
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_Rem[[i]] = result
  mods_Y1_Rem[[i]] = mod
}

results_Y1_Rem_seq = NULL
results_Y1_Rem_seq$pred = c(results_Y1_Rem[[1]]$pred$Y,results_Y1_Rem[[2]]$pred$Y, results_Y1_Rem[[3]]$pred$Y,
                            results_Y1_Rem[[4]]$pred$Y,results_Y1_Rem[[5]]$pred$Y, results_Y1_Rem[[6]]$pred$Y,
                            results_Y1_Rem[[7]]$pred$Y,results_Y1_Rem[[8]]$pred$Y,results_Y1_Rem[[9]]$pred$Y,
                            results_Y1_Rem[[10]]$pred$Y,results_Y1_Rem[[11]]$pred$Y,results_Y1_Rem[[12]]$pred$Y,
                            results_Y1_Rem[[13]]$pred$Y,results_Y1_Rem[[14]]$pred$Y)
results_Y1_Rem_seq$obs = factor(c(as.character(results_Y1_Rem[[1]]$obs),as.character(results_Y1_Rem[[2]]$obs),
                                     as.character(results_Y1_Rem[[3]]$obs),as.character(results_Y1_Rem[[4]]$obs),
                                     as.character(results_Y1_Rem[[5]]$obs),as.character(results_Y1_Rem[[6]]$obs),
                                     as.character(results_Y1_Rem[[7]]$obs),as.character(results_Y1_Rem[[8]]$obs),
                                     as.character(results_Y1_Rem[[9]]$obs),as.character(results_Y1_Rem[[10]]$obs),
                                     as.character(results_Y1_Rem[[11]]$obs),as.character(results_Y1_Rem[[12]]$obs),
                                     as.character(results_Y1_Rem[[13]]$obs),as.character(results_Y1_Rem[[14]]$obs)))
results_Y1_Rem_seq$row = c(results_Y1_Rem[[1]]$row,results_Y1_Rem[[2]]$row, results_Y1_Rem[[3]]$row,
                            results_Y1_Rem[[4]]$row,results_Y1_Rem[[5]]$row, results_Y1_Rem[[6]]$row,
                            results_Y1_Rem[[7]]$row,results_Y1_Rem[[8]]$row,results_Y1_Rem[[9]]$row,
                            results_Y1_Rem[[10]]$row,results_Y1_Rem[[11]]$row,results_Y1_Rem[[12]]$row,
                            results_Y1_Rem[[13]]$row,results_Y1_Rem[[14]]$row)

auc_Y1_Rem_seq = roc(predictor = results_Y1_Rem_seq$pred, response = results_Y1_Rem_seq$obs)$auc
roc_Y1_Rem_seq = roc(predictor = results_Y1_Rem_seq$pred, response = results_Y1_Rem_seq$obs)

plot(roc_Y1_Rem_seq)
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
set.seed(987)
ci.coords(roc_Y1_Rem_seq,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

#make ROC diagram
png("figure1A.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(roc_Y1_Rem_seq, auc.polygon=TRUE, grid=c(0.1, 0.1))
text(0.3,0.3,"AUC = 0·703\n(0·664, 0·742)\np<0·0001")
dev.off()

set.seed(987)
Y1_Rem_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_Rem_perm = permute(results_Y1_Rem_seq$obs)
  Y1_Rem_auc_null = c(Y1_Rem_auc_null, roc(predictor = results_Y1_Rem_seq$pred, response = Y1_Rem_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_Rem_auc_null >= auc_Y1_Rem_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_Rem_seq$pred, results_Y1_Rem_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_Rem_seq$pred, results_Y1_Rem_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_Rem_seq$pred, results_Y1_Rem_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_Rem_seq$pred, results_Y1_Rem_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_Rem_seq$pred, results_Y1_Rem_seq$obs=="Yes")$var))

#Test Remission model for each of the other three outcomes on EDEN data
#Y1 GAF
#
results_Y1_Rem_GAF = list()
#Add columns back in that were removed in dataset preprocessing but are required for remission model, but just as NA to allow imputation
prep_Y1_GAF_Rem = prep_Y1_GAF
setdiff(colnames(mods_Y1_Rem[[1]]$trainingData),colnames(prep_Y1_GAF_Rem))
prep_Y1_GAF_Rem$BL_MostSeriousHarm_Premeditation.Self.harm.contemplated.for.more.than.three.hours.prior.to.attempt = NA
prep_Y1_GAF_Rem$BL_MostSeriousHarm_Premeditation.Self.harm.contemplated.for.more.than.three.hours.prior.to.attempt = as.numeric(prep_Y1_GAF_Rem$BL_MostSeriousHarm_Premeditation.Self.harm.contemplated.for.more.than.three.hours.prior.to.attempt)

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  #set up leave one site out cv
  test_set = prep_Y1_GAF_Rem[ which(prep_Y1_GAF_Rem$Site == levels(prep_Y1_GAF_Rem$Site)[i]), ]

  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]

  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_GAF) 
  
  result$pred = predict(mods_Y1_Rem[[i]], test_set, type = "prob", na.action = na.pass)
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_Rem_GAF[[i]] = result
}

results_Y1_Rem_GAF_seq = NULL
results_Y1_Rem_GAF_seq$pred = c(results_Y1_Rem_GAF[[1]]$pred$Y,results_Y1_Rem_GAF[[2]]$pred$Y, results_Y1_Rem_GAF[[3]]$pred$Y,
                            results_Y1_Rem_GAF[[4]]$pred$Y,results_Y1_Rem_GAF[[5]]$pred$Y, results_Y1_Rem_GAF[[6]]$pred$Y,
                            results_Y1_Rem_GAF[[7]]$pred$Y,results_Y1_Rem_GAF[[8]]$pred$Y,results_Y1_Rem_GAF[[9]]$pred$Y,
                            results_Y1_Rem_GAF[[10]]$pred$Y,results_Y1_Rem_GAF[[11]]$pred$Y,results_Y1_Rem_GAF[[12]]$pred$Y,
                            results_Y1_Rem_GAF[[13]]$pred$Y,results_Y1_Rem_GAF[[14]]$pred$Y)
results_Y1_Rem_GAF_seq$obs = factor(c(as.character(results_Y1_Rem_GAF[[1]]$obs),as.character(results_Y1_Rem_GAF[[2]]$obs),
                                  as.character(results_Y1_Rem_GAF[[3]]$obs),as.character(results_Y1_Rem_GAF[[4]]$obs),
                                  as.character(results_Y1_Rem_GAF[[5]]$obs),as.character(results_Y1_Rem_GAF[[6]]$obs),
                                  as.character(results_Y1_Rem_GAF[[7]]$obs),as.character(results_Y1_Rem_GAF[[8]]$obs),
                                  as.character(results_Y1_Rem_GAF[[9]]$obs),as.character(results_Y1_Rem_GAF[[10]]$obs),
                                  as.character(results_Y1_Rem_GAF[[11]]$obs),as.character(results_Y1_Rem_GAF[[12]]$obs),
                                  as.character(results_Y1_Rem_GAF[[13]]$obs),as.character(results_Y1_Rem_GAF[[14]]$obs)))
results_Y1_Rem_GAF_seq$row = c(results_Y1_Rem_GAF[[1]]$row,results_Y1_Rem_GAF[[2]]$row, results_Y1_Rem_GAF[[3]]$row,
                           results_Y1_Rem_GAF[[4]]$row,results_Y1_Rem_GAF[[5]]$row, results_Y1_Rem_GAF[[6]]$row,
                           results_Y1_Rem_GAF[[7]]$row,results_Y1_Rem_GAF[[8]]$row,results_Y1_Rem_GAF[[9]]$row,
                           results_Y1_Rem_GAF[[10]]$row,results_Y1_Rem_GAF[[11]]$row,results_Y1_Rem_GAF[[12]]$row,
                           results_Y1_Rem_GAF[[13]]$row,results_Y1_Rem_GAF[[14]]$row)

auc_Y1_Rem_GAF_seq = roc(predictor = results_Y1_Rem_GAF_seq$pred, response = results_Y1_Rem_GAF_seq$obs)$auc
roc_Y1_Rem_GAF_seq = roc(predictor = results_Y1_Rem_GAF_seq$pred, response = results_Y1_Rem_GAF_seq$obs)

set.seed(987)
Y1_Rem_GAF_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_Rem_GAF_perm = permute(results_Y1_Rem_GAF_seq$obs)
  Y1_Rem_GAF_auc_null = c(Y1_Rem_GAF_auc_null, roc(predictor = results_Y1_Rem_GAF_seq$pred, response = Y1_Rem_GAF_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_Rem_GAF_auc_null >= auc_Y1_Rem_GAF_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_Rem_GAF_seq$pred, results_Y1_Rem_GAF_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_Rem_GAF_seq$pred, results_Y1_Rem_GAF_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_Rem_GAF_seq$pred, results_Y1_Rem_GAF_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_Rem_GAF_seq$pred, results_Y1_Rem_GAF_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_Rem_GAF_seq$pred, results_Y1_Rem_GAF_seq$obs=="Yes")$var))

#Y1 EET
#
results_Y1_Rem_EET = list()
#Add columns back in that were removed and required for this model
prep_Y1_EET_Rem = prep_Y1_EET
setdiff(colnames(mods_Y1_Rem[[1]]$trainingData),colnames(prep_Y1_EET_Rem))
prep_Y1_EET_Rem$BL_MostSeriousHarm_Premeditation.Self.harm.contemplated.for.more.than.three.hours.prior.to.attempt = NA
prep_Y1_EET_Rem$BL_MostSeriousHarm_Premeditation.Self.harm.contemplated.for.more.than.three.hours.prior.to.attempt = as.numeric(prep_Y1_EET_Rem$BL_MostSeriousHarm_Premeditation.Self.harm.contemplated.for.more.than.three.hours.prior.to.attempt)

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  #set up leave one site out cv
  test_set = prep_Y1_EET_Rem[ which(prep_Y1_EET_Rem$Site == levels(prep_Y1_EET_Rem$Site)[i]), ]
  
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_EET) 
  result$pred = predict(mods_Y1_Rem[[i]], test_set, type = "prob", na.action = na.pass)
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_Rem_EET[[i]] = result
}

results_Y1_Rem_EET_seq = NULL
results_Y1_Rem_EET_seq$pred = c(results_Y1_Rem_EET[[1]]$pred$Y,results_Y1_Rem_EET[[2]]$pred$Y, results_Y1_Rem_EET[[3]]$pred$Y,
                                results_Y1_Rem_EET[[4]]$pred$Y,results_Y1_Rem_EET[[5]]$pred$Y, results_Y1_Rem_EET[[6]]$pred$Y,
                                results_Y1_Rem_EET[[7]]$pred$Y,results_Y1_Rem_EET[[8]]$pred$Y,results_Y1_Rem_EET[[9]]$pred$Y,
                                results_Y1_Rem_EET[[10]]$pred$Y,results_Y1_Rem_EET[[11]]$pred$Y,results_Y1_Rem_EET[[12]]$pred$Y,
                                results_Y1_Rem_EET[[13]]$pred$Y,results_Y1_Rem_EET[[14]]$pred$Y)
results_Y1_Rem_EET_seq$obs = factor(c(as.character(results_Y1_Rem_EET[[1]]$obs),as.character(results_Y1_Rem_EET[[2]]$obs),
                                      as.character(results_Y1_Rem_EET[[3]]$obs),as.character(results_Y1_Rem_EET[[4]]$obs),
                                      as.character(results_Y1_Rem_EET[[5]]$obs),as.character(results_Y1_Rem_EET[[6]]$obs),
                                      as.character(results_Y1_Rem_EET[[7]]$obs),as.character(results_Y1_Rem_EET[[8]]$obs),
                                      as.character(results_Y1_Rem_EET[[9]]$obs),as.character(results_Y1_Rem_EET[[10]]$obs),
                                      as.character(results_Y1_Rem_EET[[11]]$obs),as.character(results_Y1_Rem_EET[[12]]$obs),
                                      as.character(results_Y1_Rem_EET[[13]]$obs),as.character(results_Y1_Rem_EET[[14]]$obs)))
results_Y1_Rem_EET_seq$row = c(results_Y1_Rem_EET[[1]]$row,results_Y1_Rem_EET[[2]]$row, results_Y1_Rem_EET[[3]]$row,
                               results_Y1_Rem_EET[[4]]$row,results_Y1_Rem_EET[[5]]$row, results_Y1_Rem_EET[[6]]$row,
                               results_Y1_Rem_EET[[7]]$row,results_Y1_Rem_EET[[8]]$row,results_Y1_Rem_EET[[9]]$row,
                               results_Y1_Rem_EET[[10]]$row,results_Y1_Rem_EET[[11]]$row,results_Y1_Rem_EET[[12]]$row,
                               results_Y1_Rem_EET[[13]]$row,results_Y1_Rem_EET[[14]]$row)

auc_Y1_Rem_EET_seq = roc(predictor = results_Y1_Rem_EET_seq$pred, response = results_Y1_Rem_EET_seq$obs)$auc
roc_Y1_Rem_EET_seq = roc(predictor = results_Y1_Rem_EET_seq$pred, response = results_Y1_Rem_EET_seq$obs)

set.seed(987)
Y1_Rem_EET_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_Rem_EET_perm = permute(results_Y1_Rem_EET_seq$obs)
  Y1_Rem_EET_auc_null = c(Y1_Rem_EET_auc_null, roc(predictor = results_Y1_Rem_EET_seq$pred, response = Y1_Rem_EET_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_Rem_EET_auc_null >= auc_Y1_Rem_EET_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_Rem_EET_seq$pred, results_Y1_Rem_EET_seq$obs=="Yes")
#95CI = 1.96*SE#95CI = 1.96*SE
(roc.area.test(results_Y1_Rem_EET_seq$pred, results_Y1_Rem_EET_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_Rem_EET_seq$pred, results_Y1_Rem_EET_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_Rem_EET_seq$pred, results_Y1_Rem_EET_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_Rem_EET_seq$pred, results_Y1_Rem_EET_seq$obs=="Yes")$var))

#Y1 EQ5D_TTO
#
results_Y1_Rem_EQ5D_TTO = list()
#Add columns back in that were removed and required for this model (none)
prep_Y1_EQ5D_TTO_Rem = prep_Y1_EQ5D_TTO
setdiff(colnames(mods_Y1_Rem[[1]]$trainingData),colnames(prep_Y1_EQ5D_TTO_Rem))

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  #set up leave one site out cv
  test_set = prep_Y1_EQ5D_TTO_Rem[ which(prep_Y1_EQ5D_TTO_Rem$Site == levels(prep_Y1_EQ5D_TTO_Rem$Site)[i]), ]
  
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_EQ5D_UK_TTO_Index_Binary) 
  result$pred = predict(mods_Y1_Rem[[i]], test_set, type = "prob", na.action = na.pass)
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_Rem_EQ5D_TTO[[i]] = result
}

results_Y1_Rem_EQ5D_TTO_seq = NULL
results_Y1_Rem_EQ5D_TTO_seq$pred = c(results_Y1_Rem_EQ5D_TTO[[1]]$pred$Y,results_Y1_Rem_EQ5D_TTO[[2]]$pred$Y, results_Y1_Rem_EQ5D_TTO[[3]]$pred$Y,
                                results_Y1_Rem_EQ5D_TTO[[4]]$pred$Y,results_Y1_Rem_EQ5D_TTO[[5]]$pred$Y, results_Y1_Rem_EQ5D_TTO[[6]]$pred$Y,
                                results_Y1_Rem_EQ5D_TTO[[7]]$pred$Y,results_Y1_Rem_EQ5D_TTO[[8]]$pred$Y,results_Y1_Rem_EQ5D_TTO[[9]]$pred$Y,
                                results_Y1_Rem_EQ5D_TTO[[10]]$pred$Y,results_Y1_Rem_EQ5D_TTO[[11]]$pred$Y,results_Y1_Rem_EQ5D_TTO[[12]]$pred$Y,
                                results_Y1_Rem_EQ5D_TTO[[13]]$pred$Y,results_Y1_Rem_EQ5D_TTO[[14]]$pred$Y)
results_Y1_Rem_EQ5D_TTO_seq$obs = factor(c(as.character(results_Y1_Rem_EQ5D_TTO[[1]]$obs),as.character(results_Y1_Rem_EQ5D_TTO[[2]]$obs),
                                      as.character(results_Y1_Rem_EQ5D_TTO[[3]]$obs),as.character(results_Y1_Rem_EQ5D_TTO[[4]]$obs),
                                      as.character(results_Y1_Rem_EQ5D_TTO[[5]]$obs),as.character(results_Y1_Rem_EQ5D_TTO[[6]]$obs),
                                      as.character(results_Y1_Rem_EQ5D_TTO[[7]]$obs),as.character(results_Y1_Rem_EQ5D_TTO[[8]]$obs),
                                      as.character(results_Y1_Rem_EQ5D_TTO[[9]]$obs),as.character(results_Y1_Rem_EQ5D_TTO[[10]]$obs),
                                      as.character(results_Y1_Rem_EQ5D_TTO[[11]]$obs),as.character(results_Y1_Rem_EQ5D_TTO[[12]]$obs),
                                      as.character(results_Y1_Rem_EQ5D_TTO[[13]]$obs),as.character(results_Y1_Rem_EQ5D_TTO[[14]]$obs)))
results_Y1_Rem_EQ5D_TTO_seq$row = c(results_Y1_Rem_EQ5D_TTO[[1]]$row,results_Y1_Rem_EQ5D_TTO[[2]]$row, results_Y1_Rem_EQ5D_TTO[[3]]$row,
                               results_Y1_Rem_EQ5D_TTO[[4]]$row,results_Y1_Rem_EQ5D_TTO[[5]]$row, results_Y1_Rem_EQ5D_TTO[[6]]$row,
                               results_Y1_Rem_EQ5D_TTO[[7]]$row,results_Y1_Rem_EQ5D_TTO[[8]]$row,results_Y1_Rem_EQ5D_TTO[[9]]$row,
                               results_Y1_Rem_EQ5D_TTO[[10]]$row,results_Y1_Rem_EQ5D_TTO[[11]]$row,results_Y1_Rem_EQ5D_TTO[[12]]$row,
                               results_Y1_Rem_EQ5D_TTO[[13]]$row,results_Y1_Rem_EQ5D_TTO[[14]]$row)

auc_Y1_Rem_EQ5D_TTO_seq = roc(predictor = results_Y1_Rem_EQ5D_TTO_seq$pred, response = results_Y1_Rem_EQ5D_TTO_seq$obs)$auc
roc_Y1_Rem_EQ5D_TTO_seq = roc(predictor = results_Y1_Rem_EQ5D_TTO_seq$pred, response = results_Y1_Rem_EQ5D_TTO_seq$obs)

set.seed(987)
Y1_Rem_EQ5D_TTO_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_Rem_EQ5D_TTO_perm = permute(results_Y1_Rem_EQ5D_TTO_seq$obs)
  Y1_Rem_EQ5D_TTO_auc_null = c(Y1_Rem_EQ5D_TTO_auc_null, roc(predictor = results_Y1_Rem_EQ5D_TTO_seq$pred, response = Y1_Rem_EQ5D_TTO_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_Rem_EQ5D_TTO_auc_null >= auc_Y1_Rem_EQ5D_TTO_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_Rem_EQ5D_TTO_seq$pred, results_Y1_Rem_EQ5D_TTO_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_Rem_EQ5D_TTO_seq$pred, results_Y1_Rem_EQ5D_TTO_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_Rem_EQ5D_TTO_seq$pred, results_Y1_Rem_EQ5D_TTO_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_Rem_EQ5D_TTO_seq$pred, results_Y1_Rem_EQ5D_TTO_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_Rem_EQ5D_TTO_seq$pred, results_Y1_Rem_EQ5D_TTO_seq$obs=="Yes")$var))

coefs_Y1_Rem = NULL
for (i in seq(1:14))
{
  coefs_Y1_Rem = c(coefs_Y1_Rem, coef(mods_Y1_Rem[[i]]$finalModel, mods_Y1_Rem[[i]]$bestTune$lambda))
}

#just get numbers
coefs_Y1_Rem_extract = NULL
for(i in seq(1:14))
{
  coefs_Y1_Rem_extract = rbind(coefs_Y1_Rem_extract, coefs_Y1_Rem[[i]][1:164])
}

#get matrix of coefficients presence (1) or absence (0)
#Presence or absence of predictors across all 14 LOSOCV models
coefs_Y1_Rem_presence = NULL
coefs_Y1_Rem_presence = coefs_Y1_Rem_extract[1:14,2:164]
coefs_Y1_Rem_presence[coefs_Y1_Rem_presence != 0] <- 1

#stability of feature selection http://jmlr.org/papers/volume18/17-514/17-514.pdf
getStability(coefs_Y1_Rem_presence)

#work out number of predictors shared across all models
table(colMeans(coefs_Y1_Rem_presence))

#get rank of coef by importance as in sports ranking
coefs_Y1_Rem_baseline_rank = NULL

for(i in seq(c(1:14)))
{
  #rank absolute value excluding the intercept for each model
  coefs_Y1_Rem_baseline_rank = rbind(coefs_Y1_Rem_baseline_rank, rank(abs(coefs_Y1_Rem_extract[i,2:164]), ties.method = "min"))
}

# rank the mean ranks of each column across all models
coefs_Y1_Rem_baseline_rank_mean = colMeans(coefs_Y1_Rem_baseline_rank)

#Invert order of rank to identify top models
coefs_Y1_Rem_baseline_order = rank(-coefs_Y1_Rem_baseline_rank_mean)

#Get the column names (not the intercept)
coef_Y1_Rem_names = dimnames(coefs_Y1_Rem[[1]])[[1]][2:164]

coefs_Y1_Rem_means = colMeans(coefs_Y1_Rem_extract)[2:164]

df_Y1_Rem_names = data.frame(coef_Y1_Rem_names, coefs_Y1_Rem_baseline_order, coefs_Y1_Rem_means, colMeans(coefs_Y1_Rem_presence))

#View the best predictor variables by absolute value in order (lower is better)
View(df_Y1_Rem_names)

#Build model on whole EDEN as training set (one SE rule)
prep_Y1_Rem_all = prep_Y1_Rem[ ,!(colnames(prep_Y1_Rem) %in% c("Site"))]
controlAll <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary, selectionFunction ="oneSE")
set.seed(987)
mod_Y1_Rem_all = train(M12_PANSS_Period_Rem ~ ., data=prep_Y1_Rem_all, method="glmnet", metric="ROC", tuneLength = 10, preProc = c("center", "scale","knnImpute"), trControl=controlAll, na.action = na.pass)
coef_Y1_Rem_all = coef(mod_Y1_Rem_all$finalModel, mod_Y1_Rem_all$bestTune$lambda)
coefs_Y1_Rem_extract_all = coef_Y1_Rem_all[1:164]
coefs_Y1_Rem_baseline_rank_all = rank(-abs(coefs_Y1_Rem_extract_all[2:164]), ties.method = "min")
df_Y1_Rem_names_all = data.frame(coef_Y1_Rem_names, coefs_Y1_Rem_baseline_rank_all,coefs_Y1_Rem_extract_all[2:164])
View(df_Y1_Rem_names_all)

#Internal-External Validation for OPUS only columns to check still works
#tried glmnet, glm, rf, linear and radial svm - all very similar. Just use GLM
#Y1 Remission
#
results_Y1_Rem_OPUS = list()
mods_Y1_Rem_OPUS = list()

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  
  #set up leave one site out cv
  test_set = prep_Y1_Rem_OPUS_Stand[ which(prep_Y1_Rem_OPUS_Stand$Site == levels(prep_Y1_Rem_OPUS_Stand$Site)[i]), ]
  train_set = prep_Y1_Rem_OPUS_Stand[ -which(prep_Y1_Rem_OPUS_Stand$Site == levels(prep_Y1_Rem_OPUS_Stand$Site)[i]), ]
  
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  train_set = train_set[ ,!(colnames(train_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_PANSS_Period_Rem) 
  
  #train model over grid of 10 by 10 lambda alpha, knn impute, standardise
  mod <- train(M12_PANSS_Period_Rem ~ ., data=train_set, method="glm", metric="ROC", tuneLength = 10, preProc = c("knnImpute"), trControl=control,na.action = na.pass)
  result$pred = predict(mod, test_set, type = "prob", na.action = na.pass)
  
  results_Y1_Rem_OPUS[[i]] = result
  mods_Y1_Rem_OPUS[[i]] = mod
}

results_Y1_Rem_OPUS_seq = NULL
results_Y1_Rem_OPUS_seq$pred = c(results_Y1_Rem_OPUS[[1]]$pred$Y,results_Y1_Rem_OPUS[[2]]$pred$Y, results_Y1_Rem_OPUS[[3]]$pred$Y,
                                 results_Y1_Rem_OPUS[[4]]$pred$Y,results_Y1_Rem_OPUS[[5]]$pred$Y, results_Y1_Rem_OPUS[[6]]$pred$Y,
                                 results_Y1_Rem_OPUS[[7]]$pred$Y,results_Y1_Rem_OPUS[[8]]$pred$Y,results_Y1_Rem_OPUS[[9]]$pred$Y,
                                 results_Y1_Rem_OPUS[[10]]$pred$Y,results_Y1_Rem_OPUS[[11]]$pred$Y,results_Y1_Rem_OPUS[[12]]$pred$Y,
                                 results_Y1_Rem_OPUS[[13]]$pred$Y,results_Y1_Rem_OPUS[[14]]$pred$Y)
results_Y1_Rem_OPUS_seq$obs = factor(c(as.character(results_Y1_Rem_OPUS[[1]]$obs),as.character(results_Y1_Rem_OPUS[[2]]$obs),
                                  as.character(results_Y1_Rem_OPUS[[3]]$obs),as.character(results_Y1_Rem_OPUS[[4]]$obs),
                                  as.character(results_Y1_Rem_OPUS[[5]]$obs),as.character(results_Y1_Rem_OPUS[[6]]$obs),
                                  as.character(results_Y1_Rem_OPUS[[7]]$obs),as.character(results_Y1_Rem_OPUS[[8]]$obs),
                                  as.character(results_Y1_Rem_OPUS[[9]]$obs),as.character(results_Y1_Rem_OPUS[[10]]$obs),
                                  as.character(results_Y1_Rem_OPUS[[11]]$obs),as.character(results_Y1_Rem_OPUS[[12]]$obs),
                                  as.character(results_Y1_Rem_OPUS[[13]]$obs),as.character(results_Y1_Rem_OPUS[[14]]$obs)))

auc_Y1_Rem_OPUS_seq = roc(predictor = results_Y1_Rem_OPUS_seq$pred, response = results_Y1_Rem_OPUS_seq$obs)$auc
roc_Y1_Rem_OPUS_seq = roc(predictor = results_Y1_Rem_OPUS_seq$pred, response = results_Y1_Rem_OPUS_seq$obs)

plot(roc_Y1_Rem_OPUS_seq)
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
coords(roc_Y1_Rem_OPUS_seq,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

set.seed(987)
Y1_Rem_auc_OPUS_null = NULL
for(i in seq (1:10001))
{
  Y1_Rem_OPUS_perm = permute(results_Y1_Rem_OPUS_seq$obs)
  Y1_Rem_auc_OPUS_null = c(Y1_Rem_auc_OPUS_null, roc(predictor = results_Y1_Rem_OPUS_seq$pred, response = Y1_Rem_OPUS_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_Rem_auc_OPUS_null >= auc_Y1_Rem_OPUS_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_Rem_OPUS_seq$pred, results_Y1_Rem_OPUS_seq$obs=="Yes")
#95CI = 1.96*SE
1.96*sqrt(roc.area.test(results_Y1_Rem_OPUS_seq$pred, results_Y1_Rem_OPUS_seq$obs=="Yes")$var)

#Build GLM model with OPUS columns on whole EDEN as training set and use this model for external validation
prep_Y1_Rem_OPUS_all = prep_Y1_Rem_OPUS_Stand[ ,!(colnames(prep_Y1_Rem_OPUS_Stand) %in% c("Site"))]
controlAll <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary)
set.seed(987)
mod_Y1_Rem_EDEN_all = train(M12_PANSS_Period_Rem ~ ., data=prep_Y1_Rem_OPUS_all, method="glm", metric="ROC", preProc = c("knnImpute"), trControl=controlAll, na.action = na.pass)
#remove training data for sharing model
mod_Y1_Rem_EDEN_all$trainingData = NULL
save(mod_Y1_Rem_EDEN_all, file = "mod_Y1_Rem_OPUS.rda")

#Internal-External Validation for CRISPFEP only columns to check it still works
#Y1 Remission
#
results_Y1_Rem_CRISPFEP = list()
mods_Y1_Rem_CRISPFEP = list()

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  
  #set up leave one site out cv
  test_set = prep_Y1_Rem_CRISPFEP_Stand[ which(prep_Y1_Rem_CRISPFEP_Stand$Site == levels(prep_Y1_Rem_CRISPFEP_Stand$Site)[i]), ]
  train_set = prep_Y1_Rem_CRISPFEP_Stand[ -which(prep_Y1_Rem_CRISPFEP_Stand$Site == levels(prep_Y1_Rem_CRISPFEP_Stand$Site)[i]), ]
  
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  train_set = train_set[ ,!(colnames(train_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_PANSS_Period_Rem) 
  
  #train model over grid of 10 by 10 lambda alpha, knn impute, standardise
  mod <- train(M12_PANSS_Period_Rem ~ ., data=train_set, method="glm", metric="ROC", tuneLength = 10, preProc = c("knnImpute"), trControl=control,na.action = na.pass)
  result$pred = predict(mod, test_set, type = "prob", na.action = na.pass)
  
  results_Y1_Rem_CRISPFEP[[i]] = result
  mods_Y1_Rem_CRISPFEP[[i]] = mod
}

results_Y1_Rem_CRISPFEP_seq = NULL
results_Y1_Rem_CRISPFEP_seq$pred = c(results_Y1_Rem_CRISPFEP[[1]]$pred$Y,results_Y1_Rem_CRISPFEP[[2]]$pred$Y, results_Y1_Rem_CRISPFEP[[3]]$pred$Y,
                                     results_Y1_Rem_CRISPFEP[[4]]$pred$Y,results_Y1_Rem_CRISPFEP[[5]]$pred$Y, results_Y1_Rem_CRISPFEP[[6]]$pred$Y,
                                     results_Y1_Rem_CRISPFEP[[7]]$pred$Y,results_Y1_Rem_CRISPFEP[[8]]$pred$Y,results_Y1_Rem_CRISPFEP[[9]]$pred$Y,
                                     results_Y1_Rem_CRISPFEP[[10]]$pred$Y,results_Y1_Rem_CRISPFEP[[11]]$pred$Y,results_Y1_Rem_CRISPFEP[[12]]$pred$Y,
                                     results_Y1_Rem_CRISPFEP[[13]]$pred$Y,results_Y1_Rem_CRISPFEP[[14]]$pred$Y)
results_Y1_Rem_CRISPFEP_seq$obs = factor(c(as.character(results_Y1_Rem_CRISPFEP[[1]]$obs),as.character(results_Y1_Rem_CRISPFEP[[2]]$obs),
                                           as.character(results_Y1_Rem_CRISPFEP[[3]]$obs),as.character(results_Y1_Rem_CRISPFEP[[4]]$obs),
                                           as.character(results_Y1_Rem_CRISPFEP[[5]]$obs),as.character(results_Y1_Rem_CRISPFEP[[6]]$obs),
                                           as.character(results_Y1_Rem_CRISPFEP[[7]]$obs),as.character(results_Y1_Rem_CRISPFEP[[8]]$obs),
                                           as.character(results_Y1_Rem_CRISPFEP[[9]]$obs),as.character(results_Y1_Rem_CRISPFEP[[10]]$obs),
                                           as.character(results_Y1_Rem_CRISPFEP[[11]]$obs),as.character(results_Y1_Rem_CRISPFEP[[12]]$obs),
                                           as.character(results_Y1_Rem_CRISPFEP[[13]]$obs),as.character(results_Y1_Rem_CRISPFEP[[14]]$obs)))

auc_Y1_Rem_CRISPFEP_seq = roc(predictor = results_Y1_Rem_CRISPFEP_seq$pred, response = results_Y1_Rem_CRISPFEP_seq$obs)$auc
roc_Y1_Rem_CRISPFEP_seq = roc(predictor = results_Y1_Rem_CRISPFEP_seq$pred, response = results_Y1_Rem_CRISPFEP_seq$obs)

plot(roc_Y1_Rem_CRISPFEP_seq)
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
coords(roc_Y1_Rem_CRISPFEP_seq,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

set.seed(987)
Y1_Rem_auc_CRISPFEP_null = NULL
for(i in seq (1:10001))
{
  Y1_Rem_CRISPFEP_perm = permute(results_Y1_Rem_CRISPFEP_seq$obs)
  Y1_Rem_auc_CRISPFEP_null = c(Y1_Rem_auc_CRISPFEP_null, roc(predictor = results_Y1_Rem_CRISPFEP_seq$pred, response = Y1_Rem_CRISPFEP_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_Rem_auc_CRISPFEP_null >= auc_Y1_Rem_CRISPFEP_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_Rem_CRISPFEP_seq$pred, results_Y1_Rem_CRISPFEP_seq$obs=="Yes")
#95CI = 1.96*SE
1.96*sqrt(roc.area.test(results_Y1_Rem_CRISPFEP_seq$pred, results_Y1_Rem_CRISPFEP_seq$obs=="Yes")$var)

#Build GLM model with CRISPFEP columns on whole EDEN as training set for external validation on Scottish data
prep_Y1_Rem_CRISPFEP_all = prep_Y1_Rem_CRISPFEP_Stand[ ,!(colnames(prep_Y1_Rem_CRISPFEP_Stand) %in% c("Site"))]
controlAll <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary)
set.seed(987)
mod_Y1_Rem_EDEN_CF_all = train(M12_PANSS_Period_Rem ~ ., data=prep_Y1_Rem_CRISPFEP_all, method="glm", metric="ROC", preProc = c("knnImpute"), trControl=controlAll, na.action = na.pass)
#remove training data for sharing model
mod_Y1_Rem_EDEN_CF_all$trainingData = NULL
save(mod_Y1_Rem_EDEN_CF_all, file = "mod_Y1_Rem.rda")

#####################################
#Y1 GAF Outcome and Models
#
#
results_Y1_GAF = list()
mods_Y1_GAF = list()

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  
  #set up leave one site out cv
  test_set = prep_Y1_GAF[ which(prep_Y1_GAF$Site == levels(prep_Y1_GAF$Site)[i]), ]
  train_set = prep_Y1_GAF[ -which(prep_Y1_GAF$Site == levels(prep_Y1_GAF$Site)[i]), ]
  
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  train_set = train_set[ ,!(colnames(train_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_GAF_Binary) 
  
  #train model over grid of 10 by 10 lambda alpha, knn impute, standardise
  mod <- train(M12_GAF_Binary ~ ., data=train_set, method="glmnet", metric="ROC", tuneLength = 10, preProc = c("center", "scale","knnImpute"), trControl=control,na.action = na.pass)
  result$pred = predict(mod, test_set, type = "prob", na.action = na.pass)
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_GAF[[i]] = result
  mods_Y1_GAF[[i]] = mod
}

results_Y1_GAF_seq = NULL
results_Y1_GAF_seq$pred = c(results_Y1_GAF[[1]]$pred$Y,results_Y1_GAF[[2]]$pred$Y, results_Y1_GAF[[3]]$pred$Y,
                            results_Y1_GAF[[4]]$pred$Y,results_Y1_GAF[[5]]$pred$Y, results_Y1_GAF[[6]]$pred$Y,
                            results_Y1_GAF[[7]]$pred$Y,results_Y1_GAF[[8]]$pred$Y,results_Y1_GAF[[9]]$pred$Y,
                            results_Y1_GAF[[10]]$pred$Y,results_Y1_GAF[[11]]$pred$Y,results_Y1_GAF[[12]]$pred$Y,
                            results_Y1_GAF[[13]]$pred$Y,results_Y1_GAF[[14]]$pred$Y)
results_Y1_GAF_seq$obs = factor(c(as.character(results_Y1_GAF[[1]]$obs),as.character(results_Y1_GAF[[2]]$obs),
                                  as.character(results_Y1_GAF[[3]]$obs),as.character(results_Y1_GAF[[4]]$obs),
                                  as.character(results_Y1_GAF[[5]]$obs),as.character(results_Y1_GAF[[6]]$obs),
                                  as.character(results_Y1_GAF[[7]]$obs),as.character(results_Y1_GAF[[8]]$obs),
                                  as.character(results_Y1_GAF[[9]]$obs),as.character(results_Y1_GAF[[10]]$obs),
                                  as.character(results_Y1_GAF[[11]]$obs),as.character(results_Y1_GAF[[12]]$obs),
                                  as.character(results_Y1_GAF[[13]]$obs),as.character(results_Y1_GAF[[14]]$obs)))
results_Y1_GAF_seq$row = c(results_Y1_GAF[[1]]$row,results_Y1_GAF[[2]]$row, results_Y1_GAF[[3]]$row,
                           results_Y1_GAF[[4]]$row,results_Y1_GAF[[5]]$row, results_Y1_GAF[[6]]$row,
                           results_Y1_GAF[[7]]$row,results_Y1_GAF[[8]]$row,results_Y1_GAF[[9]]$row,
                           results_Y1_GAF[[10]]$row,results_Y1_GAF[[11]]$row,results_Y1_GAF[[12]]$row,
                           results_Y1_GAF[[13]]$row,results_Y1_GAF[[14]]$row)

auc_Y1_GAF_seq = roc(predictor = results_Y1_GAF_seq$pred, response = results_Y1_GAF_seq$obs)$auc
roc_Y1_GAF_seq = roc(predictor = results_Y1_GAF_seq$pred, response = results_Y1_GAF_seq$obs)

plot(roc_Y1_GAF_seq)
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
set.seed(987)
ci.coords(roc_Y1_GAF_seq,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

png("figure1B.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(roc_Y1_GAF_seq, auc.polygon=TRUE, grid=c(0.1, 0.1))
text(0.3,0.3,"AUC = 0·731\n(0·697, 0·765)\np<0·0001")
dev.off()

#permutation test for the overall auc vs null distribution
set.seed(987)
Y1_GAF_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_GAF_perm = permute(results_Y1_GAF_seq$obs)
  Y1_GAF_auc_null = c(Y1_GAF_auc_null, roc(predictor = results_Y1_GAF_seq$pred, response = Y1_GAF_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_GAF_auc_null >= auc_Y1_GAF_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_GAF_seq$pred, results_Y1_GAF_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_GAF_seq$pred, results_Y1_GAF_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_GAF_seq$pred, results_Y1_GAF_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_GAF_seq$pred, results_Y1_GAF_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_GAF_seq$pred, results_Y1_GAF_seq$obs=="Yes")$var))

#Test GAF model for each of the other three outcomes on EDEN data
#Y1 EET
#
results_Y1_GAF_EET = list()
#Add columns back in that were removed and required for this model
prep_Y1_EET_GAF = prep_Y1_EET
setdiff(colnames(mods_Y1_GAF[[1]]$trainingData),colnames(prep_Y1_EET_GAF))

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  #set up leave one site out cv
  test_set = prep_Y1_EET_GAF[ which(prep_Y1_EET_GAF$Site == levels(prep_Y1_EET_GAF$Site)[i]), ]
  
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_EET) 
  result$pred = predict(mods_Y1_GAF[[i]], test_set, type = "prob", na.action = na.pass)
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_GAF_EET[[i]] = result
}

results_Y1_GAF_EET_seq = NULL
results_Y1_GAF_EET_seq$pred = c(results_Y1_GAF_EET[[1]]$pred$Y,results_Y1_GAF_EET[[2]]$pred$Y, results_Y1_GAF_EET[[3]]$pred$Y,
                                results_Y1_GAF_EET[[4]]$pred$Y,results_Y1_GAF_EET[[5]]$pred$Y, results_Y1_GAF_EET[[6]]$pred$Y,
                                results_Y1_GAF_EET[[7]]$pred$Y,results_Y1_GAF_EET[[8]]$pred$Y,results_Y1_GAF_EET[[9]]$pred$Y,
                                results_Y1_GAF_EET[[10]]$pred$Y,results_Y1_GAF_EET[[11]]$pred$Y,results_Y1_GAF_EET[[12]]$pred$Y,
                                results_Y1_GAF_EET[[13]]$pred$Y,results_Y1_GAF_EET[[14]]$pred$Y)
results_Y1_GAF_EET_seq$obs = factor(c(as.character(results_Y1_GAF_EET[[1]]$obs),as.character(results_Y1_GAF_EET[[2]]$obs),
                                      as.character(results_Y1_GAF_EET[[3]]$obs),as.character(results_Y1_GAF_EET[[4]]$obs),
                                      as.character(results_Y1_GAF_EET[[5]]$obs),as.character(results_Y1_GAF_EET[[6]]$obs),
                                      as.character(results_Y1_GAF_EET[[7]]$obs),as.character(results_Y1_GAF_EET[[8]]$obs),
                                      as.character(results_Y1_GAF_EET[[9]]$obs),as.character(results_Y1_GAF_EET[[10]]$obs),
                                      as.character(results_Y1_GAF_EET[[11]]$obs),as.character(results_Y1_GAF_EET[[12]]$obs),
                                      as.character(results_Y1_GAF_EET[[13]]$obs),as.character(results_Y1_GAF_EET[[14]]$obs)))
results_Y1_GAF_EET_seq$row = c(results_Y1_GAF_EET[[1]]$row,results_Y1_GAF_EET[[2]]$row, results_Y1_GAF_EET[[3]]$row,
                               results_Y1_GAF_EET[[4]]$row,results_Y1_GAF_EET[[5]]$row, results_Y1_GAF_EET[[6]]$row,
                               results_Y1_GAF_EET[[7]]$row,results_Y1_GAF_EET[[8]]$row,results_Y1_GAF_EET[[9]]$row,
                               results_Y1_GAF_EET[[10]]$row,results_Y1_GAF_EET[[11]]$row,results_Y1_GAF_EET[[12]]$row,
                               results_Y1_GAF_EET[[13]]$row,results_Y1_GAF_EET[[14]]$row)

auc_Y1_GAF_EET_seq = roc(predictor = results_Y1_GAF_EET_seq$pred, response = results_Y1_GAF_EET_seq$obs)$auc
roc_Y1_GAF_EET_seq = roc(predictor = results_Y1_GAF_EET_seq$pred, response = results_Y1_GAF_EET_seq$obs)

set.seed(987)
Y1_GAF_EET_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_GAF_EET_perm = permute(results_Y1_GAF_EET_seq$obs)
  Y1_GAF_EET_auc_null = c(Y1_GAF_EET_auc_null, roc(predictor = results_Y1_GAF_EET_seq$pred, response = Y1_GAF_EET_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_GAF_EET_auc_null >= auc_Y1_GAF_EET_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_GAF_EET_seq$pred, results_Y1_GAF_EET_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_GAF_EET_seq$pred, results_Y1_GAF_EET_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_GAF_EET_seq$pred, results_Y1_GAF_EET_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_GAF_EET_seq$pred, results_Y1_GAF_EET_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_GAF_EET_seq$pred, results_Y1_GAF_EET_seq$obs=="Yes")$var))

#Y1 Rem
#
results_Y1_GAF_Rem = list()
#Add columns back in that were removed and required for this model
prep_Y1_Rem_GAF = prep_Y1_Rem
setdiff(colnames(mods_Y1_GAF[[1]]$trainingData),colnames(prep_Y1_Rem_GAF))
prep_Y1_Rem_GAF$PW1_First_Contact.Police = NA
prep_Y1_Rem_GAF$PW1_First_Contact.Police = as.numeric(prep_Y1_Rem_GAF$PW1_First_Contact.Police)

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  #set up leave one site out cv
  test_set = prep_Y1_Rem_GAF[ which(prep_Y1_Rem_GAF$Site == levels(prep_Y1_Rem_GAF$Site)[i]), ]
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_PANSS_Period_Rem) 
  result$pred = predict(mods_Y1_GAF[[i]], test_set, type = "prob", na.action = na.pass)

  result$row = as.numeric(rownames(test_set))
  
  results_Y1_GAF_Rem[[i]] = result
}

results_Y1_GAF_Rem_seq = NULL
results_Y1_GAF_Rem_seq$pred = c(results_Y1_GAF_Rem[[1]]$pred$Y,results_Y1_GAF_Rem[[2]]$pred$Y, results_Y1_GAF_Rem[[3]]$pred$Y,
                                results_Y1_GAF_Rem[[4]]$pred$Y,results_Y1_GAF_Rem[[5]]$pred$Y, results_Y1_GAF_Rem[[6]]$pred$Y,
                                results_Y1_GAF_Rem[[7]]$pred$Y,results_Y1_GAF_Rem[[8]]$pred$Y,results_Y1_GAF_Rem[[9]]$pred$Y,
                                results_Y1_GAF_Rem[[10]]$pred$Y,results_Y1_GAF_Rem[[11]]$pred$Y,results_Y1_GAF_Rem[[12]]$pred$Y,
                                results_Y1_GAF_Rem[[13]]$pred$Y,results_Y1_GAF_Rem[[14]]$pred$Y)
results_Y1_GAF_Rem_seq$obs = factor(c(as.character(results_Y1_GAF_Rem[[1]]$obs),as.character(results_Y1_GAF_Rem[[2]]$obs),
                                      as.character(results_Y1_GAF_Rem[[3]]$obs),as.character(results_Y1_GAF_Rem[[4]]$obs),
                                      as.character(results_Y1_GAF_Rem[[5]]$obs),as.character(results_Y1_GAF_Rem[[6]]$obs),
                                      as.character(results_Y1_GAF_Rem[[7]]$obs),as.character(results_Y1_GAF_Rem[[8]]$obs),
                                      as.character(results_Y1_GAF_Rem[[9]]$obs),as.character(results_Y1_GAF_Rem[[10]]$obs),
                                      as.character(results_Y1_GAF_Rem[[11]]$obs),as.character(results_Y1_GAF_Rem[[12]]$obs),
                                      as.character(results_Y1_GAF_Rem[[13]]$obs),as.character(results_Y1_GAF_Rem[[14]]$obs)))
results_Y1_GAF_Rem_seq$row = c(results_Y1_GAF_Rem[[1]]$row,results_Y1_GAF_Rem[[2]]$row, results_Y1_GAF_Rem[[3]]$row,
                               results_Y1_GAF_Rem[[4]]$row,results_Y1_GAF_Rem[[5]]$row, results_Y1_GAF_Rem[[6]]$row,
                               results_Y1_GAF_Rem[[7]]$row,results_Y1_GAF_Rem[[8]]$row,results_Y1_GAF_Rem[[9]]$row,
                               results_Y1_GAF_Rem[[10]]$row,results_Y1_GAF_Rem[[11]]$row,results_Y1_GAF_Rem[[12]]$row,
                               results_Y1_GAF_Rem[[13]]$row,results_Y1_GAF_Rem[[14]]$row)

auc_Y1_GAF_Rem_seq = roc(predictor = results_Y1_GAF_Rem_seq$pred, response = results_Y1_GAF_Rem_seq$obs)$auc
roc_Y1_GAF_Rem_seq = roc(predictor = results_Y1_GAF_Rem_seq$pred, response = results_Y1_GAF_Rem_seq$obs)

set.seed(987)
Y1_GAF_Rem_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_GAF_Rem_perm = permute(results_Y1_GAF_Rem_seq$obs)
  Y1_GAF_Rem_auc_null = c(Y1_GAF_Rem_auc_null, roc(predictor = results_Y1_GAF_Rem_seq$pred, response = Y1_GAF_Rem_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_GAF_Rem_auc_null >= auc_Y1_GAF_Rem_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_GAF_Rem_seq$pred, results_Y1_GAF_Rem_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_GAF_Rem_seq$pred, results_Y1_GAF_Rem_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_GAF_Rem_seq$pred, results_Y1_GAF_Rem_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_GAF_Rem_seq$pred, results_Y1_GAF_Rem_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_GAF_Rem_seq$pred, results_Y1_GAF_Rem_seq$obs=="Yes")$var))

#Y1 EQ5D_TTO
#
results_Y1_GAF_EQ5D_TTO = list()
#Add columns back in that were removed and required for this model
prep_Y1_EQ5D_TTO_GAF = prep_Y1_EQ5D_TTO
setdiff(colnames(mods_Y1_GAF[[1]]$trainingData),colnames(prep_Y1_EQ5D_TTO_GAF))

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  #set up leave one site out cv
  test_set = prep_Y1_EQ5D_TTO_GAF[ which(prep_Y1_EQ5D_TTO_GAF$Site == levels(prep_Y1_EQ5D_TTO_GAF$Site)[i]), ]
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_EQ5D_UK_TTO_Index_Binary) 
  result$pred = predict(mods_Y1_GAF[[i]], test_set, type = "prob", na.action = na.pass)
  
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_GAF_EQ5D_TTO[[i]] = result
}

results_Y1_GAF_EQ5D_TTO_seq = NULL
results_Y1_GAF_EQ5D_TTO_seq$pred = c(results_Y1_GAF_EQ5D_TTO[[1]]$pred$Y,results_Y1_GAF_EQ5D_TTO[[2]]$pred$Y, results_Y1_GAF_EQ5D_TTO[[3]]$pred$Y,
                                results_Y1_GAF_EQ5D_TTO[[4]]$pred$Y,results_Y1_GAF_EQ5D_TTO[[5]]$pred$Y, results_Y1_GAF_EQ5D_TTO[[6]]$pred$Y,
                                results_Y1_GAF_EQ5D_TTO[[7]]$pred$Y,results_Y1_GAF_EQ5D_TTO[[8]]$pred$Y,results_Y1_GAF_EQ5D_TTO[[9]]$pred$Y,
                                results_Y1_GAF_EQ5D_TTO[[10]]$pred$Y,results_Y1_GAF_EQ5D_TTO[[11]]$pred$Y,results_Y1_GAF_EQ5D_TTO[[12]]$pred$Y,
                                results_Y1_GAF_EQ5D_TTO[[13]]$pred$Y,results_Y1_GAF_EQ5D_TTO[[14]]$pred$Y)
results_Y1_GAF_EQ5D_TTO_seq$obs = factor(c(as.character(results_Y1_GAF_EQ5D_TTO[[1]]$obs),as.character(results_Y1_GAF_EQ5D_TTO[[2]]$obs),
                                      as.character(results_Y1_GAF_EQ5D_TTO[[3]]$obs),as.character(results_Y1_GAF_EQ5D_TTO[[4]]$obs),
                                      as.character(results_Y1_GAF_EQ5D_TTO[[5]]$obs),as.character(results_Y1_GAF_EQ5D_TTO[[6]]$obs),
                                      as.character(results_Y1_GAF_EQ5D_TTO[[7]]$obs),as.character(results_Y1_GAF_EQ5D_TTO[[8]]$obs),
                                      as.character(results_Y1_GAF_EQ5D_TTO[[9]]$obs),as.character(results_Y1_GAF_EQ5D_TTO[[10]]$obs),
                                      as.character(results_Y1_GAF_EQ5D_TTO[[11]]$obs),as.character(results_Y1_GAF_EQ5D_TTO[[12]]$obs),
                                      as.character(results_Y1_GAF_EQ5D_TTO[[13]]$obs),as.character(results_Y1_GAF_EQ5D_TTO[[14]]$obs)))
results_Y1_GAF_EQ5D_TTO_seq$row = c(results_Y1_GAF_EQ5D_TTO[[1]]$row,results_Y1_GAF_EQ5D_TTO[[2]]$row, results_Y1_GAF_EQ5D_TTO[[3]]$row,
                               results_Y1_GAF_EQ5D_TTO[[4]]$row,results_Y1_GAF_EQ5D_TTO[[5]]$row, results_Y1_GAF_EQ5D_TTO[[6]]$row,
                               results_Y1_GAF_EQ5D_TTO[[7]]$row,results_Y1_GAF_EQ5D_TTO[[8]]$row,results_Y1_GAF_EQ5D_TTO[[9]]$row,
                               results_Y1_GAF_EQ5D_TTO[[10]]$row,results_Y1_GAF_EQ5D_TTO[[11]]$row,results_Y1_GAF_EQ5D_TTO[[12]]$row,
                               results_Y1_GAF_EQ5D_TTO[[13]]$row,results_Y1_GAF_EQ5D_TTO[[14]]$row)

auc_Y1_GAF_EQ5D_TTO_seq = roc(predictor = results_Y1_GAF_EQ5D_TTO_seq$pred, response = results_Y1_GAF_EQ5D_TTO_seq$obs)$auc
roc_Y1_GAF_EQ5D_TTO_seq = roc(predictor = results_Y1_GAF_EQ5D_TTO_seq$pred, response = results_Y1_GAF_EQ5D_TTO_seq$obs)

set.seed(987)
Y1_GAF_EQ5D_TTO_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_GAF_EQ5D_TTO_perm = permute(results_Y1_GAF_EQ5D_TTO_seq$obs)
  Y1_GAF_EQ5D_TTO_auc_null = c(Y1_GAF_EQ5D_TTO_auc_null, roc(predictor = results_Y1_GAF_EQ5D_TTO_seq$pred, response = Y1_GAF_EQ5D_TTO_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_GAF_EQ5D_TTO_auc_null >= auc_Y1_GAF_EQ5D_TTO_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_GAF_EQ5D_TTO_seq$pred, results_Y1_GAF_EQ5D_TTO_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_GAF_EQ5D_TTO_seq$pred, results_Y1_GAF_EQ5D_TTO_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_GAF_EQ5D_TTO_seq$pred, results_Y1_GAF_EQ5D_TTO_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_GAF_EQ5D_TTO_seq$pred, results_Y1_GAF_EQ5D_TTO_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_GAF_EQ5D_TTO_seq$pred, results_Y1_GAF_EQ5D_TTO_seq$obs=="Yes")$var))

coefs_Y1_GAF = NULL
for (i in seq(1:14))
{
  coefs_Y1_GAF = c(coefs_Y1_GAF, coef(mods_Y1_GAF[[i]]$finalModel, mods_Y1_GAF[[i]]$bestTune$lambda))
}

#just get numbers
coefs_Y1_GAF_extract = NULL
for(i in seq(1:14))
{
  coefs_Y1_GAF_extract = rbind(coefs_Y1_GAF_extract, coefs_Y1_GAF[[i]][1:164])
}

#Presence or absence of predictors across models
coefs_Y1_GAF_presence = NULL
coefs_Y1_GAF_presence = coefs_Y1_GAF_extract[1:14,2:164]
coefs_Y1_GAF_presence[coefs_Y1_GAF_presence != 0] <- 1

#stability of feature selection http://jmlr.org/papers/volume18/17-514/17-514.pdf
getStability(coefs_Y1_GAF_presence)

#work out number of predictors shared across all models
table(colMeans(coefs_Y1_GAF_presence))

#get rank of coef by importance as in sports ranking
coefs_Y1_GAF_baseline_rank = NULL

for(i in seq(c(1:14)))
{
  #rank absolute value excluding the intercept for each model
  coefs_Y1_GAF_baseline_rank = rbind(coefs_Y1_GAF_baseline_rank, rank(abs(coefs_Y1_GAF_extract[i,2:164]), ties.method = "min"))
}

# rank the mean ranks of each column across all models
coefs_Y1_GAF_baseline_rank_mean = colMeans(coefs_Y1_GAF_baseline_rank)

#Invert order of rank to identify top models
coefs_Y1_GAF_baseline_order = rank(-coefs_Y1_GAF_baseline_rank_mean)

#Get the column names (not the intercept)
coef_Y1_GAF_names = dimnames(coefs_Y1_GAF[[1]])[[1]][2:164]

coefs_Y1_GAF_means = colMeans(coefs_Y1_GAF_extract)[2:164]

df_Y1_GAF_names = data.frame(coef_Y1_GAF_names, coefs_Y1_GAF_baseline_order, coefs_Y1_GAF_means, colMeans(coefs_Y1_GAF_presence))

#View the best predictor variables by absolute value in order (lower is better)
View(df_Y1_GAF_names)

#Internal-External Validation for OPUS only columns to check still works
#tried glmnet, glm, rf, linear and radial svm - all very similar. Just use GLM
#Y1 GAF
#
results_Y1_GAF_OPUS = list()
mods_Y1_GAF_OPUS = list()

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  
  #set up leave one site out cv
  test_set = prep_Y1_GAF_OPUS_Stand[ which(prep_Y1_GAF_OPUS_Stand$Site == levels(prep_Y1_GAF_OPUS_Stand$Site)[i]), ]
  train_set = prep_Y1_GAF_OPUS_Stand[ -which(prep_Y1_GAF_OPUS_Stand$Site == levels(prep_Y1_GAF_OPUS_Stand$Site)[i]), ]
  
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  train_set = train_set[ ,!(colnames(train_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_GAF_Binary) 
  
  #train model over grid of 10 by 10 lambda alpha, knn impute, standardise
  mod <- train(M12_GAF_Binary ~ ., data=train_set, method="glm", metric="ROC", tuneLength = 10, preProc = c("knnImpute"), trControl=control,na.action = na.pass)
  result$pred = predict(mod, test_set, type = "prob", na.action = na.pass)
  
  results_Y1_GAF_OPUS[[i]] = result
  mods_Y1_GAF_OPUS[[i]] = mod
}

results_Y1_GAF_OPUS_seq = NULL
results_Y1_GAF_OPUS_seq$pred = c(results_Y1_GAF_OPUS[[1]]$pred$Y,results_Y1_GAF_OPUS[[2]]$pred$Y, results_Y1_GAF_OPUS[[3]]$pred$Y,
                                 results_Y1_GAF_OPUS[[4]]$pred$Y,results_Y1_GAF_OPUS[[5]]$pred$Y, results_Y1_GAF_OPUS[[6]]$pred$Y,
                                 results_Y1_GAF_OPUS[[7]]$pred$Y,results_Y1_GAF_OPUS[[8]]$pred$Y,results_Y1_GAF_OPUS[[9]]$pred$Y,
                                 results_Y1_GAF_OPUS[[10]]$pred$Y,results_Y1_GAF_OPUS[[11]]$pred$Y,results_Y1_GAF_OPUS[[12]]$pred$Y,
                                 results_Y1_GAF_OPUS[[13]]$pred$Y,results_Y1_GAF_OPUS[[14]]$pred$Y)
results_Y1_GAF_OPUS_seq$obs = factor(c(as.character(results_Y1_GAF_OPUS[[1]]$obs),as.character(results_Y1_GAF_OPUS[[2]]$obs),
                                       as.character(results_Y1_GAF_OPUS[[3]]$obs),as.character(results_Y1_GAF_OPUS[[4]]$obs),
                                       as.character(results_Y1_GAF_OPUS[[5]]$obs),as.character(results_Y1_GAF_OPUS[[6]]$obs),
                                       as.character(results_Y1_GAF_OPUS[[7]]$obs),as.character(results_Y1_GAF_OPUS[[8]]$obs),
                                       as.character(results_Y1_GAF_OPUS[[9]]$obs),as.character(results_Y1_GAF_OPUS[[10]]$obs),
                                       as.character(results_Y1_GAF_OPUS[[11]]$obs),as.character(results_Y1_GAF_OPUS[[12]]$obs),
                                       as.character(results_Y1_GAF_OPUS[[13]]$obs),as.character(results_Y1_GAF_OPUS[[14]]$obs)))

auc_Y1_GAF_OPUS_seq = roc(predictor = results_Y1_GAF_OPUS_seq$pred, response = results_Y1_GAF_OPUS_seq$obs)$auc
roc_Y1_GAF_OPUS_seq = roc(predictor = results_Y1_GAF_OPUS_seq$pred, response = results_Y1_GAF_OPUS_seq$obs)

plot(roc_Y1_GAF_OPUS_seq)
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
coords(roc_Y1_GAF_OPUS_seq,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

set.seed(987)
Y1_GAF_auc_OPUS_null = NULL
for(i in seq (1:10001))
{
  Y1_GAF_OPUS_perm = permute(results_Y1_GAF_OPUS_seq$obs)
  Y1_GAF_auc_OPUS_null = c(Y1_GAF_auc_OPUS_null, roc(predictor = results_Y1_GAF_OPUS_seq$pred, response = Y1_GAF_OPUS_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_GAF_auc_OPUS_null >= auc_Y1_GAF_OPUS_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_GAF_OPUS_seq$pred, results_Y1_GAF_OPUS_seq$obs=="Yes")
#95CI = 1.96*SE
1.96*sqrt(roc.area.test(results_Y1_GAF_OPUS_seq$pred, results_Y1_GAF_OPUS_seq$obs=="Yes")$var)

#Build GLM model with OPUS columns on whole EDEN as training set for external validation
prep_Y1_GAF_OPUS_all = prep_Y1_GAF_OPUS_Stand[ ,!(colnames(prep_Y1_GAF_OPUS_Stand) %in% c("Site"))]
controlAll <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary)
set.seed(987)
mod_Y1_GAF_EDEN_all = train(M12_GAF_Binary ~ ., data=prep_Y1_GAF_OPUS_all, method="glm", metric="ROC", preProc = c("knnImpute"), trControl=controlAll, na.action = na.pass)
#remove training data for sharing model
mod_Y1_GAF_EDEN_all$trainingData = NULL
save(mod_Y1_GAF_EDEN_all, file = "mod_Y1_GAF_OPUS.rda")

###############################################
#Y1 EET Outcome and Models
#
#
results_Y1_EET = list()
mods_Y1_EET = list()

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  
  #set up leave one site out cv
  test_set = prep_Y1_EET[ which(prep_Y1_EET$Site == levels(prep_Y1_EET$Site)[i]), ]
  train_set = prep_Y1_EET[ -which(prep_Y1_EET$Site == levels(prep_Y1_EET$Site)[i]), ]
  
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  train_set = train_set[ ,!(colnames(train_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_EET) 
  
  #train model over grid of 10 by 10 lambda alpha, knn impute, standardise
  mod <- train(M12_EET ~ ., data=train_set, method="glmnet", metric="ROC", tuneLength = 10, preProc = c("center", "scale","knnImpute"), trControl=control,na.action = na.pass)
  result$pred = predict(mod, test_set, type = "prob", na.action = na.pass)
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_EET[[i]] = result
  mods_Y1_EET[[i]] = mod
}

results_Y1_EET_seq = NULL
results_Y1_EET_seq$pred = c(results_Y1_EET[[1]]$pred$Y,results_Y1_EET[[2]]$pred$Y, results_Y1_EET[[3]]$pred$Y,
                            results_Y1_EET[[4]]$pred$Y,results_Y1_EET[[5]]$pred$Y, results_Y1_EET[[6]]$pred$Y,
                            results_Y1_EET[[7]]$pred$Y,results_Y1_EET[[8]]$pred$Y,results_Y1_EET[[9]]$pred$Y,
                            results_Y1_EET[[10]]$pred$Y,results_Y1_EET[[11]]$pred$Y,results_Y1_EET[[12]]$pred$Y,
                            results_Y1_EET[[13]]$pred$Y,results_Y1_EET[[14]]$pred$Y)
results_Y1_EET_seq$obs = factor(c(as.character(results_Y1_EET[[1]]$obs),as.character(results_Y1_EET[[2]]$obs),
                                  as.character(results_Y1_EET[[3]]$obs),as.character(results_Y1_EET[[4]]$obs),
                                  as.character(results_Y1_EET[[5]]$obs),as.character(results_Y1_EET[[6]]$obs),
                                  as.character(results_Y1_EET[[7]]$obs),as.character(results_Y1_EET[[8]]$obs),
                                  as.character(results_Y1_EET[[9]]$obs),as.character(results_Y1_EET[[10]]$obs),
                                  as.character(results_Y1_EET[[11]]$obs),as.character(results_Y1_EET[[12]]$obs),
                                  as.character(results_Y1_EET[[13]]$obs),as.character(results_Y1_EET[[14]]$obs)))
results_Y1_EET_seq$row = c(results_Y1_EET[[1]]$row,results_Y1_EET[[2]]$row, results_Y1_EET[[3]]$row,
                           results_Y1_EET[[4]]$row,results_Y1_EET[[5]]$row, results_Y1_EET[[6]]$row,
                           results_Y1_EET[[7]]$row,results_Y1_EET[[8]]$row,results_Y1_EET[[9]]$row,
                           results_Y1_EET[[10]]$row,results_Y1_EET[[11]]$row,results_Y1_EET[[12]]$row,
                           results_Y1_EET[[13]]$row,results_Y1_EET[[14]]$row)

auc_Y1_EET_seq = roc(predictor = results_Y1_EET_seq$pred, response = results_Y1_EET_seq$obs)$auc
roc_Y1_EET_seq = roc(predictor = results_Y1_EET_seq$pred, response = results_Y1_EET_seq$obs)

plot(roc_Y1_EET_seq)
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
set.seed(987)
ci.coords(roc_Y1_EET_seq,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

png("figure1C.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(roc_Y1_EET_seq, auc.polygon=TRUE, grid=c(0.1, 0.1))
text(0.3,0.3,"AUC = 0·736\n(0·702, 0·771)\np<0·0001")
dev.off()

#permutation test for the overall auc vs null distribution
set.seed(987)
Y1_EET_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_EET_perm = permute(results_Y1_EET_seq$obs)
  Y1_EET_auc_null = c(Y1_EET_auc_null, roc(predictor = results_Y1_EET_seq$pred, response = Y1_EET_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_EET_auc_null >= auc_Y1_EET_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_EET_seq$pred, results_Y1_EET_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_EET_seq$pred, results_Y1_EET_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_EET_seq$pred, results_Y1_EET_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_EET_seq$pred, results_Y1_EET_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_EET_seq$pred, results_Y1_EET_seq$obs=="Yes")$var))

#Test model on different outcomes
#Y1 Rem
#
results_Y1_EET_Rem = list()
#Add columns back in that were removed and required for this model
prep_Y1_Rem_EET = prep_Y1_Rem
setdiff(colnames(mods_Y1_EET[[1]]$trainingData),colnames(prep_Y1_Rem_EET))
prep_Y1_Rem_EET$PW1_First_Contact.Police = NA
prep_Y1_Rem_EET$PW1_First_Contact.Police = as.numeric(prep_Y1_Rem_EET$PW1_First_Contact.Police)

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  #set up leave one site out cv
  test_set = prep_Y1_Rem_EET[ which(prep_Y1_Rem_EET$Site == levels(prep_Y1_Rem_EET$Site)[i]), ]
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_PANSS_Period_Rem) 
  result$pred = predict(mods_Y1_EET[[i]], test_set, type = "prob", na.action = na.pass)
  
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_EET_Rem[[i]] = result
}

results_Y1_EET_Rem_seq = NULL
results_Y1_EET_Rem_seq$pred = c(results_Y1_EET_Rem[[1]]$pred$Y,results_Y1_EET_Rem[[2]]$pred$Y, results_Y1_EET_Rem[[3]]$pred$Y,
                                results_Y1_EET_Rem[[4]]$pred$Y,results_Y1_EET_Rem[[5]]$pred$Y, results_Y1_EET_Rem[[6]]$pred$Y,
                                results_Y1_EET_Rem[[7]]$pred$Y,results_Y1_EET_Rem[[8]]$pred$Y,results_Y1_EET_Rem[[9]]$pred$Y,
                                results_Y1_EET_Rem[[10]]$pred$Y,results_Y1_EET_Rem[[11]]$pred$Y,results_Y1_EET_Rem[[12]]$pred$Y,
                                results_Y1_EET_Rem[[13]]$pred$Y,results_Y1_EET_Rem[[14]]$pred$Y)
results_Y1_EET_Rem_seq$obs = factor(c(as.character(results_Y1_EET_Rem[[1]]$obs),as.character(results_Y1_EET_Rem[[2]]$obs),
                                      as.character(results_Y1_EET_Rem[[3]]$obs),as.character(results_Y1_EET_Rem[[4]]$obs),
                                      as.character(results_Y1_EET_Rem[[5]]$obs),as.character(results_Y1_EET_Rem[[6]]$obs),
                                      as.character(results_Y1_EET_Rem[[7]]$obs),as.character(results_Y1_EET_Rem[[8]]$obs),
                                      as.character(results_Y1_EET_Rem[[9]]$obs),as.character(results_Y1_EET_Rem[[10]]$obs),
                                      as.character(results_Y1_EET_Rem[[11]]$obs),as.character(results_Y1_EET_Rem[[12]]$obs),
                                      as.character(results_Y1_EET_Rem[[13]]$obs),as.character(results_Y1_EET_Rem[[14]]$obs)))
results_Y1_EET_Rem_seq$row = c(results_Y1_EET_Rem[[1]]$row,results_Y1_EET_Rem[[2]]$row, results_Y1_EET_Rem[[3]]$row,
                               results_Y1_EET_Rem[[4]]$row,results_Y1_EET_Rem[[5]]$row, results_Y1_EET_Rem[[6]]$row,
                               results_Y1_EET_Rem[[7]]$row,results_Y1_EET_Rem[[8]]$row,results_Y1_EET_Rem[[9]]$row,
                               results_Y1_EET_Rem[[10]]$row,results_Y1_EET_Rem[[11]]$row,results_Y1_EET_Rem[[12]]$row,
                               results_Y1_EET_Rem[[13]]$row,results_Y1_EET_Rem[[14]]$row)

auc_Y1_EET_Rem_seq = roc(predictor = results_Y1_EET_Rem_seq$pred, response = results_Y1_EET_Rem_seq$obs)$auc
roc_Y1_EET_Rem_seq = roc(predictor = results_Y1_EET_Rem_seq$pred, response = results_Y1_EET_Rem_seq$obs)

set.seed(987)
Y1_EET_Rem_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_EET_Rem_perm = permute(results_Y1_EET_Rem_seq$obs)
  Y1_EET_Rem_auc_null = c(Y1_EET_Rem_auc_null, roc(predictor = results_Y1_EET_Rem_seq$pred, response = Y1_EET_Rem_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_EET_Rem_auc_null >= auc_Y1_EET_Rem_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_EET_Rem_seq$pred, results_Y1_EET_Rem_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_EET_Rem_seq$pred, results_Y1_EET_Rem_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_EET_Rem_seq$pred, results_Y1_EET_Rem_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_EET_Rem_seq$pred, results_Y1_EET_Rem_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_EET_Rem_seq$pred, results_Y1_EET_Rem_seq$obs=="Yes")$var))

#Y1 GAF
#
results_Y1_EET_GAF = list()
#Add columns back in that were removed and required for this model
prep_Y1_GAF_EET = prep_Y1_GAF
setdiff(colnames(mods_Y1_EET[[1]]$trainingData),colnames(prep_Y1_GAF_EET))

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  #set up leave one site out cv
  test_set = prep_Y1_GAF_EET[ which(prep_Y1_GAF_EET$Site == levels(prep_Y1_GAF_EET$Site)[i]), ]
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_GAF) 
  result$pred = predict(mods_Y1_EET[[i]], test_set, type = "prob", na.action = na.pass)
  
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_EET_GAF[[i]] = result
}

results_Y1_EET_GAF_seq = NULL
results_Y1_EET_GAF_seq$pred = c(results_Y1_EET_GAF[[1]]$pred$Y,results_Y1_EET_GAF[[2]]$pred$Y, results_Y1_EET_GAF[[3]]$pred$Y,
                                results_Y1_EET_GAF[[4]]$pred$Y,results_Y1_EET_GAF[[5]]$pred$Y, results_Y1_EET_GAF[[6]]$pred$Y,
                                results_Y1_EET_GAF[[7]]$pred$Y,results_Y1_EET_GAF[[8]]$pred$Y,results_Y1_EET_GAF[[9]]$pred$Y,
                                results_Y1_EET_GAF[[10]]$pred$Y,results_Y1_EET_GAF[[11]]$pred$Y,results_Y1_EET_GAF[[12]]$pred$Y,
                                results_Y1_EET_GAF[[13]]$pred$Y,results_Y1_EET_GAF[[14]]$pred$Y)
results_Y1_EET_GAF_seq$obs = factor(c(as.character(results_Y1_EET_GAF[[1]]$obs),as.character(results_Y1_EET_GAF[[2]]$obs),
                                      as.character(results_Y1_EET_GAF[[3]]$obs),as.character(results_Y1_EET_GAF[[4]]$obs),
                                      as.character(results_Y1_EET_GAF[[5]]$obs),as.character(results_Y1_EET_GAF[[6]]$obs),
                                      as.character(results_Y1_EET_GAF[[7]]$obs),as.character(results_Y1_EET_GAF[[8]]$obs),
                                      as.character(results_Y1_EET_GAF[[9]]$obs),as.character(results_Y1_EET_GAF[[10]]$obs),
                                      as.character(results_Y1_EET_GAF[[11]]$obs),as.character(results_Y1_EET_GAF[[12]]$obs),
                                      as.character(results_Y1_EET_GAF[[13]]$obs),as.character(results_Y1_EET_GAF[[14]]$obs)))
results_Y1_EET_GAF_seq$row = c(results_Y1_EET_GAF[[1]]$row,results_Y1_EET_GAF[[2]]$row, results_Y1_EET_GAF[[3]]$row,
                               results_Y1_EET_GAF[[4]]$row,results_Y1_EET_GAF[[5]]$row, results_Y1_EET_GAF[[6]]$row,
                               results_Y1_EET_GAF[[7]]$row,results_Y1_EET_GAF[[8]]$row,results_Y1_EET_GAF[[9]]$row,
                               results_Y1_EET_GAF[[10]]$row,results_Y1_EET_GAF[[11]]$row,results_Y1_EET_GAF[[12]]$row,
                               results_Y1_EET_GAF[[13]]$row,results_Y1_EET_GAF[[14]]$row)

auc_Y1_EET_GAF_seq = roc(predictor = results_Y1_EET_GAF_seq$pred, response = results_Y1_EET_GAF_seq$obs)$auc
roc_Y1_EET_GAF_seq = roc(predictor = results_Y1_EET_GAF_seq$pred, response = results_Y1_EET_GAF_seq$obs)

set.seed(987)
Y1_EET_GAF_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_EET_GAF_perm = permute(results_Y1_EET_GAF_seq$obs)
  Y1_EET_GAF_auc_null = c(Y1_EET_GAF_auc_null, roc(predictor = results_Y1_EET_GAF_seq$pred, response = Y1_EET_GAF_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_EET_GAF_auc_null >= auc_Y1_EET_GAF_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_EET_GAF_seq$pred, results_Y1_EET_GAF_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_EET_GAF_seq$pred, results_Y1_EET_GAF_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_EET_GAF_seq$pred, results_Y1_EET_GAF_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_EET_GAF_seq$pred, results_Y1_EET_GAF_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_EET_GAF_seq$pred, results_Y1_EET_GAF_seq$obs=="Yes")$var))

#Y1 EQ5D_TTO
#
results_Y1_EET_EQ5D_TTO = list()
#Add columns back in that were removed and required for this model
prep_Y1_EQ5D_TTO_EET = prep_Y1_EQ5D_TTO
setdiff(colnames(mods_Y1_EET[[1]]$trainingData),colnames(prep_Y1_EQ5D_TTO_EET))

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  #set up leave one site out cv
  test_set = prep_Y1_EQ5D_TTO_EET[ which(prep_Y1_EQ5D_TTO_EET$Site == levels(prep_Y1_EQ5D_TTO_EET$Site)[i]), ]
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_EQ5D_UK_TTO_Index_Binary) 
  result$pred = predict(mods_Y1_EET[[i]], test_set, type = "prob", na.action = na.pass)
  
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_EET_EQ5D_TTO[[i]] = result
}

results_Y1_EET_EQ5D_TTO_seq = NULL
results_Y1_EET_EQ5D_TTO_seq$pred = c(results_Y1_EET_EQ5D_TTO[[1]]$pred$Y,results_Y1_EET_EQ5D_TTO[[2]]$pred$Y, results_Y1_EET_EQ5D_TTO[[3]]$pred$Y,
                                results_Y1_EET_EQ5D_TTO[[4]]$pred$Y,results_Y1_EET_EQ5D_TTO[[5]]$pred$Y, results_Y1_EET_EQ5D_TTO[[6]]$pred$Y,
                                results_Y1_EET_EQ5D_TTO[[7]]$pred$Y,results_Y1_EET_EQ5D_TTO[[8]]$pred$Y,results_Y1_EET_EQ5D_TTO[[9]]$pred$Y,
                                results_Y1_EET_EQ5D_TTO[[10]]$pred$Y,results_Y1_EET_EQ5D_TTO[[11]]$pred$Y,results_Y1_EET_EQ5D_TTO[[12]]$pred$Y,
                                results_Y1_EET_EQ5D_TTO[[13]]$pred$Y,results_Y1_EET_EQ5D_TTO[[14]]$pred$Y)
results_Y1_EET_EQ5D_TTO_seq$obs = factor(c(as.character(results_Y1_EET_EQ5D_TTO[[1]]$obs),as.character(results_Y1_EET_EQ5D_TTO[[2]]$obs),
                                      as.character(results_Y1_EET_EQ5D_TTO[[3]]$obs),as.character(results_Y1_EET_EQ5D_TTO[[4]]$obs),
                                      as.character(results_Y1_EET_EQ5D_TTO[[5]]$obs),as.character(results_Y1_EET_EQ5D_TTO[[6]]$obs),
                                      as.character(results_Y1_EET_EQ5D_TTO[[7]]$obs),as.character(results_Y1_EET_EQ5D_TTO[[8]]$obs),
                                      as.character(results_Y1_EET_EQ5D_TTO[[9]]$obs),as.character(results_Y1_EET_EQ5D_TTO[[10]]$obs),
                                      as.character(results_Y1_EET_EQ5D_TTO[[11]]$obs),as.character(results_Y1_EET_EQ5D_TTO[[12]]$obs),
                                      as.character(results_Y1_EET_EQ5D_TTO[[13]]$obs),as.character(results_Y1_EET_EQ5D_TTO[[14]]$obs)))
results_Y1_EET_EQ5D_TTO_seq$row = c(results_Y1_EET_EQ5D_TTO[[1]]$row,results_Y1_EET_EQ5D_TTO[[2]]$row, results_Y1_EET_EQ5D_TTO[[3]]$row,
                               results_Y1_EET_EQ5D_TTO[[4]]$row,results_Y1_EET_EQ5D_TTO[[5]]$row, results_Y1_EET_EQ5D_TTO[[6]]$row,
                               results_Y1_EET_EQ5D_TTO[[7]]$row,results_Y1_EET_EQ5D_TTO[[8]]$row,results_Y1_EET_EQ5D_TTO[[9]]$row,
                               results_Y1_EET_EQ5D_TTO[[10]]$row,results_Y1_EET_EQ5D_TTO[[11]]$row,results_Y1_EET_EQ5D_TTO[[12]]$row,
                               results_Y1_EET_EQ5D_TTO[[13]]$row,results_Y1_EET_EQ5D_TTO[[14]]$row)

auc_Y1_EET_EQ5D_TTO_seq = roc(predictor = results_Y1_EET_EQ5D_TTO_seq$pred, response = results_Y1_EET_EQ5D_TTO_seq$obs)$auc
roc_Y1_EET_EQ5D_TTO_seq = roc(predictor = results_Y1_EET_EQ5D_TTO_seq$pred, response = results_Y1_EET_EQ5D_TTO_seq$obs)

set.seed(987)
Y1_EET_EQ5D_TTO_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_EET_EQ5D_TTO_perm = permute(results_Y1_EET_EQ5D_TTO_seq$obs)
  Y1_EET_EQ5D_TTO_auc_null = c(Y1_EET_EQ5D_TTO_auc_null, roc(predictor = results_Y1_EET_EQ5D_TTO_seq$pred, response = Y1_EET_EQ5D_TTO_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_EET_EQ5D_TTO_auc_null >= auc_Y1_EET_EQ5D_TTO_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_EET_EQ5D_TTO_seq$pred, results_Y1_EET_EQ5D_TTO_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_EET_EQ5D_TTO_seq$pred, results_Y1_EET_EQ5D_TTO_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_EET_EQ5D_TTO_seq$pred, results_Y1_EET_EQ5D_TTO_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_EET_EQ5D_TTO_seq$pred, results_Y1_EET_EQ5D_TTO_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_EET_EQ5D_TTO_seq$pred, results_Y1_EET_EQ5D_TTO_seq$obs=="Yes")$var))

coefs_Y1_EET = NULL
for (i in seq(1:14))
{
  coefs_Y1_EET = c(coefs_Y1_EET, coef(mods_Y1_EET[[i]]$finalModel, mods_Y1_EET[[i]]$bestTune$lambda))
}

#just get numbers
coefs_Y1_EET_extract = NULL
for(i in seq(1:14))
{
  coefs_Y1_EET_extract = rbind(coefs_Y1_EET_extract, coefs_Y1_EET[[i]][1:166])
}

#Presence or absence of predictors across models
coefs_Y1_EET_presence = NULL
coefs_Y1_EET_presence = coefs_Y1_EET_extract[1:14,2:166]
coefs_Y1_EET_presence[coefs_Y1_EET_presence != 0] <- 1

#stability of feature selection http://jmlr.org/papers/volume18/17-514/17-514.pdf
getStability(coefs_Y1_EET_presence)

#work out number of predictors shared across all models
table(colMeans(coefs_Y1_EET_presence))
#number of features in each model
rowSums(coefs_Y1_EET_presence)

#get rank of coef by importance as in sports ranking
coefs_Y1_EET_baseline_rank = NULL

for(i in seq(c(1:14)))
{
  #rank absolute value excluding the intercept for each model
  coefs_Y1_EET_baseline_rank = rbind(coefs_Y1_EET_baseline_rank, rank(abs(coefs_Y1_EET_extract[i,2:166]), ties.method = "min"))
}

# rank the mean ranks of each column across all models
coefs_Y1_EET_baseline_rank_mean = colMeans(coefs_Y1_EET_baseline_rank)

#Invert order of rank to identify top models
coefs_Y1_EET_baseline_order = rank(-coefs_Y1_EET_baseline_rank_mean)

#Get the column names (not the intercept)
coef_Y1_EET_names = dimnames(coefs_Y1_EET[[1]])[[1]][2:166]

coefs_Y1_EET_means = colMeans(coefs_Y1_EET_extract)[2:166]

df_Y1_EET_names = data.frame(coef_Y1_EET_names, coefs_Y1_EET_baseline_order, coefs_Y1_EET_means, colMeans(coefs_Y1_EET_presence))

#View the best predictor variables by absolute value in order (lower is better)
View(df_Y1_EET_names)

#Internal-External Validation for CRISPFEP only columns to ensure still works
#Y1 EET
#
results_Y1_EET_CRISPFEP = list()
mods_Y1_EET_CRISPFEP = list()

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  
  #set up leave one site out cv
  test_set = prep_Y1_EET_CRISPFEP_Stand[ which(prep_Y1_EET_CRISPFEP_Stand$Site == levels(prep_Y1_EET_CRISPFEP_Stand$Site)[i]), ]
  train_set = prep_Y1_EET_CRISPFEP_Stand[ -which(prep_Y1_EET_CRISPFEP_Stand$Site == levels(prep_Y1_EET_CRISPFEP_Stand$Site)[i]), ]
  
  #Remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  train_set = train_set[ ,!(colnames(train_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_EET) 
  mod <- train(M12_EET ~ ., data=train_set, method="glm", metric="ROC", tuneLength = 10, preProc = c("knnImpute"), trControl=control,na.action = na.pass)
  result$pred = predict(mod, test_set, type = "prob", na.action = na.pass)
  
  results_Y1_EET_CRISPFEP[[i]] = result
  mods_Y1_EET_CRISPFEP[[i]] = mod
}

results_Y1_EET_CRISPFEP_seq = NULL
results_Y1_EET_CRISPFEP_seq$pred = c(results_Y1_EET_CRISPFEP[[1]]$pred$Y,results_Y1_EET_CRISPFEP[[2]]$pred$Y, results_Y1_EET_CRISPFEP[[3]]$pred$Y,
                                     results_Y1_EET_CRISPFEP[[4]]$pred$Y,results_Y1_EET_CRISPFEP[[5]]$pred$Y, results_Y1_EET_CRISPFEP[[6]]$pred$Y,
                                     results_Y1_EET_CRISPFEP[[7]]$pred$Y,results_Y1_EET_CRISPFEP[[8]]$pred$Y,results_Y1_EET_CRISPFEP[[9]]$pred$Y,
                                     results_Y1_EET_CRISPFEP[[10]]$pred$Y,results_Y1_EET_CRISPFEP[[11]]$pred$Y,results_Y1_EET_CRISPFEP[[12]]$pred$Y,
                                     results_Y1_EET_CRISPFEP[[13]]$pred$Y,results_Y1_EET_CRISPFEP[[14]]$pred$Y)
results_Y1_EET_CRISPFEP_seq$obs = factor(c(as.character(results_Y1_EET_CRISPFEP[[1]]$obs),as.character(results_Y1_EET_CRISPFEP[[2]]$obs),
                                           as.character(results_Y1_EET_CRISPFEP[[3]]$obs),as.character(results_Y1_EET_CRISPFEP[[4]]$obs),
                                           as.character(results_Y1_EET_CRISPFEP[[5]]$obs),as.character(results_Y1_EET_CRISPFEP[[6]]$obs),
                                           as.character(results_Y1_EET_CRISPFEP[[7]]$obs),as.character(results_Y1_EET_CRISPFEP[[8]]$obs),
                                           as.character(results_Y1_EET_CRISPFEP[[9]]$obs),as.character(results_Y1_EET_CRISPFEP[[10]]$obs),
                                           as.character(results_Y1_EET_CRISPFEP[[11]]$obs),as.character(results_Y1_EET_CRISPFEP[[12]]$obs),
                                           as.character(results_Y1_EET_CRISPFEP[[13]]$obs),as.character(results_Y1_EET_CRISPFEP[[14]]$obs)))

auc_Y1_EET_CRISPFEP_seq = roc(predictor = results_Y1_EET_CRISPFEP_seq$pred, response = results_Y1_EET_CRISPFEP_seq$obs)$auc
roc_Y1_EET_CRISPFEP_seq = roc(predictor = results_Y1_EET_CRISPFEP_seq$pred, response = results_Y1_EET_CRISPFEP_seq$obs)

plot(roc_Y1_EET_CRISPFEP_seq)
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
coords(roc_Y1_EET_CRISPFEP_seq,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

set.seed(987)
Y1_EET_auc_CRISPFEP_null = NULL
for(i in seq (1:10001))
{
  Y1_EET_CRISPFEP_perm = permute(results_Y1_EET_CRISPFEP_seq$obs)
  Y1_EET_auc_CRISPFEP_null = c(Y1_EET_auc_CRISPFEP_null, roc(predictor = results_Y1_EET_CRISPFEP_seq$pred, response = Y1_EET_CRISPFEP_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_EET_auc_CRISPFEP_null >= auc_Y1_EET_CRISPFEP_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_EET_CRISPFEP_seq$pred, results_Y1_EET_CRISPFEP_seq$obs=="Yes")
#95CI = 1.96*SE
1.96*sqrt(roc.area.test(results_Y1_EET_CRISPFEP_seq$pred, results_Y1_EET_CRISPFEP_seq$obs=="Yes")$var)


#Build GLM model with CRISPFEP columns on whole EDEN as training set for external validation
prep_Y1_EET_CRISPFEP_all = prep_Y1_EET_CRISPFEP_Stand[ ,!(colnames(prep_Y1_EET_CRISPFEP_Stand) %in% c("Site"))]
controlAll <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary)
set.seed(987)
mod_Y1_EET_EDEN_CF_all = train(M12_EET ~ ., data=prep_Y1_EET_CRISPFEP_all, method="glm", metric="ROC", preProc = c("knnImpute"), trControl=controlAll, na.action = na.pass)
#EETove training data for sharing model
mod_Y1_EET_EDEN_CF_all$trainingData = NULL
save(mod_Y1_EET_EDEN_CF_all, file = "mod_Y1_EET.rda")

#Internal-External Validation for OPUS only columns to check still works
#tried glmnet, glm, rf, linear and radial svm - all very similar. Just use GLM
#Y1 EET
#
results_Y1_EET_OPUS = list()
mods_Y1_EET_OPUS = list()

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  
  #set up leave one site out cv
  test_set = prep_Y1_EET_OPUS_Stand[ which(prep_Y1_EET_OPUS_Stand$Site == levels(prep_Y1_EET_OPUS_Stand$Site)[i]), ]
  train_set = prep_Y1_EET_OPUS_Stand[ -which(prep_Y1_EET_OPUS_Stand$Site == levels(prep_Y1_EET_OPUS_Stand$Site)[i]), ]
  
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  train_set = train_set[ ,!(colnames(train_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_EET) 
  
  #train model over grid of 10 by 10 lambda alpha, knn impute, standardise
  mod <- train(M12_EET ~ ., data=train_set, method="glm", metric="ROC", tuneLength = 10, preProc = c("knnImpute"), trControl=control,na.action = na.pass)
  result$pred = predict(mod, test_set, type = "prob", na.action = na.pass)
  
  results_Y1_EET_OPUS[[i]] = result
  mods_Y1_EET_OPUS[[i]] = mod
}

results_Y1_EET_OPUS_seq = NULL
results_Y1_EET_OPUS_seq$pred = c(results_Y1_EET_OPUS[[1]]$pred$Y,results_Y1_EET_OPUS[[2]]$pred$Y, results_Y1_EET_OPUS[[3]]$pred$Y,
                                 results_Y1_EET_OPUS[[4]]$pred$Y,results_Y1_EET_OPUS[[5]]$pred$Y, results_Y1_EET_OPUS[[6]]$pred$Y,
                                 results_Y1_EET_OPUS[[7]]$pred$Y,results_Y1_EET_OPUS[[8]]$pred$Y,results_Y1_EET_OPUS[[9]]$pred$Y,
                                 results_Y1_EET_OPUS[[10]]$pred$Y,results_Y1_EET_OPUS[[11]]$pred$Y,results_Y1_EET_OPUS[[12]]$pred$Y,
                                 results_Y1_EET_OPUS[[13]]$pred$Y,results_Y1_EET_OPUS[[14]]$pred$Y)
results_Y1_EET_OPUS_seq$obs = factor(c(as.character(results_Y1_EET_OPUS[[1]]$obs),as.character(results_Y1_EET_OPUS[[2]]$obs),
                                       as.character(results_Y1_EET_OPUS[[3]]$obs),as.character(results_Y1_EET_OPUS[[4]]$obs),
                                       as.character(results_Y1_EET_OPUS[[5]]$obs),as.character(results_Y1_EET_OPUS[[6]]$obs),
                                       as.character(results_Y1_EET_OPUS[[7]]$obs),as.character(results_Y1_EET_OPUS[[8]]$obs),
                                       as.character(results_Y1_EET_OPUS[[9]]$obs),as.character(results_Y1_EET_OPUS[[10]]$obs),
                                       as.character(results_Y1_EET_OPUS[[11]]$obs),as.character(results_Y1_EET_OPUS[[12]]$obs),
                                       as.character(results_Y1_EET_OPUS[[13]]$obs),as.character(results_Y1_EET_OPUS[[14]]$obs)))

auc_Y1_EET_OPUS_seq = roc(predictor = results_Y1_EET_OPUS_seq$pred, response = results_Y1_EET_OPUS_seq$obs)$auc
roc_Y1_EET_OPUS_seq = roc(predictor = results_Y1_EET_OPUS_seq$pred, response = results_Y1_EET_OPUS_seq$obs)

plot(roc_Y1_EET_OPUS_seq)
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
coords(roc_Y1_EET_OPUS_seq,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

set.seed(987)
Y1_EET_auc_OPUS_null = NULL
for(i in seq (1:10001))
{
  Y1_EET_OPUS_perm = permute(results_Y1_EET_OPUS_seq$obs)
  Y1_EET_auc_OPUS_null = c(Y1_EET_auc_OPUS_null, roc(predictor = results_Y1_EET_OPUS_seq$pred, response = Y1_EET_OPUS_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_EET_auc_OPUS_null >= auc_Y1_EET_OPUS_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_EET_OPUS_seq$pred, results_Y1_EET_OPUS_seq$obs=="Yes")
#95CI = 1.96*SE
1.96*sqrt(roc.area.test(results_Y1_EET_OPUS_seq$pred, results_Y1_EET_OPUS_seq$obs=="Yes")$var)

#Build GLM model with OPUS columns on whole EDEN as training set for external validation
prep_Y1_EET_OPUS_all = prep_Y1_EET_OPUS_Stand[ ,!(colnames(prep_Y1_EET_OPUS_Stand) %in% c("Site"))]
controlAll <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary)
set.seed(987)
mod_Y1_EET_EDEN_all = train(M12_EET ~ ., data=prep_Y1_EET_OPUS_all, method="glm", metric="ROC", preProc = c("knnImpute"), trControl=controlAll, na.action = na.pass)
#remove training data for sharing model
mod_Y1_EET_EDEN_all$trainingData = NULL
save(mod_Y1_EET_EDEN_all, file = "mod_Y1_EET_OPUS.rda")

###############################################################
#
#Y1 EQ5D_TTO Outcome and Models
#
results_Y1_EQ5D_TTO = list()
mods_Y1_EQ5D_TTO = list()

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  
  #set up leave one site out cv
  test_set = prep_Y1_EQ5D_TTO[ which(prep_Y1_EQ5D_TTO$Site == levels(prep_Y1_EQ5D_TTO$Site)[i]), ]
  train_set = prep_Y1_EQ5D_TTO[ -which(prep_Y1_EQ5D_TTO$Site == levels(prep_Y1_EQ5D_TTO$Site)[i]), ]
  
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  train_set = train_set[ ,!(colnames(train_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_EQ5D_UK_TTO_Index_Binary) 
  
  #train model over grid of 10 by 10 lambda alpha, knn impute, standardise
  mod <- train(M12_EQ5D_UK_TTO_Index_Binary ~ ., data=train_set, method="glmnet", metric="ROC", tuneLength = 10, preProc = c("center", "scale","knnImpute"), trControl=control,na.action = na.pass)
  result$pred = predict(mod, test_set, type = "prob", na.action = na.pass)
  
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_EQ5D_TTO[[i]] = result
  mods_Y1_EQ5D_TTO[[i]] = mod
}

results_Y1_EQ5D_TTO_seq = NULL
results_Y1_EQ5D_TTO_seq$pred = c(results_Y1_EQ5D_TTO[[1]]$pred$Y,results_Y1_EQ5D_TTO[[2]]$pred$Y, results_Y1_EQ5D_TTO[[3]]$pred$Y,
                            results_Y1_EQ5D_TTO[[4]]$pred$Y,results_Y1_EQ5D_TTO[[5]]$pred$Y, results_Y1_EQ5D_TTO[[6]]$pred$Y,
                            results_Y1_EQ5D_TTO[[7]]$pred$Y,results_Y1_EQ5D_TTO[[8]]$pred$Y,results_Y1_EQ5D_TTO[[9]]$pred$Y,
                            results_Y1_EQ5D_TTO[[10]]$pred$Y,results_Y1_EQ5D_TTO[[11]]$pred$Y,results_Y1_EQ5D_TTO[[12]]$pred$Y,
                            results_Y1_EQ5D_TTO[[13]]$pred$Y,results_Y1_EQ5D_TTO[[14]]$pred$Y)
results_Y1_EQ5D_TTO_seq$obs = factor(c(as.character(results_Y1_EQ5D_TTO[[1]]$obs),as.character(results_Y1_EQ5D_TTO[[2]]$obs),
                                  as.character(results_Y1_EQ5D_TTO[[3]]$obs),as.character(results_Y1_EQ5D_TTO[[4]]$obs),
                                  as.character(results_Y1_EQ5D_TTO[[5]]$obs),as.character(results_Y1_EQ5D_TTO[[6]]$obs),
                                  as.character(results_Y1_EQ5D_TTO[[7]]$obs),as.character(results_Y1_EQ5D_TTO[[8]]$obs),
                                  as.character(results_Y1_EQ5D_TTO[[9]]$obs),as.character(results_Y1_EQ5D_TTO[[10]]$obs),
                                  as.character(results_Y1_EQ5D_TTO[[11]]$obs),as.character(results_Y1_EQ5D_TTO[[12]]$obs),
                                  as.character(results_Y1_EQ5D_TTO[[13]]$obs),as.character(results_Y1_EQ5D_TTO[[14]]$obs)))

results_Y1_EQ5D_TTO_seq$row = c(results_Y1_EQ5D_TTO[[1]]$row,results_Y1_EQ5D_TTO[[2]]$row, results_Y1_EQ5D_TTO[[3]]$row,
                           results_Y1_EQ5D_TTO[[4]]$row,results_Y1_EQ5D_TTO[[5]]$row, results_Y1_EQ5D_TTO[[6]]$row,
                           results_Y1_EQ5D_TTO[[7]]$row,results_Y1_EQ5D_TTO[[8]]$row,results_Y1_EQ5D_TTO[[9]]$row,
                           results_Y1_EQ5D_TTO[[10]]$row,results_Y1_EQ5D_TTO[[11]]$row,results_Y1_EQ5D_TTO[[12]]$row,
                           results_Y1_EQ5D_TTO[[13]]$row,results_Y1_EQ5D_TTO[[14]]$row)

auc_Y1_EQ5D_TTO_seq = roc(predictor = results_Y1_EQ5D_TTO_seq$pred, response = results_Y1_EQ5D_TTO_seq$obs)$auc
roc_Y1_EQ5D_TTO_seq = roc(predictor = results_Y1_EQ5D_TTO_seq$pred, response = results_Y1_EQ5D_TTO_seq$obs)

plot(roc_Y1_EQ5D_TTO_seq)
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
set.seed(987)
ci.coords(roc_Y1_EQ5D_TTO_seq,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

png("figure1D.png",res = 300, width = 15, height = 15, units = 'cm')
plot.roc(roc_Y1_EQ5D_TTO_seq, auc.polygon=TRUE, grid=c(0.1, 0.1))
text(0.3,0.3,"AUC = 0·704\n(0·667, 0·742)\np<0·0001")
dev.off()

#permutation test for the overall auc vs null distribution
set.seed(987)
Y1_EQ5D_TTO_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_EQ5D_TTO_perm = permute(results_Y1_EQ5D_TTO_seq$obs)
  Y1_EQ5D_TTO_auc_null = c(Y1_EQ5D_TTO_auc_null, roc(predictor = results_Y1_EQ5D_TTO_seq$pred, response = Y1_EQ5D_TTO_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_EQ5D_TTO_auc_null >= auc_Y1_EQ5D_TTO_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_EQ5D_TTO_seq$pred, results_Y1_EQ5D_TTO_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_EQ5D_TTO_seq$pred, results_Y1_EQ5D_TTO_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_EQ5D_TTO_seq$pred, results_Y1_EQ5D_TTO_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_EQ5D_TTO_seq$pred, results_Y1_EQ5D_TTO_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_EQ5D_TTO_seq$pred, results_Y1_EQ5D_TTO_seq$obs=="Yes")$var))

#Test model on different outcomes
#Y1 Rem
#
results_Y1_EQ5D_TTO_Rem = list()
#Add columns back in that were removed and required for this model
prep_Y1_Rem_EQ5D_TTO = prep_Y1_Rem
setdiff(colnames(mods_Y1_EQ5D_TTO[[1]]$trainingData),colnames(prep_Y1_Rem_EQ5D_TTO))
prep_Y1_Rem_EQ5D_TTO$PW1_First_Contact.Police = NA
prep_Y1_Rem_EQ5D_TTO$PW1_First_Contact.Police = as.numeric(prep_Y1_Rem_EQ5D_TTO$PW1_First_Contact.Police)

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  #set up leave one site out cv
  test_set = prep_Y1_Rem_EQ5D_TTO[ which(prep_Y1_Rem_EQ5D_TTO$Site == levels(prep_Y1_Rem_EQ5D_TTO$Site)[i]), ]
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_PANSS_Period_Rem) 
  result$pred = predict(mods_Y1_EQ5D_TTO[[i]], test_set, type = "prob", na.action = na.pass)
  
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_EQ5D_TTO_Rem[[i]] = result
}

results_Y1_EQ5D_TTO_Rem_seq = NULL
results_Y1_EQ5D_TTO_Rem_seq$pred = c(results_Y1_EQ5D_TTO_Rem[[1]]$pred$Y,results_Y1_EQ5D_TTO_Rem[[2]]$pred$Y, results_Y1_EQ5D_TTO_Rem[[3]]$pred$Y,
                                results_Y1_EQ5D_TTO_Rem[[4]]$pred$Y,results_Y1_EQ5D_TTO_Rem[[5]]$pred$Y, results_Y1_EQ5D_TTO_Rem[[6]]$pred$Y,
                                results_Y1_EQ5D_TTO_Rem[[7]]$pred$Y,results_Y1_EQ5D_TTO_Rem[[8]]$pred$Y,results_Y1_EQ5D_TTO_Rem[[9]]$pred$Y,
                                results_Y1_EQ5D_TTO_Rem[[10]]$pred$Y,results_Y1_EQ5D_TTO_Rem[[11]]$pred$Y,results_Y1_EQ5D_TTO_Rem[[12]]$pred$Y,
                                results_Y1_EQ5D_TTO_Rem[[13]]$pred$Y,results_Y1_EQ5D_TTO_Rem[[14]]$pred$Y)
results_Y1_EQ5D_TTO_Rem_seq$obs = factor(c(as.character(results_Y1_EQ5D_TTO_Rem[[1]]$obs),as.character(results_Y1_EQ5D_TTO_Rem[[2]]$obs),
                                      as.character(results_Y1_EQ5D_TTO_Rem[[3]]$obs),as.character(results_Y1_EQ5D_TTO_Rem[[4]]$obs),
                                      as.character(results_Y1_EQ5D_TTO_Rem[[5]]$obs),as.character(results_Y1_EQ5D_TTO_Rem[[6]]$obs),
                                      as.character(results_Y1_EQ5D_TTO_Rem[[7]]$obs),as.character(results_Y1_EQ5D_TTO_Rem[[8]]$obs),
                                      as.character(results_Y1_EQ5D_TTO_Rem[[9]]$obs),as.character(results_Y1_EQ5D_TTO_Rem[[10]]$obs),
                                      as.character(results_Y1_EQ5D_TTO_Rem[[11]]$obs),as.character(results_Y1_EQ5D_TTO_Rem[[12]]$obs),
                                      as.character(results_Y1_EQ5D_TTO_Rem[[13]]$obs),as.character(results_Y1_EQ5D_TTO_Rem[[14]]$obs)))
results_Y1_EQ5D_TTO_Rem_seq$row = c(results_Y1_EQ5D_TTO_Rem[[1]]$row,results_Y1_EQ5D_TTO_Rem[[2]]$row, results_Y1_EQ5D_TTO_Rem[[3]]$row,
                               results_Y1_EQ5D_TTO_Rem[[4]]$row,results_Y1_EQ5D_TTO_Rem[[5]]$row, results_Y1_EQ5D_TTO_Rem[[6]]$row,
                               results_Y1_EQ5D_TTO_Rem[[7]]$row,results_Y1_EQ5D_TTO_Rem[[8]]$row,results_Y1_EQ5D_TTO_Rem[[9]]$row,
                               results_Y1_EQ5D_TTO_Rem[[10]]$row,results_Y1_EQ5D_TTO_Rem[[11]]$row,results_Y1_EQ5D_TTO_Rem[[12]]$row,
                               results_Y1_EQ5D_TTO_Rem[[13]]$row,results_Y1_EQ5D_TTO_Rem[[14]]$row)

auc_Y1_EQ5D_TTO_Rem_seq = roc(predictor = results_Y1_EQ5D_TTO_Rem_seq$pred, response = results_Y1_EQ5D_TTO_Rem_seq$obs)$auc
roc_Y1_EQ5D_TTO_Rem_seq = roc(predictor = results_Y1_EQ5D_TTO_Rem_seq$pred, response = results_Y1_EQ5D_TTO_Rem_seq$obs)

set.seed(987)
Y1_EQ5D_TTO_Rem_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_EQ5D_TTO_Rem_perm = permute(results_Y1_EQ5D_TTO_Rem_seq$obs)
  Y1_EQ5D_TTO_Rem_auc_null = c(Y1_EQ5D_TTO_Rem_auc_null, roc(predictor = results_Y1_EQ5D_TTO_Rem_seq$pred, response = Y1_EQ5D_TTO_Rem_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_EQ5D_TTO_Rem_auc_null >= auc_Y1_EQ5D_TTO_Rem_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_EQ5D_TTO_Rem_seq$pred, results_Y1_EQ5D_TTO_Rem_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_EQ5D_TTO_Rem_seq$pred, results_Y1_EQ5D_TTO_Rem_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_EQ5D_TTO_Rem_seq$pred, results_Y1_EQ5D_TTO_Rem_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_EQ5D_TTO_Rem_seq$pred, results_Y1_EQ5D_TTO_Rem_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_EQ5D_TTO_Rem_seq$pred, results_Y1_EQ5D_TTO_Rem_seq$obs=="Yes")$var))

#Y1 GAF
#
results_Y1_EQ5D_TTO_GAF = list()
#Add columns back in that were removed and required for this model
prep_Y1_GAF_EQ5D_TTO = prep_Y1_GAF
setdiff(colnames(mods_Y1_EQ5D_TTO[[1]]$trainingData),colnames(prep_Y1_GAF_EQ5D_TTO))
prep_Y1_GAF_EQ5D_TTO$BL_MostSeriousHarm_Premeditation.Self.harm.contemplated.for.more.than.three.hours.prior.to.attempt = NA
prep_Y1_GAF_EQ5D_TTO$BL_MostSeriousHarm_Premeditation.Self.harm.contemplated.for.more.than.three.hours.prior.to.attempt = as.numeric(prep_Y1_GAF_EQ5D_TTO$BL_MostSeriousHarm_Premeditation.Self.harm.contemplated.for.more.than.three.hours.prior.to.attempt)

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  #set up leave one site out cv
  test_set = prep_Y1_GAF_EQ5D_TTO[ which(prep_Y1_GAF_EQ5D_TTO$Site == levels(prep_Y1_GAF_EQ5D_TTO$Site)[i]), ]
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_GAF) 
  result$pred = predict(mods_Y1_EQ5D_TTO[[i]], test_set, type = "prob", na.action = na.pass)
  
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_EQ5D_TTO_GAF[[i]] = result
}

results_Y1_EQ5D_TTO_GAF_seq = NULL
results_Y1_EQ5D_TTO_GAF_seq$pred = c(results_Y1_EQ5D_TTO_GAF[[1]]$pred$Y,results_Y1_EQ5D_TTO_GAF[[2]]$pred$Y, results_Y1_EQ5D_TTO_GAF[[3]]$pred$Y,
                                     results_Y1_EQ5D_TTO_GAF[[4]]$pred$Y,results_Y1_EQ5D_TTO_GAF[[5]]$pred$Y, results_Y1_EQ5D_TTO_GAF[[6]]$pred$Y,
                                     results_Y1_EQ5D_TTO_GAF[[7]]$pred$Y,results_Y1_EQ5D_TTO_GAF[[8]]$pred$Y,results_Y1_EQ5D_TTO_GAF[[9]]$pred$Y,
                                     results_Y1_EQ5D_TTO_GAF[[10]]$pred$Y,results_Y1_EQ5D_TTO_GAF[[11]]$pred$Y,results_Y1_EQ5D_TTO_GAF[[12]]$pred$Y,
                                     results_Y1_EQ5D_TTO_GAF[[13]]$pred$Y,results_Y1_EQ5D_TTO_GAF[[14]]$pred$Y)
results_Y1_EQ5D_TTO_GAF_seq$obs = factor(c(as.character(results_Y1_EQ5D_TTO_GAF[[1]]$obs),as.character(results_Y1_EQ5D_TTO_GAF[[2]]$obs),
                                           as.character(results_Y1_EQ5D_TTO_GAF[[3]]$obs),as.character(results_Y1_EQ5D_TTO_GAF[[4]]$obs),
                                           as.character(results_Y1_EQ5D_TTO_GAF[[5]]$obs),as.character(results_Y1_EQ5D_TTO_GAF[[6]]$obs),
                                           as.character(results_Y1_EQ5D_TTO_GAF[[7]]$obs),as.character(results_Y1_EQ5D_TTO_GAF[[8]]$obs),
                                           as.character(results_Y1_EQ5D_TTO_GAF[[9]]$obs),as.character(results_Y1_EQ5D_TTO_GAF[[10]]$obs),
                                           as.character(results_Y1_EQ5D_TTO_GAF[[11]]$obs),as.character(results_Y1_EQ5D_TTO_GAF[[12]]$obs),
                                           as.character(results_Y1_EQ5D_TTO_GAF[[13]]$obs),as.character(results_Y1_EQ5D_TTO_GAF[[14]]$obs)))
results_Y1_EQ5D_TTO_GAF_seq$row = c(results_Y1_EQ5D_TTO_GAF[[1]]$row,results_Y1_EQ5D_TTO_GAF[[2]]$row, results_Y1_EQ5D_TTO_GAF[[3]]$row,
                                    results_Y1_EQ5D_TTO_GAF[[4]]$row,results_Y1_EQ5D_TTO_GAF[[5]]$row, results_Y1_EQ5D_TTO_GAF[[6]]$row,
                                    results_Y1_EQ5D_TTO_GAF[[7]]$row,results_Y1_EQ5D_TTO_GAF[[8]]$row,results_Y1_EQ5D_TTO_GAF[[9]]$row,
                                    results_Y1_EQ5D_TTO_GAF[[10]]$row,results_Y1_EQ5D_TTO_GAF[[11]]$row,results_Y1_EQ5D_TTO_GAF[[12]]$row,
                                    results_Y1_EQ5D_TTO_GAF[[13]]$row,results_Y1_EQ5D_TTO_GAF[[14]]$row)

auc_Y1_EQ5D_TTO_GAF_seq = roc(predictor = results_Y1_EQ5D_TTO_GAF_seq$pred, response = results_Y1_EQ5D_TTO_GAF_seq$obs)$auc
roc_Y1_EQ5D_TTO_GAF_seq = roc(predictor = results_Y1_EQ5D_TTO_GAF_seq$pred, response = results_Y1_EQ5D_TTO_GAF_seq$obs)

set.seed(987)
Y1_EQ5D_TTO_GAF_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_EQ5D_TTO_GAF_perm = permute(results_Y1_EQ5D_TTO_GAF_seq$obs)
  Y1_EQ5D_TTO_GAF_auc_null = c(Y1_EQ5D_TTO_GAF_auc_null, roc(predictor = results_Y1_EQ5D_TTO_GAF_seq$pred, response = Y1_EQ5D_TTO_GAF_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_EQ5D_TTO_GAF_auc_null >= auc_Y1_EQ5D_TTO_GAF_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_EQ5D_TTO_GAF_seq$pred, results_Y1_EQ5D_TTO_GAF_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_EQ5D_TTO_GAF_seq$pred, results_Y1_EQ5D_TTO_GAF_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_EQ5D_TTO_GAF_seq$pred, results_Y1_EQ5D_TTO_GAF_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_EQ5D_TTO_GAF_seq$pred, results_Y1_EQ5D_TTO_GAF_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_EQ5D_TTO_GAF_seq$pred, results_Y1_EQ5D_TTO_GAF_seq$obs=="Yes")$var))

#Y1 EET
#
results_Y1_EQ5D_TTO_EET = list()
#Add columns back in that were removed and required for this model
prep_Y1_EET_EQ5D_TTO = prep_Y1_EET
setdiff(colnames(mods_Y1_EQ5D_TTO[[1]]$trainingData),colnames(prep_Y1_EET_EQ5D_TTO))
prep_Y1_EET_EQ5D_TTO$BL_MostSeriousHarm_Premeditation.Self.harm.contemplated.for.more.than.three.hours.prior.to.attempt = NA
prep_Y1_EET_EQ5D_TTO$BL_MostSeriousHarm_Premeditation.Self.harm.contemplated.for.more.than.three.hours.prior.to.attempt = as.numeric(prep_Y1_EET_EQ5D_TTO$BL_MostSeriousHarm_Premeditation.Self.harm.contemplated.for.more.than.three.hours.prior.to.attempt)

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  #set up leave one site out cv
  test_set = prep_Y1_EET_EQ5D_TTO[ which(prep_Y1_EET_EQ5D_TTO$Site == levels(prep_Y1_EET_EQ5D_TTO$Site)[i]), ]
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_EET) 
  result$pred = predict(mods_Y1_EQ5D_TTO[[i]], test_set, type = "prob", na.action = na.pass)
  
  result$row = as.numeric(rownames(test_set))
  
  results_Y1_EQ5D_TTO_EET[[i]] = result
}

results_Y1_EQ5D_TTO_EET_seq = NULL
results_Y1_EQ5D_TTO_EET_seq$pred = c(results_Y1_EQ5D_TTO_EET[[1]]$pred$Y,results_Y1_EQ5D_TTO_EET[[2]]$pred$Y, results_Y1_EQ5D_TTO_EET[[3]]$pred$Y,
                                     results_Y1_EQ5D_TTO_EET[[4]]$pred$Y,results_Y1_EQ5D_TTO_EET[[5]]$pred$Y, results_Y1_EQ5D_TTO_EET[[6]]$pred$Y,
                                     results_Y1_EQ5D_TTO_EET[[7]]$pred$Y,results_Y1_EQ5D_TTO_EET[[8]]$pred$Y,results_Y1_EQ5D_TTO_EET[[9]]$pred$Y,
                                     results_Y1_EQ5D_TTO_EET[[10]]$pred$Y,results_Y1_EQ5D_TTO_EET[[11]]$pred$Y,results_Y1_EQ5D_TTO_EET[[12]]$pred$Y,
                                     results_Y1_EQ5D_TTO_EET[[13]]$pred$Y,results_Y1_EQ5D_TTO_EET[[14]]$pred$Y)
results_Y1_EQ5D_TTO_EET_seq$obs = factor(c(as.character(results_Y1_EQ5D_TTO_EET[[1]]$obs),as.character(results_Y1_EQ5D_TTO_EET[[2]]$obs),
                                           as.character(results_Y1_EQ5D_TTO_EET[[3]]$obs),as.character(results_Y1_EQ5D_TTO_EET[[4]]$obs),
                                           as.character(results_Y1_EQ5D_TTO_EET[[5]]$obs),as.character(results_Y1_EQ5D_TTO_EET[[6]]$obs),
                                           as.character(results_Y1_EQ5D_TTO_EET[[7]]$obs),as.character(results_Y1_EQ5D_TTO_EET[[8]]$obs),
                                           as.character(results_Y1_EQ5D_TTO_EET[[9]]$obs),as.character(results_Y1_EQ5D_TTO_EET[[10]]$obs),
                                           as.character(results_Y1_EQ5D_TTO_EET[[11]]$obs),as.character(results_Y1_EQ5D_TTO_EET[[12]]$obs),
                                           as.character(results_Y1_EQ5D_TTO_EET[[13]]$obs),as.character(results_Y1_EQ5D_TTO_EET[[14]]$obs)))
results_Y1_EQ5D_TTO_EET_seq$row = c(results_Y1_EQ5D_TTO_EET[[1]]$row,results_Y1_EQ5D_TTO_EET[[2]]$row, results_Y1_EQ5D_TTO_EET[[3]]$row,
                                    results_Y1_EQ5D_TTO_EET[[4]]$row,results_Y1_EQ5D_TTO_EET[[5]]$row, results_Y1_EQ5D_TTO_EET[[6]]$row,
                                    results_Y1_EQ5D_TTO_EET[[7]]$row,results_Y1_EQ5D_TTO_EET[[8]]$row,results_Y1_EQ5D_TTO_EET[[9]]$row,
                                    results_Y1_EQ5D_TTO_EET[[10]]$row,results_Y1_EQ5D_TTO_EET[[11]]$row,results_Y1_EQ5D_TTO_EET[[12]]$row,
                                    results_Y1_EQ5D_TTO_EET[[13]]$row,results_Y1_EQ5D_TTO_EET[[14]]$row)

auc_Y1_EQ5D_TTO_EET_seq = roc(predictor = results_Y1_EQ5D_TTO_EET_seq$pred, response = results_Y1_EQ5D_TTO_EET_seq$obs)$auc
roc_Y1_EQ5D_TTO_EET_seq = roc(predictor = results_Y1_EQ5D_TTO_EET_seq$pred, response = results_Y1_EQ5D_TTO_EET_seq$obs)

set.seed(987)
Y1_EQ5D_TTO_EET_auc_null = NULL
for(i in seq (1:10001))
{
  Y1_EQ5D_TTO_EET_perm = permute(results_Y1_EQ5D_TTO_EET_seq$obs)
  Y1_EQ5D_TTO_EET_auc_null = c(Y1_EQ5D_TTO_EET_auc_null, roc(predictor = results_Y1_EQ5D_TTO_EET_seq$pred, response = Y1_EQ5D_TTO_EET_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_EQ5D_TTO_EET_auc_null >= auc_Y1_EQ5D_TTO_EET_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_EQ5D_TTO_EET_seq$pred, results_Y1_EQ5D_TTO_EET_seq$obs=="Yes")
#95CI = 1.96*SE
(roc.area.test(results_Y1_EQ5D_TTO_EET_seq$pred, results_Y1_EQ5D_TTO_EET_seq$obs=="Yes")$area) - (1.96*sqrt(roc.area.test(results_Y1_EQ5D_TTO_EET_seq$pred, results_Y1_EQ5D_TTO_EET_seq$obs=="Yes")$var))
(roc.area.test(results_Y1_EQ5D_TTO_EET_seq$pred, results_Y1_EQ5D_TTO_EET_seq$obs=="Yes")$area) + (1.96*sqrt(roc.area.test(results_Y1_EQ5D_TTO_EET_seq$pred, results_Y1_EQ5D_TTO_EET_seq$obs=="Yes")$var))

coefs_Y1_EQ5D_TTO = NULL
for (i in seq(1:14))
{
  coefs_Y1_EQ5D_TTO = c(coefs_Y1_EQ5D_TTO, coef(mods_Y1_EQ5D_TTO[[i]]$finalModel, mods_Y1_EQ5D_TTO[[i]]$bestTune$lambda))
}

#just get numbers
coefs_Y1_EQ5D_TTO_extract = NULL
for(i in seq(1:14))
{
  coefs_Y1_EQ5D_TTO_extract = rbind(coefs_Y1_EQ5D_TTO_extract, coefs_Y1_EQ5D_TTO[[i]][1:165])
}

#Presence or absence of predictors across models
coefs_Y1_EQ5D_TTO_presence = NULL
coefs_Y1_EQ5D_TTO_presence = coefs_Y1_EQ5D_TTO_extract[1:14,2:165]
coefs_Y1_EQ5D_TTO_presence[coefs_Y1_EQ5D_TTO_presence != 0] <- 1

#stability of feature selection http://jmlr.org/papers/volume18/17-514/17-514.pdf
getStability(coefs_Y1_EQ5D_TTO_presence)

#work out number of predictors shared across all models
table(colMeans(coefs_Y1_EQ5D_TTO_presence))
#number of features in each model
rowSums(coefs_Y1_EQ5D_TTO_presence)

#get rank of coef by importance as in sports ranking
coefs_Y1_EQ5D_TTO_baseline_rank = NULL

for(i in seq(c(1:14)))
{
  #rank absolute value excluding the intercept for each model
  coefs_Y1_EQ5D_TTO_baseline_rank = rbind(coefs_Y1_EQ5D_TTO_baseline_rank, rank(abs(coefs_Y1_EQ5D_TTO_extract[i,2:165]), ties.method = "min"))
}

# rank the mean ranks of each column across all models
coefs_Y1_EQ5D_TTO_baseline_rank_mean = colMeans(coefs_Y1_EQ5D_TTO_baseline_rank)

#Invert order of rank to identify top models
coefs_Y1_EQ5D_TTO_baseline_order = rank(-coefs_Y1_EQ5D_TTO_baseline_rank_mean)

#Get the column names (not the intercept)
coef_Y1_EQ5D_TTO_names = dimnames(coefs_Y1_EQ5D_TTO[[1]])[[1]][2:165]

coefs_Y1_EQ5D_TTO_means = colMeans(coefs_Y1_EQ5D_TTO_extract)[2:165]

df_Y1_EQ5D_TTO_names = data.frame(coef_Y1_EQ5D_TTO_names, coefs_Y1_EQ5D_TTO_baseline_order, coefs_Y1_EQ5D_TTO_means, colMeans(coefs_Y1_EQ5D_TTO_presence))

#View the best predictor variables by absolute value in order (lower is better)
View(df_Y1_EQ5D_TTO_names)

#Internal-External Validation for CRISPFEP only columns to check still works
#Y1 QoL
#
results_Y1_QoL_CRISPFEP = list()
mods_Y1_QoL_CRISPFEP = list()

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  
  #set up leave one site out cv
  test_set = prep_Y1_QoL_CRISPFEP_Stand[ which(prep_Y1_QoL_CRISPFEP_Stand$Site == levels(prep_Y1_QoL_CRISPFEP_Stand$Site)[i]), ]
  train_set = prep_Y1_QoL_CRISPFEP_Stand[ -which(prep_Y1_QoL_CRISPFEP_Stand$Site == levels(prep_Y1_QoL_CRISPFEP_Stand$Site)[i]), ]
  
  #Remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  train_set = train_set[ ,!(colnames(train_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_QoL) 
  
  #train model over grid of 10 by 10 lambda alpha, knn impute, standardise
  mod <- train(M12_QoL ~ ., data=train_set, method="glm", metric="ROC", tuneLength = 10, preProc = c("knnImpute"), trControl=control,na.action = na.pass)
  result$pred = predict(mod, test_set, type = "prob", na.action = na.pass)
  
  results_Y1_QoL_CRISPFEP[[i]] = result
  mods_Y1_QoL_CRISPFEP[[i]] = mod
}

results_Y1_QoL_CRISPFEP_seq = NULL
results_Y1_QoL_CRISPFEP_seq$pred = c(results_Y1_QoL_CRISPFEP[[1]]$pred$Y,results_Y1_QoL_CRISPFEP[[2]]$pred$Y, results_Y1_QoL_CRISPFEP[[3]]$pred$Y,
                                     results_Y1_QoL_CRISPFEP[[4]]$pred$Y,results_Y1_QoL_CRISPFEP[[5]]$pred$Y, results_Y1_QoL_CRISPFEP[[6]]$pred$Y,
                                     results_Y1_QoL_CRISPFEP[[7]]$pred$Y,results_Y1_QoL_CRISPFEP[[8]]$pred$Y,results_Y1_QoL_CRISPFEP[[9]]$pred$Y,
                                     results_Y1_QoL_CRISPFEP[[10]]$pred$Y,results_Y1_QoL_CRISPFEP[[11]]$pred$Y,results_Y1_QoL_CRISPFEP[[12]]$pred$Y,
                                     results_Y1_QoL_CRISPFEP[[13]]$pred$Y,results_Y1_QoL_CRISPFEP[[14]]$pred$Y)
results_Y1_QoL_CRISPFEP_seq$obs = factor(c(as.character(results_Y1_QoL_CRISPFEP[[1]]$obs),as.character(results_Y1_QoL_CRISPFEP[[2]]$obs),
                                           as.character(results_Y1_QoL_CRISPFEP[[3]]$obs),as.character(results_Y1_QoL_CRISPFEP[[4]]$obs),
                                           as.character(results_Y1_QoL_CRISPFEP[[5]]$obs),as.character(results_Y1_QoL_CRISPFEP[[6]]$obs),
                                           as.character(results_Y1_QoL_CRISPFEP[[7]]$obs),as.character(results_Y1_QoL_CRISPFEP[[8]]$obs),
                                           as.character(results_Y1_QoL_CRISPFEP[[9]]$obs),as.character(results_Y1_QoL_CRISPFEP[[10]]$obs),
                                           as.character(results_Y1_QoL_CRISPFEP[[11]]$obs),as.character(results_Y1_QoL_CRISPFEP[[12]]$obs),
                                           as.character(results_Y1_QoL_CRISPFEP[[13]]$obs),as.character(results_Y1_QoL_CRISPFEP[[14]]$obs)))

auc_Y1_QoL_CRISPFEP_seq = roc(predictor = results_Y1_QoL_CRISPFEP_seq$pred, response = results_Y1_QoL_CRISPFEP_seq$obs)$auc
roc_Y1_QoL_CRISPFEP_seq = roc(predictor = results_Y1_QoL_CRISPFEP_seq$pred, response = results_Y1_QoL_CRISPFEP_seq$obs)

plot(roc_Y1_QoL_CRISPFEP_seq)
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
coords(roc_Y1_QoL_CRISPFEP_seq,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

set.seed(987)
Y1_QoL_auc_CRISPFEP_null = NULL
for(i in seq (1:10001))
{
  Y1_QoL_CRISPFEP_perm = permute(results_Y1_QoL_CRISPFEP_seq$obs)
  Y1_QoL_auc_CRISPFEP_null = c(Y1_QoL_auc_CRISPFEP_null, roc(predictor = results_Y1_QoL_CRISPFEP_seq$pred, response = Y1_QoL_CRISPFEP_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_QoL_auc_CRISPFEP_null >= auc_Y1_QoL_CRISPFEP_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_QoL_CRISPFEP_seq$pred, results_Y1_QoL_CRISPFEP_seq$obs=="Yes")
#95CI = 1.96*SE
1.96*sqrt(roc.area.test(results_Y1_QoL_CRISPFEP_seq$pred, results_Y1_QoL_CRISPFEP_seq$obs=="Yes")$var)

#Build GLM model with CRISPFEP columns on whole EDEN as training set for external validation
prep_Y1_QoL_CRISPFEP_all = prep_Y1_QoL_CRISPFEP_Stand[ ,!(colnames(prep_Y1_QoL_CRISPFEP_Stand) %in% c("Site"))]
controlAll <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary)
set.seed(987)
mod_Y1_QoL_EDEN_CF_all = train(M12_QoL ~ ., data=prep_Y1_QoL_CRISPFEP_all, method="glm", metric="ROC", preProc = c("knnImpute"), trControl=controlAll, na.action = na.pass)
#Remove training data for sharing model
mod_Y1_QoL_EDEN_CF_all$trainingData = NULL
save(mod_Y1_QoL_EDEN_CF_all, file = "mod_Y1_QoL.rda")

#Internal-External Validation for OPUS only columns to ensure still works
#tried glmnet, glm, rf, linear and radial svm - all very similar. Just use GLM
#Y1 QoL
#
results_Y1_QoL_OPUS = list()
mods_Y1_QoL_OPUS = list()

#14 sites
for(i in seq(1:14))
{
  set.seed(987)
  
  #set up leave one site out cv
  test_set = prep_Y1_QoL_OPUS_Stand[ which(prep_Y1_QoL_OPUS_Stand$Site == levels(prep_Y1_QoL_OPUS_Stand$Site)[i]), ]
  train_set = prep_Y1_QoL_OPUS_Stand[ -which(prep_Y1_QoL_OPUS_Stand$Site == levels(prep_Y1_QoL_OPUS_Stand$Site)[i]), ]
  
  #remove site column
  test_set = test_set[ ,!(colnames(test_set) %in% c("Site"))]
  train_set = train_set[ ,!(colnames(train_set) %in% c("Site"))]
  
  #Get the observed outcome classes for this test set
  result = data.frame(obs=test_set$M12_QoL) 
  
  #train model over grid of 10 by 10 lambda alpha, knn impute, standardise
  mod <- train(M12_QoL ~ ., data=train_set, method="glm", metric="ROC", tuneLength = 10, preProc = c("knnImpute"), trControl=control,na.action = na.pass)
  result$pred = predict(mod, test_set, type = "prob", na.action = na.pass)
  
  results_Y1_QoL_OPUS[[i]] = result
  mods_Y1_QoL_OPUS[[i]] = mod
}

results_Y1_QoL_OPUS_seq = NULL
results_Y1_QoL_OPUS_seq$pred = c(results_Y1_QoL_OPUS[[1]]$pred$Y,results_Y1_QoL_OPUS[[2]]$pred$Y, results_Y1_QoL_OPUS[[3]]$pred$Y,
                                 results_Y1_QoL_OPUS[[4]]$pred$Y,results_Y1_QoL_OPUS[[5]]$pred$Y, results_Y1_QoL_OPUS[[6]]$pred$Y,
                                 results_Y1_QoL_OPUS[[7]]$pred$Y,results_Y1_QoL_OPUS[[8]]$pred$Y,results_Y1_QoL_OPUS[[9]]$pred$Y,
                                 results_Y1_QoL_OPUS[[10]]$pred$Y,results_Y1_QoL_OPUS[[11]]$pred$Y,results_Y1_QoL_OPUS[[12]]$pred$Y,
                                 results_Y1_QoL_OPUS[[13]]$pred$Y,results_Y1_QoL_OPUS[[14]]$pred$Y)
results_Y1_QoL_OPUS_seq$obs = factor(c(as.character(results_Y1_QoL_OPUS[[1]]$obs),as.character(results_Y1_QoL_OPUS[[2]]$obs),
                                       as.character(results_Y1_QoL_OPUS[[3]]$obs),as.character(results_Y1_QoL_OPUS[[4]]$obs),
                                       as.character(results_Y1_QoL_OPUS[[5]]$obs),as.character(results_Y1_QoL_OPUS[[6]]$obs),
                                       as.character(results_Y1_QoL_OPUS[[7]]$obs),as.character(results_Y1_QoL_OPUS[[8]]$obs),
                                       as.character(results_Y1_QoL_OPUS[[9]]$obs),as.character(results_Y1_QoL_OPUS[[10]]$obs),
                                       as.character(results_Y1_QoL_OPUS[[11]]$obs),as.character(results_Y1_QoL_OPUS[[12]]$obs),
                                       as.character(results_Y1_QoL_OPUS[[13]]$obs),as.character(results_Y1_QoL_OPUS[[14]]$obs)))

auc_Y1_QoL_OPUS_seq = roc(predictor = results_Y1_QoL_OPUS_seq$pred, response = results_Y1_QoL_OPUS_seq$obs)$auc
roc_Y1_QoL_OPUS_seq = roc(predictor = results_Y1_QoL_OPUS_seq$pred, response = results_Y1_QoL_OPUS_seq$obs)

plot(roc_Y1_QoL_OPUS_seq)
#PSI = (PPV+NPV)-1
#LR+ = sens/(1-spec)
#LR- = (1-sens)/spec
coords(roc_Y1_QoL_OPUS_seq,"best", best.method = "youden", ret = c("specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv"))

set.seed(987)
Y1_QoL_auc_OPUS_null = NULL
for(i in seq (1:10001))
{
  Y1_QoL_OPUS_perm = permute(results_Y1_QoL_OPUS_seq$obs)
  Y1_QoL_auc_OPUS_null = c(Y1_QoL_auc_OPUS_null, roc(predictor = results_Y1_QoL_OPUS_seq$pred, response = Y1_QoL_OPUS_perm)$auc)
}

#get p value by taking proportion of permutated values greater or equal to the actual value
(1+sum(Y1_QoL_auc_OPUS_null >= auc_Y1_QoL_OPUS_seq))/10001

#Computes the nonparametric area under the ROC curve and its variance based on U-statistic theory (DDCP).
roc.area.test(results_Y1_QoL_OPUS_seq$pred, results_Y1_QoL_OPUS_seq$obs=="Yes")
#95CI = 1.96*SE
1.96*sqrt(roc.area.test(results_Y1_QoL_OPUS_seq$pred, results_Y1_QoL_OPUS_seq$obs=="Yes")$var)

#Build GLM model with OPUS columns on whole EDEN as training set for external validation
prep_Y1_QoL_OPUS_all = prep_Y1_QoL_OPUS_Stand[ ,!(colnames(prep_Y1_QoL_OPUS_Stand) %in% c("Site"))]
controlAll <- trainControl(method="repeatedcv", number=10, repeats=10, classProbs=TRUE, summaryFunction=twoClassSummary)
set.seed(987)
mod_Y1_QoL_EDEN_all = train(M12_QoL ~ ., data=prep_Y1_QoL_OPUS_all, method="glm", metric="ROC", preProc = c("knnImpute"), trControl=controlAll, na.action = na.pass)
#remove training data for sharing model
mod_Y1_QoL_EDEN_all$trainingData = NULL
save(mod_Y1_QoL_EDEN_all, file = "mod_Y1_QoL_OPUS.rda")

###################################################
#
#Correlation of probability outputs for each model
probabilityRem = data.frame(prob = results_Y1_Rem_seq$pred, row = results_Y1_Rem_seq$row)
for (i in seq (1:1027))
{
  if(!any(probabilityRem$row==i))
    probabilityRem[nrow(probabilityRem) + 1,] = list(NA,i)
}
probabilityRem = data.frame(prob = probabilityRem$prob, row.names = probabilityRem$row)
probabilityRem = probabilityRem[order(as.numeric(rownames(probabilityRem))),,drop=FALSE]

probabilityGAF = data.frame(prob = results_Y1_GAF_seq$pred, row = results_Y1_GAF_seq$row)
for (i in seq (1:1027))
{
  if(!any(probabilityGAF$row==i))
    probabilityGAF[nrow(probabilityGAF) + 1,] = list(NA,i)
}
probabilityGAF = data.frame(prob = probabilityGAF$prob, row.names = probabilityGAF$row)
probabilityGAF = probabilityGAF[order(as.numeric(rownames(probabilityGAF))),,drop=FALSE]

probabilityEET = data.frame(prob = results_Y1_EET_seq$pred, row = results_Y1_EET_seq$row)
for (i in seq (1:1027))
{
  if(!any(probabilityEET$row==i))
    probabilityEET[nrow(probabilityEET) + 1,] = list(NA,i)
}
probabilityEET = data.frame(prob = probabilityEET$prob, row.names = probabilityEET$row)
probabilityEET = probabilityEET[order(as.numeric(rownames(probabilityEET))),,drop=FALSE]

probabilityEQ5D_TTO = data.frame(prob = results_Y1_EQ5D_TTO_seq$pred, row = results_Y1_EQ5D_TTO_seq$row)
for (i in seq (1:1027))
{
  if(!any(probabilityEQ5D_TTO$row==i))
    probabilityEQ5D_TTO[nrow(probabilityEQ5D_TTO) + 1,] = list(NA,i)
}
probabilityEQ5D_TTO = data.frame(prob = probabilityEQ5D_TTO$prob, row.names = probabilityEQ5D_TTO$row)
probabilityEQ5D_TTO = probabilityEQ5D_TTO[order(as.numeric(rownames(probabilityEQ5D_TTO))),,drop=FALSE]

probabilityDf = data.frame(probabilityRem = probabilityRem$prob, probabilityGAF = probabilityGAF$prob,
                           probabilityEET = probabilityEET$prob, probabilityEQ5D_TTO = probabilityEQ5D_TTO$prob)

print(corr.test(probabilityDf), short = F)
