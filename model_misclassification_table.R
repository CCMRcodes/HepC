# This code computes for the misclassification table under one representative split.
# The split is selected to have the closest concordance to average in the boosting longitudinal model.
# Authors of code: Boang Liu, Xuefei Zhang

# Load results of boosting longitudinal model in result_lgt
head(result_lgt)
mean(result_lgt$CIndex)
ind_best = which.min(abs(result_lgt$CIndex-mean(result_lgt$CIndex)))
result_lgt[ind_best,]

# fit again the models on that split
library(gbm)
library(survival)
library(pROC)

fillByMedian <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
}

# Read in data in a data.frame named dataAll.
# Read in ids for 30 splits of training/test.

# longitudinal models
labName <- c('ALB','ALKRATIO','ALTRATIO','ASTRATIO','BIL','BUN','CI','CRE','GLU','HEM','PLT','K','NA', 'TOTP', 'WBC','APRI', 'ASTALT', 'BMI')
varStart <- setdiff(names(dataAll), c(paste(labName, "_MEAN", sep=""), paste(labName, "_DIFFMEAN", sep=""),
                                      paste(labName, "_DDMEAN", sep=""), paste(labName, "_TVARRAW", sep="")))

i <- ind_best
dataTrain <- dataAll[is.element(dataAll$PatientID, ids_noTreated[indexSplit_noTreated[[i]]$indexTrain]), varStart]
dataTest <- dataAll[is.element(dataAll$PatientID, ids_noTreated[indexSplit_noTreated[[i]]$indexTest]), varStart]
dataTrain[, -c(1:9)] <- sapply(dataTrain[, -c(1:9)], fillByMedian)
dataTest[, -c(1:9)] <- sapply(dataTest[, -c(1:9)], fillByMedian)

# boosting longitudinal
set.seed(90)
fit1 <- gbm(Surv(TimeToOutcome, Outcome)~., distribution = "coxph", data = dataTrain[, -1],
            n.trees = 2000, interaction.depth = 2, shrinkage = 0.03)
best.iter <- gbm.perf(fit1,method="OOB")
print(best.iter)
pred1 <- predict(fit1,newdata = dataTest,n.trees = best.iter)
CIndex1 <- survConcordance(Surv(TimeToOutcome, Outcome) ~ pred1, dataTest)
print(CIndex1$concordance)
tmp <- basehaz.gbm(dataTrain$TimeToOutcome, dataTrain$Outcome, predict(fit1,newdata = dataTrain,n.trees = best.iter), 1:5)

for (K in c(1, 3, 5)) {
  print(K)
  indKyear <- !(dataTest$Outcome == 0 & dataTest$TimeToOutcome < K)
  dataTestKyear_outcome <- (dataTest$Outcome[indKyear] == 1 & dataTest$TimeToOutcome[indKyear] < K)
  tmp2 <- 1 - exp(-exp(pred1[indKyear])*tmp[K])
  rocTestKyear <- roc(dataTestKyear_outcome, tmp2, ci = TRUE)
  print(as.numeric(rocTestKyear$auc))
  
  print(length(dataTestKyear_outcome)) 
  print(mean(dataTestKyear_outcome))
  
  rocTestKyear$direction
  co = coords(rocTestKyear,"b", best.method="closest.topleft",ret=c("thre","spec","sens","ppv","npv","tp","tn","fp","fn"))
  print(co)
  
}
##

# cox longitudinal
fit_cox_lgt <- coxph(Surv(TimeToOutcome, Outcome)~., dataTrain[, -1])
CIndex_lgt <- survConcordance(Surv(TimeToOutcome, Outcome) ~ predict(fit_cox_lgt, newdata = dataTest), dataTest)
print(CIndex_lgt$concordance)

for (K in c(1,3,5)){
  print(K)
  dataTestKyear <- dataTest[!(dataTest$Outcome == 0 & dataTest$TimeToOutcome < K),]
  dataTestKyear_TimeToOutcome <- dataTestKyear$TimeToOutcome
  dataTestKyear$TimeToOutcome <- K
  predProbTestKyear <- 1 - exp(-predict(fit_cox_lgt, newdata = dataTestKyear[, -1], type = "expected"))
  dataTestKyear_outcome <- (dataTestKyear$Outcome == 1 & dataTestKyear_TimeToOutcome < K)
  rocTestKyear <- roc(dataTestKyear_outcome, predProbTestKyear, ci = TRUE)
  print(as.numeric(rocTestKyear$auc))
  
  rocTestKyear$direction
  co = coords(rocTestKyear,"b", best.method="closest.topleft",ret=c("thre","spec","sens","ppv","npv","tp","tn","fp","fn"))
  print(co)
}
##

#cross-sectional models
labName <- c('ALB','ALKRATIO','ALTRATIO','ASTRATIO','BIL','BUN','CI','CRE','GLU','HEM','PLT','K','NA', 'TOTP', 'WBC','APRI', 'ASTALT', 'BMI')
varStart <- c(names(dataAll)[1:9], paste(labName, "_LAST", sep=""))

i <- ind_best
dataTrain <- dataAll[is.element(dataAll$PatientID, ids_noTreated[indexSplit_noTreated[[i]]$indexTrain]), varStart]
dataTest <- dataAll[is.element(dataAll$PatientID, ids_noTreated[indexSplit_noTreated[[i]]$indexTest]), varStart]
dataTrain[, -c(1:9)] <- sapply(dataTrain[, -c(1:9)], fillByMedian)
dataTest[, -c(1:9)] <- sapply(dataTest[, -c(1:9)], fillByMedian)

# boosting cross-sectional
set.seed(90)
fit1 <- gbm(Surv(TimeToOutcome, Outcome)~., distribution = "coxph", data = dataTrain[, -1],
            n.trees = 2000, interaction.depth = 2, shrinkage = 0.03)
best.iter <- gbm.perf(fit1,method="OOB")
print(best.iter)
pred1 <- predict(fit1,newdata = dataTest,n.trees = best.iter)
CIndex1 <- survConcordance(Surv(TimeToOutcome, Outcome) ~ pred1, dataTest)
print(CIndex1$concordance)
tmp <- basehaz.gbm(dataTrain$TimeToOutcome, dataTrain$Outcome, predict(fit1,newdata = dataTrain,n.trees = best.iter), 1:5)

for (K in c(1, 3, 5)) {
  print(K)
  indKyear <- !(dataTest$Outcome == 0 & dataTest$TimeToOutcome < K)
  dataTestKyear_outcome <- (dataTest$Outcome[indKyear] == 1 & dataTest$TimeToOutcome[indKyear] < K)
  tmp2 <- 1 - exp(-exp(pred1[indKyear])*tmp[K])
  rocTestKyear <- roc(dataTestKyear_outcome, tmp2, ci = TRUE)
  print(as.numeric(rocTestKyear$auc))
  
  rocTestKyear$direction
  co = coords(rocTestKyear,"b", best.method="closest.topleft",ret=c("thre","spec","sens","ppv","npv","tp","tn","fp","fn"))
  print(co)
}
##

# cox cross-sectional
fit_cox_lgt <- coxph(Surv(TimeToOutcome, Outcome)~., dataTrain[, -1])
CIndex_lgt <- survConcordance(Surv(TimeToOutcome, Outcome) ~ predict(fit_cox_lgt, newdata = dataTest), dataTest)
print(CIndex_lgt$concordance)

for (K in c(1,3,5)){
  print(K)
  dataTestKyear <- dataTest[!(dataTest$Outcome == 0 & dataTest$TimeToOutcome < K),]
  dataTestKyear_TimeToOutcome <- dataTestKyear$TimeToOutcome
  dataTestKyear$TimeToOutcome <- K
  predProbTestKyear <- 1 - exp(-predict(fit_cox_lgt, newdata = dataTestKyear[, -1], type = "expected"))
  dataTestKyear_outcome <- (dataTestKyear$Outcome == 1 & dataTestKyear_TimeToOutcome < K)
  rocTestKyear <- roc(dataTestKyear_outcome, predProbTestKyear, ci = TRUE)
  print(as.numeric(rocTestKyear$auc))
  
  rocTestKyear$direction
  co = coords(rocTestKyear,"b", best.method="closest.topleft",ret=c("thre","spec","sens","ppv","npv","tp","tn","fp","fn"))
  print(co)
}
##