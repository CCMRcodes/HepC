# This code performs cross-sectional boosting survival model with 30 random splits of training/test data.
# Authors of code: Boang Liu, Xuefei Zhang

library(gbm)
library(survival)
library(pROC)

fillByMedian <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
}

# Read in data in a data.frame named dataAll.
# Read in ids for 30 splits of training/test.

labName <- c('ALB','ALKRATIO','ALTRATIO','ASTRATIO','BIL','BUN','CI','CRE','GLU','HEM','PLT','K','NA', 'TOTP', 'WBC','APRI', 'ASTALT', 'BMI')
varStart <- c(names(dataAll)[1:9], paste(labName, "_LAST", sep=""))

# 30 splits
result <- matrix(NA, nrow = 30, ncol = 7)
colnames(result) <- c("CIndex", "1yearAUC", "1yearTestSize", "3yearAUC", "3yearTestSize", "5yearAUC", "5yearTestSize")
result <- as.data.frame(result)
result_iter <- rep(NA, 30)
result_importance <- vector("list", 30)

for (i in 1:30) {
  tm1 <- proc.time()
  print(i)
  dataTrain <- dataAll[is.element(dataAll$PatientID, ids_noTreated[indexSplit_noTreated[[i]]$indexTrain]), varStart]
  dataTest <- dataAll[is.element(dataAll$PatientID, ids_noTreated[indexSplit_noTreated[[i]]$indexTest]), varStart]
  dataTrain[, -c(1:9)] <- sapply(dataTrain[, -c(1:9)], fillByMedian)
  dataTest[, -c(1:9)] <- sapply(dataTest[, -c(1:9)], fillByMedian)
  
  set.seed(90)
  fit1 <- gbm(Surv(TimeToOutcome, Outcome)~., distribution = "coxph", data = dataTrain[, -1],
              n.trees = 2000, interaction.depth = 2, shrinkage = 0.03)
  
  best.iter <- gbm.perf(fit1,method="OOB")
  print(best.iter)
  pred1 <- predict(fit1,newdata = dataTest,n.trees = best.iter)
  result_iter[i] <- best.iter
  result_importance[[i]] <- summary(fit1, plotit = F)
  
  CIndex1 <- survConcordance(Surv(TimeToOutcome, Outcome) ~ pred1, dataTest)
  result$CIndex[i] <- CIndex1$concordance
  print(CIndex1$concordance)
  
  # K year auc for K=1,3,5
  tmpResult <- NULL
  for (K in c(1, 3, 5)) {
    indKyear <- !(dataTest$Outcome == 0 & dataTest$TimeToOutcome < K)
    dataTestKyear_outcome <- (dataTest$Outcome[indKyear] == 1 & dataTest$TimeToOutcome[indKyear] < K)
    rocTestKyear <- roc(dataTestKyear_outcome, pred1[indKyear], ci = TRUE)
    tmpResult <- c(tmpResult, rocTestKyear$auc, length(dataTestKyear_outcome))
  }
  result[i, 2:7] <- tmpResult
  print(tmpResult[c(1,2)])
  tm2 <- proc.time()
  print(tm2 - tm1)
}
print(sapply(result, mean))
print(sapply(result, sd))

result_cs <- result
result_iter_cs <- result_iter
result_importance_cs <- result_importance
# Save results