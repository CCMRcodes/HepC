# This code performs cross-sectional and longitudinal Cox models with 30 random splits of training/test data.
# Authors of code: Boang Liu, Xuefei Zhang

library(survival)
library(pROC)

# Read in data in a data.frame named data_no_treatment.
# Read in ids for 30 splits of training/test.

fillByMedian <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
}
demo_var_list <- c('PatientID','Outcome','TimeToOutcome','GENDERM','RACEHISPANIC','RACEMISSING','RACEOTHER','RACEWHITE','AGE')
lab_var_list <- c('ALB','ALKRATIO','ALTRATIO','ASTRATIO','BIL','BUN','CI','CRE','GLU','HEM','PLT',
              'K','NA', 'TOTP', 'WBC','APRI', 'ASTALT', 'BMI')

data_no_treatment[,paste0(lab_var_list, "_TVARRAW")] <- NULL

result_baseline <- as.data.frame(matrix(0, nrow = 30, ncol = 10))
colnames(result_baseline) <- c("CIndex", "1yearAUC", "1yearTestSize", "1yearOutcomeProp", 
                               "3yearAUC", "3yearTestSize", "3yearOutcomeProp",
                               "5yearAUC", "5yearTestSize", "5yearOutcomeProp")
result_lgt <- as.data.frame(matrix(0, nrow = 30, ncol = 10))
colnames(result_lgt) <- c("CIndex", "1yearAUC", "1yearTestSize", "1yearOutcomeProp", 
                          "3yearAUC", "3yearTestSize", "3yearOutcomeProp",
                          "5yearAUC", "5yearTestSize", "5yearOutcomeProp")
set.seed(100)
for (i in 1:30){
  dataTrain <- data_no_treatment[is.element(data_no_treatment$PatientID, ids_noTreated[indexSplit_noTreated[[i]]$indexTrain]),]
  dataTest <- data_no_treatment[is.element(data_no_treatment$PatientID, ids_noTreated[indexSplit_noTreated[[i]]$indexTest]),]
  dataTrain[, -(which(names(dataTrain) %in% demo_var_list))] <- sapply(dataTrain[, -(which(names(dataTrain) %in% demo_var_list))], fillByMedian)
  dataTest[, -(which(names(dataTest) %in% demo_var_list))] <- sapply(dataTest[, -(which(names(dataTest) %in% demo_var_list))], fillByMedian)
  # fit baseline model
  baseline_var_names = c(demo_var_list,paste(lab_var_list, '_LAST', sep = ""))
  dataTrainBaseline <- dataTrain[,baseline_var_names]
  dataTestBaseline <- dataTest[, baseline_var_names]
  fit_cox_baseline <- coxph(Surv(TimeToOutcome, Outcome)~., dataTrainBaseline[, -1])
  CIndex_baseline <- survConcordance(Surv(TimeToOutcome, Outcome) ~ predict(fit_cox_baseline, newdata = dataTestBaseline[,-1]),dataTestBaseline[,-1])
  result_baseline$CIndex[i] <- CIndex_baseline$concordance
  tmp_result <- NULL
  for (K in c(1,3,5)){
    dataTestBaselineKyear <- dataTestBaseline[!(dataTestBaseline$Outcome == 0 & dataTestBaseline$TimeToOutcome < K),]
    dataTestBaselineKyear_TimeToOutcome <- dataTestBaselineKyear$TimeToOutcome
    dataTestBaselineKyear$TimeToOutcome <- K
    predProbTestBaselineKyear <- 1 - exp(-predict(fit_cox_baseline, newdata = dataTestBaselineKyear[, -1], type = "expected"))
    dataTestBaselineKyear_outcome <- (dataTestBaselineKyear$Outcome == 1 & dataTestBaselineKyear_TimeToOutcome < K)
    rocTestBaselineKyear <- roc(dataTestBaselineKyear_outcome, predProbTestBaselineKyear, ci = TRUE)
    tmp_result <- c(tmp_result, rocTestBaselineKyear$auc, nrow(dataTestBaselineKyear), mean(dataTestBaselineKyear_outcome))
  }
  result_baseline[i, 2:10] <- tmp_result
  print(result_baseline[i,])

  # fit lgt model
  lgt_to_remove <- c('_MEAN','_DIFFMEAN','_DDMEAN','_TVARRAW')
  var_to_remove <- outer(lab_var_list, lgt_to_remove, FUN = 'paste0')
  fit_cox_lgt <- coxph(Surv(TimeToOutcome, Outcome)~., dataTrain[, -c(1, which(names(dataTrain) %in% var_to_remove))])
  CIndex_lgt <- survConcordance(Surv(TimeToOutcome, Outcome) ~ predict(fit_cox_lgt, newdata = dataTest[,-c(1, which(names(dataTest) %in% var_to_remove))]),
                                dataTest[,-c(1, which(names(dataTest) %in% var_to_remove))])
  result_lgt$CIndex[i] <- CIndex_lgt$concordance
  tmp_result <- NULL
  for (K in c(1,3,5)){
    dataTestKyear <- dataTest[!(dataTest$Outcome == 0 & dataTest$TimeToOutcome < K),]
    dataTestKyear_TimeToOutcome <- dataTestKyear$TimeToOutcome
    dataTestKyear$TimeToOutcome <- K
    predProbTestKyear <- 1 - exp(-predict(fit_cox_lgt, newdata = dataTestKyear[, -1], type = "expected"))
    dataTestKyear_outcome <- (dataTestKyear$Outcome == 1 & dataTestKyear_TimeToOutcome < K)
    rocTestKyear <- roc(dataTestKyear_outcome, predProbTestKyear, ci = TRUE)
    tmp_result <- c(tmp_result, rocTestKyear$auc, nrow(dataTestKyear), mean(dataTestKyear_outcome))
  }
  result_lgt[i, 2:10] <- tmp_result
  print(result_lgt[i,])
}

# Save results
round(sapply(result_baseline, mean),3)
round(sapply(result_lgt, mean),3)
round(sapply(result_baseline, sd),3)
round(sapply(result_lgt, sd),3)

round(sapply(result_baseline[,c(1,2,5,8)], mean) - 1.96/sqrt(30) * sapply(result_baseline[,c(1,2,5,8)], sd),3)
round(sapply(result_lgt[,c(1,2,5,8)], mean) + 1.96/sqrt(30) * sapply(result_lgt[,c(1,2,5,8)], sd),3)

########################################
ttest_concordance <- t.test(result_baseline$CIndex, result_lgt$CIndex, paired = T)
ttest_concordance$p.value
ttest_1yearAUC <- t.test(result_baseline$`1yearAUC`, result_lgt$`1yearAUC`, paired = T)
ttest_1yearAUC$p.value
ttest_3yearAUC <- t.test(result_baseline$`3yearAUC`, result_lgt$`3yearAUC`, paired = T)
ttest_3yearAUC$p.value
ttest_5yearAUC <- t.test(result_baseline$`5yearAUC`, result_lgt$`5yearAUC`, paired = T)
ttest_5yearAUC$p.value
