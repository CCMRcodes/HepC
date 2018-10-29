# This code performs variable selection for longitudinal Cox model with 30 random splits of training/test data.
# Authors of code: Boang Liu, Xuefei Zhang

library(glmnet)
library(survival)
library(pROC)

fillByMedian <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
}

# Read in data in a data.frame named dataAll.
# Read in ids for 30 splits of training/test.

alpha <- 0.5 #elastic net parameter
labName <- c('ALB','ALKRATIO','ALTRATIO','ASTRATIO','BIL','BUN','CI','CRE','GLU','HEM','PLT','K','NA', 'TOTP', 'WBC','APRI', 'ASTALT', 'BMI')
varStart <- setdiff(names(dataAll), c(paste(labName, "_MEAN", sep=""), paste(labName, "_DIFFMEAN", sep=""),
                                      paste(labName, "_DDMEAN", sep=""), paste(labName, "_TVARRAW", sep="")))

# 30 splits
cv_fit_all <- vector("list", 30)
CIndex_all <- rep(NA, 30)
pvalues_all <- vector("list", 30)
varFinal_all <- vector("list", 30)

for (i in 1:30) {
  tm1 <- proc.time()
  print(i)
  dataTrain <- dataAll[is.element(dataAll$PatientID, ids_noTreated[indexSplit_noTreated[[i]]$indexTrain]), varStart]
  dataTest <- dataAll[is.element(dataAll$PatientID, ids_noTreated[indexSplit_noTreated[[i]]$indexTest]), varStart]
  dataTrain[, -c(1:9)] <- sapply(dataTrain[, -c(1:9)], fillByMedian)
  dataTest[, -c(1:9)] <- sapply(dataTest[, -c(1:9)], fillByMedian)
  
  set.seed(115)
  cv_fit <- cv.glmnet(x=as.matrix(dataTrain[,-c(1,2,3)]),y=Surv(dataTrain$TimeToOutcome, dataTrain$Outcome),
                      family="cox",alpha=alpha,nfolds=10)
  
  # fit unpenalized cox model with selected variables from elastic net
  nameTmp <- rownames(coef(cv_fit,s="lambda.1se"))[which(coef(cv_fit,s="lambda.1se")!=0)]
  varSelect <- c("Outcome","TimeToOutcome",nameTmp)
  print(length(nameTmp))
  
  set.seed(150)
  fit_cox_select <- coxph(Surv(TimeToOutcome, Outcome)~., dataTrain[, varSelect])
  CIndexSelect <- survConcordance(Surv(TimeToOutcome, Outcome) ~ predict(fit_cox_select,newdata = dataTest[,nameTmp]),dataTest[,varSelect])
  print(CIndexSelect$concordance)
  pvalues_cox <- summary(fit_cox_select)$coefficients[,5]
  
  #save cv_fit, concordance, pvalues, variables
  cv_fit_all[[i]] <- cv_fit
  CIndex_all[i] <- CIndexSelect$concordance
  pvalues_all[[i]] <- pvalues_cox
  varFinal_all[[i]] <- names(pvalues_cox)
  
  tm2 <- proc.time()
  print(tm2-tm1)
}

##################### evaluation
head(varStart, 20)
length(varStart) - 3

mean(CIndex_all)
sd(CIndex_all)

mean(sapply(varFinal_all, length))
sd(sapply(varFinal_all, length))
summary(sapply(varFinal_all, length))

mean(sapply(pvalues_all, function(x) mean(x < 0.1)))
sd(sapply(pvalues_all, function(x) mean(x < 0.1)))

varIntersect <- Reduce(intersect, varFinal_all)
length(varIntersect)
varIntersect
length(Reduce(union,varFinal_all))

varAll <- unlist(varFinal_all)
tableVar <- sort(table(varAll),decreasing=T)
sum(tableVar >= 20)
print(tableVar)