# This code produces partial dependence plots for boosting cross-sectional model.
# Authors of code: Boang Liu, Xuefei Zhang

library(gbm)
library(survival)
library(pROC)

fillByMedian <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
}

# Read in data in a data.frame named dataAll.

labName <- c('ALB','ALKRATIO','ALTRATIO','ASTRATIO','BIL','BUN','CI','CRE','GLU','HEM','PLT','K','NA', 'TOTP', 'WBC','APRI', 'ASTALT', 'BMI')

varStart <- c(names(dataAll)[1:9], paste(labName, "_LAST", sep=""))
dataTrain <- dataAll[, varStart]
dataTrain[, -c(1:9)] <- sapply(dataTrain[, -c(1:9)], fillByMedian)

############################ partial dependene plots
# partial dependence for cross-sectionl model predictors
# Read in fitted results

x = dataTrain[, -c(1,2,3)]
rm(dataAll, dataTrain)
var_imp_all <- summary(fit1, plotit = F)
imp_names = as.character(var_imp_all[, 1])
impdata = x[imp_names]
length_v = rep(0, length(imp_names))
uniqueV = list()

for(i in 1:length(length_v)){
  uniqueV[[i]] = unique(impdata[,i])
  length_v[i] = length(uniqueV[[i]])
  if(length_v[i]>50){
    uniqueV[[i]] = quantile(impdata[,i], seq(0.01, 0.99, 0.02))
    length_v[i] = length(uniqueV[[i]])
  }
}

best.iter <- gbm.perf(fit1,method="OOB")

est = list()
for(i in 1:length(length_v)){
  tm1 <- proc.time()
  print(i)
  est[[i]] = rep(0, length_v[i])
  for(j in 1:length_v[i]){
    newdata_j = x
    newdata_j[imp_names[i]] = uniqueV[[i]][j]
    est[[i]][j] = mean(predict(fit1, newdata=newdata_j, n.trees=best.iter)) #predict on the log hazard scale
    gc()
  }
  
  #plot
  plot(uniqueV[[i]], est[[i]], xlab = imp_names[i], ylab = "Log Hazard")
  
  tm2 <- proc.time()
  print(tm2 - tm1)
}

# Save results