# This code produces variable importance plots for boosting cross-sectional and longitudinal models.
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

# cross-sectional model fit and variable importance
varStart <- c(names(dataAll)[1:9], paste(labName, "_LAST", sep=""))
dataTrain <- dataAll[, varStart]
dataTrain[, -c(1:9)] <- sapply(dataTrain[, -c(1:9)], fillByMedian)

tm1 <- proc.time()
set.seed(90)
fit1 <- gbm(Surv(TimeToOutcome, Outcome)~., distribution = "coxph", data = dataTrain[, -1],
            n.trees = 2000, interaction.depth = 2, shrinkage = 0.03)
tm2 <- proc.time()
print(tm2 - tm1)

tmp <- summary(fit1, plotit = F)
tmp <- tmp[nrow(tmp):1,]
var_importance <- tmp[, 2]
names(var_importance) <- tmp[, 1]

# plot
barplot(var_importance,space=0.5,horiz=T,cex.names=0.75,axis.lty=0.5,las=1,mgp=c(3,0.7,0),xlim=c(0,50))

# longitudinal model
varStart <- setdiff(names(dataAll), c(paste(labName, "_MEAN", sep=""), paste(labName, "_DIFFMEAN", sep=""),
                                      paste(labName, "_DDMEAN", sep=""), paste(labName, "_TVARRAW", sep="")))
dataTrain <- dataAll[, varStart]
dataTrain[, -c(1:9)] <- sapply(dataTrain[, -c(1:9)], fillByMedian)

tm1 <- proc.time()
set.seed(90)
fit1 <- gbm(Surv(TimeToOutcome, Outcome)~., distribution = "coxph", data = dataTrain[, -1],
            n.trees = 2000, interaction.depth = 2, shrinkage = 0.03)
tm2 <- proc.time()
print(tm2 - tm1)

tmp <- summary(fit1, plotit = F)
tmp <- tmp[30:1,]
var_importance <- tmp[, 2]
names(var_importance) <- tmp[, 1]

# plot
barplot(var_importance,space=0.5,horiz=T,cex.names=0.75,axis.lty=0.5,las=1,mgp=c(3,0.7,0),xlim=c(0,25))