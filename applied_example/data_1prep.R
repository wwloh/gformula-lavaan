rm(list=ls())
library("data.table")

## prep data ##################################################################
# How to access the public-use data files: https://doi.org/10.3886/ICPSR34598.v1
# The following code is based on data files being in a subfolder "ICPSR_34598"
load("ICPSR_34598/DS0001/34598-0001-Data.rda")
load("ICPSR_34598/DS0002/34598-0002-Data.rda")
load("ICPSR_34598/DS0003/34598-0003-Data.rda")
load("ICPSR_34598/DS0004/34598-0004-Data.rda")
WAVE0 <- data.table(da34598.0001)
WAVE0[,RESPID := NULL]
setkey(WAVE0)
WAVE1 <- data.table(da34598.0002)
WAVE1[,RESPID := NULL]
setkey(WAVE1)
WAVE2 <- data.table(da34598.0003)
WAVE2[,RESPID := NULL]
setkey(WAVE2)
WAVE3 <- data.table(da34598.0004)
WAVE3[,RESPID := NULL]
setkey(WAVE3)
rm(da34598.0001,da34598.0002,da34598.0003,da34598.0004)

rawData <- Reduce(merge,list(WAVE0,WAVE1,WAVE2,WAVE3))
setkey(rawData)
rm(WAVE0,WAVE1,WAVE2,WAVE3)

# Past week how often: felt lonely
rawData[, LONELY__1 := ifelse(W2G1D=="(1) Not at all uncomfortable", 0, 1)]
rawData[, LONELY__2 := ifelse(W3G1D=="(1) Not at all uncomfortable", 0, 1)]
table(rawData[, list(LONELY__1,LONELY__2)])

# parental support
rawData[, FAMILY__1 := rowMeans(cbind(
  rawData[, lapply(.SD, as.integer), 
          .SDcols=paste0("W2C4",LETTERS[1:6])],
  5-rawData[, lapply(.SD, as.integer), 
            .SDcols=paste0("W2C4","G")]), na.rm=TRUE)-1.0]
rawData[, FAMILY__2 := rowMeans(cbind(
  rawData[, lapply(.SD, as.integer), 
          .SDcols=paste0("W3C4",LETTERS[1:6])],
  5-rawData[, lapply(.SD, as.integer), 
            .SDcols=paste0("W3C4","G")]), na.rm=TRUE)-1.0]
summary(rawData[, list(FAMILY__1,FAMILY__2)])

# friend support
rawData[, FRIEND__1 := rowMeans(rawData[, lapply(.SD, as.integer), 
                                        .SDcols=paste0("W2C1",LETTERS[1:5])], na.rm=TRUE)]
rawData[, FRIEND__2 := rowMeans(rawData[, lapply(.SD, as.integer), 
                                        .SDcols=paste0("W3C1",LETTERS[1:5])], na.rm=TRUE)]
rawData[, FRIEND__3 := rowMeans(rawData[, lapply(.SD, as.integer), 
                                        .SDcols=paste0("W4C1",LETTERS[1:5])], na.rm=TRUE)]
summary(rawData[, list(FRIEND__1,FRIEND__2,FRIEND__3)])

# W1A1_3: Participant's year of birth
# W1A2: Participant's gender
# W1A3: Participant's race
# W2B13: Family receives income support
summary(rawData[, list(W1A1_3,W1A2,W1A3,W2B13)])

Data <- rawData
rm(rawData)

Data[, id := as.integer(CASEID)]
Data[, GENDER := ifelse(W1A2=="(2) Female", 1, 0)]
Data[, RACE := ifelse(W1A3=="(1) Black or African American", 1, 0)]
Data[, AGE := (W1A1_3==1980)*1L]
Data[, INCOME := ifelse(W2B13=="(1) Yes", 0, 1)]

Data.id <- Data[,list(id,GENDER,RACE,AGE,INCOME,
                      LONELY__1,LONELY__2,
                      FAMILY__1,FAMILY__2,
                      FRIEND__1,FRIEND__2,FRIEND__3)]
rm(Data)
setkey(Data.id)

# remove missing data
Data.id <- Data.id[!apply(apply(Data.id,1,is.na),2,any)]
setkey(Data.id)
nrow(Data.id)

summary(Data.id)

# compare treatments over time
round(table(Data.id[,list(LONELY__1,LONELY__2)])/nrow(Data.id),2)

# mean-center continuous variables
Data.id[, FAMILY__1 := scale(FAMILY__1, center=TRUE, scale=FALSE)]
Data.id[, FAMILY__2 := scale(FAMILY__2, center=TRUE, scale=FALSE)]
Data.id[, FRIEND__1 := scale(FRIEND__1, center=TRUE, scale=FALSE)]
Data.id[, FRIEND__2 := scale(FRIEND__2, center=TRUE, scale=FALSE)]
Data.id[, FRIEND__3 := scale(FRIEND__3, center=TRUE, scale=FALSE)]

summary(Data.id)

library("xtable")
print(xtable(Data.id[(1:7)*100], digits=c(rep(0,8),rep(2,5))),
      include.rownames=FALSE)

# relabel time-varying variable names with time point indicator
trtT <- 2
Xname <- "LONELY__"
Cnames <- c("GENDER","RACE","AGE","INCOME","FAMILY__1","FRIEND__1")
Lnames <- "FAMILY__"
Yname <- "FRIEND__"
finalY <- "FRIEND__3"
