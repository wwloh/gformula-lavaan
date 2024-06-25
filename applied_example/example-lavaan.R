# load data and organize for analysis
source("data_1prep.R")
library("lavaan")
library("MASS")
source("../help_funs-g_formula.R")

Data <- Data.id
rm(Data.id)

# create interaction terms
Data[, LONELY__2_x_FRIEND__2 := LONELY__2*FRIEND__2]
Data[, LONELY__2_x_FAMILY__2 := LONELY__2*FAMILY__2]
Data[, LONELY__1_x_FRIEND__1 := LONELY__1*FRIEND__1]
Data[, LONELY__1_x_FAMILY__1 := LONELY__1*FAMILY__1]
setkey(Data)




# this is line 20
# model for outcome at wave 3
model_Y3 <- '
	FRIEND__3~LONELY__2+LONELY__2_x_FRIEND__2+LONELY__2_x_FAMILY__2+LONELY__1+FRIEND__2+FRIEND__1+FAMILY__2+FAMILY__1+GENDER+RACE+AGE+INCOME'

# model for time-varying variables at wave 2
model_Y2L2 <- '
	# conditional means
	FAMILY__2~LONELY__1+LONELY__1_x_FRIEND__1+LONELY__1_x_FAMILY__1+FRIEND__1+FAMILY__1+GENDER+RACE+AGE+INCOME
	FRIEND__2~LONELY__1+LONELY__1_x_FRIEND__1+LONELY__1_x_FAMILY__1+FRIEND__1+FAMILY__1+GENDER+RACE+AGE+INCOME
	# conditional covariance
  FAMILY__2~~FRIEND__2'

FITmodel <- c(model_Y3,model_Y2L2)
fit_joint <- lavaan::sem(FITmodel,data=Data,
                         meanstructure=TRUE, fixed.x=TRUE, se="none")
# extract parameter estimates from lavaan fitted models
est_joint <- lavaan::parameterEstimates(fit_joint)
subset(est_joint, subset=lhs=="FAMILY__2" & rhs=="FRIEND__2")
#          lhs op       rhs   est
# 31 FAMILY__2 ~~ FRIEND__2 0.112
subset(est_joint, subset=lhs=="FRIEND__3" & rhs=="LONELY__1")
#         lhs op       rhs    est
# 4 FRIEND__3  ~ LONELY__1 -0.214






# this is line 50
# create a list with vectors of names of post-treatment variables
post_treatment_variables <- list(
  c("FAMILY__2","FRIEND__2"), # wave 2
  "FRIEND__3") # wave 3
# plug in parameter estimates to construct g-formula
GFmodel <- FITTED_MODEL.lavaan(dep.var=post_treatment_variables,
                               est.lavaan=est_joint)












# this is line 70
res <- FitMSM(Data, 
              trt.name="LONELY__",
              trt.T=2,
              baseline.names=c("GENDER","RACE","AGE","INCOME","FAMILY__1","FRIEND__1"), 
              post_treatment_variables=post_treatment_variables, 
              regfit=GFmodel,
              n.MC=100)
round(res,2)
# LONELY__1 LONELY__2 
# -0.22     -0.18



















# this is line 100
# model for outcome at wave 3
# interaction terms X2_x_L2 and X2_x_Y2 were manually created in the data
model_Y3 <- '
	Y3 ~ X2 + L2 + Y2 + X2_x_L2 + X2_x_Y2 + X1 + L1 + Y1 + C
'
# model for time-varying variables at wave 2
# interaction terms X1_x_L1 and X1_x_L1 were manually created in the data
model_Y2L2 <- '
	# conditional means
	L2 ~ X1 + L1 + Y1 + X1_x_L1 + X1_x_Y1 + C
	Y2 ~ X1 + L1 + Y1 + X1_x_L1 + X1_x_Y1 + C
	# conditional covariance to be freely estimated
  L2 ~~ Y2
'
FITmodel <- c(model_Y3,model_Y2L2)
fit_joint <- lavaan::sem(FITmodel,data=Data,
                         meanstructure=TRUE, fixed.x=TRUE, se="none")
# create a list with vectors of names of post-treatment variables
post_treatment_variables <- list(
  c("L2","Y2"), # wave 2
  "Y3") # wave 3
# extract parameter estimates from lavaan fitted models
est_joint <- lavaan::parameterEstimates(fit_joint)
# plug in parameter estimates from steps A1 and A2 to construct g-formula
GFmodel <- FITTED_MODEL.lavaan(dep.var=post_treatment_variables,
                               est.lavaan=est_joint)



# this is line 130
res <- FitMSM(Data, 
              trt.name="X",
              trt.T=2,
              baseline.names=c("C","L1","Y1"), 
              post_treatment_variables=post_treatment_variables, 
              regfit=GFmodel,
              n.MC=100)
