# doubly robust estimation
rm(list=ls())
library("lavaan")
library("data.table")
library("MASS")
source("help_funs-g_formula.R")

# true parameter values
coefs <- c(0.1,0.2)
names(coefs) <- paste0("psi",c(1:2))

# sample size
N <- 1000

# generate data
# single Y observed at final time point
# single confounder L
OneData <- function(N, non.linear) {
  # unmeasured confounder
  U <- rnorm(N)
  # baseline confounder
  C <- rexp(N, rate=2)
  ## t=1
  L1 <- 2*U + rnorm(N)
  X1.ast <- 0.1*L1 + 0.1*C
  X1 <- rbinom(N,1,exp(X1.ast)/(1+exp(X1.ast)))
  rm(X1.ast)
  ## t=2
  L2 <- 1*L1 + 1*C + (X1-0.5) + 2*U + rnorm(N)
  X2.ast <- 0.3*L2 + 0.1*C
  X2 <- rbinom(N,1,exp(X2.ast)/(1+exp(X2.ast)))
  rm(X2.ast)
  ## t=3
  if (non.linear==FALSE) {
    Y <- 1*L1 + 1*C
  } else {
    ## arbitrary functions of the baseline covariates
    Y <- abs(L1)*log(C) + (U<0)*((C-0.5)^2) + (U>=0)*(L1<0)*1.0
  }
  a1 <- coefs["psi1"]+0.1
  Y <- Y + coefs["psi2"]*X2 + a1*X1 + (-0.1)*L2 + 2*U + rnorm(N)
  return(data.frame(C,L1,X1,L2,X2,Y))
}

OneSim <- function(non.linear) {
  Data <- OneData(N, non.linear)
  
  # vector of names of baseline pre-treatment covariates
  baseline.names <- c("C","L1")
  # list with vectors of names of post-treatment variables for sampling jointly
  post_treatment_variables <- list(
    "L2",
    "Y") # last entry is the end-of-study outcome
  
  # separate models defined by the treatment instances ========================
  FITmodels <- est_joint <- NULL
  ## models with X1 only
  FITmodels[[1]] <- '
    ## t=2
    L2 ~ L1 + b.l2_x1*X1 + C
  '
  ## models with X1 and X2
  FITmodels[[2]] <- '
    ## t=3
    Y ~ de2*X2 + de1*X1 + L2 + L1 + C
  '
  
  # single model without weights ==============================================
  FITmodel <- unlist(FITmodels)
  fit_joint <- lavaan::sem(FITmodel,data=Data,meanstructure=TRUE)
  est_joint <- lavaan::parameterEstimates(fit_joint)
  
  # plug in parameter estimates to construct g-formula
  GFmodel <- FITTED_MODEL.lavaan(dep.var=post_treatment_variables,
                                 est.lavaan=est_joint)
  
  res <- FitMSM(data=Data, 
                trt.name="X",
                trt.T=2,
                baseline.names=baseline.names, 
                post_treatment_variables=post_treatment_variables, 
                regfit=GFmodel,
                simple=TRUE,
                n.MC=100)
  
  res.unwt <- res
  rm(res)
  
  # compare to raw coefficient estimates
  poc.est <- subset(est_joint, label %in% c("de1","de2"))[,c("label","est")]
  res.poc <- poc.est[,"est"]
  names(res.poc) <- poc.est[,"label"]
  rm(poc.est)
  
  rm(FITmodel,fit_joint,est_joint,GFmodel)
  
  ## doubly robust variant ====================================================
  # fit PS models for treatments
  PS.models <- list("X1"=as.formula(X1 ~ L1 + C),
                    "X2"=as.formula(X2 ~ X1 + L2 + L1 + C))
  PS.fit <- lapply(PS.models, function(m) 
    glm(m, family=binomial('logit'), data=Data))
  PS.pred <- lapply(PS.fit,predict,type="response")
  # IPW for each treatment instance
  IPW <- lapply(names(PS.models), function(x) {
    Xt <- Data[,x]
    phat <- PS.pred[[x]]
    return( Xt/phat + (1-Xt)/(1-phat) )
  })
  names(IPW) <- names(PS.models)
  # cumulative product of weights
  IPW.cumprod <- t(apply(do.call(cbind,IPW),1,cumprod))
  colnames(IPW.cumprod) <- paste0("IPW..",names(IPW))
  ## standardize to sum to one
  IPW.cumprod <- apply(IPW.cumprod,2,function(ipw) ipw/sum(ipw))
  Data <- cbind(Data,IPW.cumprod)
  rm(PS.models,PS.fit,PS.pred,IPW,IPW.cumprod)
  
  est_joint <- NULL
  for (x in 1:length(FITmodels)) {
    fit_X <- lavaan::sem(FITmodels[[x]],data=Data,meanstructure=TRUE,
                         sampling.weights=paste0("IPW..X",x))
    est_X <- lavaan::parameterEstimates(fit_X)
    # parameter estimates specific to each post-treatment variable,
    # including intercept and (co)variances
    dv.x <- unique(subset(est_X, op=="~")[,c("lhs")])
    est_joint[[x]] <- subset(est_X, lhs %in% dv.x)
  }
  est_joint <- do.call(rbind,est_joint)
  
  # plug in parameter estimates to construct parametric g-formula
  GFmodel <- FITTED_MODEL.lavaan(dep.var=post_treatment_variables,
                                 est.lavaan=est_joint)

  res <- FitMSM(data=Data, 
                trt.name="X",
                trt.T=2,
                baseline.names=baseline.names, 
                post_treatment_variables=post_treatment_variables, 
                regfit=GFmodel,
                simple=TRUE,
                n.MC=100)
  
  res.dr <- res
  rm(res)
  
  return(unlist(list("msm"=list("unwt"=res.unwt,
                                "dr"=res.dr),
                     "poc"=list(res.poc))))
}

set.seed(30322)
simres <- lapply(c(FALSE,TRUE), function(x) 
  replicate(n=5000,OneSim(non.linear=x)))
save.image("g_formula-t3-lavaan-dr.Rdata")

# load and tidy results #######################################################
rm(list=ls())
load("g_formula-t3-lavaan-dr.Rdata")
round(do.call(rbind,lapply(simres, rowMeans)),2)

res.sum <- lapply(simres, function(res) {
  res <- t(res)
  t(sapply(c("poc","unwt","dr"), function(est.type) {
    res.type <- res[,grepl(pattern=est.type,colnames(res))]
    res.type <- res.type[,order(colnames(res.type))]
    colnames(res.type) <- paste0("X",1:2)
    bias.type <- t(t(res.type)-coefs)  # subtract by row 
    c("mean"=colMeans(res.type),
      "bias"=colMeans(bias.type),
      "ese"=apply(res.type,2,sd),
      "rmse"=sqrt(colMeans(bias.type^2))
    )
  }))
})
lapply(res.sum,round,digits=2)

library("xtable")
lapply(res.sum,xtable)
