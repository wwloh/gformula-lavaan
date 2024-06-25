# load data and organize for analysis
source("data_1prep.R")
library("lavaan")
library("MASS")
source("../help_funs-g_formula.R")

for (seed in 1:1000) {
## set seed to differentiate the jobs #########################################
set.seed(30322+seed) # ensures different seed for each bootstrap

OneEst <- function(Data, n.MC=100) {
  
  # sample size
  N <- nrow(Data)
  # vector of names of baseline pre-treatment covariates
  baseline.names <- Cnames
  # list with vectors of names of post-treatment variables for sampling jointly
  post_treatment_variables <- c(
    lapply(2:trtT,function(tt) paste0(c(Lnames,Yname),tt)),
    finalY) # last entry is the end-of-study outcome
  
  # construct regression models to be fited ###################################
  # predictors of each treatment instance 
  Xt.preds <- lapply(1:trtT, function(tt) {
    # historical values
    # previous treatment(s)
    hist.idx <- tt-(1:1)
    hist.idx <- hist.idx[hist.idx>0]
    if(length(hist.idx)==0) {
      hist.X <- NULL
    } else {
      hist.X <- paste0(Xname,hist.idx)
    }
    # time-varying covariates
    hist.idx <- c(tt,hist.idx)
    hist.Y <- paste0(Yname,hist.idx) # previous outcomes
    hist.L <- as.vector(sapply(Lnames, function(l) paste0(l,hist.idx)))
    # remove duplicates already in baseline covariates
    hist.Y <- hist.Y[!(hist.Y %in% Cnames)]
    hist.L <- hist.L[!(hist.L %in% Cnames)]
    Xtt.preds <- c(hist.X,hist.Y,hist.L,Cnames)
    # keep only those in the dataset
    Xtt.preds <- Xtt.preds[Xtt.preds %in% colnames(Data)]
    # drop any that are singular
    Xtt.preds <- Xtt.preds[apply(Data[,..Xtt.preds],2,var,na.rm=TRUE)>0]
    return( Xtt.preds )
  })
  names(Xt.preds) <- paste0(Xname,1:trtT)
  
  # separate models defined by the treatment instances
  FITmodels <- est_joint <- NULL
  FITmodels <- lapply(1:trtT, function(tt) {
    Lt <- post_treatment_variables[[tt]]
    # regression models
    reg.tt <- sapply(Lt, function(lt) {
      lt.xt <- names(Xt.preds)[[tt]] # treatment at time t
      lt.pred <- Xt.preds[[tt]] # predictors of treatment at time t
      ## interactions between treatment at time t and time-varying variables
      lt.ints <- unlist(sapply(c(Yname,Lnames), function(lt.int) {
        lt.int.term <- grep(lt.int,lt.pred,value=TRUE)
        if (length(lt.int.term)>0) {
          # just the latest instance
          lt.int.term <- sort(lt.int.term,decreasing=TRUE)[1]
          return(paste0(paste0(lt.xt,"_x_"),lt.int.term))
        }
      }))
      names(lt.ints) <- NULL
      lt.formula <- paste0(
        paste(c(lt,"~"),collapse=""), # time-varying covariate at time t+1
        paste(c(lt.xt,lt.ints,lt.pred),collapse="+")
      )
      return(lt.formula)
    })
    if (length(Lt)>1) {
      cov.tt <- NULL
      for (i in 1:length(Lt)) {
        for (j in (i+1):length(Lt)) {
          if (j>i & j<=length(Lt)) {
            cov.tt <- c(cov.tt,paste0(Lt[i],"~~",Lt[j]))
          }
        }
      }
      return(c(reg.tt,cov.tt))
    } else {
      return(reg.tt)
    }
  })
  
  # relevant interaction terms in dataset
  for (fmi in 1:length(FITmodels)) {
    fm <- FITmodels[[fmi]]
    int.names <- grep("_x_",strsplit(fm,split="[+]")[[1]],value=TRUE)
    for (int.name in int.names) {
      if (!(int.name %in% colnames(Data))) {
        # create the product
        int.name.terms <- strsplit(int.name,"_x_")[[1]]
        int.term <- drop(apply(Data[,..int.name.terms],1,prod))
        Data[, eval(int.name) := int.term]
      }
    }
  }
  setkey(Data)
  Data <- data.frame(Data)
  
  # single model without weights ##############################################
  FITmodel <- unlist(FITmodels)
  fit_joint <- lavaan::sem(FITmodel,data=Data,
                           meanstructure=TRUE, fixed.x=TRUE, se="none")
  est_joint <- lavaan::parameterEstimates(fit_joint)
  
  # plug in parameter estimates to construct g-formula
  GFmodel <- FITTED_MODEL.lavaan(dep.var=post_treatment_variables,
                                 est.lavaan=est_joint)
  
  res <- FitMSM(data=Data, 
                trt.name=Xname,
                trt.T=trtT,
                baseline.names=baseline.names, 
                post_treatment_variables=post_treatment_variables, 
                regfit=GFmodel,
                simple=TRUE,
                n.MC=n.MC)
  
  res.unwt <- res
  rm(res)
  
  res_lavaan.unwt <- est_joint[,c("lhs","op","rhs","est")]
  
  rm(FITmodel,fit_joint,est_joint,GFmodel)
  
  ## doubly robust variant ####################################################
  # fit PS models for treatments
  PS.models <- lapply(1:trtT, function(tt) {
    return( as.formula(paste0(
      paste(c(Xname,tt,"~"),collapse=""), # treatment at time t
      paste(Xt.preds[[tt]], collapse="+")
    )) )
  })
  names(PS.models) <- names(Xt.preds)
  PS.fit <- lapply(PS.models, function(m) 
    glm(m, family=binomial('logit'), data=Data))
  PS.pred <- lapply(PS.fit,predict,type="response")
  # IPW for each treatment instance
  IPW <- lapply(names(PS.models), function(x) {
    # observed treatment
    Xt <- Data[,x]
    names(Xt) <- NULL
    # predicted PS
    Pt <- PS.pred[[x]]
    phat <- rep(NA,N)
    phat[as.integer(names(Pt))] <- Pt
    rm(Pt)
    # inverse weights
    ipw <- Xt/phat + (1-Xt)/(1-phat)
    # set unit weight for missing propensity score
    ipw[is.na(phat)] <- 1.0
    # set zero weights for missing observed treatments
    ipw[is.na(Xt)] <- 0.0
    return( ipw )
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
                         sampling.weights=paste0("IPW..",names(Xt.preds))[x])
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
                trt.name=Xname,
                trt.T=trtT,
                baseline.names=baseline.names, 
                post_treatment_variables=post_treatment_variables, 
                regfit=GFmodel,
                simple=TRUE,
                n.MC=n.MC)
  
  res.dr <- res
  rm(res)
  
  res_lavaan.dr <- est_joint[,c("lhs","op","rhs","est")]
  
  return(list(
    "est"=unlist(list("unwt"=res.unwt,
                      "dr"=res.dr)),
    "fit_lavaan"=list("unwt"=res_lavaan.unwt,
                      "dr"=res_lavaan.dr))
  )
}

ptm <- proc.time()[3]
if (seed==1) {
  msm.est <- OneEst(Data=Data.id, n.MC=100)
} else {
   # resample with replacement for bootstrap
  ids <- unique(Data.id$id)
  boot.i <- sort(sample(x=ids,size=length(ids),replace=TRUE))
  Data.boot <- do.call(rbind,lapply(1:length(boot.i), function(ib) {
    D.ib <- Data.id[id==boot.i[ib]]
    D.ib[,id := ib]
    return(D.ib)
  }))
  setkey(Data.boot)
  msm.est <- OneEst(Data=Data.boot, n.MC=100)
}
round((proc.time()[3]-ptm)/60)
# 70 mins for 100 MC draws

myfilename <- "data-2g_formula-flint/bootstrap_results"
save(msm.est,file=paste0(myfilename,"-seed_",seed,".Rdata"))

}
