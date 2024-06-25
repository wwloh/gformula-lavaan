# extract model parameter estimates from fitted lavaan model
FITTED_MODEL.lavaan <- function(dep.var, est.lavaan) {
  GFmodel <- NULL
  for (v.names in dep.var) {
    v.ests <- subset(est.lavaan, lhs %in% v.names & 
                       op %in% c("~", "~~", "~1", "|"))[,c(
                         "lhs","op","rhs","est")]
    # initialize results
    ## mean vector
    v.est.mean <- sapply(v.names,function(x) NULL)
    ## covariance matrix
    v.est.cov <- matrix(NA,nrow=length(v.names),ncol=length(v.names),
                        dimnames=list(v.names,v.names))
    ## type of variable
    v.type <- unlist(sapply(v.names,function(x) "gaussian"))
    # each variable in turn
    for (v.name in v.names) {
      v.est <- subset(v.ests, lhs==v.name)
      # regression coefficients
      reg.coef <- rep(NA,sum(v.est$op=="~")+1)
      if ("|" %in% v.est$op) {
        v.type[v.name] <- "binomial" # binary variable
        ## intercept
        reg.coef[1] <- -1*subset(v.est, op=="|")[,"est"]
      } else {
        ## intercept
        reg.coef[1] <- subset(v.est, op=="~1")[,"est"]
      }
      names(reg.coef)[1] <- "(Intercept)"
      if (length(reg.coef)>1) {
        ## slopes
        reg.coef[-1] <- subset(v.est, op=="~")[,"est"]
        names(reg.coef)[-1] <- subset(v.est, op=="~")[,"rhs"]
      }
      v.est.mean[[v.name]] <- reg.coef
      # residual (co)variance(s)
      s2 <- subset(v.est, op=="~~" & lhs==v.name)
      ## fill in symmetric covariance matrix
      v.est.cov[v.name,s2$rhs] <- v.est.cov[s2$rhs,v.name] <- s2$est
    }
    GFmodel.v <- list("reg"=v.est.mean,"resid.var"=v.est.cov,"type"=v.type)
    GFmodel <- c(GFmodel,list(GFmodel.v))
    rm(GFmodel.v)
  }
  names(GFmodel) <- lapply(dep.var,paste,collapse=",")
  return(GFmodel)
}


# create duplicated data with all treatment combinations for a given dataset
Dupdata <- function(data, # observed data
                    trt.name, # name of treatment variable
                    trt.T, # number of treatment time points
                    baseline.names # pre-treatment covariates
                    ) {
  
  # helper function for a single individual
  AllXs <- function(trt.name,t) {
    out <- expand.grid(lapply(1:t, function(i) 0:1))
    colnames(out) <- paste0(trt.name,1:t)
    out <- cbind("x"=1:nrow(out),out)
    return(out)
  }
  
  N <- nrow(data) # sample size
  D <- data.table("id"=1:N,"(Intercept)"=1,data)
  setkey(D)
  to_keep <- c("id","(Intercept)",baseline.names)
  dat <- D[,..to_keep]
  setkey(dat)
  alevels <- dat[,as.data.table(AllXs(trt.name,trt.T)),by=id]
  setkey(alevels)
  dat <- merge(dat,alevels,all.x=TRUE)
  setkey(dat)
  rm(alevels,D)
  return(dat)
}

# make random draws of post-treatment variables
OneMCdraw <- function(dupdata, # duplicated data from Dupdata function
                      trt.name, trt.T, 
                      post_treatment_variables, # post-treatment variables
                      regfit # fitted regression models
                      ) {
  # final outcome
  Y.name <- rev(post_treatment_variables)[[1]]
  
  datMC <- data.table(dupdata)
  setkey(datMC)
  for (v.names in post_treatment_variables) {
    fitted_pars <- regfit[[paste(v.names,collapse=",")]]
    # predicted counterfactuals
    drawn.mean <- sapply(v.names,function(x) NULL)
    for (v.name in v.names) {
      reg.coef <- fitted_pars$reg[[v.name]]
      coef_names <- names(reg.coef)
      # check for interaction terms
      int.terms <- grep(pattern="_x_",x=coef_names,value=TRUE)
      if (length(int.terms)>0) {
        for (ints in int.terms) {
          ints.v <- strsplit(x=ints,split="_x_")[[1]]
          datMC[, ints := Reduce(`*`, .SD), .SDcols = ints.v]
          setnames(datMC,"ints",ints)
          setkey(datMC)
        }
      }
      drawn.mean[[v.name]] <- drop(as.matrix(
        datMC[, ..coef_names]) %*% reg.coef)
      rm(reg.coef,coef_names)
    }
    drawn.mean <- do.call(cbind,drawn.mean)
    if (all(v.names==Y.name)) {
      # set to predicted outcome (random draw is unnecessary)
      drawn.vals <- drop(drawn.mean)
      if (fitted_pars$type=="binomial") {
        # predicted probability
        drawn.vals <- pnorm(drawn.vals)
      }
    } else {
      # random draw from normal distribution
      if (nrow(fitted_pars$resid.var)>1) {
        drawn.vals <- as.list(data.frame(t(apply(drawn.mean, 1, function(x) 
          MASS::mvrnorm(n=1, mu=x, Sigma=fitted_pars$resid.var)))))
      } else {
        drawn.vals <- list(rnorm(n=nrow(drawn.mean),
                                 mean=drop(drawn.mean),
                                 sd=drop(fitted_pars$resid.var)))
      }
      names(drawn.vals) <- v.names
    }
    datMC[, eval(v.names) := drawn.vals]
    setkey(datMC)
    rm(fitted_pars,drawn.mean,drawn.vals)
  }
  
  # return only treatment values and outcome to save memory
  return.cols <- c(paste0(trt.name,1:trt.T),Y.name)
  return(datMC[!is.na(get(Y.name)),..return.cols])
}

FitMSM <- function(data, trt.name, trt.T,
                   baseline.names, post_treatment_variables, regfit, 
                   simple=TRUE, n.MC=100) {
  
  dd <- Dupdata(data,trt.name,trt.T,baseline.names)
  
  dat.msm <- rbindlist(replicate(
    n=n.MC, 
    expr=OneMCdraw(dupdata=dd,trt.name,trt.T,post_treatment_variables,regfit), 
    simplify=FALSE))
  
  # fit MSM
  Y.name <- rev(colnames(dat.msm))[1]
  msm.form <- as.formula(
    paste0(Y.name,"~",
           paste(paste0(trt.name,1:trt.T),
                 collapse=ifelse(simple==TRUE,"+","*"))))
  
  return(coef(lm(msm.form,data=dat.msm))[-1])
}