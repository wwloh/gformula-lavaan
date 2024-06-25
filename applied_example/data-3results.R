rm(list=ls())
library("data.table")

subfolder <- "data-2g_formula-flint/"
myfiles <- list.files(subfolder)
myfiles <- myfiles[grep(pattern=".Rdata",myfiles)]
simres <- lavfit <- NULL
for (ll in myfiles) {
  load(paste0(subfolder,ll))
  simres <- c(simres,list(msm.est$est))
  lavfit <- c(lavfit,list(msm.est$fit_lavaan))
  rm(msm.est)
}
sim_res <- do.call(rbind,simres); rm(simres)

msm.est <- sim_res[1,] # observed estimate
sim_res <- sim_res[-1,]
nrow(sim_res)

res <- lapply(c("unwt","dr"), function(x) {
  obs <- msm.est[grepl(pattern=x,names(msm.est))]
  boots <- sim_res[,grepl(pattern=x,colnames(sim_res))]
  res.x <- cbind(obs,
                 "se"=apply(boots,2,sd),
                 t(apply(boots,2,quantile,probs=c(.025,.975)))
                 )
  rownames(res.x) <- paste0("$psi_",1:length(obs),"$")
  return(res.x)
})
names(res) <- c("unwt","dr")
lapply(res,round, digits=2)

library("xtable")
lapply(res,xtable)

# fitted models
lavres <- lapply(c("unwt","dr"), function(x) {
  par.est.x <- lapply(lavfit,"[[",x)
  obs <- data.table(par.est.x[[1]])
  setkey(obs)
  boots <- data.table(do.call(rbind,par.est.x[-1]))
  setkey(boots)
  boots.ci <- boots[, as.list(
    c("se"=sd(est),quantile(est, probs=c(.025,.975)))),
    by=c("lhs","op","rhs")]
  res.x <- merge(obs,boots.ci)
  res.x <- cbind("meth"=x,res.x)
  setkey(res.x)
  return( res.x[op=="~" | (op=="~~" & lhs!=rhs)] )
})
write.csv(rbindlist(lavres),file="data-3results-lavaan_fits.csv", row.names=FALSE)

