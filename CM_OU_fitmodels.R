#!/usr/bin/env Rscript
#--get-user-env

#Model Simulation & Fitting
#Dec 5

#Get basic environment
library(ape); library(phytools); library(MASS); library(scales); library(parallel)
#This file contains all the the equations and functions needed to fit the optimization
load("~/Functions_Mar2018.RData")

#read in the data
read.delim("./yourdata.txt", header=F) -> simdata
##Just a step to deal with my data not having a header row and having a first column of names rather than data
###what we want is a phenotypes-by-samples matrix
simdata[,1] -> names; names -> row.names(simdata)
simdata[,-1] -> simdata

#This is a step to parallelize. I ran this on the cluster with each R file taking 100 sites, replacing "NUMBER" with 1 to sites/100
MAXVAL <- 100*NUMBER
MINVAL <- MAXVAL-99
if (MAXVAL > nrow(simdata)) (MAXVAL <- nrow(simdata))

#Get set of sites to optimize now
simdata[MINVAL:MAXVAL,] -> simdata

#Fit function for each site
###When some parameters are too close to 0, terms within the covariance matrix are reduced to 0 which prevents optimization from working
###This constant keeps parameters above that level
min1 <- 10^-3
###Run in parallel using lApply
detectCores() -> cores
c1 <- makeCluster(cores)
clusterExport(c1, ls())
clusterEvalQ(c1, library(mvtnorm))
##Fit for 1 to the number of sites being run here
result <- parLapply(c1, 1:100, fit_multi, type="TYPE")
stopCluster(c1)

#save output
##This will need to be updated based upon the models you actually fit
##Here it includes outputs for drift (BM), stabilizing selection, and models with shifts on the northern lineage (ou2n), southern lineage (ou2s), anubis-Guinea lineage (ou2ag), hamadryas (ou2h), yellow (ou2y), and an identity matrix where all species are independent (id)
as.data.frame(matrix(ncol=43, nrow=100)) -> params
colnames(params) <- c("x0",'bm_delta','bm_tau', 'bm_lik','theta',"delta","alpha","tau","lik", 'ou2n_th1', 'ou2n_th2', 'ou2n_delta', 'ou2n_alpha', 'ou2n_tau','ou2n_lik','ou2s_th1', 'ou2s_th2', 'ou2s_delta', 'ou2s_alpha', 'ou2s_tau','ou2s_lik','ou2ag_th1', 'ou2ag_th2', 'ou2ag_delta', 'ou2ag_alpha', 'ou2ag_tau','ou2ag_lik','ou2h_th1', 'ou2h_th2', 'ou2h_delta', 'ou2h_alpha', 'ou2h_tau','ou2h_lik','ou2y_th1', 'ou2y_th2', 'ou2y_delta', 'ou2y_alpha', 'ou2y_tau','ou2y_lik','id_theta','id_delta','id_tau', 'id_lik')
for (i in 1:100) {
  result[[i]] -> params[i,]
}
row.names(params) <- names[MINVAL:MAXVAL]

write.table(params, "res_real_NUMBER.txt", row.names=T, col.names=F, sep="\t", quote=F)

#If type='BM', colnames = c('x0','delta','tau', 'lik')
#If type='OU', colnames = c('theta','delta','alpha','tau', 'lik')
#for DS models, colnames = c('theta1', 'theta2','delta','alpha','tau', 'lik')
