#Likelihoods and functions for OU model optimization
library(ape); library(phytools); library(MASS); library(scales)

#Setup: get covariance matrix and individuals per species
{#read in the phylogenetic tree, convert to covariance matrix, and set up the study design (individuals per species)
{read.tree("~/CM/Newick_NoChacma.txt") -> tree
  vcv.phylo(tree) -> covar_matrix
  covar_matrix[c(1,4,3,2,5,6),c(1,4,3,2,5,6)] -> covar_matrix #I needed to reorder the matrix to match my other files
  covar_matrix <- covar_matrix/max(covar_matrix) #just scales the matrix, this doesn't matter and just relates to delta
  rm(tree)}

###Get number of individuals per species & produce the z-matrix (which we will later multiply by the covariance matrix)
{read.delim('~/CM/CM44_nindiv.txt', header=F, sep=' ') -> indiv_spec
  sum(indiv_spec) -> nindiv
  for (i in 2:6) {indiv_spec[i] <- sum(indiv_spec[(i-1):i])}; rm(i)
  nspec <- length(indiv_spec)
  zmat <- matrix(0,nrow=nindiv, ncol=nspec)
  for (k in seq(1,nspec)) {
    if (k == max(nspec)) {
      mini <- indiv_spec[k-1]+1
      maxi <- nindiv
    }
    if (k < max(nspec) & k > 1) {
      mini <- indiv_spec[k-1]+1
      maxi <- indiv_spec[k]
    }
    if (k == 1) {
      mini <- 1
      maxi <- indiv_spec[k]
    }
    zmat[c(seq(unlist(mini),unlist(maxi))),k] <- 1
  }
  rm(k); rm(maxi); rm(mini)}}

#fit_multi: function for fitting using parallel
##This is the function that will fit all models for a site and return a list of parameter estimates and likelihoods
###uses functions: BM_all, OU_all, OU2_all, indep_All
###Requires: simdata, min

fit_multi <- function(number, type) {
  #number is the site number in the data table file
  #type is the type of model being run: BM for brownian motion, OU for stabilizing selection, or the lineage name for directional selection ('north', 'south', 'anugui', 'anubis', 'guinea', 'hamadryas', 'kinda', 'yellow')
  n <- as.numeric(number[1])
  #Get data
  data1 = simdata[n,]
  #Guess at the ancestral state using the mean across all individuals
  prevU = mean(unlist(data1), na.rm=T)
  
  #Based upon "type", fit the appropriate BM or OU model
  
  if (type[1]="BM") {
    #Fit BM model
    ##Make temporary files for the parameter estimates (model_pars) and -log-likelihood (Model_Lik)
    ###We'll use these objects to store the best fitting model
    Model_Lik <- Inf; model_pars <- c(NA, NA, NA)
    ##Prepare 'output' to be filled with optimization output
    output <- NULL
    
    #Grid search through starting parameters
    #Optimize BM model with those starting parameters
    #Replace parameters & likelihood if (i) the model fits better than the previous best model (rounded to 6 decimal places) and (ii) the optimization found a definite solution (-log likehood != -Inf)
    for (testS in seq(0,10,1)) { #set of values for sigma2 to test
      for (testT in seq(0,30,5)) { #set of values of tau2 to test
        output <- optim(par=c(prevU, testS, testT), data1=data1, fn=BM_all, method="BFGS")
        output$value -> lik2
        if (lik2 <= ceiling(Model_Lik*100000)/100000 & output$value > -Inf) {
          Model_Lik <- output$value; output$par -> model_pars[1:3] }
      }
    }
    #Store best parameters for the site in a single table
    as.data.frame(matrix(ncol=4,nrow=1)) -> pars
    model_pars -> pars[1,1:3]; Model_Lik -> pars[1,4]
    }
  
  else { 
    if (type[1]="OU") {
    #Fit OU model for stabilizing selection
    ##Make temporary files for the parameter estimates (model_pars) and -log-likelihood (Model_Lik)
    ###We'll use these objects to store the best fitting model
    Model_Lik <- Inf; model_pars <- c(NA, NA, NA, NA)
    ##Prepare 'output' to be filled with optimization output
    output <- NULL
    
    #Grid search through starting parameters
    #Optimize OU model with those starting parameters
    #Replace parameters & likelihood if (i) the model fits better than the previous best model (rounded to 6 decimal places) and (ii) the optimization found a definite solution (-log likehood != -Inf)
    for (testS in seq(0,10,1)) { #set of values for sigma2 to test
      for (testT in seq(0,30,5)) { #set of values of tau2 to test
        for (testA in c(.1,2,4,10)) { #set of values for alpha to test
          if (testA < min1) {testA <- min1} #mostly unnecessary, but depending on the phylogeny you may need to set a minimum level for testA
          output <- optim(par=c(prevU, testS, testA, testT), data1=data1, fn=OU_all, method="BFGS")
          output$value -> lik2
          if (lik2 <= ceiling(Model_Lik*100000)/100000 & output$value > -Inf) {
            Model_Lik <- output$value; output$par -> model_pars[1:4] }
        }
      }
    }
    #Store best parameters for the site in a single table
    as.data.frame(matrix(ncol=5,nrow=1)) -> pars
    model_pars -> pars[1,1:4]; Model_Lik -> pars[1,5]
    }
    
    else {
      #Fit OU model for directional selection
      ##Make temporary files for the parameter estimates (model_pars) and -log-likelihood (Model_Lik)
      ###We'll use these objects to store the best fitting model
      Model_Lik <- Inf; model_pars <- c(NA, NA, NA, NA, NA)
      ##Prepare 'output' to be filled with optimization output
      output <- NULL
      
      #Grid search through starting parameters
      #Optimize OU model with those starting parameters
      #Replace parameters & likelihood if (i) the model fits better than the previous best model (rounded to 6 decimal places) and (ii) the optimization found a definite solution (-log likehood != -Inf)
      for (testS in seq(0,10,1)) { #set of values for sigma2 to test
        for (testT in seq(0,30,5)) { #set of values of tau2 to test
          for (testA in c(.01,.1,1,4,10)) { #set of values for alpha to test
            if (testA < min1) {testA <- min1} #mostly unnecessary, but depending on the phylogeny you may need to set a minimum level for testA
            output <- optim(par=c(prevU, prevU, prevS, testA, prevT), data1=data1, lineage=type[1], fn=OU2_all, method="BFGS")
              #PrevU is repeated twice in the parameters so it can be optimized for the lineage with ancestral DNA methylation and the lineage with a shift due to positive selection
            output$value -> lik2
            if (lik2 <= ceiling(Model_Lik*100000)/100000 & output$value > -Inf) {
              Model_Lik <- output$value; output$par -> model_pars[1:5] }
          }
        }
      }
      #Store best parameters for the site in a single table
      as.data.frame(matrix(ncol=6,nrow=1)) -> pars
      model_pars -> pars[1,1:5]; Model_Lik -> pars[1,6]
  } }
  
  return(pars)
}

#BM functions: these two functions are used for the likelihood of a BM model. 
BMlik <- function(delta, tau2, x0, info) {
  Edata <- matrix(x0,ncol=1, nrow=nindiv)
  #Scale the covariance matrix by the diffusion rate
  covar_matrix * delta -> VCV_species
  #Change from species-level covariance matrix to include individuals
  zmat %*% VCV_species %*% t(zmat) -> varcovar
  diag(varcovar) <- diag(varcovar) + tau2
  for (k in nindiv:1) { if ( is.na(info[k]) == T) {as.data.frame(matrix(Edata[-k,],ncol=1)) -> Edata; varcovar[-k,-k] -> varcovar; info[,-k] -> info }}
  #return the -loglikelihood
  return(-(dmvnorm(unlist(info),unlist(Edata),varcovar,log=TRUE)))
}
BM_all <- function(par, data1) {
  delta <- as.numeric(par[2])
  tau2 <- as.numeric(par[3])
  x0 <- as.numeric(par[1])
  #safeguard on delta and tau2 to prevent degeneracy in the covariance matrix
  if ( delta>0.01 & tau2>0) {
    BMlik(delta, tau2, x0, info=data1)
  }
  #return 10,000 if safeguards are not met
  else { 10000 } 
}

#OU functions: these two functions are used for the likelihood of a OU model for stabilizing selection. 
OUlik <- function(delta, alpha, tau2, theta, info) {
  Edata <- matrix(theta, ncol=1, nrow=nindiv)
  ##Restructure and scale the covariance matrix
  (exp(2*alpha*covar_matrix)-1) -> cvm2
  cvm2 <- cvm2/max(cvm2)
  cvm2 * delta -> VCV_species
  #Change from species-level covariance matrix to include individuals
  zmat %*% VCV_species %*% t(zmat) -> varcovar
  diag(varcovar) <- diag(varcovar) + tau2
  for (k in nindiv:1) { if ( is.na(info[k]) == T) {as.data.frame(matrix(Edata[-k,],ncol=1)) -> Edata; varcovar[-k,-k] -> varcovar; info[,-k] -> info }}
  #return the -loglikelihood
  return(-(dmvnorm(unlist(info),unlist(Edata),varcovar,log=TRUE)))
}
OU_all <- function(par, data1) {
  delta <- as.numeric(par[2])
  alpha <- as.numeric(par[3])
  theta <- as.numeric(par[1])
  tau2 <- as.numeric(par[4])
  #safeguard on delta, alpha, and tau2 to prevent degeneracy in the covariance matrix
  if ( delta>0.001 & tau2>0 & alpha>0.001 ) {
    OUlik(delta, alpha, tau2, theta, info=data1)
  }
  #return 10,000 if safeguards are not met
  else { 10000 } 
}

#OU2 functions: these two functions are used for the likelihood of an OU model. 
##OU2_all includes safeguards to keep parameters away from 0 where the covariance matrix becomes degenerate
##OU2lik just codes the likelihood for the model
OU2_all <- function(par, data1, lineage) {
  delta <- as.numeric(par[3])
  alpha <- as.numeric(par[4])
  theta1 <- as.numeric(par[1])
  theta2 <- as.numeric(par[2])
  tau2 <- as.numeric(par[5])
  #safeguard on theta values to prevent them from taking meaningless values far outside the observed data
  #safeguard on delta, tau2, and alpha to prevent degeneracy in the covariance matrix
  if ( delta > 0.001 & tau2>0 & alpha>0.001 & min(c(theta1,theta2)) >= min(data1, na.rm=T) & max(c(theta1,theta2)) <= max(data1, na.rm=T)) {
    OU2lik(delta, alpha, tau2, theta1,theta2,lineage, info=data1)
  }
  #return 10,000 if safeguards are not met
  else { 10000 } 
}
OU2lik <- function(delta, alpha, tau2, theta1, theta2, lin, info) {
  #paint entire tree with theta1, then apply theta2 to the relevant lineage which we're testing selection upon
  Edata <- matrix(theta1, ncol=1, nrow=nindiv)
  if (lin == "north") {Edata[6:34,1] <- theta1*exp(-alpha*(1-0.4955752))+theta2*(1-exp(-alpha*0.4955752))}
  if (lin == "south") {Edata[35:44,1] <- theta1*exp(-alpha*(1-0.4955752))+theta2*(1-exp(-alpha*0.4955752))}
  if (lin == "anugui") {Edata[6:20,1] <- theta1*exp(-alpha*(1-0.5663717))+theta2*(1-exp(-alpha*0.5663717))}
  if (lin == "anu") {Edata[6:14,1] <- theta1*exp(-alpha*(1-0.6371681))+theta2*(1-exp(-alpha*0.6371681))}
  if (lin == "gui") {Edata[15:20,1] <- theta1*exp(-alpha*(1-0.6371681))+theta2*(1-exp(-alpha*0.6371681))}
  if (lin == "ham") {Edata[21:34,1] <- theta1*exp(-alpha*(1-0.5663717))+theta2*(1-exp(-alpha*0.5663717))}
  if (lin == "yel") {Edata[39:44,1] <- theta1*exp(-alpha*(1-0.5752212))+theta2*(1-exp(-alpha*0.5752212))}
  if (lin == "kin") {Edata[35:38,1] <- theta1*exp(-alpha*(1-0.5752212))+theta2*(1-exp(-alpha*0.5752212))}
  #0.4955752 is the shared distance between the N/S split
  
  #Adjust the covariance matrix to reflect selection
  (exp(2*alpha*covar_matrix)-1) -> cvm2
  cvm2 <- cvm2/max(cvm2)
  cvm2 * delta -> VCV_species
  #Change from species-level covariance matrix to include individuals
  zmat %*% VCV_species %*% t(zmat) -> varcovar
  diag(varcovar) <- diag(varcovar) + tau2
  for (k in nindiv:1) { if ( is.na(info[k]) == T) {as.data.frame(matrix(Edata[-k,],ncol=1)) -> Edata; varcovar[-k,-k] -> varcovar; info[,-k] -> info }}
  #return the -loglikelihood
  return(-(dmvnorm(unlist(info),unlist(Edata),varcovar,log=TRUE)))
}

save.image('~/Functions_Mar2018.RData')