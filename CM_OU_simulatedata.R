#Simulate data for one phenotype based on set parameter values
library(ape); library(phytools); library(MASS); library(scales)

##Modify Block 7-14 & line 77

#set paramaters we simulate based off of
{alpha=2 #the strength of selection, 0 represents brownian motion
delta=2 #the diffusion factor or rate of neutral change
#delta = sigma2/2alpha for OU models; sigma2 for BM models
#note that for these factors, we assume the mutational input follows a normal distribution
tau2=.3 #the within-species variance, constrained to be constant across all species
diff=0 #the difference between the branch with selection and branch without selection. 0 for stabilizing selection, not relevant for Brownian motion
#we treat the ancestral state/theta as 0 -- the value of this parameter doesn't affect parameter fitting, just shifts everything up or down
}

#read in the phylogenetic tree, convert to covariance matrix, and set up the study design (individuals per species)
{read.tree("./ExNewick.txt") -> tree
vcv.phylo(tree) -> covar_matrix
covar_matrix[c(1,4,3,2,5,6),c(1,4,3,2,5,6)] -> covar_matrix #I needed to reorder the matrix to match my other files
covar_matrix <- covar_matrix/max(covar_matrix) #just scales the matrix, this doesn't matter and just relates to delta
rm(tree)}

###Get number of individuals per species & produce the z-matrix (which we will later multiply by the covariance matrix)
{read.delim('./Ex_nindiv.txt', header=F, sep=' ') -> indiv_spec
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
rm(k); rm(maxi); rm(mini)}

#Begin procedure to simulate data
##create table to put data into
simdata <- as.data.frame(matrix(ncol=max(indiv_spec), nrow=1))
sim_pars <- as.data.frame(matrix(ncol=6, nrow=1))
##Set optimum value
{x0 <- 0 - diff/2
x1 <- x0 + diff}
##If genetic drift (alpha=0)
if (alpha==0) {
  #draw species means from a multivariate normal
  species_means <- mvrnorm(mu=rep(x0, nspec), Sigma = delta * covar_matrix)
  #for each species, draw individuals for that species from a normal distribution, and put them into the simdata table
  for (k in 1:nspec) {
    if (k == 1) {
      N <- as.numeric(indiv_spec[k])
      simdata[1,1:N] <- rnorm(n=N, mean=species_means[k], sd=sqrt(tau2))
      rm(N)
    }
    if (k > 1) {
      N <- as.numeric(indiv_spec[k]-indiv_spec[k-1])
      n1 <- as.numeric(indiv_spec[k]); n0 <- as.numeric(indiv_spec[k-1])+1
      simdata[1,n0:n1] <- rnorm(n = N, mean = species_means[k], sd = sqrt(tau2))
      rm(N); rm(n1); rm(n0)
    }}
}
##If there is selection (alpha>0)
if (alpha > 0) {
  ##Restructure covariance matrix to reflect the effects of selection
  (exp(2*alpha*covar_matrix)-1) -> cvm2;   cvm2 <- cvm2/max(cvm2);   cvm2 * delta -> spec_covar
  #define means based upon theta, accounting for a shift between lineages
  x1 <- x0 + diff; means <- rep(x0,nrow(covar_matrix)); 
#This is one of the few lines that needs to be individually adjusted for which species selection has acted on. 
  means[c(2,3,4)] <- x1
  #draw species means from a multivariate normal
  species_means <- mvrnorm(mu=means, Sigma = spec_covar)
  #for each species, draw individuals for that species from a normal distribution, and put them into the simdata table
  for (k in 1:nspec) {
    if (k == 1) {
      N <- as.numeric(indiv_spec[k])
      simdata[1,1:N] <- rnorm(n=N, mean=species_means[k], sd=sqrt(tau2)); rm(N)
    }
    if (k > 1) {
      N <- as.numeric(indiv_spec[k]-indiv_spec[k-1])
      n1 <- as.numeric(indiv_spec[k]); n0 <- as.numeric(indiv_spec[k-1])+1
      simdata[1,n0:n1] <- rnorm(n = N, mean = species_means[k], sd = sqrt(tau2))
      rm(N); rm(n1); rm(n0)
    }
  }
  ##Alternate formulation that does the same thing (is equivalent), but creates different intermediaries
  #exp(-2*alpha*(1-covar_matrix)) -> cvm3; cvm3 <- (cvm3-min(cvm3)); cvm3 <- cvm3/max(cvm3)
}
##Save parameters to accompany the data
sim_pars[1,] <- c(x0,x1,delta,alpha,tau2,diff)


##To simulate multiple sites, put this withing a for-loop shell
###Create simdata and sim_pars files with the number of rows for the number of sites. 
###simdata <- as.data.frame(matrix(ncol=max(indiv_spec), nrow=100))
###sim_pars <- as.data.frame(matrix(ncol=6, nrow=100))
###for (i in 1:100) {save ouptuts in lines 61, 67, 85, 90, and 98 to i rather than 1}
##can simulate sets over multiple parameter values by drawing or setting parameter within the for loop as well
##output file
