# PapioMethylation

This repository contains code to simulate and analyze data using Ornstein-Uhlenbeck models for phenotypic evolution. 

1) CM_OU_simulatedata.R takes a phylogeny (Newick format: ExNewick.txt), number of individuals per species (Ex_nindiv.txt), and OU model parameters (lines 7-14, 77) to simulate the phenotype of a trait. It is easily modifiable (lines 101-107) to simulate multiple traits, potentially under multiple starting parameters. 

2) CM_OU_functions.R loads into R some basic information about the phylogeny and study design (again, based on ExNewick.txt and Ex_nindiv.txt), and defines functions which will be used to fit the OU models. These functions are saved in a Rdata file called OU_Functions.RData. 

3) Loading your data and OU_Functions.RData, CM_OU_fitmodels.R fits OU models to data which match your study design and phylogeny. You will need to specify which type(s) of OU models you are running (line 33), assign the appropriate column names to the output (line 40), and name the output file (line 46). 


These tools require the R packages ape, phytools, MASS, scales, and parallel. 