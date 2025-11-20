#### Demonstrate blessing of dimensionality via estimation of
#### covariance of fixed set of variables. 
#### Incorporate more variables.

library(tidyverse)
library(dplyr)
library(infinitefactor)

## Specify the rawData ##

source("extras/cleanData.R")
source("extras/AutomaticSparsity.R")
source("R/FABLE_code.R")

file.dir = "." #provide directory where raw data is saved

rawData = readr::read_csv(file.dir)

#### Log likelihood function 

LL_gauss<-function(Y, Psi) {
  
  Y = as.matrix(Y)
  N = nrow(Y)
  p = ncol(Y)
  
  eigen_Psi = eigen(Psi)
  
  #Psi_det = prod(eigen_Psi$values)
  Psi_inv = eigen_Psi$vectors %*% diag(1 / eigen_Psi$values) %*% t(eigen_Psi$vectors)
  
  term1 = (-0.5 * N) * sum(log(eigen_Psi$values))
  term2 = (-0.5) * sum(diag(Psi_inv %*% (t(Y) %*% Y)))
  
  return(term1 + term2)
  
}

Y = cleanData(rawData, needs.processing = TRUE, var.prop = 0.10)

N = nrow(Y)
p = ncol(Y)

## Investigate covariance of top nVars variables by M3Drop/gold standard network/variance.

## First rank the variances of "important" variables.

allColVars = apply(Y, 2, var)
allColVarsRank = data.table::frank(-allColVars)

nVars = 100 #number of relevant variables, of which submatrix required

# Obtain column position of each relevant gene in original Y matrix

# allVarIndices = rep(0, p)
# for(j in 1:p) {
#   
#   allVarIndices[j] = which(colnames(Y) == M3Drop_genes$Gene[j])
#   
# }

#relGeneNames = M3Drop_genes$Gene[1:nVars] #top nVars genes selected using effect size estimate
#relevant_ind = allVarIndices[1:nVars]
relevant_ind = which(allColVarsRank <= nVars)

test_length = 50
train_length = N - test_length

## Begin blessing of dimensionality exploration ##

root.dir = "~/Research/FABLE/FABLE-GitHub/extras/GeneDataApplication/GSE109125_output_updated/nVars=100/"

#p_train_seq = as.integer(exp(seq(log(100), log(3000), length.out = 10)))
p_train_seq = c(seq(100, 1000, by = 100), 2000, 4000)
write.csv(p_train_seq, file = paste0(root.dir, "p_train_seq.csv"))
#p_train_seq = c(ncol(Y) - nVars)

v1 = 1
R = 10

test_loglik_stor = matrix(0, nrow = R, ncol = length(p_train_seq))
test_baseline_stor = rep(0, R)

test_loglik_stor_ROTATE = matrix(0, nrow = R, ncol = length(p_train_seq))
test_baseline_stor_ROTATE = rep(0, R)

test_loglik_stor_MGSP = matrix(0, nrow = R, ncol = length(p_train_seq))
test_baseline_stor_MGSP = rep(0, R)

r = 1
lambda0 = 10

#print(paste0("Rerun Application with ROTATE lambda0 = ", lambda0))

for(r in 1:R)
{
  
  #### Vary training set and take submatrix of resulting estimate using FABLE.
  
  set.seed(2001+r)
  
  test_ind = sort(sample(1:N, size = test_length, replace = FALSE))
  train_ind = c(1:N)[-test_ind]
  
  Y_test = Y[test_ind, relevant_ind]
  
  ## First fit baseline / comparison model with only the p_cov = 100 variables.
  
  Y_train_baseline = Y[train_ind, relevant_ind]
  
  train_mod_baseline = FABLEPostmean(as.matrix(Y_train_baseline),
                                     gamma0 = 1,
                                     delta0sq = 1)
  train_mod_baseline_ROTATE = autoSparse(as.matrix(Y_train_baseline), lambda0)
  train_mod_baseline_MGSP = infinitefactor::linearMGSP(as.matrix(Y_train_baseline),
                                                       nrun = 3000,
                                                       burn = 1000,
                                                       output = "covMean")
  
  test_loglik_baseline = LL_gauss(Y_test, train_mod_baseline)
  test_loglik_baseline_ROTATE = LL_gauss(Y_test, train_mod_baseline_ROTATE)
  test_loglik_baseline_MGSP = LL_gauss(Y_test, train_mod_baseline_MGSP$covMean)
  
  #print(c(test_loglik_baseline, test_loglik_baseline_ROTATE, test_loglik_baseline_MGSP))
  
  test_baseline_stor[r] = test_loglik_baseline / 10^3
  test_baseline_stor_ROTATE[r] = test_loglik_baseline_ROTATE / 10^3
  test_baseline_stor_MGSP[r] = test_loglik_baseline_MGSP / 10^3
  
  #### SAVE
  
  testLogLikBaselineSaveObject = c(test_baseline_stor[r],
                                   test_baseline_stor_ROTATE[r],
                                   test_baseline_stor_MGSP[r])
  testLogLikBaselineSaveObject = as.matrix(t(testLogLikBaselineSaveObject))
  colnames(testLogLikBaselineSaveObject) = c("FABLE", "ROTATE", "MGSP")
  write.csv(testLogLikBaselineSaveObject, file = paste0(root.dir, "baselineLogLik_nVar=", nVars, "_r=", r, ".csv"))
  
  for(v1 in 1:length(p_train_seq))
  {
    
    p_train = p_train_seq[v1]
    
    ## Now increase number of variables
    
    #train_var_ind = allVarIndices[1:(nVars + p_train)]
    train_var_ind = c(relevant_ind, which((allColVarsRank > nVars) & (allColVarsRank <= nVars + p_train)))
    
    ## Train model on extra variables + original variables
    
    Y_train_extra_var = Y[train_ind, train_var_ind]
    
    cov_hat_extra_var = FABLEPostmean(as.matrix(Y_train_extra_var),
                                      gamma0 = 1,
                                      delta0sq = 1) 
    cov_hat_extra_var_ROTATE = autoSparse(as.matrix(Y_train_extra_var), 10)
    cov_hat_extra_var_MGSP = infinitefactor::linearMGSP(as.matrix(Y_train_extra_var),
                                                        nrun = 3000,
                                                        burn = 1000,
                                                        output = "covMean")
    
    ## Subset nVars * nVars submatrix.
    
    cov_hat_relevant_var = cov_hat_extra_var[1:nVars, 1:nVars]
    cov_hat_relevant_var_ROTATE = cov_hat_extra_var_ROTATE[1:nVars, 1:nVars]
    cov_hat_relevant_var_MGSP = cov_hat_extra_var_MGSP$covMean[1:nVars, 1:nVars]
    
    test_loglik_stor[r,v1] = LL_gauss(Y_test, cov_hat_relevant_var) / 10^3
    test_loglik_stor_ROTATE[r,v1] = LL_gauss(Y_test, cov_hat_relevant_var_ROTATE) / 10^3
    test_loglik_stor_MGSP[r,v1] = LL_gauss(Y_test, cov_hat_relevant_var_MGSP) / 10^3
    
    #print(c(test_loglik_stor[r,v1], test_loglik_stor_ROTATE[r,v1], test_loglik_stor_MGSP[r,v1]))
    
    #### SAVE
    
    testLogLikSaveObject = c(test_loglik_stor[r,v1],
                             test_loglik_stor_ROTATE[r,v1],
                             test_loglik_stor_MGSP[r,v1])
    testLogLikSaveObject = as.matrix(t(testLogLikSaveObject))
    colnames(testLogLikSaveObject) = c("FABLE", paste0("ROTATE_lambda0=", lambda0) , "MGSP")
    write.csv(testLogLikSaveObject, file = paste0(root.dir, "LogLik_extraVar=", p_train, "_nVar=", nVars, "_r=", r, ".csv"))
    
    print(paste("Dimension Index: ", v1, sep = ""))
    
  }
  
  print(paste("Replicate: ", r, sep = ""))
  
}






