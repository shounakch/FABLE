
#### Train-test split ####

library(tidyverse)
library(dplyr)
library(infinitefactor)

## Specify the rawData ##

source("extras/GeneDataApplication/cleanData.R")
source("extras/AutomaticSparsity.R")
source("R/FABLE_code.R")

file.dir = "." #provide directory where raw data is saved

rawData <- readr::read_csv(file.dir)

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

#top 4% variables by variance to improve speed
Y = cleanData(rawData, needs.processing = TRUE, var.prop = 0.04) 

N = nrow(Y) # = 205
p = ncol(Y) # = 2120
print(c(N,p))

root.dir = "~/Research/FABLE/FABLE-GitHub/extras/GeneDataApplication/GSE109125_output_updated/trainTest/"

trainSeq = c(110, 130, 150, 170, 190)
write.csv(trainSeq, file = paste0(root.dir, "trainSamples_seq.csv"))

## Begin train/test split exploration ##

R = 10
r = 1

for(r in 1:R)
{
  
  #### Vary training set
  
  set.seed(2001+r)
  
  for(trainInd in 1:length(trainSeq)) {
    
    # set.seed(2001+r)
    
    nTrain = trainSeq[trainInd]
    nTest = N - nTrain
    
    train_ind = sample(1:N, size = nTrain, replace = FALSE)
    test_ind = c(1:N)[-train_ind]
    
    YTrain = Y[train_ind,]
    YTest = Y[test_ind,]
    
    # train each model
    
    t1 = proc.time()
    
    FABLEFit = FABLEPostmean(as.matrix(YTrain),
                             gamma0 = 1,
                             delta0sq = 1)
    
    # t3 = proc.time()
    
    MGSPFit = infinitefactor::linearMGSP(as.matrix(YTrain),
                                         nrun = 3000,
                                         burn = 1000,
                                         output = "covMean")
    
    # t4 = proc.time()
    
    ROTATEFit_lambda01 = autoSparse(as.matrix(YTrain), 1.0)
    ROTATEFit_lambda05 = autoSparse(as.matrix(YTrain), 5.0)
    ROTATEFit_lambda010 = autoSparse(as.matrix(YTrain), 10.0)
    
    # t5 = proc.time()
    
    # evaluate out of sample log likelihoods
    
    FABLEOOSLL = LL_gauss(YTest, FABLEFit)
    MGSPOOSLL = LL_gauss(YTest, MGSPFit$covMean)
    ROTATEOOSLL_lambda01 = LL_gauss(YTest, ROTATEFit_lambda01)
    ROTATEOOSLL_lambda05 = LL_gauss(YTest, ROTATEFit_lambda05)
    ROTATEOOSLL_lambda010 = LL_gauss(YTest, ROTATEFit_lambda010)
    
    t2 = proc.time()
    
    # save the values
    
    allOOSLL = c(FABLEOOSLL,
                 MGSPOOSLL,
                 ROTATEOOSLL_lambda01,
                 ROTATEOOSLL_lambda05,
                 ROTATEOOSLL_lambda010) / 10^3
    
    allOOSLL = as.matrix(t(allOOSLL))
    
    colnames(allOOSLL) = c("FABLE",
                           "MGSP",
                           paste0("ROTATE_lambda0=", 1), 
                           paste0("ROTATE_lambda0=", 5),
                           paste0("ROTATE_lambda0=", 10))
    
    write.csv(allOOSLL, 
              file = paste0(root.dir, "LogLik_trainInd=", trainInd, "_r=", r, ".csv"))
    
    # print("Time taken = ", as.numeric((t2-t1)[3]))
    
  }
  
  print(paste("Replicate: ", r, sep = ""))
  
}





