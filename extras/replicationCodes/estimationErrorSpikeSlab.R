##### Install FABLE package first #####

Rcpp::sourceCpp("src/updated-FABLE-functions.cpp") #replace by library(FABLE)
source("extras/AutomaticSparsity.R")
library(cvCovEst)

#### Below function provides estimation errors ####
#### Data generating scheme is spike-and-slab ####
#### n = 500, 1000; p = 1000, 5000; 
#### k = 10; R = 50; pi0 = 0.5, 0.85; lambdasd = 0.5
#### save.dir is the directory under which result should be saved

provideEstimationErrorsCase1 <- function(n, 
                                         p, 
                                         k, 
                                         repIndex, 
                                         pi0,
                                         lambdasd,
                                         save.dir = NA) {
  
  
  #### Generate true Lambda and Sigma ####
  
  set.seed(1)
  
  Lambda = matrix(rnorm(p*k, mean = 0, sd = lambdasd), nrow = p, ncol = k)
  BinMat = matrix(rbinom(p*k, 1, 1-pi0), nrow = p, ncol = k) #pi0 = P(zero)
  Lambda = Lambda * BinMat
  
  Sigma0 = runif(p, 0.5, 5)
  
  Psi0 = Matrix::tcrossprod(Lambda) + diag(Sigma0)
  normPsi0 = norm(Psi0, type = "2")
  
  #### MGSP and FABLE hyperparameters ####
  
  MGSPnrun = 3000
  MGSPburn = 1000
  
  gamma0 = 1
  delta0sq = 1
  
  #### Begin generating and fitting ####
  
  r = repIndex
  
  print(paste0("Replicate: ", r))
  
  set.seed(2001+r)
  
  M = matrix(rnorm(n*k), nrow = n, ncol = k)
  E = matrix(rnorm(n*p), nrow = n, ncol = p)
  E = sweep(E, 2, sqrt(Sigma0), "*")
  
  Y = (M %*% t(Lambda)) + E
  
  ## Fit model
  
  #### MGSP ####
  
  MGSPSamplingOutput = infinitefactor::linearMGSP(X = Y,
                                                  nrun = MGSPnrun,
                                                  burn = MGSPburn,
                                                  output = c("covMean", "sigSamples"))
  
  avgESS = mean(coda::effectiveSize(t(MGSPSamplingOutput$sigmaSamps)))
  write.csv(avgESS, file = paste0(dir.name, "ESS_MGSP_burn=", MGSPburn, "_nrun=", MGSPnrun, "_rep=", r, ".csv"))
  
  MGSPError = norm(MGSPSamplingOutput$covMean - Psi0, type = "2") / normPsi0
  print(paste0("MGSP Error = ", MGSPError))
  
  #### FABLE ####
  
  svdY = svd(Y)
  U_Y = svdY$u
  V_Y = svdY$v
  svalsY = svdY$d
  
  kMax = min(which(cumsum(svalsY) / sum(svalsY) >= 0.95))
  
  FABLEOutput = CPPFABLEPostMean(Y, gamma0, delta0sq, U_Y, V_Y, svalsY, kMax)
  
  FABLEError = norm(FABLEOutput$CovPostMean - Psi0, type = "2") / normPsi0
  print(paste0("FABLE Error = ", FABLEError))
  
  #### ROTATE WITH lambda0 as tuning parameter gradually increased
  
  ROTATEOutput = autoSparse_paper(Y, 20)
  ROTATEError = norm(ROTATEOutput - Psi0, type = "2") / normPsi0
  print(paste0("ROTATE Error = ", ROTATEError))
  
  #### Other competitors
  
  # thresholdingEst : hard thresholding of Bickel & Levina, 2008
  # scadEst : SCAD penalty (Fan and Li, 2001)
  # linearShrinkLWEst : Ledoit and Wolf (2004)
  # spcovEst : Bien and Tibshirani (2011)
  
  # need to scale data first
  
  YScaled = scale(Y, center = TRUE, scale = TRUE)
  
  vectorSds = apply(Y, 2, sd)
  scaleMatrix = vectorSds %*% t(vectorSds)
  
  # thresholdingEst
  
  thresholdingResults <- cvCovEst(
    dat = YScaled,
    estimators = c(thresholdingEst),
    estimator_params = list(
      thresholdingEst = list(gamma = seq(0.05, 0.5, length.out = 10))),
    cv_loss = cvMatrixFrobeniusLoss,
    cv_scheme = "v_fold",
    v_folds = 5
  )
  thresholdingResultsEstimate = thresholdingResults$estimate * scaleMatrix
  
  thresholdingError = norm(thresholdingResultsEstimate - Psi0, type = "2") / normPsi0
  print(paste0("Hard Thresholding Error = ", thresholdingError))
  
  # scadEst
  
  scadResults <- cvCovEst(
    dat = YScaled,
    estimators = c(scadEst),
    estimator_params = list(
      scadEst        = list(lambda = seq(0.01, 0.5, length.out = 10))),
    cv_loss = cvMatrixFrobeniusLoss,
    cv_scheme = "v_fold",
    v_folds = 5
  )
  scadResultsEstimate = scadResults$estimate * scaleMatrix
  
  scadError = norm(scadResultsEstimate - Psi0, type = "2") / normPsi0
  print(paste0("SCAD Error = ", scadError))
  
  # linearShrinkLWEst
  
  linearShrinkLWEstResults <- linearShrinkLWEst(YScaled)
  LSLWEResultsEstimate = linearShrinkLWEstResults * scaleMatrix
  
  linearShrinkLWEstError = norm(LSLWEResultsEstimate - Psi0, type = "2") / normPsi0
  print(paste0("LSLWE Error = ", linearShrinkLWEstError))
  
  #### Return output ####
  
  allErrors <- list("MGSP" = MGSPError,
                    "FABLE" = FABLEError,
                    "ROTATE" = ROTATEError,
                    "HT" = thresholdingError,
                    "SCAD" = scadError,
                    "LW" = linearShrinkLWEstError)
  
  allErrors <- matrix(unlist(allErrors), 
                      nrow = 1, 
                      dimnames = list(NULL, names(allErrors)))
  
  if(!is.na(save.dir)) {
    
    write.csv(allErrors, save.dir)
    
  }
  
  return(allErrors)
  
}

#### Run the function ####

R = 50
n = 500
p = 1000
k = 10
pi0 = 0.5
lambdasd = 0.5

for(r in 1:R) {
  
  provideEstimationErrors(n = n, 
                          p, 
                          k, 
                          repIndex = r, 
                          pi0,
                          lambdasd,
                          save.dir = NA)
  
}

