source("extras/AutomaticSparsity.R")
#Rcpp::sourceCpp("src/updated-FABLE-functions.cpp")
library(FABLE)
library(infinitefactor)
library(cvCovEst)

#n = no. of samples, p = no. of dimensions, R = no. of replicates
#function to replicate runtime results for given seedValue
#run for 20 seeds 2001 + r, 1 \leq r \leq 20.

runTimeResults<-function(n, p, seedValue) {
  
  set.seed(1)
  
  pi0 = 0.5
  k = 10
  lambdasd = 0.5
  
  Lambda0 = matrix(rnorm(p*k, 0, lambdasd), nrow = p, ncol = k)
  BinMat0 = matrix(rbinom(p*k, 1, 1 - pi0), nrow = p, ncol = k)
  Lambda0 = Lambda0 * BinMat0
  
  Sigma0 = runif(p, 0.5, 5)
  
  timeStor = matrix(0, nrow = 1, ncol = 7)
  timeStor = as.data.frame(timeStor)
  colnames(timeStor) = c("FABLESamples", 
                         "FABLEPostMean", 
                         "MGSP", 
                         "ROTATE", 
                         "HT", 
                         "SCAD", 
                         "LW")
  
  set.seed(seedValue)
  
  E = matrix(rnorm(n*p), nrow = n, ncol = p)
  E = sweep(E, 2, sqrt(Sigma0), "*")
  M = matrix(rnorm(n*k), nrow = n, ncol = k)
  
  Y = (M %*% t(Lambda0)) + E
  
  #### FIT FABLE ####
  
  # ---------- SVD BEGINS ------------
  
  tSVD1 = proc.time()
  
  svdY = svd(Y)
  U_Y = svdY$u
  V_Y = svdY$v
  svalsY = svdY$d
  kMax = min(which(cumsum(svalsY) / sum(svalsY) >= 0.5))
  
  tSVD2 = proc.time()
  
  # --------- SVD ENDS --------------
  
  # --------- SAMPLING BEGINS ---------------
  
  tFABLESample1 = proc.time()
  
  kEst = CPPRankEstimator(Y,
                          U_Y,
                          V_Y,
                          svalsY,
                          kMax)
  
  FABLEHypPars = FABLEHyperParameters(Y,
                                      U_Y,
                                      V_Y,
                                      svalsY,
                                      kEst)
  
  CovCorrectMatrix = CPPcov_correct_matrix(FABLEHypPars$SigmaSqEstimate,
                                           FABLEHypPars$G)
  
  varInflation = (mean(CovCorrectMatrix))^2
  
  FABLESampler = CPPFABLESampler(Y,
                                 1.0,
                                 1.0,
                                 1000,
                                 U_Y,
                                 V_Y,
                                 svalsY,
                                 kEst,
                                 varInflation)
  
  tFABLESample2 = proc.time()
  
  # ------------- SAMPLING ENDS -----------------
  
  # ------------- PSEUDO-POSTERIOR MEAN -----------
  
  tFABLEPostMean1 = proc.time()
  
  FABLEPostMean = CPPFABLEPostMean(Y, 
                                   1.0, 
                                   1.0,
                                   U_Y,
                                   V_Y,
                                   svalsY,
                                   kMax)
  
  tFABLEPostMean2 = proc.time()
  
  # ------------ PSEUDO-POSTERIOR MEAN ENDS ----------
  
  # ----------- SAVE FABLE RESULTS ----------
  
  timeStor[1,1] = (tFABLESample2 - tFABLESample1)[3] + (tSVD2 - tSVD1)[3]
  timeStor[1,2] = (tFABLEPostMean2 - tFABLEPostMean1)[3] + (tSVD2 - tSVD1)[3]
  
  #### NOW MOVE ON TO OTHER APPROACHES #########
  
  # ----------- MGSP ------------------
  
  tMGSP1 = proc.time()
  
  MGSPSample = infinitefactor::linearMGSP(X = Y,
                                          nrun = 3000,
                                          burn = 1000,
                                          output = c("factSamples", "sigSamples"))
  
  tMGSP2 = proc.time()
  
  timeStor[1,3] = (tMGSP2 - tMGSP1)[3]
  
  # ------------ ROTATE ----------
  
  tROTATE1 = proc.time()
  
  ROTATEFit = autoSparse(Y, 10)
  
  tROTATE2 = proc.time()
  
  timeStor[1,4] = (tROTATE2 - tROTATE1)[3]
  
  # ------ Standardize Y for latter approaches ------
  
  YScaled = scale(Y, center = TRUE, scale = TRUE)
  
  vectorSds = apply(Y, 2, sd)
  scaleMatrix = vectorSds %*% t(vectorSds)
  
  # ------------ HT ----------
  
  tHT1 = proc.time()
  
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
  
  tHT2 = proc.time()
  
  timeStor[1,5] = (tHT2 - tHT1)[3]
  
  # ------------ SCAD ----------
  
  tSCAD1 = proc.time()
  
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
  
  tSCAD2 = proc.time()
  
  timeStor[1,6] = (tSCAD2 - tSCAD1)[3]
  
  # ------------ LW ----------
  
  tLW1 = proc.time()
  
  linearShrinkLWEstResults <- linearShrinkLWEst(YScaled)
  
  LSLWEResultsEstimate = linearShrinkLWEstResults * scaleMatrix
  
  tLW2 = proc.time()
  
  timeStor[1,7] = (tLW2 - tLW1)[3]
  
  return(timeStor)
  
}








