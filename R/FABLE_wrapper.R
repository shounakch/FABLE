FABLEPosteriorSampler <- function(Y,
                                  gamma0 = 1,
                                  delta0sq = 1,
                                  maxProp = 0.95,
                                  MC = 1000) {
  
  tFABLESample1 = proc.time()
  
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  svdY = svd(Y)
  U_Y = svdY$u
  V_Y = svdY$v
  svalsY = svdY$d
  kMax = min(which(cumsum(svalsY) / sum(svalsY) >= maxProp))
  
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
  
  CovCorrectMatrix = cov_correct_matrix(FABLEHypPars$SigmaSqEstimate, 
                                        FABLEHypPars$G)
  
  varInflation = (sum(CovCorrectMatrix) / (p*(p+1)/2))^2
  
  FABLESamples = CPPFABLESampler(Y, 
                                 gamma0, 
                                 delta0sq, 
                                 MC,
                                 U_Y,
                                 V_Y,
                                 svalsY,
                                 kEst,
                                 varInflation)
  
  tFABLESample2 = proc.time()
  tSample = (tFABLESample2 - tFABLESample1)[3]
  
  OutputList = list("CCFABLESamples" = FABLESamples,
                    "FABLEHyperParameters" = FABLEHypPars,
                    "svdY" = svdY,
                    "estRank" = kEst,
                    "varInflation" = varInflation,
                    "runTime" = tSample)
  
  return(OutputList)
  
}

FABLEPosteriorMean <- function(Y,
                               gamma0 = 1,
                               delta0sq = 1,
                               maxProp = 0.5) {
  
  tFABLEPostMean1 = proc.time()
  
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  svdY = svd(Y)
  U_Y = svdY$u
  V_Y = svdY$v
  svalsY = svdY$d
  kMax = min(which(cumsum(svalsY) / sum(svalsY) >= maxProp))
  
  # kEst = RankEstimator(Y, 
  #                      U_Y,
  #                      V_Y,
  #                      svalsY,
  #                      kMax)
  
  CovEstMod = CPPFABLEPostMean(Y, gamma0, delta0sq, U_Y, V_Y, svalsY, kMax)
  CovEst = CovEstMod$CovPostMean
  kEst = CovEstMod$kEst
  
  # FABLEHypPars = FABLEHyperParameters(Y,
  #                                     U_Y,
  #                                     V_Y,
  #                                     svalsY,
  #                                     kEst,
  #                                     gamma0,
  #                                     delta0sq)
  
  # Part1 = FABLEHypPars$G
  # Part2 = as.numeric(FABLEHypPars$gammaDeltasq / (FABLEHypPars$gamman - 2))
  # CovEst = Part1 + diag(Part2)
  
  tFABLEPostMean2 = proc.time()
  tPostMean = (tFABLEPostMean2 - tFABLEPostMean1)[3]
  
  OutputList = list("FABLEPostMean" = CovEst,
                    "svdY" = svdY,
                    "estRank" = kEst,
                    "runTime" = tPostMean)
  
  return(OutputList)
  
}