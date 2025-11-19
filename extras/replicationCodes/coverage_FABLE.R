#Rcpp::sourceCpp("src/NEW-FABLE-functions.cpp")
Rcpp::sourceCpp("src/updated-FABLE-functions.cpp")

n = 1000
p = 1000
pi0 = 0.5
alpha = 0.05
lambdasd = 0.5
#relevantIndices = c(1:100)


dir.name = "." #set directory to save
dir.create(dir.name)

set.seed(1) #set the seed here

relevantIndices = sample(1:p, size = 100, replace = FALSE)
pSub = length(relevantIndices)

k = 10

R = 100
covStor = matrix(0, nrow = R, ncol = pSub * (pSub+1)/2)
widthStor = rep(0, R)

Lambda = matrix(rnorm(p*k, mean = 0, sd = lambdasd), nrow = p, ncol = k)
BinMat = matrix(rbinom(p*k, 1, 1-pi0), nrow = p, ncol = k) #pi0 = P(zero)
Lambda = Lambda * BinMat

Sigma0 = runif(p, 0.5, 5)

gamma0 = 1
delta0sq = 1
MC = 1000
Psi0 = Matrix::tcrossprod(Lambda) + diag(Sigma0)

r = 1

for(r in 1:R) {
  
  print(paste0("Replicate: ", r))
  
  set.seed(2001+r)
  
  M = matrix(rnorm(n*k), nrow = n, ncol = k)
  E = matrix(rnorm(n*p), nrow = n, ncol = p)
  E = sweep(E, 2, sqrt(Sigma0), "*")
  
  Y = (M %*% t(Lambda)) + E
  
  #Y = sweep(Y, 2, colMeans(Y_uncen))
  
  svdmod = svd(Y)
  U_Y = svdmod$u
  V_Y = svdmod$v
  svalsY = svdmod$d
  
  ###### TESTING POSTERIOR SAMPLING CODE #########
  
  ##### GET MATRIX OF COVERAGE CORRECTION COEFFICIENTS #####
  
  kEst = CPPRankEstimator(Y, U_Y, V_Y, svalsY, 50)
  
  FABLEHypPars = FABLEHyperParameters(Y, U_Y, V_Y, svalsY, kEst)
  covCorrectEntries = CPPcov_correct_matrix(FABLEHypPars$SigmaSqEstimate,
                                            FABLEHypPars$G)
  
  # varInfSolveFunction <- function(eta) {
  # 
  #   zAlpha = qnorm(1 - (alpha/2))
  #   
  #   diagIndices = sapply(1:p, function(x) {return(x*(x+1)/2)})
  #   
  #   LHSTerm1 = sum((2 * pnorm(eta * zAlpha / covCorrectEntries[-diagIndices])) - 1)
  #   coefLHSTerm2 = sqrt((1 + (4 * eta^2 * (covCorrectEntries[diagIndices]^2 - 1))) /
  #                         ((2 * covCorrectEntries[diagIndices]^2) - 1)^2)
  #   LHSTerm2 = sum((2 * pnorm(zAlpha * coefLHSTerm2)) - 1)
  #   
  #   LHS = (LHSTerm1 + LHSTerm2) / (p*(p+1)/2)
  #   RHS = 1 - alpha
  #   output = LHS - RHS
  # 
  #   return(output)
  # 
  # }
  
  #NonLinEqSolv = nleqslv::nleqslv(x = mean(covCorrectEntries), varInfSolveFunction)
  #varInflation = quantile(covCorrectEntries, 0.50)^2
  #varInflation = median(covCorrectEntries)
  #varInflation = NonLinEqSolv$x
  
  ## Fit model
  
  varInflation = mean(covCorrectEntries)^2
  
  CPPSamplingOutput = CPPFABLESampler(Y, 
                                      gamma0, 
                                      delta0sq, 
                                      MC,
                                      U_Y,
                                      V_Y,
                                      svalsY,
                                      kEst,
                                      varInflation)
  
  ##Check coverage of 100*100 submatrix
  
  t11 = proc.time()
  
  CPPPostProcess = CCFABLEPostProcessingSubmatrix(CPPSamplingOutput,
                                                  alpha,
                                                  relevantIndices)
  
  # CPPPostProcess = CCFABLEPostProcessingSubmatrix_Optimized(CPPSamplingOutput,
  #                                                 alpha,
  #                                                 relevantIndices)
  
  t12 = proc.time()
  
  print((t12-t11))
  
  truePsi0 = Psi0[relevantIndices, relevantIndices]
  lowPsi0 = CPPPostProcess$LowerQuantileMatrix
  highPsi0 = CPPPostProcess$UpperQuantileMatrix
  
  truePsi0Vec = truePsi0[upper.tri(truePsi0, diag = TRUE)]
  lowPsi0Vec = lowPsi0[upper.tri(lowPsi0, diag = TRUE)]
  highPsi0Vec = highPsi0[upper.tri(highPsi0, diag = TRUE)]
  
  covStor[r,] = as.numeric((lowPsi0Vec <= truePsi0Vec) & (truePsi0Vec <= highPsi0Vec))
  widthStor[r] = mean(highPsi0Vec - lowPsi0Vec)
  
  print(c(mean(covStor[r,]), mean(widthStor[r])))
  
  ## Save information
  
  write.csv(covStor[r,], file = paste0(dir.name, "/", "coverage_rep=", r,  "_n=", n, "_p=", p, "_lambdasd=", lambdasd, "_k=", k, "_pi0=", pi0, ".csv"))
  write.csv(widthStor[r], file = paste0(dir.name, "/", "width_rep=", r,  "_n=", n, "_p=", p, "_lambdasd=", lambdasd, "_k=", k, "_pi0=", pi0, ".csv"))
  write.csv(mean(covStor[r,]), file = paste0(dir.name, "/", "avg_coverage_rep=", r,  "_n=", n, "_p=", p, "_lambdasd=", lambdasd, "_k=", k, "_pi0=", pi0, ".csv"))
  
}

print(c(mean(covStor), mean(widthStor)))