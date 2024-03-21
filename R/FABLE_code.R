##Symmetrizes an upper triangular matrix.
SymmetrizeMatrix <- function(A) {
  
  A1 = A + t(A)
  diag(A1) = diag(A1) / 2
  
  return(A1)
  
}

## Evaluates the likelihood function for y_{ij} ~ N(\eta_i' \lambda_j, \sigma_j^2).
LogLikelihoodEvaluator <- function(Ymat, M, Lambda, SigmaSq)
{
  
  n = nrow(Ymat)
  
  SigmaVec = sqrt(SigmaSq)
  
  DiffMat = (Ymat - (M %*% t(Lambda)))
  DiffMat = sweep(DiffMat, 2, SigmaVec, "/")
  Term2 = -0.5 * sum(DiffMat^2)
  Term1 = -n * sum(log(SigmaVec))
  
  Term = Term1 + Term2
  
  return(Term/(nrow(Ymat) * ncol(Ymat)))
  
}

# incorporated bisection search in new RankEstimator function
# keep testing this function
RankEstimator <- function(Ymat, svdmod, kMax) {
  
  #### Use the idea in `Determining the number of factors in high-dimensional generalized latent factor models` (Chen & Li, 2022)
  #### kMaxProportion is the proportion of variation explained (using singular values).
  #### svd of full matrix passed as svdmod.
  #### Use bisection search to identify minimizer k.
  
  #### Ymat must be pre-processed to have zero row mean before.
  
  nwhole = nrow(Ymat)
  pwhole = ncol(Ymat)
  
  #### Obtain SVD of entire matrix
  
  #svdY = svd(Ymat)
  svdY = svdmod
  svalsY = svdY$d
  #kMax = min(which(cumsum(svalsY) / sum(svalsY) >= 0.95)) #upper bound on latent factor dimension
  U_Y = svdY$u
  V_Y = svdY$v
  
  ## Define penalized criterion as a function of k
  
  CriterionFunction <- function(kInd) {
    
    # Extract the n*k left singular matrix U
    
    U_K = as.matrix(U_Y[,1:kInd])
    V_K = as.matrix(V_Y[,1:kInd])
    
    if(kInd == 1)
    {
      
      D_K = as.matrix(svalsY[1])
      
    }else
    {
      
      D_K = diag(svalsY[1:kInd])
      
    }
    
    # Estimate \sigma_j^2 and \tau^2
    
    #t1 = proc.time()
    
    # Most time taking step below!
    # sigsq_hat_diag = colSums(((diag(nwhole) - (U_K %*% t(U_K))) %*% Ymat)^2) / nwhole
    # tausq_est = mean( (colSums(((U_K %*% t(U_K)) %*% Ymat)^2) / nwhole) / (kInd * sigsq_hat_diag))
    
    #Alternative (MUCH FASTER - include in main code as well.)
    
    #t2 = proc.time()
    
    UDVt = U_K %*% D_K %*% t(V_K)
    sigsq_hat_diag = colSums((Ymat-UDVt)^2) / nwhole
    tausq_est = mean( (colSums(UDVt^2) / nwhole) / (kInd * sigsq_hat_diag))
    
    #t3 = proc.time()
    
    # Estimate factor loadings \Lambda
    
    YtU = as.matrix(sweep(V_K, 2, svalsY[1:kInd], "*"))
    LambdaEst = (sqrt(nwhole) / (nwhole + tausq_est)) * (YtU)
    
    #t4 = proc.time()
    
    # Estimate latent factor matrix M = \sqrt{n} * U_K
    
    FactorEst = sqrt(nwhole) * U_K
    
    # Obtain log-likelihood of whole data set
    
    LogLikDataSet = LogLikelihoodEvaluator(Ymat = Ymat, M = FactorEst, Lambda = LambdaEst, SigmaSq = sigsq_hat_diag)
    
    #t5 = proc.time()
    
    # Store this value
    # Loglikelihood value is already scaled by n*p
    
    CriterionStor = (-2 * LogLikDataSet) + 
      (kInd * max(nwhole, pwhole) * log(min(nwhole, pwhole)) / (nwhole * pwhole)) 
    
    return(as.numeric(CriterionStor))
    
  }
  
  BisectionRecursion <- function(lower, upper) {
    
    if(upper - lower > 5) {
      
      midpoint = as.integer(0.5 * (lower + upper))
      
      #print(midpoint)
      
      # print(c(lower, midpoint, upper))
      
      CriterionLower = CriterionFunction(lower)
      # CriterionUpper = CriterionFunction(upper)
      CriterionMidPoint = CriterionFunction(midpoint)
      
      if(CriterionLower < CriterionMidPoint) {
        
        #print(BisectionRecursion(lower, midpoint))
        return(BisectionRecursion(lower, midpoint))
        
      }else {
        
        ThreeFourthPoint = as.integer(0.5 * (midpoint + upper))
        CriterionThreeFourth = CriterionFunction(ThreeFourthPoint)
        
        #print(ThreeFourthPoint)
        
        if(CriterionMidPoint <= CriterionThreeFourth) {
          
          #print(BisectionRecursion(lower, midpoint))
          return(BisectionRecursion(lower, midpoint))
          
        }else {
          
          #print(BisectionRecursion(midpoint, upper))
          return(BisectionRecursion(midpoint, upper))
          
        }
        
      }
      
    }else {
      
      return(c(lower, upper))
      
    }
    
  }
  
  BisectionInterval = BisectionRecursion(1, kMax)
  #print(BisectionInterval)
  
  kEvalSeq = seq(BisectionInterval[1], BisectionInterval[2], by = 1)
  CriterionValues = rep(0, length(kEvalSeq))
  
  for(j in 1:length(kEvalSeq)) { 
    
    # Store this value
    
    CriterionValues[j] = CriterionFunction(kEvalSeq[j])
    
  }
  
  kMinimizer = kEvalSeq[which.min(CriterionValues)]
  
  return(kMinimizer)
  
}

#new version - should be faster. Test time and errors.
cov_correct_matrix <- function(sigsq_hat, llprime_hat) {
  
  p = length(sigsq_hat)
  diag_llprime = as.numeric(diag(llprime_hat))
  
  #B = matrix(0, nrow = p, ncol = p)
  
  num_matrix = (diag_llprime %*% t(diag_llprime)) + (llprime_hat^2)
  den_matrix = (diag_llprime %*% t(sigsq_hat)) + (sigsq_hat %*% t(diag_llprime))
  
  B = num_matrix / den_matrix
  diag(B) = diag(B) / 2.0
  
  B = sqrt(1.0 + B)
  
  return(B)
  
}

##Coverage-Corrected FABLE pseudo-posterior samples of the covariance matrix.
##Time-optimized code.
##Not space optimized - can be very bad choice for large p!
##Check out FABLE and PostProcessing functions below, more efficient.
CCFABLE_DirectSampler <- function(Y, gamma0 = 1, delta0sq = 1, MC = 1000) {
  
  n = dim(Y)[1]
  p = dim(Y)[2]
  
  svdmod = svd(Y)
  
  ####### CHOOSE k ##########
  
  kMax = min(which(cumsum(svdmod$d) / sum(svdmod$d) >= 0.95))
  
  k = RankEstimator(Y, svdmod, kMax)
  U = svdmod$u[,1:k]
  V = svdmod$v[,1:k]
  
  if(k == 1)
  {
    
    D = as.matrix(svdmod$d[1])
    
  }else
  {
    
    D = diag(svdmod$d[1:k])
    
  }
  
  ##### CHOOSE \tau^2 #####
  
  UDVt = U %*% D %*% t(V)
  sigsq_hat_diag = colSums((Y-UDVt)^2) / n
  tausq_est = mean( (colSums(UDVt^2) / n) / (k * sigsq_hat_diag))
  
  ##### Obtain the hyperparameters #####
  
  YtU = as.matrix(sweep(V, 2, svdmod$d[1:k], "*"))
  G0 = (sqrt(n) / (n + (1/tausq_est))) * YtU
  G = Matrix::tcrossprod(G0)
  
  CCMatrix = cov_correct_matrix(sigsq_hat_diag, G)
  
  gamma_n = gamma0 + n                              #\gamma_n = \gamma_0 + n
  
  gamma_n_deltasq = rep(gamma0*delta0sq, p) + 
    as.numeric(apply(Y^2, 2, sum)) -
    ((n / (n + (1 / tausq_est))) * as.numeric(apply((t(YtU))^2, 2, sum)))
  
  ##### Now proceed to store samples of \Psi = \Lambda \Lambda' + \Sigma #####
  ##### WARNING: Massive RAM consumption for large n, p ... #####
  
  CovSampleStor = array(0, dim = c(MC, p, p))
  a_n = 1 / sqrt(n + (1 / tausq_est))
  
  t1 = proc.time()
  
  for(m in 1:MC) {
    
    ## Sample \sigma_j^2 \sim IG(\gamma_n / 2, \gamma_n \delta_j^2 / 2), j=1,\ldots,p. ##
    
    sigmaSqSample = 1 / rgamma(n = p, shape = gamma_n / 2,
                               rate = gamma_n_deltasq / 2)
    
    SigmaSample = diag(sigmaSqSample)
    
    ## Sample low-rank part L_C = \Lambda \Lambda' ##
    
    ZSample = matrix(rnorm(p*k), nrow = p, ncol = k)
    SigmaHalfZ = sweep(ZSample, 1, sqrt(sigmaSqSample), "*")
    SigmaHalfZG0t = SigmaHalfZ %*% t(G0)
    
    Part2Matrix = (a_n * (SigmaHalfZG0t + t(SigmaHalfZG0t))) + 
      (a_n^2 * ((ZSample %*% t(ZSample)) * (sqrt(sigmaSqSample) %*% t(sqrt(sigmaSqSample)))))
    
    LLtSample = G + (CCMatrix * Part2Matrix)
    
    ## Store sample of \Psi ##
    
    PsiSample = LLtSample + SigmaSample
    CovSampleStor[m,,] = PsiSample
    
  }
  
  t2 = proc.time()
  
  return(CovSampleStor)
  
}

#function to provide pseudo-posterior mean without carrying out any sampling.
FABLEPostmean <- function(Y, gamma0 = 1, delta0sq = 1, kMax) {
  
  n = dim(Y)[1]
  p = dim(Y)[2]
  
  svdmod = svd(Y)
  
  # t1 = proc.time()
  
  ####### CHOOSE k ##########
  
  k = RankEstimator(Y, svdmod, kMax)
  
  U = svdmod$u[,1:k]
  V = svdmod$v[,1:k]
  
  if(k == 1)
  {
    
    D = as.matrix(svdmod$d[1])
    
  }else
  {
    
    D = diag(svdmod$d[1:k])
    
  }
  
  ##### CHOOSE \tau^2 #####
  
  UDVt = U %*% D %*% t(V) #consider using tcrossprod to make faster
  YtU = as.matrix(sweep(V, 2, svdmod$d[1:k], "*"))
  sigsq_hat_diag = colSums((Y-UDVt)^2) / n
  tausq_est = mean( (colSums((t(YtU))^2) / n) / (k * sigsq_hat_diag))
  
  ##### Obtain the hyperparameters #####
  
  #YtU = as.matrix(sweep(V, 2, svdmod$d[1:k], "*"))
  G0 = (sqrt(n) / (n + (1/tausq_est))) * YtU
  #G = G0 %*% t(G0)
  G = Matrix::tcrossprod(G0)
  
  # CCMatrix = cov_correct_matrix(sigsq_hat_diag, G)
  
  gamma_n = gamma0 + n                              #\gamma_n = \gamma_0 + n
  
  gamma_n_deltasq = rep(gamma0*delta0sq, p) + 
    as.numeric(apply(Y^2, 2, sum)) -
    ((n / (n + (1 / tausq_est))) * as.numeric(apply((t(YtU))^2, 2, sum)))
  
  LLtEstimate = G
  SigmaEstimate = diag(gamma_n_deltasq / (gamma_n - 2))
  
  CovPostMean = LLtEstimate + SigmaEstimate
  
  # t2 = proc.time()
  
  return(CovPostMean)
  
}

#### STUFF TO DO:
## 1. Write wrapper function for implementing FABLE. 
## 2. Write everything in RCpp + Armadillo: 
## 2.1. SVD - either exact or approximate 
## 2.2. Bisection search would be faster. (done)
## 2.3. MCMC sampling (done)

#### NEW STUFF TO DO (MARCH 5, 2024):

### Having FABLE sampler is important for space constraint reasons.
### Obtain \Lambda, \Sigma samples using regular FABLE.
### Then correct for coverage.
### Simple operations provide credible interval for \Psi_{uv}.

##Vanilla FABLE sampler.
##Useful for space constraint! Post-process to get CC-FABLE intervals.
##Write RCpp code equivalent.
FABLESampler <- function(Y, gamma0 = 1, delta0sq = 1, MC = 1000, kMax) {
  
  n = dim(Y)[1]
  p = dim(Y)[2]
  
  svdmod = svd(Y)
  
  ####### CHOOSE k ##########
  
  k = RankEstimator(Y, svdmod, kMax)
  U = svdmod$u[,1:k]
  V = svdmod$v[,1:k]
  
  if(k == 1)
  {
    
    D = as.matrix(svdmod$d[1])
    
  }else
  {
    
    D = diag(svdmod$d[1:k])
    
  }
  
  ##### CHOOSE \tau^2 #####
  
  UDVt = U %*% D %*% t(V)
  sigsq_hat_diag = colSums((Y-UDVt)^2) / n
  tausq_est = mean( (colSums(UDVt^2) / n) / (k * sigsq_hat_diag))
  
  ##### Obtain the hyperparameters #####
  
  YtU = as.matrix(sweep(V, 2, svdmod$d[1:k], "*"))
  G0 = (sqrt(n) / (n + (1/tausq_est))) * YtU
  G = G0 %*% t(G0)
  
  CCMatrix = cov_correct_matrix(sigsq_hat_diag, G)
  
  gamma_n = gamma0 + n                              #\gamma_n = \gamma_0 + n
  
  gamma_n_deltasq = rep(gamma0*delta0sq, p) + 
    as.numeric(apply(Y^2, 2, sum)) -
    ((n / (n + (1 / tausq_est))) * as.numeric(apply((t(YtU))^2, 2, sum)))
  
  ##### Now proceed to store samples of \Psi = \Lambda \Lambda' + \Sigma #####
  ##### WARNING: Massive RAM consumption for large n, p ... #####
  
  SigmaSampleStor = matrix(1, nrow = MC, ncol = p)
  LambdaSampleStor = array(0, dim = c(MC, p, k))
  
  a_n = 1 / sqrt(n + (1 / tausq_est))
  
  #t1 = proc.time()
  
  for(m in 1:MC) {
    
    ## Sample \sigma_j^2 \sim IG(\gamma_n / 2, \gamma_n \delta_j^2 / 2), j=1,\ldots,p. ##
    
    sigmaSqSample = 1 / rgamma(n = p, shape = gamma_n / 2,
                               rate = gamma_n_deltasq / 2)
    
    SigmaSampleStor[m,] = sigmaSqSample
    
    ## Sample \Lambda ##
    
    ZSample = matrix(rnorm(p*k), nrow = p, ncol = k)
    LambdaSampleStor[m,,] = G0 + (a_n * sweep(ZSample, 1, sqrt(sigmaSqSample), "*"))
    
    ## Sample low-rank part L_C = \Lambda \Lambda' ##
    
    # ZSample = matrix(rnorm(p*k), nrow = p, ncol = k)
    # SigmaHalfZ = sweep(ZSample, 1, sqrt(sigmaSqSample), "*")
    # SigmaHalfZG0t = SigmaHalfZ %*% t(G0)
    # 
    # Part2Matrix = (a_n * (SigmaHalfZG0t + t(SigmaHalfZG0t))) + 
    #   (a_n^2 * ((ZSample %*% t(ZSample)) * (sqrt(sigmaSqSample) %*% t(sqrt(sigmaSqSample)))))
    # 
    # LLtSample = G + (CCMatrix * Part2Matrix)
    # 
    # ## Store sample of \Psi ##
    # 
    # PsiSample = LLtSample + SigmaSample
    # CovSampleStor[m,,] = PsiSample
    
  }
  
  #t2 = proc.time()
  
  output = list("LambdaSamples" = LambdaSampleStor,
                "SigmaSqSamples" = SigmaSampleStor,
                "G" = G,
                "SigSqEstimate" = gamma_n_deltasq / (gamma_n - 2))
  
  return(output)
  
}

##Post-process FABLE draws to obtain entrywise pseudo-posterior
##mean, lower quantile, and upper quantiles of the covariance matrix.
##RCpp code equivalent will probably be much faster. (YES!)
CCFABLEPostProcessing <- function(FABLEOutput,
                                  CovCorrectMatrix,
                                  alpha = 0.05) {
  
  LambdaSamples = FABLEOutput$LambdaSamples
  SigmaSqSamples = FABLEOutput$SigmaSqSamples
  G = FABLEOutput$G
  
  nMC = nrow(SigmaSqSamples)
  p = ncol(SigmaSqSamples)
  
  CovMatPostMean = matrix(0, nrow = p, ncol = p)
  CovMatLower = matrix(0, nrow = p, ncol = p)
  CovMatUpper = matrix(0, nrow = p, ncol = p)
  
  ## Define the sample extractor function
  
  # SampleExtractor <- function(m, ind1, ind2) {
  #   
  #   lowRankTerm = sum(LambdaSamples[m,ind1,] * LambdaSamples[m,ind2,])
  #   return(as.numeric(lowRankTerm))
  #   
  # }
  
  for(j1 in 1:p) {
    
    t1 = proc.time()
    
    for(j2 in j1:p) {
      
      ## Extract the MC samples
      
      #t1 = proc.time()
      
      Lambj1Samples = t(LambdaSamples[,j1,])
      Lambj2Samples = t(LambdaSamples[,j2,])
      
      #t2 = proc.time()
      
      Lambj1Lambj2Samples = as.numeric(apply(Lambj1Samples * Lambj2Samples, 2, sum))
      
      #t3 = proc.time()
      
      Lambj1Lambj2SamplesCorrected = rep(G[j1,j2], nMC) + 
        (CovCorrectMatrix[j1,j2] * (Lambj1Lambj2Samples - rep(G[j1,j2], nMC)))
      
      #t4 = proc.time()
      
      coefj1j2 = ifelse(j1 == j2, 1, 0)
      
      Covj1j2Samples = Lambj1Lambj2SamplesCorrected + (coefj1j2 * SigmaSqSamples[,j1])
      
      #t5 = proc.time()
      
      ## Extract mean, low, high.
      
      # CovMatPostMean[j1,j2] = sum(Covj1j2Samples) / nMC
      # 
      # t51 = proc.time()
      # 
      # CovMatPostMean[j1,j2] = mean(Covj1j2Samples)
      # 
      # t511 = proc.time()
      # 
      # #print(c((t511 - t51)[3], (t51-t5)[3]))
      # 
      # CovMatLower[j1,j2] = quantile(Covj1j2Samples, alpha/2)
      # 
      # t52 = proc.time()
      # 
      # CovMatUpper[j1,j2] = quantile(Covj1j2Samples, 1 - (alpha/2))
      # 
      # t6 = proc.time()
      
    }
    
    t2 = proc.time()
    
    print(j1)
    
  }
  
  CovMatPostMean = SymmetrizeMatrix(CovMatPostMean)
  CovMatLower = SymmetrizeMatrix(CovMatLower)
  CovMatUpper = SymmetrizeMatrix(CovMatUpper)
  
  output = list("PostMeanMatrix" = CovMatPostMean,
                "LowerQuantileMatrix" = CovMatLower,
                "UpperQuantileMatrix" = CovMatUpper)
  
  return(output)
  
}
