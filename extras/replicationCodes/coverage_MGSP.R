Rcpp::sourceCpp("extras/MGSP-functions.cpp")

n = 1000
p = 1000
pi0 = 0.5
alpha = 0.05
lambdasd = 0.5

dir.name = "." #set directory name here
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
Psi0 = Matrix::tcrossprod(Lambda) + diag(Sigma0)

r = 1

for(r in 1:R) {
  
  print(paste0("Coverage for n = ", n, ", p = ", p))
  print(paste0("Replicate: ", r))
  
  set.seed(2001+r)
  
  M = matrix(rnorm(n*k), nrow = n, ncol = k)
  E = matrix(rnorm(n*p), nrow = n, ncol = p)
  E = sweep(E, 2, sqrt(Sigma0), "*")
  
  Y = (M %*% t(Lambda)) + E
  
  ## Fit model
  
  #t1 = proc.time()
  MGSPSamplingOutput = infinitefactor::linearMGSP(X = Y,
                                                  nrun = 3000,
                                                  burn = 1000,
                                                  output = c("factSamples",
                                                             "sigSamples"))
  #t2 = proc.time()
  ##Check coverage of 100*100 submatrix
  
  #t11 = proc.time()
  
  VY = apply(Y, 2, var) #correction factor
  
  #t3 = proc.time()
  
  MGSPPostProcess = MGSPPostProcessingSubmatrix_Optimized(MGSPSamplingOutput$lambdaSamps,
                                                          MGSPSamplingOutput$sigmaSamps,
                                                          alpha,
                                                          relevantIndices,
                                                          VY)
  
  #t4 = proc.time()  
  
  #t12 = proc.time()
  
  truePsi0 = Psi0[relevantIndices, relevantIndices]
  lowPsi0 = MGSPPostProcess$LowerQuantileMatrix
  highPsi0 = MGSPPostProcess$UpperQuantileMatrix
  
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