#### OLD CODE ####

#Earlier version, much slower than new version using bisection search
# RankEstimatorOld <- function(Ymat, kMaxProportion, svdmod){
# 
#   #### Use the idea in `Determining the number of factors in high-dimensional generalized latent factor models` (Chen & Li, 2022)
#   #### kMaxProportion is the proportion of variation explained (using singular values).
#   #### svd of full matrix passed as svdmod.
#   #### Use bisection search to identify minimizer k.
# 
#   #### Ymat must be pre-processed to have zero row mean before.
# 
#   nwhole = nrow(Ymat)
#   pwhole = ncol(Ymat)
# 
#   #### Obtain SVD of entire matrix
# 
#   #svdY = svd(Ymat)
#   svdY = svdmod
#   svalsY = svdY$d
#   kMax = min(which(cumsum(svalsY) / sum(svalsY) >= kMaxProportion)) #upper bound on latent factor dimension
#   U_Y = svdY$u
#   V_Y = svdY$v
# 
#   ## Now obtain penalized criterion value for each k
# 
#   CriterionStor = rep(0, kMax)
#   for(kInd in 1:kMax){
# 
#     # Extract the n*k left singular matrix U
# 
#     U_K = as.matrix(U_Y[,1:kInd])
#     V_K = as.matrix(V_Y[,1:kInd])
# 
#     if(kInd == 1)
#     {
# 
#       D_K = as.matrix(svalsY[1])
# 
#     }else
#     {
# 
#       D_K = diag(svalsY[1:kInd])
# 
#     }
# 
#     # Estimate \sigma_j^2 and \tau^2
# 
#     #t1 = proc.time()
# 
#     # Most time taking step below!
#     # sigsq_hat_diag = colSums(((diag(nwhole) - (U_K %*% t(U_K))) %*% Ymat)^2) / nwhole
#     # tausq_est = mean( (colSums(((U_K %*% t(U_K)) %*% Ymat)^2) / nwhole) / (kInd * sigsq_hat_diag))
# 
#     #Alternative (MUCH FASTER - include in main code as well.)
# 
#     #t2 = proc.time()
# 
#     UDVt = U_K %*% D_K %*% t(V_K)
#     sigsq_hat_diag = colSums((Ymat-UDVt)^2) / nwhole
#     tausq_est = mean( (colSums(UDVt^2) / nwhole) / (kInd * sigsq_hat_diag))
# 
#     #t3 = proc.time()
# 
#     # Estimate factor loadings \Lambda
# 
#     LambdaEst = (sqrt(nwhole) / (nwhole + tausq_est)) * (t(Ymat) %*% U_K)
# 
#     #t4 = proc.time()
# 
#     # Estimate latent factor matrix M = \sqrt{n} * U_K
# 
#     FactorEst = sqrt(nwhole) * U_K
# 
#     # Obtain log-likelihood of whole data set
# 
#     LogLikDataSet = LogLikelihoodEvaluator(Ymat = Ymat, M = FactorEst, Lambda = LambdaEst, SigmaSq = sigsq_hat_diag)
# 
#     #t5 = proc.time()
# 
#     # Store this value
# 
#     CriterionStor[kInd] = (-2 * LogLikDataSet) +
#       (kInd * max(nwhole, pwhole) * log(min(nwhole, pwhole)) /
#          (nwhole * pwhole))
# 
#     # Print progress
# 
#     # if(kInd %% 10 == 0)
#     # {
#     #
#     #   print(kInd)
#     #
#     # }
# 
#     #print(kInd)
# 
#   }
# 
#   return(CriterionStor)
# 
# }

#old version - slow.
# cov_correct_matrix_old <- function(sigsq_hat, llprime_hat) {
# 
#   p = length(sigsq_hat)
# 
#   B = matrix(0, nrow = p, ncol = p)
# 
#   for(u in 1:p)
#   {
# 
#     for(v in 1:p)
#     {
# 
#       den_uv = (sigsq_hat[u] * llprime_hat[v,v]) + (sigsq_hat[v] * llprime_hat[u,u])
#       num_uv = (llprime_hat[u,u] * llprime_hat[v,v]) + (llprime_hat[u,v])^2
# 
#       if(den_uv == 0)
#       {
# 
#         B[u,v] = 1
# 
#       }else
#       {
# 
#         if(u == v)
#         {
# 
#           B[u,v] = sqrt(1 + (llprime_hat[u,u] / (2 * sigsq_hat[u])))
# 
#         }else
#         {
# 
#           B[u,v] = sqrt(1 + (num_uv / den_uv))
# 
#         }
# 
#       }
# 
#     }
# 
#   }
# 
#   return(B)
# 
# }

#old code
# FABLE_vanilla<-function(Y, gamma0, delta0sq, MC, kMaxProportion) {
#   
#   n = dim(Y)[1]
#   p = dim(Y)[2]
#   
#   ## Perform SVD of Y ##
#   
#   svdmod = svd(Y)
#   
#   ####### CHOOSE k ##########
#   
#   #CritValues = RankEstimator(Y, svdmod)
#   k = RankEstimator(Y, svdmod)
#   
#   # normsvec = svdmod$d
#   # 
#   # min_fn<-function(j)
#   # {
#   #   
#   #   return(normsvec[j] - normsvec[j+1]) 
#   #   
#   # }
#   # 
#   # min_fn_vals = rep(0, n-1)
#   # 
#   # for(j in 1:(n-1))
#   # {
#   #   
#   #   min_fn_vals[j] = min_fn(j)
#   #   
#   # }
#   # 
#   # k = which.max(min_fn_vals)
#   
#   # print(paste("The chosen value of k: ", k, sep = ""))
#   # print(paste("The estimated SNR is: ", max(min_fn_vals) / normsvec[1], sep = ""))
#   
#   ## Obtain spectral estimate of latent factors from SVD ##
#   
#   U = svdmod$u[,1:k]
#   Mhat = sqrt(n) * U
#   
#   ############## FABLE SAMPLING ALGORITHM ####################
#   
#   ####Drawing samples for factor loadings
#   
#   p_new = p
#   
#   MCMC = MC
#   
#   lambda_stor_MCMC = array(0, dim = c(MCMC, k, p_new))
#   sigsq_stor_MCMC = matrix(0, nrow = MCMC, ncol = p_new)
#   
#   ## Choosing \tau^2 by Empirical Bayes ##
#   
#   llphat_mat = t(U) %*% Y
#   
#   llprime_hat = (t(llphat_mat) %*% llphat_mat) / n
#   
#   sigsq_hat_diag = colSums(((diag(n) - (U %*% t(U))) %*% Y)^2) / n
#   tausq_est = mean( (colSums(((U %*% t(U)) %*% Y)^2) / n) / (k * sigsq_hat_diag))
#   
#   tausq = rep(tausq_est, p)
#   
#   #### BEGIN PSEUDO-POSTERIOR SAMPLING ####
#   
#   for(j in 1:p_new)
#   {
#     
#     #### Sample (\lambda_j, \sigma_j^2) by our method using NIG prior ####
#     
#     L_n = diag(k) / (n + (1/tausq[j]))
#     A_n = (n + (1/tausq[j])) * diag(k)
#     mu_n = (t(Mhat) %*% Y[,j]) / (n + (1/tausq[j]))
#     
#     sig_rate = as.numeric(t(mu_n) %*% A_n %*% mu_n)
#     
#     for(u in 1:MCMC)
#     {
#       
#       sigsq_samp = 1/rgamma(n = 1, shape = (gamma0/2) + (n/2),
#                             rate = (gamma0*delta0sq/2) +
#                               0.5*(sum(Y[,j]^2) - sig_rate))
#       
#       lambda_stor_MCMC[u,,j] = rnorm(n = k, 
#                                      mean = mu_n,
#                                      sd = sqrt(sigsq_samp) * (n + (1/tausq[j]))^(-1/2))
#       
#       sigsq_stor_MCMC[u,j] = sigsq_samp
#       
#     }
#     
#   }
#   
#   # ##Get Lambda Lambda' samples
#   # 
#   # llprime_MCMC = array(0, dim = c(p_new, p_new, MCMC))
#   # 
#   # for(m in 1:MCMC)
#   # {
#   #   
#   #   llprime_MCMC[,,m] = t(lambda_stor_MCMC[m,,]) %*% lambda_stor_MCMC[m,,]
#   #   
#   # }
#   # 
#   # mean_llprime_MCMC = apply(llprime_MCMC, c(1,2), mean)
#   # 
#   # #### Samples for whole matrix \Psi = \Lambda\Lambda' + \Sigma
#   # 
#   # Psi_MCMC = array(0, dim = c(p_new, p_new, MCMC))
#   # 
#   # for(m in 1:MCMC)
#   # {
#   #   
#   #   Psi_MCMC[,,m] = llprime_MCMC[,,m] + diag(sigsq_stor_MCMC[m,1:p_new])
#   #   
#   # }
#   # 
#   # Psi_mean_MCMC = apply(Psi_MCMC, c(1,2), mean)
#   
#   hyperpar = list("n" = n,
#                   "p" = p,
#                   "MC" = MC,
#                   "Y" = Y,
#                   "SVD_Y" = svdmod,
#                   "tausq_est" = tausq_est,
#                   "k_est" = k,
#                   "llprime_hat" = llprime_hat,
#                   "sigsq_hat" = sigsq_hat_diag)
#   
#   output = list("Lambda_samples" = lambda_stor_MCMC,
#                 "Sigmasq_samples" = sigsq_stor_MCMC,
#                 "hyperpar" = hyperpar)
#   
#   return(output)
#   
# }
# 
# FABLE_corrected<-function(FABLEmod) {
#   
#   Y = FABLEmod$hyperpar$Y
#   n = dim(Y)[1]
#   p = dim(Y)[2]
#   MC = FABLEmod$hyperpar$MC
#   
#   llprime_hat = FABLEmod$hyperpar$llprime_hat
#   sigsq_hat = FABLEmod$hyperpar$sigsq_hat
#   
#   cov_correct_mat = cov_correct_matrix(sigsq_hat, llprime_hat)
#   
#   ##### Obtain coverage-corrected samples #####
#   
#   cov_correct_llprime_samples = array(0, dim = c(p, p, MC))
#   cov_correct_sigmasq_samples = matrix(0, nrow = MC, ncol = p)
#   Lambda_samples = FABLEmod$Lambda_samples
#   Sigmasq_samples = FABLEmod$Sigmasq_samples
#   
#   for(m in 1:MC)
#   {
#     
#     cov_correct_llprime_samples[,,m] = llprime_hat + 
#       (cov_correct_mat * (Lambda_samples[,,m] - llprime_hat))
#     
#   }
#   
#   cov_correct_sigmasq_samples = Sigmasq_samples
#   
#   output = list("FABLEmod" = FABLEmod,
#                 "cov_correct_llprime_samples" = cov_correct_llprime_samples,
#                 "cov_correct_sigmasq_samples" = cov_correct_sigmasq_samples)
#   
#   return(output)
#   
# }

# FABLE_postmean <- function(Y, gamma0, delta0sq)
# {
#   
#   n = dim(Y)[1]
#   p = dim(Y)[2]
#   
#   svdmod = svd(Y)
#   
#   ####### CHOOSE k ##########
#   
#   k = RankEstimator(Y, svdmod)
#   
#   U = svdmod$u[,1:k]
#   V = svdmod$v[,1:k]
#   
#   if(k == 1)
#   {
#     
#     D = as.matrix(svdmod$d[1])
#     
#   }else
#   {
#     
#     D = diag(svdmod$d[1:kInd])
#     
#   }
#   
#   ##### CHOOSE \tau^2 #####
#   
#   UDVt = U %*% D %*% t(V)
#   sigsq_hat_diag = colSums((Y-UDVt)^2) / n
#   tausq_est = mean( (colSums(UDVt^2) / n) / (k * sigsq_hat_diag))
#   
#   ##### Obtain the hyperparameters #####
#   
#   YtU = as.matrix(sweep(V, 2, svdmod$d[1:k], "*"))
#   G0 = (sqrt(n) / (n + (1/tausq_est))) * YtU
#   G = G0 %*% t(G0)
#   
#   CCMatrix = cov_correct_matrix(sigsq_hat_diag, G)
#   
#   gamma_n = gamma0 + n                              #\gamma_n = \gamma_0 + n
#   
#   gamma_n_deltasq = rep(gamma0*delta0sq, p) + 
#     as.numeric(apply(Y^2, 2, sum)) -
#     ((n / (n + (1 / tausq_est))) * as.numeric(apply((t(YtU))^2, 2, sum)))
#   
#   LLtEstimate = G
#   SigmaEstimate = diag(gamma_n_deltasq / (gamma_n - 2))
#   
#   CovPostMean = LLtEstimate + SigmaEstimate
#   
#   return(CovPostMean)
#   
#   # ##### Estimate \Lambda0 \Lambda0' #####
#   # 
#   # UtY = t(U) %*% Y
#   # 
#   # FABLE_est_LLprime = (n / (n + (1/tausq_est))^2) * (t(UtY) %*% UtY)
#   # 
#   # ##### Estimate \Sigma_0 #####
#   # 
#   # gamma_n = gamma0 + n
#   # part1 = rep((gamma0 * delta0sq / gamma_n), p)
#   # part21 = colSums(Y^2) - ((n / (n + (1/tausq_est))) * colSums((UtY)^2))
#   # part2 = part21 / gamma_n
#   # 
#   # delta_vec = part1 + part2
#   # 
#   # FABLE_est_Sigma = (gamma_n / (gamma_n - 2)) * delta_vec
#   # 
#   # ##### Estimate \Psi_0 = \Lambda_0 \Lambda_0' + \Sigma #####
#   # 
#   # # FABLE_est_Psi = FABLE_est_LLprime + diag(FABLE_est_Sigma)
#   # 
#   # #### Return output
#   # 
#   # # output = list("PostMean" = FABLE_est_Psi,
#   # #               "tauSq" = tausq_est,
#   # #               "U" = U)
#   # #
#   # # return(output)
#   
# }