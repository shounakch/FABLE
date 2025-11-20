source("extras/AutomaticSparsity.R")
#Rcpp::sourceCpp("src/updated-FABLE-functions.cpp")
library(FABLE)
library(infinitefactor)
library(cvCovEst)

#n = no. of samples, p = no. of dimensions, R = no. of replicates
#function to replicate runtime results
runTimeResults<-function(n, p, seedValue) {
  
  set.seed(seedValue)
  
  pi0 = 0.5
  k = 10
  lambdasd = 0.5
  
  Lambda0 = matrix(rnorm(p*k, 0, lambdasd), nrow = p, ncol = k)
  BinMat0 = matrix(rbinom(p*k, 1, 1 - pi0), nrow = p, ncol = k)
  Lambda0 = Lambda0 * BinMat0
  
  Sigma0 = runif(p, 0.5, 5)
  
  timeStor = matrix(0, nrow = 1, ncol = 8)
  timeStor = as.data.frame(timeStor)
  colnames(timeStor) = c("FABLESamples", 
                         "FABLEPostMean", 
                         "MGSP", 
                         "ROTATE", 
                         "HT", 
                         "POET", 
                         "SCAD", 
                         "LW")
  
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
                                   kMax,
                                   1)
  
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
  
  # ------------ POET ----------
  
  tPOET1 = proc.time()
  
  poetResults <- cvCovEst(
    dat = YScaled,
    estimators = c(poetEst),   
    estimator_params = list(
      poetEst        = list(k = 1:10,
                            lambda = seq(0.01, 0.5, length.out = 10))),
    cv_loss = cvMatrixFrobeniusLoss,
    cv_scheme = "v_fold",
    v_folds = 5
  )
  
  poetResultsEstimate = poetResults$estimate * scaleMatrix
  
  tPOET2 = proc.time()
  
  timeStor[1,6] = (tPOET2 - tPOET1)[3]
  
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
  
  timeStor[1,7] = (tSCAD2 - tSCAD1)[3]
  
  # ------------ LW ----------
  
  tLW1 = proc.time()
  
  linearShrinkLWEstResults <- linearShrinkLWEst(YScaled)
  
  LSLWEResultsEstimate = linearShrinkLWEstResults * scaleMatrix
  
  tLW2 = proc.time()
  
  timeStor[1,8] = (tLW2 - tLW1)[3]
  
  #for(r in 1:R) {
  
  # set.seed(2001 + r)
  # 
  # E = matrix(rnorm(n*p), nrow = n, ncol = p)
  # E = sweep(E, 2, sqrt(Sigma0), "*")
  # M = matrix(rnorm(n*k), nrow = n, ncol = k)
  # 
  # Y = (M %*% t(Lambda0)) + E
  # 
  # #### FIT FABLE ####
  # 
  # # ---------- SVD BEGINS ------------
  # 
  # tSVD1 = proc.time()
  # 
  # svdY = svd(Y)
  # U_Y = svdY$u
  # V_Y = svdY$v
  # svalsY = svdY$d
  # kMax = min(which(cumsum(svalsY) / sum(svalsY) >= 0.5))
  # 
  # tSVD2 = proc.time()
  # 
  # # --------- SVD ENDS --------------
  # 
  # # --------- SAMPLING BEGINS ---------------
  # 
  # tFABLESample1 = proc.time()
  # 
  # kEst = CPPRankEstimator(Y, 
  #                         U_Y,
  #                         V_Y,
  #                         svalsY,
  #                         kMax)
  # 
  # FABLEHypPars = FABLEHyperParameters(Y,
  #                                     U_Y,
  #                                     V_Y,
  #                                     svalsY,
  #                                     kEst)
  # 
  # CovCorrectMatrix = CPPcov_correct_matrix(FABLEHypPars$SigmaSqEstimate, 
  #                                          FABLEHypPars$G)
  # 
  # varInflation = (mean(CovCorrectMatrix))^2
  # 
  # FABLESampler = CPPFABLESampler(Y, 
  #                                1.0, 
  #                                1.0, 
  #                                1000,
  #                                U_Y,
  #                                V_Y,
  #                                svalsY,
  #                                kEst,
  #                                varInflation)
  # 
  # tFABLESample2 = proc.time()
  # 
  # # ------------- SAMPLING ENDS -----------------
  # 
  # # ------------- PSEUDO-POSTERIOR MEAN -----------
  # 
  # tFABLEPostMean1 = proc.time()
  # 
  # FABLEPostMean = CPPFABLEPostMean(Y, 
  #                                  1.0, 
  #                                  1.0,
  #                                  U_Y,
  #                                  V_Y,
  #                                  svalsY,
  #                                  kMax)
  # 
  # tFABLEPostMean2 = proc.time()
  # 
  # # ------------ PSEUDO-POSTERIOR MEAN ENDS ----------
  # 
  # # ----------- SAVE FABLE RESULTS ----------
  # 
  # timeStor[r,1] = (tFABLESample2 - tFABLESample1)[3] + (tSVD2 - tSVD1)[3]
  # timeStor[r,2] = (tFABLEPostMean2 - tFABLEPostMean1)[3] + (tSVD2 - tSVD1)[3]
  # 
  # #### NOW MOVE ON TO OTHER APPROACHES #########
  # 
  # # ----------- MGSP ------------------
  # 
  # tMGSP1 = proc.time()
  # 
  # MGSPSample = infinitefactor::linearMGSP(X = Y,
  #                                         nrun = 3000,
  #                                         burn = 1000,
  #                                         output = c("factSamples", "sigSamples"))
  # 
  # tMGSP2 = proc.time()
  # 
  # timeStor[r,3] = (tMGSP2 - tMGSP1)[3]
  # 
  # # ------------ ROTATE ----------
  # 
  # tROTATE1 = proc.time()
  # 
  # ROTATEFit = autoSparse(Y, 10)
  # 
  # tROTATE2 = proc.time()
  # 
  # timeStor[r,4] = (tROTATE2 - tROTATE1)[3]
  
  return(timeStor)
  
}

#all runtime results

# pSeq = seq(500, 5000, by = 500)
# n = 500
# R = 5
# 
# allRunTimeResults = matrix(0, nrow = length(pSeq), ncol = 8)
# allRunTimeResults = as.data.frame(allRunTimeResults)
# colnames(allRunTimeResults) = c("FABLESamples", 
#                                 "FABLEPostMean", 
#                                 "MGSP", 
#                                 "ROTATE", 
#                                 "HT", 
#                                 "POET", 
#                                 "SCAD", 
#                                 "LW")
# #dir.name = "~/Library/CloudStorage/GoogleDrive-shounak.chattopadhyay@gmail.com/My Drive/Research_Duke/FABLE/Code/FABLE-GitHub/extras/runTimeResults/"
# dir.name = "~/Library/CloudStorage/OneDrive-UniversityofVirginia/Research/FABLE/FABLE-GitHub/extras/runTimeResults"
# 
# for(j in 1:length(pSeq)) {
#   
#   p = pSeq[j]
#   
#   print(paste0("Dimension: ", p))
#   
#   runTimeResults_p = matrix(0, nrow = R, ncol = 8)
#   
#   for(r in 1:R) {
#     
#     print(paste0("Replicate: ", r))
#     
#     runtime_rep_r_dim_p = runTimeResults(n, p, 2001+r)
#     write.csv(runtime_rep_r_dim_p, paste0(dir.name, "/runtime_newMethods_rep=", r, "_p=", p, ".csv"))
#     runTimeResults_p[r,] = as.numeric(runtime_rep_r_dim_p)
#     
#   }
#   
#   allRunTimeResults[j,] = as.numeric(colMeans(runTimeResults_p))
#   
# }

# -------------- PLOT THE RESULTS ------------------

# Plot 1: FABLE Sampling vs MGSP

library(ggplot2)
library(latex2exp)

n = 500
R = 5

pSeqNew = seq(500, 4000, by = 500)

dir.name = "~/Library/CloudStorage/OneDrive-UniversityofVirginia/Research/FABLE/FABLE-GitHub/extras/runTimeResults"

# read FABLE and MGSP sampling runtimes from saved files

FABLESamplingTimes = matrix(0, nrow = R, ncol = length(pSeqNew))
MGSPSamplingTimes = matrix(0, nrow = R, ncol = length(pSeqNew))

for(r in 1:R) {
  
  for(colInd in 1:length(pSeqNew)) {
    
    source.file = read.csv(paste0(dir.name, "/runtime_rep=", r, "_p=", pSeqNew[colInd], ".csv"))
    
    FABLESamplingTimes[r,colInd] = source.file$FABLESamples
    MGSPSamplingTimes[r,colInd] = source.file$MGSP
    
  }
  
}

FABLESamplingMean = colMeans(FABLESamplingTimes)
FABLESamplingLower = apply(FABLESamplingTimes, 2, quantile, 0.001)
FABLESamplingUpper = apply(FABLESamplingTimes, 2, quantile, 0.999)

MGSPSamplingMean = colMeans(MGSPSamplingTimes)
MGSPSamplingLower = apply(MGSPSamplingTimes, 2, quantile, 0.001)
MGSPSamplingUpper = apply(MGSPSamplingTimes, 2, quantile, 0.999)

plot_df1 = data.frame("pSeq" = rep(pSeqNew, 2),
                      "runTimes" = c(log10(FABLESamplingMean), 
                                     log10(MGSPSamplingMean)),
                      "Method" = c(rep("FABLE", length(pSeqNew)), 
                                   rep("MGSP", length(pSeqNew))),
                      "lower" = c(log10(FABLESamplingLower), 
                                  log10(MGSPSamplingLower)),
                      "upper" = c(log10(FABLESamplingUpper), 
                                  log10(MGSPSamplingUpper)))

# CB-friendly colors for the 2 methods
cb_colors_plot1 <- c("FABLE" = "black", "MGSP" = "#0072B2")

# plot1 <- ggplot(plot_df1, aes(x = pSeq, y = runTimes, col = Method, fill = Method)) +
#   geom_line(size = 1.5) +  # solid lines
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +  # shaded, no borders
#   scale_color_manual(values = cb_colors_plot1) +
#   scale_fill_manual(values = cb_colors_plot1) +
#   theme_classic(base_size = 24) +  # globally larger font
#   xlab("Number of dimensions") +
#   ylab(TeX("Log$_{10}$ of runtime in seconds")) +
#   ggtitle("Runtime comparison for posterior sampling") +
#   theme(plot.title = element_text(hjust = 0.5))

plot1 <- ggplot(plot_df1, aes(x = pSeq, y = runTimes, col = Method, fill = Method)) +
  geom_line(linewidth = 1.5) +  # solid lines
  geom_ribbon(aes(ymin = lower, ymax = upper, group = Method), alpha = 0.2, color = NA) +  # fixed ribbon
  scale_color_manual(values = cb_colors_plot1) +
  scale_fill_manual(values = cb_colors_plot1) +
  theme_classic(base_size = 24) +  # white background, large fonts
  xlab("Number of dimensions") +
  ylab(TeX("Log$_{10}$ of runtime in seconds")) +
  ggtitle("Runtime comparison for posterior sampling") +
  theme(plot.title = element_text(hjust = 0.5))


print(plot1)

# Plot 2: FABLE Estimation vs ROTATE

FABLEEstTimes = matrix(0, nrow = R, ncol = length(pSeqNew))
ROTATEEstTimes = matrix(0, nrow = R, ncol = length(pSeqNew))
HTEstTimes = matrix(0, nrow = R, ncol = length(pSeqNew))
SCADEstTimes = matrix(0, nrow = R, ncol = length(pSeqNew))
LWEstTimes = matrix(0, nrow = R, ncol = length(pSeqNew))

for(r in 1:R) {
  
  for(colInd in 1:length(pSeqNew)) {
    
    source.file = read.csv(paste0(dir.name, "/runtime_newMethods_rep=", r, "_p=", pSeqNew[colInd], ".csv"))
    
    FABLEEstTimes[r,colInd] = source.file$FABLEPostMean
    ROTATEEstTimes[r,colInd] = source.file$ROTATE
    HTEstTimes[r,colInd] = source.file$HT
    SCADEstTimes[r,colInd] = source.file$SCAD
    LWEstTimes[r,colInd] = source.file$LW
    
  }
  
}

colQuantile <- function(datSet, p) {return(apply(datSet, 2, quantile, p))}

plot_df2 = data.frame("pSeq" = rep(pSeqNew, 5),
                      "runTimes" = c(log10(colMeans(FABLEEstTimes)), 
                                     log10(colMeans(ROTATEEstTimes)),
                                     log10(colMeans(HTEstTimes)),
                                     log10(colMeans(SCADEstTimes)),
                                     log10(colMeans(LWEstTimes))),
                      "Method" = c(rep("FABLE", length(pSeqNew)), 
                                   rep("ROTATE", length(pSeqNew)),
                                   rep("HT", length(pSeqNew)),
                                   rep("SCAD", length(pSeqNew)),
                                   rep("LW", length(pSeqNew))),
                      "lower" = c(log10(colQuantile(FABLEEstTimes, 0.001)), 
                                  log10(colQuantile(ROTATEEstTimes, 0.001)),
                                  log10(colQuantile(HTEstTimes, 0.001)),
                                  log10(colQuantile(SCADEstTimes, 0.001)),
                                  log10(colQuantile(LWEstTimes, 0.001))),
                      "upper" = c(log10(colQuantile(FABLEEstTimes, 0.999)), 
                                  log10(colQuantile(ROTATEEstTimes, 0.999)),
                                  log10(colQuantile(HTEstTimes, 0.999)),
                                  log10(colQuantile(SCADEstTimes, 0.999)),
                                  log10(colQuantile(LWEstTimes, 0.999))))

# Colors already defined for the 5 methods
cb_colors_plot2 <- c(
  "FABLE" = "black",
  "ROTATE" = "#009E73",
  "HT" = "#E69F00",
  "SCAD" = "#0072B2",
  "LW" = "#D55E00"
)

plot2 <- ggplot(plot_df2, aes(x = pSeq, y = runTimes, col = Method, fill = Method)) +
  geom_line(linewidth = 1.5) +  # all solid lines
  geom_ribbon(aes(ymin = lower, ymax = upper, group = Method), alpha = 0.2, color = NA) +  # fixed ribbon
  scale_color_manual(values = cb_colors_plot2) +
  scale_fill_manual(values = cb_colors_plot2) +
  theme_classic(base_size = 24) +  # white background, globally larger fonts
  xlab("Number of dimensions") +
  ylab(TeX("Log$_{10}$ of runtime in seconds")) +
  ggtitle("Runtime comparison for point estimation") +
  theme(plot.title = element_text(hjust = 0.5))

print(plot2)

## Combine plot1 and plot2

library(patchwork)

# Determine combined y limits
y_range <- range(c(plot_df1$lower, plot_df1$upper, plot_df2$lower, plot_df2$upper))

# Apply same y-axis limits and centered legends
plot1 <- plot1 + ylim(y_range) +
  theme(legend.position = "top",
        legend.justification = "center")

plot2 <- plot2 + ylim(y_range) +
  theme(legend.position = "top",
        legend.justification = "center")

# Combine side by side
library(patchwork)
combined <- plot1 + plot2 + plot_layout(ncol = 2)

print(combined)

ggsave(
  "extras/runTimeResults/Plots/runtimePlotsShaded.pdf",
  plot = combined,
  width = 18,       # increase width
  height = 8,       # increase height
  units = "in"
)








