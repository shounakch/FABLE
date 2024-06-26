# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rGamma <- function(n, a, b) {
    .Call(`_FABLE_rGamma`, n, a, b)
}

SymmetrizeMatrix <- function(A) {
    .Call(`_FABLE_SymmetrizeMatrix`, A)
}

CPPsvd <- function(X, flag) {
    .Call(`_FABLE_CPPsvd`, X, flag)
}

LogLikelihoodEval <- function(Y, M, Lambda, SigmaSq) {
    .Call(`_FABLE_LogLikelihoodEval`, Y, M, Lambda, SigmaSq)
}

CriterionFunction <- function(kInd, Y, U_Y, V_Y, svalsY) {
    .Call(`_FABLE_CriterionFunction`, kInd, Y, U_Y, V_Y, svalsY)
}

BisectionRecursion <- function(lower, upper, Y, U_Y, V_Y, svalsY) {
    .Call(`_FABLE_BisectionRecursion`, lower, upper, Y, U_Y, V_Y, svalsY)
}

RankEstimator <- function(Y, U_Y, V_Y, svalsY, kMax) {
    .Call(`_FABLE_RankEstimator`, Y, U_Y, V_Y, svalsY, kMax)
}

FABLEHyperParameters <- function(Y, U_Y, V_Y, svalsY, kEst, gamma0, delta0sq) {
    .Call(`_FABLE_FABLEHyperParameters`, Y, U_Y, V_Y, svalsY, kEst, gamma0, delta0sq)
}

cov_correct_matrix <- function(sigsq_hat, llprime_hat) {
    .Call(`_FABLE_cov_correct_matrix`, sigsq_hat, llprime_hat)
}

FABLEPostMean <- function(Y, gamma0, delta0sq, U_Y, V_Y, svalsY, kMax) {
    .Call(`_FABLE_FABLEPostMean`, Y, gamma0, delta0sq, U_Y, V_Y, svalsY, kMax)
}

FABLESampler <- function(Y, gamma0, delta0sq, MC, U_Y, V_Y, svalsY, kEst, tausq_est, gammaDeltasq, G0, varInflation) {
    .Call(`_FABLE_FABLESampler`, Y, gamma0, delta0sq, MC, U_Y, V_Y, svalsY, kEst, tausq_est, gammaDeltasq, G0, varInflation)
}

CCFABLEPostProcessing <- function(FABLEOutput, alpha) {
    .Call(`_FABLE_CCFABLEPostProcessing`, FABLEOutput, alpha)
}

CCFABLEPostProcessingSubmatrix <- function(FABLEOutput, alpha, SelectedIndices) {
    .Call(`_FABLE_CCFABLEPostProcessingSubmatrix`, FABLEOutput, alpha, SelectedIndices)
}

