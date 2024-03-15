#include <math.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
////[[Rcpp::plugins("cpp11")]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// // [[Rcpp::export]]
// double random_gamma(double a) {
//   return R::rgamma(a, 1.0);
// }

//Generates n G(a,b) random draws, a = shape, b = rate. E(X) = a/b.
// [[Rcpp::export]]
arma::vec rGamma(int n, double a, double b) {
  
  arma::vec v(n, fill::ones);
  
  v = randg(n, distr_param(a, 1.0/b));
  
  return v;
  
}

// [[Rcpp::export]]
arma::mat SymmetrizeMatrix(arma::mat A) {
  
  arma::mat A1(A.n_rows, A.n_cols, fill::zeros);
  A1 = A + A.t();
  A1.diag() = A1.diag() / 2.0;
  
  return A1;
  
}

// [[Rcpp::export]]
Rcpp::List CPPsvd(arma::mat X,
                  int flag){
  
  int n = X.n_rows;
  int p = X.n_cols;
  int r = std::min(n,p);
  
  arma::mat U(n, n, fill::zeros);
  arma::mat V(p, r, fill::zeros);
  arma::vec d(r, fill::zeros);
  
  if(flag == 1){
    
    svd(U, d, V, X, "std");
    
  } else if(flag == 2){
    
    svd(U, d, V, X, "dc");
    
  } else if(flag == 3){
    
    svd_econ(U, d, V, X, "both", "std");
    
  } else {
    
    svd_econ(U, d, V, X, "both", "dc");
    
  }
  
  return Rcpp::List::create(Rcpp::Named("U") = U,
                            Rcpp::Named("d") = d,
                            Rcpp::Named("V") = V);
  
}

// // [[Rcpp::export]]
// arma::mat MultEachCol(arma::mat A,
//                       arma::vec b) {
//   
//   //Multiply each column of A by corresponding entry of b.
//   //A.n_col must equal length of b.
//   
//   arma::mat C(A.n_rows, A.n_cols, fill::zeros);
//   
//   for(int j=0; j < A.n_cols; ++j){
//     
//     C.col(j) = b(j) * A.col(j);
//     
//   }
//   
//   return C;
//   
// }

// [[Rcpp::export]]
double CPPLogLikelihoodEval(arma::mat Y,
                            arma::mat M,
                            arma::mat Lambda,
                            arma::vec SigmaSq) {
  
  //Evaluates the likelihood function for y_{ij} ~ N(\eta_i' \lambda_j, \sigma_j^2).
  
  int n = Y.n_rows;
  int p = Y.n_cols;
  double Term1 = 0;
  double Term2 = 0;
  double Term = 0;
  
  arma::vec SigmaVec = sqrt(SigmaSq);
  
  arma::mat DiffMat1(n, p, fill::zeros);
  DiffMat1 = Y - (M * Lambda.t());
  arma::mat DiffMat2(n, p, fill::zeros);
  // DiffMat2 = MultEachCol(DiffMat1, 1.0 / SigmaVec);
  DiffMat2.each_row() /= SigmaVec.t();
  
  Term1 = -n * accu(log(SigmaVec));
  Term2 = -0.5 * accu(square(DiffMat2));
  
  Term = Term1 + Term2;
  Term = Term / (n*p);
  
  return Term;
  
}

// [[Rcpp::export]]
double CPPCriterionFunction(int kInd,
                            arma::mat Y,
                            arma::mat U_Y,
                            arma::mat V_Y,
                            arma::vec svalsY) {
  
  // Extract the n*k left singular matrix U
  
  int n = Y.n_rows;
  int p = Y.n_cols;
  int r = svalsY.n_elem;
  
  arma::mat U_K(n, kInd, fill::zeros);
  arma::mat V_K(p, kInd, fill::zeros);
  arma::mat D_K(kInd, kInd, fill::zeros);
  
  U_K = U_Y.cols(0, kInd-1);
  V_K = V_Y.cols(0, kInd-1);
  
  D_K = diagmat(svalsY.subvec(0, kInd-1));
  
  // Estimate \sigma_j^2 and \tau^2
  
  // #t1 = proc.time()
  //     
  // # Most time taking step below!
  // # sigsq_hat_diag = colSums(((diag(nwhole) - (U_K %*% t(U_K))) %*% Ymat)^2) / nwhole
  // # tausq_est = mean( (colSums(((U_K %*% t(U_K)) %*% Ymat)^2) / nwhole) / (kInd * sigsq_hat_diag))
  //     
  // #Alternative (MUCH FASTER - include in main code as well.)
  //     
  // #t2 = proc.time()
  
  arma::mat UDVt(n, p, fill::zeros);
  
  UDVt = U_K * D_K * V_K.t();
  
  arma::vec sigsq_hat_diag(p, fill::zeros);
  sigsq_hat_diag = (sum(square(Y - UDVt), 0)).t() / n; //problematic statement
  
  double tausq_est = 1.0;
  
  tausq_est = (mean(sum(square(UDVt), 0).t() / sigsq_hat_diag)) / (n * kInd);
  
  // tausq_est = mean( (colSums(UDVt^2) / nwhole) / (kInd * sigsq_hat_diag))
  
  // Estimate factor loadings \Lambda
  
  arma::mat LambdaEst(p, kInd, fill::zeros);  
  
  arma::mat YtU(p, kInd, fill::zeros);
  YtU = V_K;
  YtU.each_row() %= svalsY.subvec(0, kInd-1).t();
  // LambdaEst = (sqrt(n) / (n + tausq_est)) * (Y.t() * U_K); //can be made faster
  LambdaEst = (sqrt(n) / (n + tausq_est)) * (YtU);
  
  // t4 = proc.time()
  
  // Estimate latent factor matrix M = \sqrt{n} * U_K
  
  arma::mat FactorEst(n, kInd, fill::zeros);    
  
  FactorEst = sqrt(n) * U_K;
  
  // Obtain log-likelihood of whole data set
  
  double LogLikDataSet = 0;    
  
  LogLikDataSet = CPPLogLikelihoodEval(Y, FactorEst, LambdaEst, sigsq_hat_diag);
  
  // t5 = proc.time()
  
  // Store this value
  
  double CriterionStor = 0;
  double maxint = std::max(n, p);
  double minint = std::min(n, p);
  
  CriterionStor = (-2 * LogLikDataSet) + (kInd * maxint * log(minint) / (n*p)); //logL already scaled by n*p
  
  return CriterionStor;
  
}

// [[Rcpp::export]]
arma::vec CPPBisectionRecursion(int lower, 
                                int upper,
                                arma::mat Y,
                                arma::mat U_Y,
                                arma::mat V_Y,
                                arma::vec svalsY) {
  
  if(upper - lower > 5) {
    
    int midpoint = (lower + upper) / 2;
    
    // Rcpp::Rcout << lower << midpoint << upper << std::endl;
    
    double CriterionLower = 0;
    // double CriterionUpper = 0;
    double CriterionMidPoint = 0;
    
    CriterionLower = CPPCriterionFunction(lower, Y, U_Y, V_Y, svalsY);
    // CriterionUpper = CPPCriterionFunction(upper, Y, U_Y, V_Y, svalsY);
    CriterionMidPoint = CPPCriterionFunction(midpoint, Y, U_Y, V_Y, svalsY);
    
    if(CriterionLower < CriterionMidPoint) {
      
      return CPPBisectionRecursion(lower, midpoint, Y, U_Y, V_Y, svalsY);
      
    }else {
      
      int ThreeFourthPoint = (midpoint + upper) / 2;
      
      double CriterionThreeFourth = 0;
      
      CriterionThreeFourth = CPPCriterionFunction(ThreeFourthPoint, Y, U_Y, V_Y, svalsY);
      
      if(CriterionMidPoint <= CriterionThreeFourth) {
        
        return CPPBisectionRecursion(lower, midpoint, Y, U_Y, V_Y, svalsY);
        
      }else {
        
        return CPPBisectionRecursion(midpoint, upper, Y, U_Y, V_Y, svalsY);
        
      }
      
    }
    
  }else {
    
    arma::vec outputvec(2, fill::zeros);
    outputvec(0) = lower;
    outputvec(1) = upper;
    
    return outputvec;
    
  }
  
}

// [[Rcpp::export]]
int CPPRankEstimator(arma::mat Y, 
                     arma::mat U_Y,
                     arma::mat V_Y,
                     arma::vec svalsY,
                     int kMax) {
  
  // Use the idea in `Determining the number of factors in high-dimensional generalized latent factor models` (Chen & Li, 2022)
  // U_Y is n*r matrix of left singular vectors.
  // V_Y is p*r matrix of right singular vectors.
  // svalsY is the r*1 vector of singular values.
  // kMax is the upper bound on search space.
  // Use bisection search to identify minimizer k.
  // r = min(n, p), n = nrow(Y), p = ncol(Y).
  
  // Y must be pre-processed to have zero row mean before.
  
  int n = Y.n_rows;
  int p = Y.n_cols;
  
  // Obtain SVD of entire matrix
  
  // arma::vec svalsY(r, fill::zeros);
  // svalsY = svdmod["d"];
  // arma::mat U_Y(n, r, fill::zeros);
  // arma::mat V_Y(p, r, fill::zeros);
  // U_Y = svdmod["u"];
  // V_Y = svdmod["v"];
  
  // Define penalized criterion as a function of k
  
  arma::vec BisectionInterval(2, fill::zeros);
  
  BisectionInterval = CPPBisectionRecursion(1, kMax, Y, U_Y, V_Y, svalsY);
  
  int l_int = BisectionInterval(1) - BisectionInterval(0) + 1 + 0.1; //add extra 0.1 to ensure floor function works, numerical overflow issues
  
  // Rcpp::Rcout << BisectionInterval << std::endl;
  // Rcpp::Rcout << "Fine till here - CPPRankEstimator" << l_int << std::endl;
  
  arma::vec kEvalSeq(l_int, fill::zeros); 
  kEvalSeq = regspace(BisectionInterval(0), 1, BisectionInterval(1));
  arma::vec CriterionValues(kEvalSeq.n_elem, fill::zeros);
  
  for(int j = 0; j < kEvalSeq.n_elem; ++j){
    
    int kEval = kEvalSeq(j) + 0.1; //add 0.1 to avoid overflow issues
    CriterionValues(j) = CPPCriterionFunction(kEval, Y, U_Y, V_Y, svalsY);
    
  }
  
  //Rcpp::Rcout << "Sequence of k" << kEvalSeq << std::endl;
  //Rcpp::Rcout << "Criterion" << CriterionValues << std::endl;
  
  int kMinimizer = 1;
  
  kMinimizer = kEvalSeq(CriterionValues.index_min());    
  
  return kMinimizer;
  
}

// [[Rcpp::export]]
arma::mat CPPcov_correct_matrix(arma::vec sigsq_hat, 
                                arma::mat llprime_hat) {
  
  int p = sigsq_hat.n_elem;
  arma::vec diag_llprime(p, fill::zeros);
  diag_llprime = llprime_hat.diag();
  
  arma::mat B(p, p, fill::zeros);
  arma::mat num_matrix(p, p, fill::zeros);
  arma::mat den_matrix(p, p, fill::zeros);
  
  num_matrix = (diag_llprime * diag_llprime.t()) + square(llprime_hat);
  den_matrix = (diag_llprime * sigsq_hat.t()) + (sigsq_hat * diag_llprime.t());
  
  B = num_matrix / den_matrix;
  B.diag() = B.diag() / 2.0;
  
  arma::mat oneMat(p, p, fill::ones);
  
  B = sqrt(oneMat + B);
  
  return B;
  
}

//Requires input of SVD.
// [[Rcpp::export]]
arma::mat CPPFABLEPostMean(arma::mat Y, 
                           double gamma0, 
                           double delta0sq,
                           arma::mat U_Y,
                           arma::mat V_Y,
                           arma::vec svalsY,
                           int kMax) {
  
  int n = Y.n_rows;
  int p = Y.n_cols;
  
  // svdmod = svd(Y)
  
  // CHOOSE k 
  
  int k = 1;
  k = CPPRankEstimator(Y, U_Y, V_Y, svalsY, kMax);
  
  arma::mat U_K(n, k, fill::zeros);
  arma::mat V_K(p, k, fill::zeros);
  arma::mat D_K(k, k, fill::zeros);
  
  U_K = U_Y.cols(0, k-1);
  V_K = V_Y.cols(0, k-1);
  
  D_K = diagmat(svalsY.subvec(0, k-1));
  
  // Estimate \tau^2 and \sigma_j^2
  
  arma::mat UDVt(n, p, fill::zeros);
  
  UDVt = U_K * D_K * V_K.t();
  
  arma::vec sigsq_hat_diag(p, fill::zeros);
  sigsq_hat_diag = (sum(square(Y - UDVt), 0)).t() / n; //problematic statement
  
  double tausq_est = 1.0;
  
  tausq_est = (mean(sum(square(UDVt), 0).t() / sigsq_hat_diag)) / (n * k);
  
  // Obtain the hyperparameters
  
  arma::mat YtU(p, k, fill::zeros);
  YtU = V_K;
  YtU.each_row() %= svalsY.subvec(0, k-1).t();
  
  arma::mat G0(p, k, fill::zeros);
  arma::mat G(p, p, fill::zeros);
  
  G0 = (sqrt(n) / (n + (1/tausq_est))) * YtU;
  G = G0 * G0.t();
  
  double gamma_n = gamma0 + n;    // \gamma_n = \gamma_0 + n
  
  arma::vec gamma_n_deltasq(p, fill::ones);
  arma::vec oneVec(p, fill::ones);
  double impCoef = n / (n + (1 / tausq_est));
  
  gamma_n_deltasq = ((gamma0 * delta0sq) * oneVec) + sum(square(Y), 0).t() - (impCoef * sum(square(YtU), 1));  
  
  // gamma_n_deltasq = rep(gamma0*delta0sq, p) + 
  //   as.numeric(apply(Y^2, 2, sum)) -
  //   ((n / (n + (1 / tausq_est))) * as.numeric(apply((t(YtU))^2, 2, sum)))
  
  arma::mat SigmaEstimate(p, p, fill::zeros);
  arma::mat CovPostMean(p, p, fill::zeros);
  SigmaEstimate = diagmat(gamma_n_deltasq / (gamma_n - 2));
  
  CovPostMean = G + SigmaEstimate;
  
  return CovPostMean;
  
}

// [[Rcpp::export]]
Rcpp::List CPPFABLESampler(arma::mat Y, 
                           double gamma0, 
                           double delta0sq, 
                           int MC,
                           arma::mat U_Y,
                           arma::mat V_Y,
                           arma::vec svalsY,
                           int kMax) {
  
  int n = Y.n_rows;
  int p = Y.n_cols;
  
  // svdmod = svd(Y)
  
  // CHOOSE k 
  
  int k = 1;
  k = CPPRankEstimator(Y, U_Y, V_Y, svalsY, kMax);
  
  arma::mat U_K(n, k, fill::zeros);
  arma::mat V_K(p, k, fill::zeros);
  arma::mat D_K(k, k, fill::zeros);
  
  U_K = U_Y.cols(0, k-1);
  V_K = V_Y.cols(0, k-1);
  
  D_K = diagmat(svalsY.subvec(0, k-1));
  
  // Estimate \tau^2 and \sigma_j^2
  
  arma::mat UDVt(n, p, fill::zeros);
  
  UDVt = U_K * D_K * V_K.t();
  
  arma::vec sigsq_hat_diag(p, fill::zeros);
  sigsq_hat_diag = (sum(square(Y - UDVt), 0)).t() / n; //problematic statement
  
  double tausq_est = 1.0;
  
  tausq_est = (mean(sum(square(UDVt), 0).t() / sigsq_hat_diag)) / (n * k);
  
  // Obtain the hyperparameters
  
  arma::mat YtU(p, k, fill::zeros);
  YtU = V_K;
  YtU.each_row() %= svalsY.subvec(0, k-1).t();
  
  arma::mat G0(p, k, fill::zeros);
  arma::mat G(p, p, fill::zeros);
  
  G0 = (sqrt(n) / (n + (1/tausq_est))) * YtU;
  G = G0 * G0.t();
  
  double gamma_n = gamma0 + n;    // \gamma_n = \gamma_0 + n
  
  arma::vec gamma_n_deltasq(p, fill::ones);
  arma::vec oneVec(p, fill::ones);
  double impCoef = n / (n + (1 / tausq_est));
  
  gamma_n_deltasq = ((gamma0 * delta0sq) * oneVec) + sum(square(Y), 0).t() - (impCoef * sum(square(YtU), 1));  
  
  arma::mat SigmaSampleStor(MC, p, fill::ones);
  arma::mat LambdaSampleStor(MC, k*p, fill::zeros); //Return column-vectorized version of samples of Lambda
  
  double a_n = 1 / sqrt(n + (1 / tausq_est));
  
  for(int m = 0; m < MC; ++m) {
    
    // Sample \sigma_j^2, j = 1, \ldots, p.
    
    arma::vec sigmaSqSample(p, fill::ones);
    sigmaSqSample = (0.5 * gamma_n_deltasq) / rGamma(p, 0.5 * gamma_n, 1.0);
    
    SigmaSampleStor.row(m) = sigmaSqSample.t();
    
    // Sample \Lambda
    
    arma::mat ZSample(p, k, fill::randn);
    arma::mat LambdaSample(p, k, fill::zeros);
    ZSample.each_col() %= sqrt(sigmaSqSample);
    
    LambdaSample = G0 + (a_n * ZSample);
    LambdaSampleStor.row(m) = vectorise(LambdaSample.t()).t();
    
  }
  
  return Rcpp::List::create(Rcpp::Named("LambdaSamples") = LambdaSampleStor,
                            Rcpp::Named("SigmaSqSamples") = SigmaSampleStor,
                            Rcpp::Named("G") = G,
                            Rcpp::Named("SigmaSqEstimate") = gamma_n_deltasq / (gamma_n - 2),
                            Rcpp::Named("EstimatedRank") = k);
  
}

// [[Rcpp::export]]
Rcpp::List CPPCCFABLEPostProcessing(Rcpp::List FABLEOutput,
                                    arma::mat CovCorrectMatrix,
                                    double alpha) {
  
  arma::mat LambdaSamples = FABLEOutput["LambdaSamples"];
  arma::mat SigmaSqSamples = FABLEOutput["SigmaSqSamples"];
  arma::mat G = FABLEOutput["G"];
  
  int nMC = SigmaSqSamples.n_rows;
  int p = SigmaSqSamples.n_cols;
  int k = FABLEOutput["EstimatedRank"];
  
  arma::mat CovMatPostMean(p, p, fill::zeros);
  arma::mat CovMatLower(p, p, fill::zeros);
  arma::mat CovMatUpper(p, p, fill::zeros);
  
  // Define the sample extractor function
  
  arma::vec oneVec(nMC, fill::ones);
  
  for(int j1 = 0; j1 < p; ++j1) {
    
    if(j1 % 100 == 0) {
      
      Rcpp::Rcout << "j1: " << j1 << std::endl;
      
    }
    
    arma::mat Lambdaj1Samples(nMC, k, fill::zeros);
    
    Lambdaj1Samples = LambdaSamples.cols(j1*k, (j1*k)+k-1);
    
    for(int j2 = j1; j2 < p; ++j2) {
      
      arma::mat Lambdaj2Samples(nMC, k, fill::zeros);
      
      Lambdaj2Samples = LambdaSamples.cols(j2*k, (j2*k)+k-1);
      
      // Extract samples for L[j1,k2], L = \Lambda \Lambda'.
      
      arma::vec Lj1j2Samples(nMC, fill::zeros);
      
      Lj1j2Samples = sum(Lambdaj1Samples % Lambdaj2Samples, 1);
      
      arma::vec Gj1j2Vec(nMC, fill::ones);
      
      Gj1j2Vec = G(j1,j2) * oneVec;
      
      // Rcpp::Rcout << G(j1,j2) << std::endl;
      
      arma::vec Lj1j2SamplesCorrected(nMC, fill::zeros);
      
      Lj1j2SamplesCorrected = Gj1j2Vec + (CovCorrectMatrix(j1,j2) * (Lj1j2Samples - Gj1j2Vec));
      
      arma::vec Covj1j2Samples(nMC, fill::zeros);
      
      if(j2 == j1) {
        
        Covj1j2Samples = Lj1j2SamplesCorrected + SigmaSqSamples.col(j1);
        
      }else {
        
        Covj1j2Samples = Lj1j2SamplesCorrected;
        
      }
      
      // Extract summary statistics for Cov[j1,j2]
      
      CovMatPostMean(j1,j2) = accu(Covj1j2Samples) / nMC;
      arma::vec LowerQuantile = {alpha / 2.0};
      arma::vec UpperQuantile = {1 - (alpha / 2.0)};
      
      CovMatLower(j1,j2) = as_scalar(quantile(Covj1j2Samples, LowerQuantile));
      CovMatUpper(j1,j2) = as_scalar(quantile(Covj1j2Samples, UpperQuantile));
      
    }
    
    if(j1 % 500 == 0) {
      
      Rcpp::checkUserInterrupt();
      
    }
    
  }
  
  CovMatPostMean = SymmetrizeMatrix(CovMatPostMean);
  CovMatLower = SymmetrizeMatrix(CovMatLower);
  CovMatUpper = SymmetrizeMatrix(CovMatUpper);
  
  return Rcpp::List::create(Rcpp::Named("PostMeanMatrix") = CovMatPostMean,
                            Rcpp::Named("LowerQuantileMatrix") = CovMatLower,
                            Rcpp::Named("UpperQuantileMatrix") = CovMatUpper);
  
}

//Need to check this function carefully.
// [[Rcpp::export]]
Rcpp::List CCFABLEPostProcessingSubmatrix(Rcpp::List FABLEOutput,
                                          arma::mat CovCorrectMatrix,
                                          double alpha,
                                          arma::uvec SelectedIndices) {
  
  SelectedIndices = SelectedIndices - 1; //Somehow need to account for different index starts in C++ and R.
  
  arma::mat LambdaSamples = FABLEOutput["LambdaSamples"];
  arma::mat SigmaSqSamples = FABLEOutput["SigmaSqSamples"];
  arma::mat G = FABLEOutput["G"];
  
  int nMC = SigmaSqSamples.n_rows;
  int p = SigmaSqSamples.n_cols;
  int k = FABLEOutput["EstimatedRank"];
  
  int SubLength = SelectedIndices.n_elem;
  arma::mat CovMatPostMean(SubLength, SubLength, fill::zeros);
  arma::mat CovMatLower(SubLength, SubLength, fill::zeros);
  arma::mat CovMatUpper(SubLength, SubLength, fill::zeros);
  
  // Define the sample extractor function
  
  arma::vec oneVec(nMC, fill::ones);
  
  for(int j1 = 0; j1 < SubLength; ++j1) {
    
    if((j1+1) % 100 == 0) {
      
      Rcpp::Rcout << "j1: " << j1 << std::endl;
      
    }
    
    int FirstIndex = SelectedIndices(j1);
    arma::mat Lambdaj1Samples(nMC, k, fill::zeros);
    
    Lambdaj1Samples = LambdaSamples.cols(FirstIndex*k, (FirstIndex*k)+k-1);
    
    for(int j2 = j1; j2 < SubLength; ++j2) {
      
      int SecondIndex = SelectedIndices(j2);
      arma::mat Lambdaj2Samples(nMC, k, fill::zeros);
      
      Lambdaj2Samples = LambdaSamples.cols(SecondIndex*k, (SecondIndex*k)+k-1);
      
      // Extract samples for L[j1,k2], L = \Lambda \Lambda'.
      
      arma::vec Lj1j2Samples(nMC, fill::zeros);
      
      Lj1j2Samples = sum(Lambdaj1Samples % Lambdaj2Samples, 1);
      
      // Correct the samples for coverage
      
      arma::vec Gj1j2Vec(nMC, fill::ones);
      
      Gj1j2Vec = G(FirstIndex,SecondIndex) * oneVec;
      
      arma::vec Lj1j2SamplesCorrected(nMC, fill::zeros);
      
      Lj1j2SamplesCorrected = Gj1j2Vec + (CovCorrectMatrix(FirstIndex,SecondIndex) * (Lj1j2Samples - Gj1j2Vec));
      
      arma::vec Covj1j2Samples(nMC, fill::zeros);
      
      if(FirstIndex == SecondIndex) {
        
        Covj1j2Samples = Lj1j2SamplesCorrected + SigmaSqSamples.col(FirstIndex);
        
      }else {
        
        Covj1j2Samples = Lj1j2SamplesCorrected;
        
      }
      
      // Store summaries
      
      CovMatPostMean(j1,j2) = accu(Covj1j2Samples) / nMC;
      arma::vec LowerQuantile = {alpha / 2.0};
      arma::vec UpperQuantile = {1 - (alpha / 2.0)};
      
      CovMatLower(j1,j2) = as_scalar(quantile(Covj1j2Samples, LowerQuantile));
      CovMatUpper(j1,j2) = as_scalar(quantile(Covj1j2Samples, UpperQuantile));
      
    }
    
    if(j1 % 100 == 0) {
      
      Rcpp::checkUserInterrupt();
      
    }
    
  }
  
  CovMatPostMean = SymmetrizeMatrix(CovMatPostMean);
  CovMatLower = SymmetrizeMatrix(CovMatLower);
  CovMatUpper = SymmetrizeMatrix(CovMatUpper);
  
  // Timer
  
  //clock.stop("naptimes");
  // for(int i = 0; i < res.size(); ++i) {
  //   
  //   res[i] = 
  //   
  // }
  
  return Rcpp::List::create(Rcpp::Named("PostMeanMatrix") = CovMatPostMean,
                            Rcpp::Named("LowerQuantileMatrix") = CovMatLower,
                            Rcpp::Named("UpperQuantileMatrix") = CovMatUpper);
  
  // return res;
  
}

///////// OLD RCpp functions below (developing)

// // [[Rcpp::export]]
// arma::cube CPPCCFABLESampler(arma::mat Y, 
//                           double gamma0, 
//                           double delta0sq, 
//                           int MC,
//                           arma::mat U_Y,
//                           arma::mat V_Y,
//                           arma::vec svalsY,
//                           int kMax) {
//   
//   int n = Y.n_rows;
//   int p = Y.n_cols;
//   int r = std::min(n,p);
//   
//   // CHOOSE k 
//   
//   int k = 1;
//   k = CPPRankEstimator(Y, U_Y, V_Y, svalsY, kMax, r);
//   
//   arma::mat U(n, k, fill::zeros);
//   arma::mat V(p, k, fill::zeros);
//   arma::mat D(k, k, fill::zeros);
//   
//   U = U_Y.cols(0,k-1);
//   V = V_Y.cols(0,k-1);
//   D = diagmat(svalsY.subvec(0,k-1));
//   
// // CHOOSE \tau^2
// 
// arma::mat UDVt(n, p, fill::zeros);
//   
//   UDVt = U * D * V.t();
//   
//   arma::vec sigsq_hat_diag(p, fill::zeros);
//   sigsq_hat_diag = sum(pow(Y - UDVt, 2.0), 0) / n;
//   
//   double tausq_est = 1.0;
//   
//   tausq_est = (mean(sum(pow(UDVt, 2.0), 0) / sigsq_hat_diag)) / (n * k);
//     
// // Obtain the hyperparameters
//     
//     arma::mat YtU(p, k, fill::zeros);
//     arma::mat G0(p, k, fill::zeros);
//     arma::mat G(p, p, fill::zeros);
//     
//     YtU = MultEachCol(V, svalsY.subvec(0,k-1));
//     G0 = (sqrt(n) / (n + (1.0/tausq_est))) * YtU;
//     G = G0 * G0.t();
//     
//     arma::mat CCMatrix(p, p, fill::zeros);
//     
//     CCMatrix = CPPcov_correct_matrix(sigsq_hat_diag, G);
//     
//     double gamma_n = gamma0 + n;
//       
//     arma::vec gamma_n_deltasq(p, fill::zeros);
//     arma::vec one_vec(p, fill::ones);
//     double gdelCoef = n / (n + (1.0 / tausq_est));
//     //Done till here.
//     
//     gamma_n_deltasq = ((gamma0*delta0sq)*one_vec) + sum(pow(Y, 2.0), 0) - (gdelCoef * sum(pow(YtU.t(), 2.0), 0));
//   
//     // gamma_n_deltasq = rep(gamma0*delta0sq, p) + 
//     //   as.numeric(apply(Y^2, 2, sum)) -
//     //   ((n / (n + (1 / tausq_est))) * as.numeric(apply((t(YtU))^2, 2, sum)))
//       
// // Now proceed to store samples of \Psi = \Lambda \Lambda' + \Sigma #####
// // WARNING: Massive RAM consumption for large n, p ... #####
// 
//     arma::cube CovSampleStor(MC, p, p, fill::zeros);
//       
//     //CovSampleStor = array(0, dim = c(MC, p, p))
//     double a_n = 1 / sqrt(n + (1 / tausq_est));
//     arma::vec sigmaSqSample(p, fill::ones);
//     arma::mat ZSample(p, k, fill::zeros);
//     
//     for(int m=0; m < MC; ++m){
//       
//       // Sample \sigma_j^2 \sim IG(\gamma_n/2, \gamma_n\delta_j^2/2)
//       
//       for(j = 0; j < p; ++j){
//         
//         
//         
//       }
//       
//     }
//         
//     //t1 = proc.time()
//         
//         for(m in 1:MC) {
//           
// ## Sample \sigma_j^2 \sim IG(\gamma_n / 2, \gamma_n \delta_j^2 / 2), j=1,\ldots,p. ##
//           
//           sigmaSqSample = 1 / rgamma(n = p, shape = gamma_n / 2,
//                                      rate = gamma_n_deltasq / 2)
//             
//             SigmaSample = diag(sigmaSqSample)
//             
// ## Sample low-rank part L_C = \Lambda \Lambda' ##
//             
//             ZSample = matrix(rnorm(p*k), nrow = p, ncol = k)
//               SigmaHalfZ = sweep(ZSample, 1, sqrt(sigmaSqSample), "*")
//               SigmaHalfZG0t = SigmaHalfZ %*% t(G0)
//               
//               Part2Matrix = (a_n * (SigmaHalfZG0t + t(SigmaHalfZG0t))) + 
//                 (a_n^2 * ((ZSample %*% t(ZSample)) * (sqrt(sigmaSqSample) %*% t(sqrt(sigmaSqSample)))))
//               
//               LLtSample = G + (CCMatrix * Part2Matrix)
//               
// ## Store sample of \Psi ##
//               
//               PsiSample = LLtSample + SigmaSample
//               CovSampleStor[m,,] = PsiSample
//               
//         }
//         
//         t2 = proc.time()
//           
//           return(CovSampleStor)
//           
// }

