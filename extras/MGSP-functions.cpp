#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
arma::mat SymmetrizeMatrix(arma::mat A) {
  
  arma::mat A1(A.n_rows, A.n_cols, fill::zeros);
  A1 = A + A.t();
  A1.diag() = A1.diag() / 2.0;
  
  return A1;
  
}

//Need VY to adjust for the fact that the code in infinitefactor::linearMGSP is wrong.
// [[Rcpp::export]]
Rcpp::List MGSPPostProcessingSubmatrix(Rcpp::List LambdaSamples,
                                       arma::mat SigmaSqSamples,
                                       double alpha,
                                       arma::uvec SelectedIndices,
                                       arma::vec VY) {
  
  SelectedIndices = SelectedIndices - 1; //Somehow need to account for different index starts in C++ and R.
  
  int nMC = SigmaSqSamples.n_cols;
  int p = SigmaSqSamples.n_rows;
  
  int SubLength = SelectedIndices.n_elem;
  arma::mat CovMatPostMean(SubLength, SubLength, fill::zeros);
  arma::mat CovMatLower(SubLength, SubLength, fill::zeros);
  arma::mat CovMatUpper(SubLength, SubLength, fill::zeros);
  
  arma::mat sqrtVY = sqrt(VY);
  
  // Define the sample extractor function
  
  arma::vec oneVec(nMC, fill::ones);
  
  for(int j1 = 0; j1 < SubLength; ++j1) {
    
    // if((j1+1) % 10 == 0) {
    //   
    //   Rcpp::Rcout << "j1: " << j1 << std::endl;
    //   
    // }
    
    Rcpp::Rcout << "j1: " << j1 << std::endl;
    
    int FirstIndex = SelectedIndices(j1);
    
    for(int j2 = j1; j2 < SubLength; ++j2) {
      
      int SecondIndex = SelectedIndices(j2);
      
      arma::vec Lambdaj1j2Samples(nMC, fill::zeros);
      arma::vec Covj1j2Samples(nMC, fill::zeros);
      
      for(int m = 0; m < nMC; ++m) {
        
        arma::mat LambdaIndiv = LambdaSamples[m];
        
        Lambdaj1j2Samples(m) = accu(LambdaIndiv.row(FirstIndex).t() % LambdaIndiv.row(SecondIndex).t());
        
      }
      
      if(FirstIndex == SecondIndex) {
        
        Covj1j2Samples = Lambdaj1j2Samples + SigmaSqSamples.row(FirstIndex).t();
        
      }else {
        
        Covj1j2Samples = Lambdaj1j2Samples;
        
      }
      
      Covj1j2Samples = sqrtVY(j1) * sqrtVY(j2) * Covj1j2Samples;
      
      // Store summaries
      
      CovMatPostMean(j1,j2) = accu(Covj1j2Samples) / nMC;
      arma::vec LowerQuantile = {alpha / 2.0};
      arma::vec UpperQuantile = {1 - (alpha / 2.0)};
      
      CovMatLower(j1,j2) = as_scalar(quantile(Covj1j2Samples, LowerQuantile));
      CovMatUpper(j1,j2) = as_scalar(quantile(Covj1j2Samples, UpperQuantile));
      
    }
    
    Rcpp::checkUserInterrupt();
    
    // if(j1 % 10 == 0) {
    //   
    //   Rcpp::checkUserInterrupt();
    //   
    // }
    
  }
  
  CovMatPostMean = SymmetrizeMatrix(CovMatPostMean);
  CovMatLower = SymmetrizeMatrix(CovMatLower);
  CovMatUpper = SymmetrizeMatrix(CovMatUpper);
  
  return Rcpp::List::create(Rcpp::Named("PostMeanMatrix") = CovMatPostMean,
                            Rcpp::Named("LowerQuantileMatrix") = CovMatLower,
                            Rcpp::Named("UpperQuantileMatrix") = CovMatUpper);
  
}

// [[Rcpp::export]]
Rcpp::List MGSPPostProcessingSubmatrix_Optimized(Rcpp::List LambdaSamples,
                                                 arma::mat SigmaSqSamples,
                                                 double alpha,
                                                 arma::uvec SelectedIndices,
                                                 arma::vec VY) {
  
  SelectedIndices = SelectedIndices - 1;  // R -> C++ indexing
  int nMC = LambdaSamples.size();
  int p = SigmaSqSamples.n_rows;
  int SubLength = SelectedIndices.n_elem;
  
  arma::mat CovMatPostMean(SubLength, SubLength, fill::zeros);
  arma::mat CovMatLower(SubLength, SubLength, fill::zeros);
  arma::mat CovMatUpper(SubLength, SubLength, fill::zeros);
  
  arma::vec sqrtVY = sqrt(VY);
  
  // Find max k over all samples
  int kMax = 0;
  std::vector<int> kVec(nMC);
  for(int m = 0; m < nMC; ++m) {
    int k_m = Rcpp::as<arma::mat>(LambdaSamples[m]).n_cols;
    kVec[m] = k_m;
    if(k_m > kMax) kMax = k_m;
  }
  
  // Pre-allocate cube with zero-padding: p x kMax x nMC
  arma::cube LambdaCube(p, kMax, nMC, arma::fill::zeros);
  for(int m = 0; m < nMC; ++m) {
    arma::mat LambdaMat = Rcpp::as<arma::mat>(LambdaSamples[m]);
    LambdaCube.slice(m).cols(0, kVec[m]-1) = LambdaMat;  // copy actual columns
  }
  
  for(int j1 = 0; j1 < SubLength; ++j1) {
    int FirstIndex = SelectedIndices(j1);
    
    for(int j2 = j1; j2 < SubLength; ++j2) {
      int SecondIndex = SelectedIndices(j2);
      
      arma::vec Lambdaj1j2Samples(nMC, fill::zeros);
      arma::vec Covj1j2Samples(nMC, fill::zeros);
      
      for(int m = 0; m < nMC; ++m) {
        int k_m = kVec[m];  // actual rank for this sample
        arma::rowvec row1 = LambdaCube.slice(m).row(FirstIndex).cols(0, k_m-1);
        arma::rowvec row2 = LambdaCube.slice(m).row(SecondIndex).cols(0, k_m-1);
        Lambdaj1j2Samples(m) = arma::dot(row1, row2);
      }
      
      if(FirstIndex == SecondIndex) {
        Covj1j2Samples = Lambdaj1j2Samples + SigmaSqSamples.row(FirstIndex).t();
      } else {
        Covj1j2Samples = Lambdaj1j2Samples;
      }
      
      Covj1j2Samples = sqrtVY(j1) * sqrtVY(j2) * Covj1j2Samples;
      
      // Store summaries
      CovMatPostMean(j1,j2) = accu(Covj1j2Samples) / nMC;
      arma::vec LowerQuantile = {alpha / 2.0};
      arma::vec UpperQuantile = {1 - (alpha / 2.0)};
      CovMatLower(j1,j2) = as_scalar(quantile(Covj1j2Samples, LowerQuantile));
      CovMatUpper(j1,j2) = as_scalar(quantile(Covj1j2Samples, UpperQuantile));
    }
    
    Rcpp::checkUserInterrupt();
  }
  
  // Symmetrize
  CovMatPostMean = SymmetrizeMatrix(CovMatPostMean);
  CovMatLower    = SymmetrizeMatrix(CovMatLower);
  CovMatUpper    = SymmetrizeMatrix(CovMatUpper);
  
  return Rcpp::List::create(
    Rcpp::Named("PostMeanMatrix") = CovMatPostMean,
    Rcpp::Named("LowerQuantileMatrix") = CovMatLower,
    Rcpp::Named("UpperQuantileMatrix") = CovMatUpper
  );
}

