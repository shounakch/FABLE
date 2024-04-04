\name{PseudoPosteriorSampler}
\alias{PseudoPosteriorSampler}
\title{FABLE Pseudo-Posterior Sampling}

\description{
  This function implements the CC-FABLE pseudo-posterior sampling, providing MC samples of the factor loadings matrix and the error variance.
}

\usage{
PseudoPosteriorSampler(Y,
                       gamma0 = 1,
                       delta0sq = 1,
                       maxProp = 0.5,
                       MC = 1000)
}

\arguments{
  \item{Y}{n by p matrix of observed data.}
  \item{gamma0}{Variance inverse-gamma shape hyperparameter.}
  \item{delta0sq}{Variance inverse-gamma mean hyperparameter.}
  \item{maxProp}{Proportion of cumulative singular value explained. Used to find upper bound for rank estimation.}
  \item{MC}{Number of Monte Carlo samples.}
}
  

\details{
  A list is returned with the following components.
}

\value{
  \item{CCFABLESamples}{Samples of factor loadings matrix, error variances, and other related information.}
  \item{FABLEHyperParameters}{List containing FABLE hyperparameters.}
  \item{svdY}{List containing components of singular value decomposition (SVD) of Y.}
  \item{estRank}{Estimated rank of the factor model.}
  \item{varInflation}{Estimated variance inflation factor in CC-FABLE.}
  \item{runTime}{Runtime for the computation.}
}

\references{
  [to be added]
}

\author{Shounak Chattopadhyay}

\examples{
[to be added]
}
