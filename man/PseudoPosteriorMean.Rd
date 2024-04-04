\name{PseudoPosteriorMean}
\alias{PseudoPosteriorMean}
\title{FABLE Pseudo-Posterior Mean}

\description{
  This function provides the analytical CC-FABLE pseudo-posterior mean.
}

\usage{
PseudoPosteriorMean(Y,
                    gamma0 = 1,
                    delta0sq = 1,
                    maxProp = 0.5)
}

\arguments{
  \item{Y}{n by p matrix of observed data.}
  \item{gamma0}{Variance inverse-gamma shape hyperparameter.}
  \item{delta0sq}{Variance inverse-gamma mean hyperparameter.}
  \item{maxProp}{Proportion of cumulative singular value explained. Used to find upper bound for rank estimation.}
}
  

\details{
  A list is returned with the following components.
}

\value{
  \item{FABLEPostMean}{Pseudo-posterior mean from FABLE/CC-FABLE.}
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
