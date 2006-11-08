\name{gmm-resp}
\alias{GMM.resp}
\title{Gaussian Mixture Models}
\description{Compute component responsibilities and posterior
  probabilities}
\usage{
GMM.resp( x, gmm )
}
\arguments{
  \item{x}{[N x p] matrix of N samples in p dimensions}
  \item{gmm}{Gaussian mixture model for K components with dimensionality p}
}
\value{
  \item{resp}{[N x K] matrix of component responsibilties}
  \item{LL}{[N x 1] matrix of log-likelihood values for each sample}
}
\author{Artem Sokolov \email{Artem.Sokolov@gmail.com}}
\keyword{models}

