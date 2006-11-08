\name{gmm-gendata}
\alias{GMM.genData}
\title{Gaussian Mixture Models}
\description{Generate data from a provided model.}
\usage{
GMM.genData( N, gmm )
}
\arguments{
  \item{N}{Number of samples to generate}
  \item{gmm}{The model of dimensionality p}
}
\details{
  Mixture coefficients are used to decide which component each datapoint
  is generated from. Component mean and covariance matrix is then
  employed to generate the actual sample.
}
\value{
  [N x p] matrix of generated data samples
}
\author{Artem Sokolov \email{Artem.Sokolov@gmail.com}}
\keyword{models}
