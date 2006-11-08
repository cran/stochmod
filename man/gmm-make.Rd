\name{gmm-make}
\alias{GMM.make}
\title{Gaussian Mixture Models}
\description{Generate an instance of GMM parameters.}
\usage{
GMM.make( K, p, mu=rep(0, p), sigma=diag(p) )
}
\arguments{
  \item{K}{Number of components}
  \item{p}{Dimensionality}
  \item{mu}{General center of the data being modeled}
  \item{sigma}{General orientation of the data being modeled}
}
\details{
  Component centers are generated from a multivariate Normal
  distribution with the provided mean mu and covariance matrix sigma.
  Component covariance matrices are all set to the provided sigma.
  Prior distribution is uniform across all states.
}
\value{
  A Gaussian Mixture Model defined by:
  \item{mu}{[K x p] matrix of component means}
  \item{sigma}{[K x p x p] array of component covariance matrices}
  \item{pi}{[K x 1] vector of mixture coefficients}
}
\author{Artem Sokolov \email{Artem.Sokolov@gmail.com}}
\keyword{models}
