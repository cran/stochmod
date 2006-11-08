\name{hmm-make}
\alias{HMM.make}
\title{Hidden Markov Models}
\description{Generate an instance of HMM parameters.}
\usage{
HMM.make( K, p, mu=rep(0, p), sigma=diag(p) )
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
  Rows of the transition matrix are random multinomial distributions.
}
\value{
  A Hidden Markov Model defined by:
  \item{mu}{[K x p] matrix of component means}
  \item{sigma}{[K x p x p] array of component covariance matrices}
  \item{pi}{[K x 1] vector of mixture coefficients}
  \item{A}{[K x K] state transition matrix with element (i,j) referring
    to transition from state i to state j}
}
\author{Artem Sokolov \email{Artem.Sokolov@gmail.com}}
\keyword{models}
