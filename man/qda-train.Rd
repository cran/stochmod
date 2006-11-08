\name{qda-train}
\alias{QDA.train}
\title{Quadratic Discriminant Analysis}
\description{Training method for QDA.}
\usage{
QDA.train( x, y, cov.reg = 0.0 )
}
\arguments{
  \item{x}{N x p data matrix of N samples in p dimensions}
  \item{y}{N x 1 vector of labels}
  \item{cov.reg}{Covariance matrix regularization (towards identity),
    value must be in [0, 1]}
}
\details{
  Models each class as a single (multivariate) Gaussian. Relaxes the
  common covariance matrix constraint of LDA. Computes the
  distribution parameters, the Bayesian class priors, and the
  discriminant functions. QDA is insensitive to temporal structure of
  the data and, therefore, only needs to work with a single observation
  sequence. This, in turn, requires a label for each sample.
}
\value{
  An QDA classifier defined by:
  \item{labels}{Vector of unique class labels}
  \item{priors}{K x 1 vector of priors, estimated as fraction of points
    from each class}
  \item{means}{K x p matrix of means approximated from the data}
  \item{covmats}{K x p x p array of covariance matrices esimated from
    the data}
  \item{icovmats}{K x p x p array of inverse covariance matrices}
  \item{bias}{K x 1 vector of bias terms for discriminant function computations}
}
\author{Artem Sokolov \email{Artem.Sokolov@gmail.com}}
\keyword{models}
