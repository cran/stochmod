\name{lda-train}
\alias{LDA.train}
\title{Linear Discriminant Analysis}
\description{Training method for LDA.}
\usage{
LDA.train( x, y, cov.reg = 0.0 )
}
\arguments{
  \item{x}{N x p data matrix of N samples in p dimensions}
  \item{y}{N x 1 vector of labels}
  \item{cov.reg}{Covariance matrix regularization (towards identity),
    value must be in [0, 1]}
}
\details{
  Models each class as a single (multivariate) Gaussian and uses a
  single common covariance matrix across all classes. Computes the
  distribution parameters, the Bayesian class priors, and the
  discriminant functions. LDA is insensitive to temporal structure of
  the data and, therefore, only needs to work with a single observation
  sequence. This, in turn, requires a label for each sample.
}
\value{
  An LDA classifier defined by:
  \item{labels}{Vector of unique class labels}
  \item{priors}{K x 1 vector of priors, estimated as fraction of points
    from each class}
  \item{means}{K x p matrix of means approximated from the data}
  \item{covmat}{The common p x p covariance matrix}
  \item{weights}{K x (p+1) matrix of weights and the bias term for each
    of the K discriminant functions}
}
\author{Artem Sokolov \email{Artem.Sokolov@gmail.com}}
\keyword{models}
