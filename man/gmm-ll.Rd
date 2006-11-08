\name{gmm-ll}
\alias{GMM.ll}
\title{Gaussian Mixture Models}
\description{Computes log-likelihood of new observation sequences given
  a previously trained model}
\usage{
GMM.ll( x, gmm )
}
\arguments{
  \item{x}{A matrix of a list of matrices containing one or more
    observation sequences}
  \item{gmm}{A Gaussian mixture model computed by, e.g., GMM.learn}
}
\value{
  If x is a single matrix, the function returns the log-likelihood of
  that observation sequence as a single value. If x is a list of
  observation sequences, the log-likelihood values are computed for each
  sequence and returned together in a list.
}
\author{Artem Sokolov \email{Artem.Sokolov@gmail.com}}
\keyword{models}

