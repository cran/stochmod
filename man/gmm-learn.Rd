\name{gmm-learn}
\alias{GMM.learn}
\title{Gaussian Mixture Models}
\description{Expectation Maximization algorithm for Gaussian Mixture
  Models}
\usage{
GMM.learn( xL, K, vL = NULL, gmm.init = NULL, cov.reg = 0.0, tol =
1e-03, LLstop = Inf, min.iter = 3, max.iter = Inf )
}
\arguments{
  \item{xL}{Either a matrix or a list of matrices containing
    training observation sequences, with one sample per row}
  \item{K}{Desired number of components}
  \item{vL}{Either a matrix or a list of matrices containing
    validation observation sequences, with one sample per row}
  \item{gmm.init}{Optional initial model, can be partially specified}
  \item{cov.reg}{Covariance matrix regularization (towards identity),
    value must be in [0, 1]}
  \item{tol}{Stopping criterion: relative tolerance on the
    log-likelihood}
  \item{LLstop}{Stopping criterion: hard bound on the log-likelihood
    value}
  \item{min.iter}{At least this number of EM iterations is preformed
    before validation and tolerance stopping criteria are triggered}
  \item{max.iter}{Stoppint criterion: maximum number of iterations}
}
\details{
  Learns a maximum likelihood GMM given the data
}
\value{
  A Gaussian Mixture Model defined by:
  \item{mu}{[K x p] matrix of component means}
  \item{sigma}{[K x p x p] array of component covariance matrices}
  \item{pi}{[K x 1] vector of mixture coefficients}
}
\author{Artem Sokolov (\email{Artem.Sokolov@gmail.com})}
\keyword{models}

  
