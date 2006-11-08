\name{hmm-learn}
\alias{HMM.learn}
\title{Hidden Markov Models}
\description{Maximum Likelihood learning via EM}
\usage{
HMM.learn( xL, K, vL=NULL, hmm.init=NULL, tol=1e-03, LLstop=Inf,
min.iter=3, max.iter=Inf )
}
\arguments{
  \item{xL}{Either a matrix or a list of matrices containing
    training observation sequences, with one sample per row}
  \item{K}{Desired number of states}
  \item{vL}{Either a matrix or a list of matrices containing
    validation observation sequences, with one sample per row}
  \item{hmm.init}{Optional initial model, can be partially specified}
  \item{tol}{Stopping criterion: relative tolerance on the
    log-likelihood}
  \item{LLstop}{Stopping criterion: hard bound on the log-likelihood value}
  \item{min.iter}{At least this number of EM iterations is preformed
    before validation and tolerance stopping criteria are triggered}
  \item{max.iter}{Stoppint criterion: maximum number of iterations}
}
\details{
  Learns a maximum likelihood HMM given the data
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
