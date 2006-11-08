\name{hmm-obs}
\alias{HMM.obs}
\title{Hidden Markov Models}
\description{Computes observation probabilities}
\usage{
HMM.obs( x, hmm )
}
\arguments{
  \item{x}{[N x p] matrix of N samples in p dimensions}
  \item{hmm}{An HMM}
}
\details{
  Computes observation probability for each (sample, state) pair
}
\value{
  \item{B}{[N x K] matrix of probabilities}
  \item{outliers}{indices of points that had zero probability for all states}
}
\author{Artem Sokolov \email{Artem.Sokolov@gmail.com}}
\keyword{models}
