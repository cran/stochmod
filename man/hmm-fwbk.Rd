\name{hmm-fwbk}
\alias{HMM.fwbk}
\title{Hidden Markov Models}
\description{Forward-Backward procedure}
\usage{
HMM.fwbk( x, hmm )
}
\arguments{
  \item{x}{[N x p] matrix of N samples in p dimensions}
  \item{hmm}{An HMM}
}
\details{
  The E-step of the EM algorithm for Hidden Markov Models
}
\value{
  \item{alpha}{[N x K] matrix of forward variable values}
  \item{beta}{[N x K] matrix of backward variable values}
  \item{c}{[N x 1] vector of scaling coefficients}
  \item{gamma}{[N x K] matrix of transition probabilities from state i: gamma[t,i] = P( Qt = i | O,hmmp )}
  \item{xi}{[N-1 x K x K] array of transition probabilties from state i to state j: xi[t,i,j] = P( Qt = i, Qt+1 = j | O,hmmp )}
}
\author{Artem Sokolov \email{Artem.Sokolov@gmail.com}}
\keyword{models}
