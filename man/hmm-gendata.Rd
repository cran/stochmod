\name{hmm-gendata}
\alias{HMM.genData}
\title{Hidden Markov Models}
\description{Generate data from a provided model.}
\usage{
HMM.genData( N, hmm )
}
\arguments{
  \item{N}{Number of samples to generate}
  \item{hmm}{The model of dimensionality p}
}
\details{
  Hidden state value is first generated from prior distribution (first
  sample) or transition matrix (second sample and on). The state mean
  and covariance matrix are then used to generate the actual sample.
}
\value{
  [N x p] matrix of generated data samples
}
\author{Artem Sokolov \email{Artem.Sokolov@gmail.com}}
\keyword{models}
