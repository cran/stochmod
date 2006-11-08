\name{lda-test}
\alias{LDA.test}
\title{Linear Discriminant Analysis}
\description{Test method for LDA.}
\usage{
LDA.test( x, clsf )
}
\arguments{
  \item{x}{N x p data matrix of N samples in p dimensions}
  \item{clsf}{LDA classifier object produced by LDA.train}
}
\details{
  Classifies new test data using LDA discriminant functions computed
  during training.
}
\value{
  N x 1 vector of predicted labels
}
\author{Artem Sokolov \email{Artem.Sokolov@gmail.com}}
\seealso{\link{LDA.train}}
\keyword{models}
