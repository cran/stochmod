\name{qda-test}
\alias{QDA.test}
\title{Quadratic Discriminant Analysis}
\description{Test method for QDA.}
\usage{
QDA.test( x, clsf )
}
\arguments{
  \item{x}{N x p data matrix of N samples in p dimensions}
  \item{clsf}{QDA classifier object produced by QDA.train}
}
\details{
  Classifies new test data using QDA discriminant functions computed
  during training.
}
\value{
  N x 1 vector of predicted labels
}
\author{Artem Sokolov \email{Artem.Sokolov@gmail.com}}
\seealso{\link{QDA.train}}
\keyword{models}
