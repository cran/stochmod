\name{feat-lagdata}
\alias{FEAT.lagData}
\title{Data Lagging}
\description{Time-lagged embedding of a time series}
\usage{
FEAT.lagData( x, k=1 )
}
\arguments{
  \item{x}{Either a matrix or a list of matrices containing
    observation sequences, with one sample per row}
  \item{k}{Number of lags}
}
\details{
  Time-delayed or time-lagged embedding appends a time-series to a
  lagged version of itself, thereby increasing dimensionality from p to
  k*p and reducing the number of samples from n to n-k.
}
\value{
  [n-k by k*p] matrix of lagged time series
}
\author{Artem Sokolov (\email{Artem.Sokolov@gmail.com})}
\keyword{models}

  
