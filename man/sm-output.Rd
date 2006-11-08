\name{sm-output}
\alias{d.cat}
\alias{d.print}
\title{Output of Debug Information}
\description{Versions of cat and print that depend on debug.output
  option}
\usage{
d.cat( lvl, ... )
d.print( lvl, ... )
}
\arguments{
  \item{lvl}{When debug.output option is above (or equal to) this level,
    output is displayed}
  \item{...}{Arguments that get passed to the actual output function}
}
\value{
  'd.cat' returns the same value as cat
  'd.print' return the same value as print
}
\author{Artem Sokolov (\email{Artem.Sokolov@gmail.com})}
\keyword{IO}
