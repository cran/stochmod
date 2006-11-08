\name{sm-options}
\alias{SM.setOption}
\alias{SM.getOption}
\alias{SM.getOptions}
\title{Manipulation of Stochastic Modeling Options}
\description{These functions provide a means to get/set options that
  affect behavior of many modeling routines}
\usage{
SM.setOption( name, value = TRUE )
SM.getOption( name )
SM.getOptions()
}
\arguments{
  \item{name}{
    Argument name. Currently supported options are:
    \item{debug.output:}{Amount of debug output. (integer)}
    \item{pure.R:}{whether exclusively R code is to be used. (logical)}
  }
  \item{value}{Agument value}
}
\value{
  'SM.setOption' returns the previous argument value.
  'SM.getOption' returns the current value.
  'SM.getOptions' returns a list of current option assignments.
}
\author{Artem Sokolov \email{Artem.Sokolov@gmail.com}}
\keyword{utilities}

  
