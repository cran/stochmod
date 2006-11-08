## Load hooks for stochmod package
##
## by Artem Sokolov

sm.options <- new.env()

.onLoad <- function( libname, pkgname )
  {
    assign( "debug.output", 1, env=sm.options )
    assign( "pure.R", FALSE, env=sm.options )
  }
