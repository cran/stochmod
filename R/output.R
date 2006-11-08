## output.R - output routines
##
## by Artem Sokolov

## Version of cat that displays only if debug.ouput option is above (or equal to) a certain level
d.cat <- function( lvl, ... )
  {
    if( SM.getOption("debug.output") >= lvl )
      cat( ... )
  }

## Version of print that displays only if debug.ouput option is above (or equal to) a certain level
d.print <- function( lvl, ... )
  {
    if( SM.getOption("debug.output") >= lvl )
      print( ... )
  }
