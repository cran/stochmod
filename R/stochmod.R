SM.getOptions <- function()
  {
    res <- list()
    for( op in ls( env=sm.options ) )
      res[[op]] <- SM.getOption( op )
    res
  }

SM.getOption <- function( name )
  {
    tryCatch(
             get( name, env=sm.options, inherits=FALSE ),
             error=function(e)
             stop( paste( "There is no such option:", name ) )
             )
  }

SM.setOption <- function( name, value=TRUE )
  {
    oldval <- tryCatch( SM.getOption( name ), error=function(e) NULL )
    assign( name, value, env=sm.options )
    invisible(oldval)
  }
