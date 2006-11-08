## lag.R - routines for lagging data
##
## by Artem Sokolov

## x - data matrix or list of data matrices
## k - number of lags (>= 0)
FEAT.lagData <- function( x, k=1 )
  {
    if( is.list(x) )
      return( lapply( x, sys.function(), k ) )

    if( k < 1 ) return( x );
    if( k >= nrow(x) ) stop( "Impossible to apply this many lags" )
    
    n <- nrow(x)
    p <- ncol(x)
    res <- matrix( 0, n-k, (k+1)*p )
    for( i in 0:k )
      res[,(i*p+1):((i+1)*p)] <- x[(i+1):(n-k+i),]
    res
  }
