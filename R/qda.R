## qda.R - Quadratic Discriminant Analysis Classifier
##
## by Artem Sokolov

## For all functions:
##
## (N) - Number of samples
## (K) - Number of classes
## (p) - Dimensionality

## Trains a QDA classifier
## x - N x p data matrix of N samples in p dimensions
## y - N x 1 vector of labels
## cov.reg - covariance matrix regularization (towards identity) parameter, must be in [0,1]
##
## Output:
## labels   - vector of class labels
## priors   - K x 1 vector of priors, estimated as fraction of
##		points from each class
## means    - K x p matrix of means approximated from the data
## covmats  - [K x p x p] array of covariance matrices
## icovmats - [K x p x p] array of inverse covariance matrices
## bias	    - K x 1 vector of bias terms for discriminant function computations
QDA.train <- function( x, y, cov.reg=0.0 )
  {
    x <- as.matrix(x)
    clsf <- list()

    ## Argument verification
    if( cov.reg < 0 ) cov.reg <- 0
    if( cov.reg > 1 ) cov.reg <- 1
    
    ## Retrieve labels
    clsf$labels <- unique(y)
    K <- length( clsf$labels )
    N <- nrow(x)
    p <- ncol(x)

    if( SM.getOption( "pure.R" ) == FALSE )
      {
        ## Reindex the labels
        y.i <- rep( 0, length(y));
        for( i in 1:length(y) ) y.i[i] <- which( y[i] == clsf$labels )

        res <- .C( "QDAtrain", as.numeric(x), as.integer( y.i ), as.integer( N ), as.integer( p ),
                  as.integer( K ), as.numeric( cov.reg ),
                  as.integer( SM.getOption( "debug.output" ) ), priors = numeric( K ),
                  means = numeric( K*p ), covmats = numeric( K*p*p ), icovmats = numeric( K*p*p ),
                  bias = numeric( K ), status = integer(1), info = integer(1) )
        
        if( res$info != 0 )
          {
            msg <- paste( "The covariance matrix for class", res$info, "is singular\n" )
            stop( msg )
          }
        else
          {
            res <- list(
                        labels = clsf$labels,
                        priors = res$priors,
                        means = matrix( res$mean, K, p ),
                        covmats = array( res$covmats, c( K, p, p ) ),
                        icovmats = array( res$icovmats, c( K, p, p ) ),
                        bias = res$bias
                        )

            test.na <- function( x ) {any(is.na(c(x)))}
            if( any(unlist(lapply( res, test.na ))) == TRUE )
              stop( "Computation in QDA.train are numerically unstable\n" )
            
            return( res )
          }
      }
    else
      {
        d.cat( 1, "Currently running R version of the QDA.train code\n" )
        
        ## Approximate priors, means, bias
        ##  and covariance matrices
        clsf$priors <- rep( 0, K )
        clsf$means <- matrix( 0, K, p )
        clsf$covmats <- array( 0, c(K,p,p) )
        clsf$icovmats <- array( 0, c(K,p,p) )
        for( k in 1:K )
          {
            mask <- y == clsf$labels[k]
            z <- x[ mask, ,drop=FALSE ]
            clsf$priors[k] <- nrow(z) / N
            clsf$means[k,] <- apply( z, 2, mean )
            clsf$covmats[k,,] <- cov.reg * diag(p) + (1-cov.reg) * cov(z)
            clsf$bias[k] <- log(clsf$priors[k]) - 0.5 * sum( log( svd(clsf$covmats[k,,])$d ) )
            clsf$icovmats[k,,] <- solve(clsf$covmats[k,,])
          }
        
        clsf
      }
  }

## Applies previously trained QDA classifier to classify rows of x
##
## x - N x p data matrix of N samples in p dimensions
## clsf - classifier produced by QDA.train()
##
## Output:
## column vector of labels for each row of x
QDA.test <- function( x, clsf )
  {
    x <- as.matrix(x)

    N <- nrow( x )
    p <- ncol( x )
    K <- length(clsf$bias)

    if( SM.getOption( "pure.R" ) == FALSE )
      {
        if( is.null( clsf$means ) == TRUE )
          stop( "Invalid means field in the classifier" )
        if( is.null( clsf$icovmats ) == TRUE )
          stop( "Invalid icovmats field in the classifier" )
        if( is.null( clsf$bias ) == TRUE )
          stop( "Invalid bias field in the classifier" )
          
        res <- .C( "QDAtest",
                  as.numeric(x), as.numeric( clsf$means ), as.numeric( clsf$icovmats ),
                  as.numeric( clsf$bias ),
                  as.integer( N ), as.integer( p ),
                  as.integer( K ), i = integer( N ) )

        i <- res$i
        }
      else
      {
        deltas <- matrix( 0, nrow(x), K )
        for( i in 1:nrow(x) )
          for( k in 1:K )
            {
              x.c <- x[i,] - clsf$means[k,]
              deltas[i,k] <- -0.5 * t(x.c) %*% clsf$icovmats[k,,] %*% x.c + clsf$bias[k]
            }
        
        i <- apply( deltas, 1, which.max )
      }

    clsf$labels[i]
  }

