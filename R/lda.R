## lda.R - Linear Discriminant Analysis Classifier
##
## by Artem Sokolov
##

## For all functions:
##
## (N) - Number of samples
## (K) - Number of classes
## (p) - Dimensionality

## Trains an LDA classifier
## x - N x p data matrix of N samples in p dimensions
## y - N x 1 vector of labels
## cov.reg - covariance matrix regularization (towards identity) parameter, must be in [0,1]
##
## Output:
## labels  - vector of class labels
## priors  - K x 1 vector of priors, estimated as fraction of
##		points from each class
## means   - K x p matrix of means approximated from the data
## covmat  - the common p x p covariance matrix
## weights - [K x (p+1)] matrix of weights and the bias term
##		for each of the K classes
LDA.train <- function( x, y, cov.reg=0.0 )
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
    x.c <- x

    if( SM.getOption( "pure.R" ) == FALSE )
      {
        ## Reindex the labels
        y.i <- rep( 0, length(y));
        for( i in 1:length(y) ) y.i[i] <- which( y[i] == clsf$labels )

        res <- .C( "LDAtrain", as.numeric(x), as.integer(y.i), as.integer(N), as.integer(p), as.integer(K),
                  as.numeric( cov.reg ),
                  as.integer( SM.getOption( "debug.output" ) ), priors=numeric(K), means=numeric(K*p),
                  covmat=numeric(p*p), weights=numeric(K*(p+1)), status=integer(1),
                  DUP=FALSE )

        if( res$status == 1 )
          stop( "The common covariance matrix is singular" )
        else
          return( list(
                       labels = clsf$labels,
                       priors = res$priors,
                       means = matrix( res$means, K, p ),
                       covmat = matrix( res$covmat, p, p ),
                       weights = matrix( res$weights, K, p+1 )
                       )
                 )
      }
    else
      {
        d.cat( 1, "Currently running R version of the LDA.train code\n" );
        
        ## Approximate priors, means,
        ##  and the common covariance matrix
        clsf$priors <- rep( 0, K )
        clsf$means <- matrix( 0, K, p )
        for( k in 1:K )
          {
            mask <- y==clsf$labels[k]
            z <- x[ mask, ,drop=FALSE ]
            N.k <- nrow(z)
            clsf$priors[k] <- N.k / N
            clsf$means[k,] <- apply( z, 2, mean )

            x.c[mask,] <- x[mask,] - matrix( clsf$means[k,], N.k, p, byrow=TRUE)
          }

        clsf$covmat <- (t(x.c)%*%x.c) / (N-K)
        if( cov.reg > 0 )
          clsf$covmat <- cov.reg * diag(p) + (1-cov.reg) * clsf$covmat
        icov <- solve( clsf$covmat )

        ## Compute weights and bias terms
        clsf$weights <- matrix( 0, K, p+1 )
        for( k in 1:K )
          {
            mu <- clsf$means[k,]		# Column vector, because drop=TRUE implicitly
            clsf$weights[k,1:p] <- icov %*% mu
            clsf$weights[k,p+1] <- -0.5 * (t(mu) %*% icov %*% mu) + log(clsf$priors[k])
          }

        clsf
      }
  }

## Applies previously trained LDA classifier to classify rows of x
##
## x - N x p data matrix of N samples in p dimensions
## clsf - classifier produced by LDA.train()
##
## Output:
## column vector of labels for each row of x
LDA.test <- function( x, clsf )
  {
    x <- as.matrix(x)

    ## Append a column of 1's
    x <- cbind( x, rep( 1, nrow(x) ) )

    ## Compute deltas
    deltas <- clsf$weights %*% t(x)

    ## Determine labels
    i <- apply( deltas, 2, which.max )
    clsf$labels[i]
  }
