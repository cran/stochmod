## gmm.R - Gaussian Mixture Models
##
## by Artem Sokolov

## For all functions:
##
## (N) - Number of samples
## (K) - Number of components
## (p) - Dimensionality
##
## Gaussian mixture model parameters gmm must contain:
##	mu    - [K x p] matrix of component means
##	sigma - [K x p x p] array of covariance matrices for each component
##	pi    - [K x 1] vector of mixture coefficients (has to sum up to 1)


## Generates a gmm instance
##
## Inputs:
##	K - number of components
##	p - dimensionality
##	mu, sigma - general orientation of the model
##
## Output:
##	the model
GMM.make <- function( K, p, mu=rep(0, p), sigma=diag(p) )
  {
    gmm <- list()
    gmm$pi <- rep( 1/K, K )
    gmm$mu <- rmvnorm( K, mu, sigma )
    gmm$sigma <- array( 0, c(K,p,p) )
    for( i in 1:K )
      gmm$sigma[i,,] <- sigma
    gmm
  }

## Generates data from a provided GMM
##
## Inputs:
##	N - number of samples
##	gmm - the model (with dimensionality p)
##
## Output:
##	[N x p] data matrix
GMM.genData <- function( N, gmm )
  {
    K <- nrow( gmm$mu )
    p <- ncol( gmm$mu )
    x <- matrix( 0, N, p )
    j <- sample( 1:K, N, TRUE, gmm$pi )

    for( i in 1:K )
      {
        k <- which( j == i )
        n <- length( k )
        if( n > 0 )
          x[ k, ] <- rmvnorm( n, gmm$mu[i,], gmm$sigma[i,,] )
      }

    x
  }

## Computes component responsibilities and posterior probabilities
##
## Inputs:
## 	x - [N x p] data matrix
## 	gmm - the model
##
## Outputs:
##	resp - [N x K] matrix of responsibilities
##	LL   - [N x 1] vector of log likelihoods
GMM.resp <- function( x, gmm )
  {
    if( is.list(x) == TRUE )
      stop( "GMM.resp: x must be a single observation sequence" )
    
    x <- as.matrix(x)
    gmm$pi <- gmm$pi/sum(gmm$pi)

    N <- nrow(x)
    p <- ncol(x)
    K <- length(gmm$pi)

    if( SM.getOption( "pure.R" ) == FALSE )
      {
        res <- .C( "GMMrespR", as.numeric(x), as.numeric(gmm$mu), as.numeric(gmm$sigma), as.numeric(gmm$pi),
                  as.integer(N), as.integer(p), as.integer(K), resp=numeric(N*K), LL=numeric(N),
                  status=integer(1), info=integer(N), DUP=FALSE )

        if( res$status == 1 )
          stop( paste("GMM.resp: Covariance matrix is singular for component", res$info[1]) )
        
        else
          return( list(
                       resp=matrix(res$resp, N, K),
                       LL=matrix(res$LL, N, 1),
                       outliers=which(res$info == 1)
                       ) )
      }
    else
      {
        ## Compute posterior probabilties using Bayes rule
        post.prob <- matrix( 0, N, K )
        for( i in 1:K )
          post.prob[,i] <- tryCatch(gmm$pi[i] * dmvnorm(x, gmm$mu[i,,drop=FALSE], as.matrix(gmm$sigma[i,,])),
                                    error=function(e)
                                    stop(paste("GMM.resp: Covariance matrix is singular for component", i))
                                    )

        ## Compute responsibilities by scaling posteriors
        p.x <- apply( post.prob, 1, sum )
        resp <- apply( post.prob, 2, "/", p.x )

        ## Return responsibilities and log likelihood of the posteriors
        return( list( resp=as.matrix(resp), LL=as.matrix(log(p.x)), outliers=which(p.x==0) ) )
      }
  }

## Computes log-likelihood of new observation sequences given a previously trained model
##
## Inputs:
##	x - Either a single observation sequence matrix or a list of multiple observation sequences
##	gmm - the model
##
## Output:
##	If x is a single observation sequence, returns log-likelihood of that sequence
##	If x is a list of multiple observation sequences, returns a list with each element being log-likelihood
##		of the corresponding sequence in x
GMM.ll <- function( x, gmm )
  {
    if( is.list( x ) == TRUE )
      lapply( x, sys.function(), gmm )
    else
      sum( GMM.resp( x, gmm )$LL )
  }

## Expectation Maximization for multiple observation sequences
##
## (Ni) - Number of samples in training sequence i
## (Mi) - Number of samples in validation sequence i
##
## Inputs:
## xL - list of [Ni x p] data matrices
## K - Number of components
## vL - list of [Mi x p] validation data matrices
## gmm.init - optional initialization point (can be partially specified)
## cov.reg - covariance matrix regularization (towards identity) parameter, must be in [0,1]
##
## Stopping Criteria:
## tol		- relative tolerance for convergence
## LLstop	- stop iterating when LL exceeds this
## max.iter	- cap on the number of iterations
##
## Outputs:
## gmm - Gaussian Mixture model
GMM.learn <- function( xL, K, vL=NULL, gmm.init=NULL, cov.reg=0.0, tol=1e-03, LLstop=Inf, min.iter=3, max.iter=Inf )
  {
    ## Argument type verification
    if( is.list(xL) == FALSE )
      xL <- list( xL )
    if( is.null(vL) == FALSE && is.list(vL) == FALSE )
      vL <- list( vL )

    ## Further argument verification
    if( cov.reg < 0 ) cov.reg <- 0
    if( cov.reg > 1 ) cov.reg <- 1
    
    nseq <- length(xL)
    nval <- length(vL)

    d.cat( 1, "Learning a Gaussian Mixture Model with", K,"component(s) for", nseq, "observation sequences\n" )
    
    ## Concatenate all training sequences together
    x <- NULL
    for( i in 1:nseq )
      x <- rbind( x, as.matrix(xL[[i]]) )
    N <- nrow(x)

    ## Concatenate all validation sequences together
    v <- NULL
    for( i in seq( 1, nval, length=nval ) )
      v <- rbind( v, as.matrix(vL[[i]]) )
    M <- nrow(v)
        
    p <- ncol( x )

    ## Special case
    if( K == 1 )
      return( list( pi=1, mu = matrix( apply( x, 2, mean ), nrow=1 ), sigma = array( cov(x), c(1,p,p) ) ) )

    ## Initialize GMM
    gmm <- GMM.make( K, p, apply( x, 2, mean ), cov(x) )
    if( is.null( gmm.init$mu ) == FALSE )	gmm$mu <- gmm.init$mu
    if( is.null( gmm.init$sigma ) == FALSE )	gmm$sigma <- gmm.init$sigma
    if( is.null( gmm.init$pi ) == FALSE )	gmm$pi <- gmm.init$pi

    ########################
    ## Attempt to run C code
    if( SM.getOption("pure.R") == FALSE )
      {
        ## There is a weird issue that occurs when C code is the first thing ran after R startup
        ## Apparently random number generation using rnorm() in Rmath.h always returns the same number
        ## The problem seems to be alleviated when the first call to the function is made from R
        ## For this purpose:
        rnorm( 1, 0, 1 )

        ## Attempt to run C code
        if( is.finite( max.iter ) == FALSE ) max.iter <- -1
        res <- .C( "GMMlearnR",
                  as.numeric(x), as.numeric(v), as.integer(N), as.integer(M), as.integer(K), as.integer(p),
                  as.numeric(gmm$mu), as.numeric(gmm$sigma), as.numeric(gmm$pi),
                  as.numeric( cov.reg ),
                  as.numeric(tol), as.numeric(LLstop),
                  as.integer(min.iter), as.integer(max.iter),
                  as.integer(SM.getOption("debug.output")),
                  pi = numeric( K ), mu = numeric( K*p ), sigma = numeric( K*p*p ), status = integer( 1 ),
                  NAOK=TRUE, DUP=FALSE )

        if( res$status == 1 )
          stop( "Covariance matrix is singular for one of the clusters\n" )
        else if( res$status == 2 )
          stop( "Some points were not assigned to a cluster during the E-Step\n" )
        else if( res$status == 3 )
          stop( "No points were assigned to one of the clusters\n" )
        else
          return( list(
                       pi = res$pi,
                       mu = matrix( res$mu, K, p ),
                       sigma = array( res$sigma, c( K, p, p ) )
                       ) )
      }

    #################################
    ## Revert to the R Implementation
    d.cat( 1, "Currently running R version of the GMM.learn code\n" );

    gmm.prev <- list()

    LLval <- 0
    LLval.prev <- 0
    LLbase <- 0
    LLprev <- 0
    iter <- 0
    while( TRUE )
      {
        ## Validation log-likelihood stopping criterion
        if( nval > 0 )
          {
            LLval <- sum( GMM.resp( v, gmm )$LL )
            d.cat( 3, "Validation log-likelihood: ", LLval, "\n" )
            if( iter <= min.iter || LLval >= LLval.prev )
              {
                LLval.prev <- LLval
                gmm.prev <- gmm
              }
            else
              {
                d.cat( 1, "Final log likelihod of training data after", iter-1, "iterations:", LLprev, "\n" )
                d.cat( 1, "Validation log likelihood:", LLval.prev, "\n" )
                return( gmm.prev )
              }
          }
            
        ##
        ## E STEP - compute responsibilities
        ##
        pp <- GMM.resp( x, gmm )

        if( length( pp$outliers ) > 0 )
          stop( "Some points were not assigned to a cluster during the E-Step\n" )
        
        resp <- pp$resp
        LL <- sum( pp$LL )

        d.cat( if( (iter%%10)==0 ) 2 else 3,
              "Log likelihood at iteration", iter, ":", LL, "\n" )
        
        ## Verify stopping criteria
        if( iter >= max.iter ) break
        if( LL > LLstop ) break
        if( is.infinite(LL) ) break
        if( LL == LLprev ) break
        
        if( iter <= (min.iter-1) ) LLbase <- LL
        else if( (LL-LLbase) < (1+tol)*(LLprev-LLbase) ) break

        ##
        ## M STEP - Compute the new parameters
        ##
        for( i in 1:K )
          {
            gamma.i <- matrix( resp[,i], N, p )
            norm.fac <- sum( resp[,i] )

            if( norm.fac == 0 )
              stop( paste( "No points were assigned to cluster", i, "during the E-Step\n" ) )
            
            gmm$mu[i,] <- apply( (gamma.i*x), 2, sum ) / norm.fac
            x.shifted <- t( t( x ) - gmm$mu[i,] )
            gmm$sigma[i,,] <- (t( x.shifted * gamma.i ) %*% x.shifted) / norm.fac

            if( cov.reg > 0 )
              gmm$sigma[i,,] <- cov.reg * diag(p) + (1-cov.reg) * gmm$sigma[i,,]
            
            gmm$pi[i] <- norm.fac / N
          }
        gmm$pi <- gmm$pi / sum(gmm$pi)

        LLprev <- LL
        iter <- iter+1
      }

    d.cat( 1, "Final log likelihod after", iter, "iterations:", LL, "\n" )
    gmm
  }

