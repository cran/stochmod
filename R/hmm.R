## hmm.R - Hidden Markov models with continuous observation distributions
##
## by Artem Sokolov
##

## For all functions:
##
## K - Number of states
## N - Length of observation sequence
## p - dimensionality
##
## HMM parameters hmm must contain:
##
## pi - [K x 1] vector of prior probabilities: pi[i] = P( Q1 = i )
## A  - [K x K] matrix of state transition probabilities: A[i,j] = P( Qt = j | Qt-1 = i )
## mu - [K x p] matrix of means for distributions associated with each state
## sigma - [K x p x p] array of covariance matrices associated with each state
##


## Generates a hmm instance
##
## Inputs:
##	K - number of components
##	p - dimensionality
##	mu, sigma - general orientation of the model
##
## Output:
##	the model
HMM.make <- function( K, p, mu=rep(0, p), sigma=diag(p) )
  {
    hmm <- GMM.make( K, p, mu, sigma )
    hmm$A <- matrix( runif( K*K ), K, K )
    for( i in 1:K )
      hmm$A[i,] <- hmm$A[i,] / sum( hmm$A[i,] )
    hmm
  }


## Generates data from a provided HMM
##
## Inputs:
##	N - number of samples
##	hmm - the model (with dimensionality p)
##
## Output:
##	[N x p] data matrix
HMM.genData <- function( N, hmm )
  {
    K <- nrow( hmm$mu )
    p <- ncol( hmm$mu )
    x <- matrix( 0, N, p )
    j <- sample( 1:K, 1, prob=hmm$pi )
    for( i in 1:N )
      {
        x[i,] <- rmvnorm( 1, hmm$mu[j,], hmm$sigma[j,,] )
        j <- sample( 1:K, 1, prob=hmm$A[j,] )
      }
    x
  }


## Computes observation probabilities
##
##
## Inputs:
## x - [N x p] observation sequence
## hmm - HMM model
##
## Output:
## B - [N x K] matrix of observation probabilities
## outliers - indices of points that weren't assigned to any state
HMM.obs <- function( x, hmm )
  {
    ## Retrieve dimensionality
    N <- nrow(x)
    p <- ncol(x)
    K <- length( hmm$pi )

    if( SM.getOption("pure.R") == FALSE )
      {
    
        ## Attempt to run C code
        res <- .C( "HMMobsR", as.numeric(x), as.numeric(hmm$pi), as.numeric(hmm$A), as.numeric(hmm$mu),
                  as.numeric(hmm$sigma), as.integer(N), as.integer(K), as.integer(p),
                  B=numeric(N*K), status=integer(1), info=integer(1), DUP=FALSE )
        if( res$status == 1 )
          stop( paste( "HMM.obs: Covariance matrix is singular for state", res$info ) )
        else
          B <- matrix( res$B, N, K )
      }
    else
      {

        ## Compute probabilities using provided parameters
        B <- matrix( 0, N, K )
        for( i in 1:K )
          B[,i] <- tryCatch( dmvnorm( x, hmm$mu[i,], as.matrix( hmm$sigma[i,,] ) ),
                            error=function(e)
                            stop(paste("HMM.obs: Covariance matrix is singular for state", i))
                            )
      }
    
    list( B = B, outliers = which( apply( B, 1, sum ) == 0 ) )
  }

## Forward-Backward procedure
##
## Inputs:
## x - [N x p] observation sequence
## hmm - HMM parameters
##
## Output:
## alpha - [N x K] matrix of forward variable values
## beta - [N x K] matrix of backward variable values
## c - [N x 1] vector of scaling coefficients
## gamma - [N x K] matrix of transitions from state i: gamma[t,i] = P( Qt = i | O,hmm )
## xi - [N-1 x K x K] matrix of transitions from state i to state j: xi[t,i,j] = P( Qt = i, Qt+1 = j | O,hmm )
HMM.fwbk <- function( x, hmm )
  {
    ## Retrieve dimensionality
    N <- nrow(x)
    p <- ncol(x)
    K <- length( hmm$pi )

    if( SM.getOption("pure.R") == FALSE )
      {
        ## Attempt to run C code
        res <- .C( "HMMfwbkR", as.numeric(x), as.numeric(hmm$pi), as.numeric(hmm$A), as.numeric(hmm$mu),
                  as.numeric(hmm$sigma), as.integer( N ), as.integer( K ), as.integer( p ),
                  alpha=numeric(N*K), beta=numeric(N*K), c=numeric(N), gamma=numeric(N*K),
                  xi=numeric( (N-1)*K*K ), status=integer(1), DUP=FALSE )

        if( res$status == 2 )
          {
            stop( "Some samples were not accounted for by any state" );
          }
        else if( res$status == 1 )
          {
            stop( "Covariance matrix for one of the states is singular" );
          }
        else
          list(
               alpha = matrix( res$alpha, N, K ),
               beta = matrix( res$beta, N, K ),
               c = res$c,
               gamma = matrix( res$gamma, N, K ),
               xi = array( res$xi, c( N-1, K, K ) )
               )
      }
    else
      {
        ## Create variables
        res <- list();
        res$alpha <- matrix( 0, N, K )
        res$beta <- matrix( 0, N, K )
        res$c <- rep( 0, N )
        res$gamma <- matrix( 0, N, K )
        res$xi <- array( 0, c(N-1,K,K) )
        
        ## Compute observation probabilities
        obs <- HMM.obs( x, hmm )
        if( length( obs$outliers ) > 0 )
          stop( "Some samples were not accounted for by any state" );
        B <- obs$B
        
        ## Compute alpha values for the first observation
        res$alpha[1,] <- hmm$pi * B[1,]
        res$c[1] <- 1 / sum( res$alpha[1,] )
        res$alpha[1,] <- res$alpha[1,] * res$c[1]

        ## Recursively compute remaining alpha values
        for( i in 2:N )
          {
            res$alpha[i,] <- (res$alpha[i-1,] %*% hmm$A) * B[i,]
            res$c[i] <- 1 / sum( res$alpha[i,] )
            res$alpha[i,] <- res$alpha[i,] * res$c[i]
          }

        ## Compute beta values for the last observation
        res$beta[N,] <- res$c[N]	# really res$c[N] * 1

        ## Recursively compute remaining beta values
        for( i in (N-1):1)
          res$beta[i,] <- res$c[i] * (hmm$A %*% (res$beta[i+1,] * B[i+1,]))

        ## Compute and normalize gamma values
        res$gamma <- res$alpha * res$beta
        for( t in 1:N )
          res$gamma[t,] <- res$gamma[t,] / sum(res$gamma[t,])
        
        ## Compute xi values
        for( t in 1:N-1 )
          for( i in 1:K )
            res$xi[t,i,] = res$alpha[t,i] * hmm$A[i,] * (B[t+1,] * res$beta[t+1,])
                                        # Normally would divide by res$xi[t,,],
                                        # but scaling takes care of it
        res
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
HMM.ll <- function( x, hmm )
  {
    if( is.list( x ) == TRUE )
      lapply( x, sys.function(), hmm )
    else
      -sum(log(HMM.fwbk(x, hmm)$c))
  }



## Constructs, trains and returns an HMM given multiple observation sequences
##
## (Ni) - Number of samples in training sequence i
## (Mi) - Number of samples in validation sequence i
##
## Inputs:
## xL - list of [Ni x p] training data matrices
## K - Number of states
## vL - list of [Mi x p] validation data matrices
## hmm.init - initial parameter estimates
##
## Stopping Criteria:
## tol		- relative tolerance for convergence
## LLstop	- stop iterating when LL exceeds this
## max.iter	- cap on the number of iterations
##
## Output:
## HMM parameters hmm
HMM.learn <- function( xL, K, vL=NULL, hmm.init=NULL , tol=1e-03, LLstop=Inf, min.iter=3, max.iter=Inf )
  {
    ## Argument verification
    if( is.list( xL ) == FALSE )
      xL <- list( xL )
    if( is.null( vL ) == FALSE && is.list( vL ) == FALSE )
      vL <- list( vL )
    
    nseq <- length( xL )
    nval <- length( vL )

    d.cat( 1, "Learning a Hidden Markov Model with", K, "state(s) for", nseq, "observation sequences\n" )

    ## Concatenate all training sequences together
    Ni <- unlist(lapply( xL, nrow ))
    x <- NULL
    for( i in 1:nseq ) x <- rbind( x, xL[[i]] )

    ## Concatenate all validation sequences together
    Mi <- rep( 0, nval )
    v <- NULL
    if( nval > 0 )
      {
        Mi <- unlist(lapply( vL, nrow ))
        for( i in 1:nval )
          v <- rbind( v, vL[[i]] )
      }
    
    N <- nrow(x)
    p <- ncol(x)

    ## Special case
    if( K == 1 )
      return( list( pi = 1, mu = matrix( apply( x, 2, mean ), nrow=1 ), sigma = array( cov(x), c(1,p,p) ),
                   A = matrix( 1, 1, 1) ) )
    
    ## Initialize HMM
    hmm <- HMM.make( K, p, apply( x, 2, mean ), cov(x) )
    if( is.null( hmm.init$mu ) == FALSE )	hmm$mu <- hmm.init$mu
    if( is.null( hmm.init$sigma ) == FALSE )	hmm$sigma <- hmm.init$sigma
    if( is.null( hmm.init$pi ) == FALSE )	hmm$pi <- hmm.init$pi
    if( is.null( hmm.init$A ) == FALSE )	hmm$A <- hmm.init$A
    
    if( SM.getOption("pure.R") == FALSE )
      {
        ## There is a weird issue that occurs when C code is the first thing ran after R startup
        ## Apparently random number generation using rnorm() in Rmath.h always returns the same number
        ## The problem seems to be alleviated when the first call to the function is made from R
        ## For this purpose:
        rnorm( 1, 0, 1 )
        
        ## Attempt to run C code
        if( is.infinite( max.iter ) ) max.iter <- -1
        res <- .C( "HMMlearnR",
                  as.numeric(x), as.numeric(v), as.integer(nseq), as.integer(nval),
                  as.integer(Ni), as.integer(Mi), as.integer(K), as.integer(p),
                  as.numeric(hmm$pi), as.numeric(hmm$A), as.numeric(hmm$mu), as.numeric(hmm$sigma),
                  as.numeric(tol), as.numeric(LLstop), as.integer(min.iter), as.integer(max.iter),
                  as.integer(SM.getOption("debug.output")),
                  pi=numeric(K), A=numeric(K*K), mu=numeric(K*p), sigma=numeric(K*p*p),
                  status = integer(1), NAOK=TRUE, DUP=FALSE )

        if( res$status == 1 )
          stop( "Covariance matrix is singular for one of the states" )
        if( res$status == 2 )
          stop( "Some samples were not accounted for by any state" )
        if( res$status == 3 )
          stop( "Computations in the C code became unstable" )
        else
          return( list(
                       pi = res$pi,
                       A = matrix( res$A, K, K ),
                       mu = matrix( res$mu, K, p ),
                       sigma = array( res$sigma, c( K, p, p ) )
                       ) )

      }

    d.cat( 1, "Currently running R version of the HMM.learn code\n" )
        
    hmm.prev <- list()

    compute.valLL <- function()
      {
        LL <- 0
        for( i in 1:nval )
          LL <- LL - sum(log(HMM.fwbk(vL[[i]], hmm)$c))
        LL
      }
    
    ## Training
    LLval <- 0
    LLval.prev <- 0
    LLprev <- 0
    LL <- 0
    LLbase <- 0
    iter <- 0
    while( TRUE )
      {
        ## Validation log-likelihood stopping criterion
        if( nval > 0)
          {
            LLval <- compute.valLL()
            d.cat( 3, "Validation log-likelihood: ", LLval, "\n" )
            if( iter <= min.iter || LLval >= LLval.prev )
              {
                LLval.prev <- LLval
                hmm.prev <- hmm
              }
            else
              {
                d.cat( 1, "Final log likelihood of training data after", iter-1, "iterations:", LLprev, "\n" )
                d.cat( 1, "Final validation log likelihood:", LLval.prev, "\n" )
                return( hmm.prev )
              }
          }

        ## Reinit cumulant variables for the current iteration
        LL <- 0
        tot.gamma1 <- rep( 0, K )
        tot.gamma <- matrix( 0, 0, K )
        tot.gammasum <- rep( 0, K )
        tot.xi <- matrix( 0, K, K )

        ## Iterate over each observation sequence
        for( si in 1:nseq )
          {
            ## Apply forward-backward algorithm
            res.fwbk <- HMM.fwbk( xL[[si]], hmm )
            LL <- LL - sum(log(res.fwbk$c))

            ## Compute contributions
            tot.gamma1 <- tot.gamma1 + res.fwbk$gamma[1,]
            tot.gamma <- rbind( tot.gamma, res.fwbk$gamma )
            tot.gammasum <- tot.gammasum + apply( res.fwbk$gamma, 2, sum )
            tot.xi <- tot.xi + apply( res.fwbk$xi, c(2,3), sum )
          }

        d.cat( if( (iter%%10)==0 ) 2 else 3, "Log likelihood at iteration", iter, ":", LL, "\n" )

        ## Verify stopping criteria
        if( iter > max.iter ) break
        if( LL > LLstop ) break
        if( is.infinite(LL) ) break
        if( LL == LLprev ) break

        if( iter <= (min.iter-1) ) LLbase <- LL
        else if( (LL-LLbase) < (1+tol)*(LLprev-LLbase) ) break

        ## Update means and transition matrix
        hmm$mu <- (t(tot.gamma) %*% x) / matrix( tot.gammasum, K, p )
        for( i in 1:K ) hmm$A[i,] <- tot.xi[i,] / sum(tot.xi[i,])

        ## Update priors
        hmm$pi <- tot.gamma1 / sum(tot.gamma1)

        ## Update covariance matrices
        for( l in 1:K )
          {
            hmm$sigma[l,,] <- matrix( 0, p, p )
            for( t in 1:N )
              {
                x.c <- x[t,] - hmm$mu[l,]
                hmm$sigma[l,,] <- hmm$sigma[l,,] + tot.gamma[t,l] * (x.c %*% t(x.c))
              }
            hmm$sigma[l,,] <- hmm$sigma[l,,] / tot.gammasum[l]
          }

        LLprev <- LL
        iter <- iter+1
      }

    d.cat( 1, "Final log likelihood after", iter, "iterations:", LL, "\n" )
    
    hmm
  }
