##' @useDynLib IndependenceTests
NULL

#' Nonparametric tests of independence between random vectors
#' 
#' @description
#' This function can be used for the following two problems:
#' 
#' 1) testing mutual independence between some 
#' numerical random vectors
#' 
#' 2) testing for serial independence of a multivariate
#' stationary quantitative time series.
#' 
#'  The proposed test does not assume
#' continuous marginals. It is valid for any probability distribution. It
#' is also invariant with respect to the affine general linear group of transformations on
#' the vectors. This test is based on a
#' characterization of mutual independence defined from probabilities
#' of half-spaces in a combinatorial formula of Möbius. As such,
#' it is a natural generalization of tests of independence between
#' univariate random variables using the empirical distribution
#' function. Without the assumption that each vector is
#' one-dimensional with a continuous cumulative distribution
#' function, any test of independence can not be distribution free.
#' The critical values of the proposed test are thus computed with
#' the bootstrap which was shown to be consistent in this context.
#' @usage dependogram(X, vecd.or.p, N = 10, B = 2000, alpha = 0.05,
#'  display = TRUE, graphics = TRUE, nbclus = 1)
#' 
#' @param X Data.frame or matrix with observations corresponding to rows and variables to columns.
#' @param vecd.or.p For the mutual independence problem 1), a vector 
#' giving the sizes of each subvector. For the serial independence
#' problem 2), an integer indicating the number of consecutive observations.
#' @param N Integer. Number of points of the discretization to obtain
#' directions on the sphere in order to evaluate the value of the test statistic.
#' @param B Integer. Number of bootstrap samples. Note that `B` can
#' be slightly modified if `nbclus` > 1
#' @param alpha Double. Global significance level of the test.
#' @param display Logical. TRUE to display values of the \eqn{A}-dependence statistics.
#' @param graphics Logical. TRUE to plot the dependogram.
#' @param nbclus Integer. Number of nodes in the cluster. Used only for parallel computations.
#' @author Bilodeau M., Lafaye de Micheaux P.
#' @encoding UTF-8
#' @references Beran R., Bilodeau M., Lafaye de Micheaux P. (2007). Nonparametric tests of independence between random vectors, Journal of Multivariate Analysis, 98, 1805-1824.
#' @md
#' @export
#' @return
#'  A list with the following components:
#'  
#'In the mutual independence case:
#'  
#'  \item{norm.RnA}{Supremum norm (Kolmogorov). Test statistic is
#'    \eqn{\|R_{n,A}\|}{||Rna||} and is computed
#'    from the \enc{Möbius}{Mobius} independence half space processes \eqn{R_{n,A}}{RnA}.}
#'\item{Rn}{Maximum value of `norm.RnA` over all subsets \eqn{A} of variables.}
#'\item{rA}{Critical value of the bootstrap distribution of the test
#'  statistic \eqn{\|R_{n,A}\|}{||Rna||}.}
#'\item{r}{Critical value of the bootstrap distribution of the test
#'  statistic \eqn{R_n}{Rn}.}
#'\item{RnAsstar}{Matrix of size \eqn{(2 ^ p - p - 1)\times B}{(2 ^ p -
#'                                                                p - 1) x B} which contains, for each of the \eqn{B} bootstrap samples, the statistics \code{norm.RnA} for all
#'  subsets \eqn{A} of variables.}
#'
#'In the serial case:
#'  
#'  \item{norm.SnA}{Supremum norm (Kolmogorov). Test statistic  is
#'    \eqn{\|S_{n,A}\|}{||Sna||} and is computed
#'    from the \enc{Möbius}{Mobius} independence half space processes \eqn{S_{n,A}}{SnA}.}
#'\item{Sn}{Maximum value of `norm.SnA` over all subsets \eqn{A} of variables.}
#'\item{sA}{Critical value of the bootstrap distribution of the test
#'  statistic \eqn{\|S_{n,A}\|}{||Sna||}.}
#'\item{s}{Critical value of the bootstrap distribution of the test
#'  statistic \eqn{S_n}{Sn}.}
#'\item{SnAsstar}{Matrix of size \eqn{(2 ^ {p - 1} - 1)\times B}{(2 ^ {p
#'  - 1} - 1) x B} which
#'  contains, for each of the \eqn{B} bootstrap samples, the statistics \code{norm.SnA} for all
#'  subsets \eqn{A} of variables.}
#'
#' @examples
#' # NOTE: In real applications, B should be set to at least 1000.
#'
#' # Example 4.1: Test of mutual independence between four discrete Poisson
#' #variables. The pair (X1,X2) is independent of the pair (X3,X4), with
#' #each pair having a correlation of 3/4.
#' # NOTE: with B=1000, this one took 65s with nbclus=1 and 15s with nbclus=7 on my computer.
#' n <- 100
#' W1 <- rpois(n, 1)
#' W3 <- rpois(n, 1)
#' W4 <- rpois(n, 1)
#' W6 <- rpois(n, 1)
#' W2 <- rpois(n, 3)
#' W5 <- rpois(n, 3)
#' X1 <- W1 + W2
#' X2 <- W2 + W3
#' X3 <- W4 + W5
#' X4 <- W5 + W6
#' X <- cbind(X1, X2, X3, X4)
#' dependogram(X, vecd.or.p = c(1, 1, 1, 1), N = 10, B = 20, alpha = 0.05,
#'            display = TRUE, graphics = TRUE)
#'
#' # Example 4.2: Test of mutual independence between three bivariate
#' #  vectors. The block-diagonal structure of the covariance matrix is
#' #such that only the second and third subvectors are dependent.
#' # NOTE: with B=2000, this one took 3.8h with nbclus=1 on my computer.
#' n <- 50
#' mu <- rep(0,6)
#' Psi <- matrix(c(1, 0, 0, 0, 0, 0,
#'                 0, 1, 0, 0, 0, 0,
#'                 0, 0, 1, 0,.4,.5,
#'                 0, 0, 0, 1,.1,.2,
#'                 0, 0,.4,.1, 1, 0,
#'                 0, 0,.5,.2, 0, 1), nrow = 6, byrow = TRUE)
#' X <- mvrnorm(n, mu, Psi)
#' \donttest{
#'   dependogram(X, vecd.or.p = c(2, 2, 2), N = 10, B = 20, alpha = 0.05,
#'               display = TRUE, graphics = TRUE)
#' }
#' 
#' # Example 4.3: Test of mutual independence between 4 dependent binary 
#' # variables which are 2-independent (pairwise) and also 3-independent
#' # (any 3 of the 4 variables are mutually independent).
#' \dontrun{
#' n <- 100
#' W <- sample(x = 1:8, size = n, TRUE)
#' X1 <- W %in% c(1, 2, 3, 5)
#' X2 <- W %in% c(1, 2, 4, 6)
#' X3 <- W %in% c(1, 3, 4, 7)
#' X4 <- W %in% c(2, 3, 4, 8)
#' X <- cbind(X1, X2, X3, X4)
#' dependogram(X, vecd.or.p = c(1, 1, 1, 1), N = 10, B = 20, alpha = 0.05,
#'             display = TRUE, graphics = TRUE)
#' # Example 4.4: Test of serial independence of binary sequences of zeros
#' # and ones. The sequence W is an i.i.d. sequence. The sequence Y is
#' # dependent at lag 3.
#' n <- 100 ; lag <- 3
#' W <- rbinom(n, 1, 0.8)
#' Y <- W[1:(n - lag)] * W[(1 + lag):n]
#' dependogram(W, vecd.or.p = 4, N = 10, B = 20, alpha = 0.05, display =
#'               TRUE, graphics = TRUE)
#' dependogram(Y, vecd.or.p = 4, N = 10, B = 20, alpha = 0.05, display =
#'               TRUE, graphics = TRUE)
#' }
#' # Example 4.5: Test of serial independence of sequences of directional
#' # data on the 2-dimensional sphere. The sequence W is an
#' # i.i.d. sequence. The sequence Y is dependent at lag 1.
#' # NOTE: with B=2000, this one took 7.9h with nbclus=1 on my computer.
#' n <- 75 ; lag <- 1
#' U <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
#' W <- U[1:(n - lag),] + sqrt(2) * U[(1 + lag):n,]
#' Y <- W / apply(W, MARGIN = 1, FUN = function(x) {sqrt(x[1] ^ 2 + x[2] ^ 2)})
#' 
#' \donttest{
#'   dependogram(Y, vecd.or.p = 3, N = 10, B = 20, alpha = 0.05, display =
#'                 TRUE, graphics = TRUE)
#' }
#' 
#' # This one always gives the same value of the test statistic:
#' \donttest{
#' x <- rnorm(100)
#' dependogram(X = cbind(x, x), vecd.or.p = c(1, 1), N = 2, B = 2, alpha =
#'               0.05, display = FALSE, graphics = FALSE, nbclus = 1)$Rn
#' # This is correct because this is equivalent to computing:
#' I <- 1:100
#' n <- 100
#' sqrt(n) * max(I/n - I ^ 2 / n ^ 2)
#' 
#' }

dependogram <- function(X,vecd.or.p,N=10,B=2000,alpha=0.05,display=TRUE,graphics=TRUE,nbclus=1) {
  
  
  if (nbclus>1) {
    
    #    suppressWarnings(parallel.pkg.present <- require(parallel))
    parallel.pkg.present <- "package:parallel" %in% search()
    Rmpi.pkg.present <- "package:Rmpi" %in% search()
    if (all(!c(parallel.pkg.present,Rmpi.pkg.present))) stop("Either package parallel or Rmpi should be installed!")
    #    suppressWarnings(rsprng.pkg.present <- require(rsprng))
    #    if (!rsprng.pkg.present) stop("Package rsprng is not installed!")
    cluster.type <- if (parallel.pkg.present) "PSOCK" else "MPI" # We prefer to use "PSOCK" (i.e. parallel) because it's easier.
    
  }
  
  
  X <- as.matrix(X)
  
  # if length(vecd.or.p)>1 then serial case, otherwise non serial case
  
  if (length(vecd.or.p) > 1) {
    # we perform the serial case
    seriel <- 0
    vecd <- vecd.or.p
    p <- length(vecd)
    taille <- 2^p-p-1
    #
    
    
    
    RnAs <- rep(0,2^p-p-1)
    Rn <- 0
    
    # We start the cluster
    if (nbclus > 1) {
      
      B <- round(B/nbclus)*nbclus
      
      
      cl <- parallel::makeCluster(nbclus, type = cluster.type) 
      
      
      myfunc <- function(B,p) {
        RnAsstar <- matrix(0,nrow=(2^p-p-1),ncol=B)
        Rnstar <- rep(0,B)
        require(IndependenceTests)
        # We call the C function dependogram
        .C("dependogramC",
           as.integer(N),
           as.integer(vecd),
           as.integer(length(vecd)),
           as.integer(p),
           as.numeric(X),
           as.integer(nrow(X)),
           as.integer(ncol(X)),
           B=as.integer(B),
           as.numeric(alpha),
           RnAs=as.numeric(RnAs),
           RnAsstar=as.numeric(RnAsstar),
           Rn=as.numeric(Rn),
           Rnstar=as.numeric(Rnstar),
           as.integer(seriel),PACKAGE="IndependenceTests")
        
      }
      
      out2 <- parallel::clusterCall(cl, myfunc, round(B/nbclus),p)
      
      
      
      # We stop the cluster
      parallel::stopCluster(cl)
      
      out <- list(RnAs=c(),RnAsstar=c(),Rn=c(),Rnstar=c())
      
      
      for (clus in 1:nbclus) {
        
        out$Rnstar <- c(out$Rnstar,out2[[clus]]$Rnstar)
        out$RnAsstar <- c(out$RnAsstar,out2[[clus]]$RnAsstar)
        
      }
      
      out$Rn <- out2[[clus]]$Rn
      out$RnAs <- out2[[clus]]$RnAs
      
      
    } else {
      
      RnAsstar <- matrix(0,nrow=(2^p-p-1),ncol=B)
      Rnstar <- rep(0,B)
      
      # We call the C function dependogram
      out <- .C("dependogramC",
                as.integer(N),
                as.integer(vecd),
                as.integer(length(vecd)),
                as.integer(p),
                as.numeric(X),
                as.integer(nrow(X)),
                as.integer(ncol(X)),
                as.integer(B),
                as.numeric(alpha),
                RnAs=as.numeric(RnAs),
                RnAsstar=as.numeric(RnAsstar),
                Rn=as.numeric(Rn),
                Rnstar=as.numeric(Rnstar),
                as.integer(seriel),PACKAGE="IndependenceTests")
      
    }
    
    
    
    
    
    
    if (display) {
      # require(combinat)
      RES <- as.list(1:(p-1))
      for (cardA in 2:p) {RES[[cardA]] <- as.matrix(combn(p,cardA))}
      nb <- 0
      for (cardA in 2:p) {
        for (j in 1:(choose(p,cardA))) {
          nb <- nb+1
          
          cat(c(nb,": A=",RES[[cardA]][,j],": ||RnA||=",round(out$RnAs[nb],3),"\n"))
          
          
        }
      }
    }
    
    
    # Sort the elements on each line
    RnAsstar <- matrix(out$RnAsstar,nrow=(2^p-p-1),ncol=B,byrow=FALSE)
    matseuils <- t(apply(RnAsstar,FUN=sort,MARGIN=1))
    
    
    beta <- (1-alpha)^(1/taille)
    # Contains the beta-quantiles of the stats RnA for each one of the 2^p-p-1 set A
    AllThresholds <- matseuils[,round(beta*B)]
    
    if (graphics) {
      # We draw a vertical bar for each A with height ||RnA||
      plot(out$RnAs,type="h",ylim=c(0,max(c(max(AllThresholds),max(out$RnAs),max(out$Rnstar)))+0.1),xlim=c(0,2^p-p),main="Dependogram",xlab="Subsets",ylab="||RnA||")
    }
    
    Rnstar <- sort(out$Rnstar)
    GlobalThreshold <- Rnstar[round((1-alpha)*B)]
    # abline(h=GlobalThreshold,col="red")
    
    
    if (graphics) { 
      # We put a star for each beta-quantile of ||R_A||
      points((1:(2^p-p-1)),AllThresholds,pch="*")
    }
    
    res <- list(norm.RnA=out$RnAs,Rn=out$Rn,rA=AllThresholds,r=GlobalThreshold,RnAsstar=RnAsstar)
    return(res)
    
    
  }
  
  
  
  if (length(vecd.or.p) == 1) {
    # We perform the serial case
    seriel <- 1
    p <- vecd.or.p
    vecd <- rep(ncol(X),p)
    taille <- 2^(p-1)-1
    
    
    
    
    SnAs <- rep(0,taille)
    Sn <- 0
    
    
    # We start the cluster
    if (nbclus > 1) {
      
      B <- round(B/nbclus)*nbclus
      
      SnAsstar <- matrix(0,nrow=taille,ncol=B)
      Snstar <- rep(0,B)
      
      cl <- parallel::makeCluster(nbclus, type = cluster.type) 
      
      
      myfunc <- function(B,p) {
        SnAsstar <- matrix(0,nrow=(2^p-p-1),ncol=B)
        Snstar <- rep(0,B)
        require(IndependenceTests)
        # We call the C function dependogram
        .C("dependogramC",
           as.integer(N),
           as.integer(vecd),
           as.integer(length(vecd)),
           as.integer(p),
           as.numeric(X),
           as.integer(nrow(X)),
           as.integer(ncol(X)),
           B=as.integer(B),
           as.numeric(alpha),
           SnAs=as.numeric(SnAs),
           SnAsstar=as.numeric(SnAsstar),
           Sn=as.numeric(Sn),
           Snstar=as.numeric(Snstar),
           as.integer(seriel),PACKAGE="IndependenceTests")
      }
      
      out2 <- parallel::clusterCall(cl, myfunc, round(B/nbclus),p)
      
      
      
      # We stop the cluster
      parallel::stopCluster(cl)
      
      out <- list(SnAs=c(),SnAsstar=c(),Sn=c(),Snstar=c())
      
      
      for (clus in 1:nbclus) {
        
        out$Snstar <- c(out$Snstar,out2[[clus]]$Snstar)
        out$SnAsstar <- c(out$SnAsstar,out2[[clus]]$SnAsstar)
        
      }
      
      out$Sn <- out2[[clus]]$Sn
      out$SnAs <- out2[[clus]]$SnAs
      
      
    } else {
      
      SnAsstar <- matrix(0,nrow=taille,ncol=B)
      Snstar <- rep(0,B)
      # We call the C function dependogram
      out <- .C("dependogramC",
                as.integer(N),
                as.integer(vecd),
                as.integer(length(vecd)),
                as.integer(p),
                as.numeric(X),
                as.integer(nrow(X)),
                as.integer(ncol(X)),
                as.integer(B),
                as.numeric(alpha),
                SnAs=as.numeric(SnAs),
                SnAsstar=as.numeric(SnAsstar),
                Sn=as.numeric(Sn),
                Snstar=as.numeric(Snstar),
                as.integer(seriel),PACKAGE="IndependenceTests")
      
    }
    
    
    
    
    
    
    
    if (display) {
      
      # require(combinat)
      RES <- as.list(1:(p-1))
      for (cardA in 2:p) {RES[[cardA]] <- as.matrix(rbind(rep(1,choose(p-1,cardA-1)),as.matrix(combn(p-1,cardA-1)+1)))}
      nb <- 0
      for (cardA in 2:p) {
        for (j in 1:(choose(p-1,cardA-1))) {
          nb <- nb+1
          
          cat(c(nb, ": A=",RES[[cardA]][,j],": ||SnA||=",round(out$SnAs[nb],3),"\n"))
          
        }
      }
    }
    
    
    SnAsstar <- matrix(out$SnAsstar,nrow=taille,ncol=B,byrow=FALSE)
    Sn <- max(out$SnAs)
    
    
    beta <- (1-alpha)^(1/taille)
    
    
    # The beta-quantile of S_n,A is computed by grouping all the values 
    # S_n,A^* with |A|=k as in the article
    # There are choose(p-1,|A|-1) sets A of size |A| containing 1
    # We therefore must take in the matrix SnAsstar groups of choose(p-1,|A|-1) rows, for |A|=2 to p.
    # There will be p-1 groups.
    # For each one of these groups, we create a vector vecA by taking all the elements in the group, 
    # then we compute AllThresholds[|A|]<-vecA[round(beta*B*choose(p-1,|A|-1))] for |A|=2 to p
    # This vector AllThresholds will thus contain the heights of the horizontal bars to put on the dependogram (one horizontal bar 
    # for each |A|)
    
    
    AllThresholds <- rep(0,p-1)
    begin <- 1
    end <- 0
    for (cardA in 2:p) {
      end <- end+choose(p-1,cardA-1)
      vecA <- as.vector(SnAsstar[begin:end,])
      vecA <- sort(vecA)
      AllThresholds[cardA-1] <- vecA[round(beta*B*choose(p-1,cardA-1))]
      begin <- end+1
    }
    
    if (graphics) {
      # We draw a vertical bar for each A of height ||SnA||
      plot(out$SnAs,type="h",ylim=c(0,max(c(max(AllThresholds),max(out$SnAs),max(out$Snstar)))+0.1),xlim=c(0,2^(p-1)),main="Dependogram",xlab="Subsets",ylab="||SnA||")
    }
    
    Snstar <- sort(out$Snstar)
    GlobalThreshold <- Snstar[round((1-alpha)*B)]
    # abline(h=GlobalThreshold,col="red")
    
    
    # It remains to put the horizontal lines of the threshold values for each |A|
    
    
    if (graphics) {
      begin <- 1
      end <- 0
      for (cardA in 2:p) {
        end <- end+choose(p-1,cardA-1)
        
        segments(begin-0.5,AllThresholds[cardA-1],end+0.5,AllThresholds[cardA-1],lty=4)
        
        begin <- end+1
      }
    }
    
    
    
    
    
    res <- list(norm.SnA=out$SnAs,Sn=out$Sn,sA=AllThresholds,s=GlobalThreshold,SnAsstar=SnAsstar)
    return(res)
  }
  
  
  
}



