#' Computation of the multidimensional distance covariance statistic for mutual independence using characteristic functions.
#'
#' @description
#' Computation of the multidimensional distance covariance statistic for mutual independence using
#' characteristic functions. Compute the eigenvalues associated with the empirical covariance of the
#' limiting Gaussian procces. Compute the \eqn{p}-value associated with the test statistic, using the Imhof procedure.
#' 
#' @usage mdcov(X, vecd, a = 1, weight.choice = 1, N = 200, cubature = FALSE, K =
#' 100, epsrel = 10 ^ -6, norming = TRUE, thresh.eigen = 10 ^ -8, estim.a =
#' FALSE, Cpp = TRUE, pval.comp = TRUE)
#' @param X Data.frame or matrix with observations corresponding to rows and variables to columns.
#' @param vecd a vector giving the sizes of each subvector.
#' @param a real number parameter for the weight function.
#' @param weight.choice Integer value in 1, 2, 3, 4, 5 corresponding to the choice in our paper.
#' @param N Number of Monte-Carlo samples.
#' @param cubature Logical. If \code{FALSE}, a Monte-Carlo approach is used. If \code{TRUE}, a cubature approach
#' is used.
#' @param K Number of eigenvalues to compute.
#' @param epsrel relative accuracy requested for the Imhof procedure.
#' @param norming Logical. Should we normalize the test statistic with \eqn{H_n}.
#' @param thresh.eigen We will not compute eigenvalues (involved in the limiting distribution) below that threshold.
#' @param estim.a Logical. Should we automatically estimate the value of \eqn{a}.
#' @param Cpp Logical. If \code{TRUE} computations will be done using a fast C code. 
#' The use of \code{FALSE} is only useful to compare the results with the one given by the C code.
#' @param pval.comp Logical. If \code{FALSE} do not compute the p-values and lambdas.
#' 
#' @export
#' 
#' @author Lafaye de Micheaux P.
#' @references Fan Y., Lafaye de Micheaux P., Penev S. and Salopek D. (2017). Multivariate nonparametric test of independence, Journal of Multivariate Analysis, 153, 189-210.
#' @examples
#' \donttest{
#' a <- 1
#' # 4.1 Dependence among four discrete variables
#' set.seed(1)
#' n <- 100
#'w1 <- rpois(n, 1)
#'w3 <- rpois(n, 1)
#'w4 <- rpois(n, 1)
#'w6 <- rpois(n, 1)
#'w2 <- rpois(n, 3)
#'w5 <- rpois(n, 3)
#'x1 <- w1 + w2
#'x2 <- w2 + w3
#'x3 <- w4 + w5
#'x4 <- w5 + w6
#'
#'X <- cbind(x1, x2, x3, x4)
#'mdcov(X, vecd = c(2, 2), a, weight.choice = 1, N = 100, cubature = TRUE)
#'mdcov(X, vecd = c(1, 1, 1, 1), a, weight.choice = 1, N = 100, cubature = TRUE)
#'
#'X <- cbind(x1, x2)
#'mdcov(X, vecd = c(1, 1), a, weight.choice = 1, N = 100, cubature = TRUE)
#'
#'X <- cbind(x3, x4)
#'mdcov(X, vecd = c(1, 1), a, weight.choice = 1, N = 100, cubature = TRUE)
#'
#'# 4.2 Dependence between three bivariate vectors
#'set.seed(2)
#'n <- 200
#'Sigma <- matrix(c(
#'  1 , 0 ,  0 ,  0 ,  0 ,  0 ,
#'  0 , 1 ,  0 ,  0 ,  0 ,  0 ,
#'  0 , 0 ,  1 ,  0 , .4 , .5 ,
#'  0 , 0 ,  0 ,  1 , .1 , .2 ,
#'  0 , 0 , .4 , .1 ,  1 ,  0 ,
#'  0 , 0 , .5 , .2 ,  0 ,  1 ) ,
#'  nrow = 6 , ncol = 6)   
#' W <- MASS::mvrnorm(n = n, mu = rep(0,6), Sigma = Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#'mdcov(W, vecd = c(2, 2, 2), a, weight.choice = 1, N = 100, cubature = TRUE, epsrel = 10 ^ -7)
#'
#'# X^{(1)} with X^{(2)}^2
#'mdcov(W[,1:4], vecd = c(2, 2), a, weight.choice = 1, N = 100, cubature = TRUE)
#'
#'# X^{(2)} with X^{(3)}^2
#'mdcov(W[,2:6], vecd = c(2, 2), a, weight.choice = 1, N = 100, cubature = TRUE)
#'
#'# X^{(1)} with X^{(3)}^2
#'mdcov(W[,c(1:2, 4:6)], vecd = c(2, 2), a, weight.choice = 1, N = 100, cubature = TRUE)
#'
#'# 4.3 Four-dependent variables which are 2-independent and 3-independent
#'set.seed(3)
#'n <- 300
#'W <- sample(1:8, n, replace = TRUE)
#'X1 <- W %in% c(1, 2, 3, 5)
#'X2 <- W %in% c(1, 2, 4, 6)
#'X3 <- W %in% c(1, 3, 4, 7)
#'X4 <- W %in% c(2, 3, 4, 8)
#'X <- cbind(X1, X2, X3, X4)
#'# pairwise independence
#'mdcov(X[,c(1, 2)], vecd = c(1, 1), a, weight.choice = 1, cubature = TRUE)
#'mdcov(X[,c(1, 3)], vecd = c(1, 1), a, weight.choice = 1, N = 100, cubature = TRUE)
#'mdcov(X[,c(1, 4)], vecd = c(1, 1), a, weight.choice = 1, cubature = TRUE)
#'mdcov(X[,c(2, 3)], vecd = c(1, 1), a, weight.choice = 1, cubature = TRUE)
#'mdcov(X[,c(2, 4)], vecd = c(1, 1), a, weight.choice = 1, cubature = TRUE)
#'mdcov(X[,c(3, 4)], vecd = c(1, 1), a, weight.choice = 1, cubature = TRUE)
#'# 3-independence
#'mdcov(X[,c(1, 2, 3)], vecd = c(1, 1, 1), a, weight.choice =
#'        1, cubature = TRUE)
#'mdcov(X[,c(1, 2, 4)], vecd = c(1, 1, 1), a, weight.choice =
#'        1, cubature = TRUE)
#'mdcov(X[,c(1, 3, 4)], vecd = c(1, 1, 1), a, weight.choice =
#'        1, cubature = TRUE)
#'mdcov(X[,c(2, 3, 4)], vecd = c(1, 1, 1), a, weight.choice =
#'        1, cubature = TRUE)
#'# 4-dependence
#'mdcov(X, vecd = c(1, 1, 1, 1), a, weight.choice = 1, cubature = TRUE)
#'}
#' @return
#' A list with the following components:
#' \item{mdcov}{the value of the statistic \eqn{nT_n(w)}{n * Tn(w)} (this value has been normed if \code{norming = TRUE})}
#' \item{Hn}{the denominator of \eqn{nT_n(w)}{nTn(w)}, namely \eqn{H_n}}
#' \item{pvalue}{the \eqn{p}-value of the test}
#' \item{lambdas}{the vector of eigenvalues computed (they have not been divided by their sum)}

mdcov <- function(X, vecd, a = 1, weight.choice = 1, N = 200, cubature = FALSE,
                  K = 100, epsrel = 10 ^ -6, norming = TRUE,
                  thresh.eigen = 10 ^ -8, estim.a = FALSE, Cpp = TRUE, pval.comp = TRUE) {
    version <- 2 # Case (b) in our Theorem 1.
    if (!(weight.choice %in% 1:5)) stop("'weight.choice' should be in 1, 2, 3, 4, 5.")
    if ((weight.choice == 1) & (a <= 0)) stop("'a' should be > 0 for 'weight.choice' = 1.")
    if ((weight.choice == 2) & (a <= 0)) stop("'a' should be > 0 for 'weight.choice' = 2.")
    if ((weight.choice == 3) & ((a <= 0) | a >= 2)) stop("You should take 0 < a < 2 for 'weight.choice' = 3.")
    if (weight.choice == 3) version <- 1 # Case (a) in our Theorem 1.
    if ((weight.choice == 5) & (a <= 0)) stop("'a' should be > 0 for 'weight.choice' = 5.")

    choosea <- function(avec, X, vecd = c(1, 1), weight.choice = 1) { # Value of the test statistic as a function of 'a'.
        n <- nrow(X)
        p <- ncol(X)
        la <- length(avec)
        res <- rep(NA, la)
        for (i in 1:la) {
            res[i] <- .C("Dcov2Cnormed", as.double(X), as.integer(vecd), as.double(avec[i]), as.integer(n), as.integer(p), as.integer(weight.choice), res = as.double(0.0), denom = 0.0, PACKAGE = "IndependenceTests")$res
        }
        return(n * res)
    }

    if (estim.a) a <- optim(0.1, choosea, gr = NULL, X = X, vecd = vecd, lower = 0, upper = 100, method = "L-BFGS-B", control = list(fnscale = -1), weight.choice = weight.choice)$par

    n <- nrow(X)
    q <- ncol(X)
    p <- length(vecd)
    if (Cpp) { # Fast version 1 that uses a C code.
        if (version == 1) res <- .C("Dcov1Cnormed", as.double(X), as.integer(vecd), as.double(a), as.integer(n), as.integer(p), as.integer(weight.choice), res = as.double(0.0), denom = 0.0, PACKAGE = "IndependenceTests")
        if (version == 2) res <- .C("Dcov2Cnormed", as.double(X), as.integer(vecd), as.double(a), as.integer(n), as.integer(p), as.integer(weight.choice), res = as.double(0.0), denom = 0.0, PACKAGE = "IndependenceTests")
        denom <- res$denom
        res <- res$res
    } else { # Slow version 1 that DOES NOT use the C code. Kept for historical reasons.
        if (version == 1) { # the one using the gamma's and beta's
            gammajl <- function(j, l, X, a, vecd, weight.choice) {
                ind.begin.bloc.l <- if (l == 1) 1 else sum(vecd[1:(l - 1)]) + 1
                prod <- 1
                if (weight.choice == 1) {
                    for (k in 1:vecd[l]) {
                        prod <- prod * exp(-(a * X[j, ind.begin.bloc.l + k - 1]) ^ 2 / 2.0)
                    }
                }
                if (weight.choice == 2) {
                    for (k in 1:vecd[l]) {
                        prod <- prod / (1 + (a * X[j, ind.begin.bloc.l + k - 1]) ^ 2)
                    }
                }
                if (weight.choice == 3) {
                    somme <- 0
                    for (k in 1:vecd[l]) {
                        somme <- somme + (X[j, ind.begin.bloc.l + k - 1]) ^ 2
                    }
                    prod <- somme ^ (a / 2)
                }
                if (weight.choice == 4) {
                    for (k in 1:vecd[l]) {
                        prod <- prod * (-2 * abs(X[j, ind.begin.bloc.l + k - 1]) + abs(X[j, ind.begin.bloc.l + k - 1] - 2 * a) + abs(X[j, ind.begin.bloc.l + k - 1] + 2 * a)) / (4 * abs(a))
                    }
                }
                if (weight.choice == 5) {
                    for (k in 1:vecd[l]) {
                        tmp <- X[j, ind.begin.bloc.l + k - 1]
                        prod <- prod * (1 - (a * tmp) ^ 2) * exp(-(a * tmp) ^ 2 / 2)
                    }
                }
                return(1 - prod)
            }
            gammajjprimel <- function(j, jprime, l, X, a, vecd, weight.choice) {
                ind.begin.bloc.l <- if (l == 1) 1 else sum(vecd[1:(l - 1)]) + 1
                prod <- 1
                if (weight.choice == 1) {
                    for (k in 1:vecd[l]) {
                        prod <- prod * exp(-(a * (X[j, ind.begin.bloc.l + k - 1] - X[jprime, ind.begin.bloc.l + k - 1])) ^ 2 / 2.0)
                    }
                }
                if (weight.choice == 2) {
                    for (k in 1:vecd[l]) {
                        prod <- prod / (1 + (a * (X[j, ind.begin.bloc.l + k - 1] - X[jprime, ind.begin.bloc.l + k - 1])) ^ 2)
                    }
                }
                if (weight.choice == 3) {
                    somme <- 0
                    for (k in 1:vecd[l]) {
                        somme <- somme + (X[j, ind.begin.bloc.l + k - 1] - X[jprime, ind.begin.bloc.l + k - 1]) ^ 2
                    }
                  prod <- somme ^ (a / 2)
                }
                if (weight.choice == 4) {
                    for (k in 1:vecd[l]) {
                        prod <- prod * (-2 * abs(X[j, ind.begin.bloc.l + k - 1] - X[jprime, ind.begin.bloc.l + k - 1]) + abs(X[j, ind.begin.bloc.l + k - 1] - X[jprime, ind.begin.bloc.l + k - 1] - 2 * a) + abs(X[j, ind.begin.bloc.l + k - 1] - X[jprime, ind.begin.bloc.l + k - 1] + 2 * a)) / (4 * abs(a))
                    }
                }
                if (weight.choice == 5) {
                    for (k in 1:vecd[l]) {
                        tmp <- X[j, ind.begin.bloc.l + k - 1] - X[jprime, ind.begin.bloc.l + k - 1]
                        prod <- prod * (1 - (a * tmp) ^ 2) * exp(-(a * tmp) ^ 2 / 2)
                    }
                }
                return(1 - prod)
            }
            betajjprimel <- function(j, jprime, l, X, a, vecd, weight.choice) {
                return(gammajl(j, l, X, a, vecd, weight.choice) + gammajl(jprime, l, X, a, vecd, weight.choice) - gammajjprimel(j, jprime, l, X, a, vecd, weight.choice))
            }
            denom1 <- function(p, n, X, a, vecd, weight.choice) {
                glvec <- blvec <- rep(0, p)
                for (l in 1:p) {
                    for (j in 1:n) {
                        for (jprime in 1:n) {
                            blvec[l] <- blvec[l] + betajjprimel(j, jprime, l, X, a, vecd, weight.choice)        
                        }
                        glvec[l] <- glvec[l] + gammajl(j, l, X, a, vecd, weight.choice)
                    }
                }
                blvec <- blvec  / n ^ 2
                glvec <- glvec / n
                
                res <- 0.0
                Ip <- 1:p
                allBsets <- lapply(Ip, function(x) combn(p, x))
                for (indB1 in 1:p) {
                    for (indB2 in 1:choose(p, indB1)) {
                        B <- allBsets[[indB1]][, indB2]
                        for (indBprime1 in 1:p) {
                            for (indBprime2 in 1:choose(p, indBprime1)) {
                                Bprime <- allBsets[[indBprime1]][, indBprime2]
                                BuBprime <- union(B, Bprime)
                                cardBuBprime <- length(BuBprime)
                                if (cardBuBprime == 0) stop("BUG")
                                if (cardBuBprime != 1) {
                                    BnBprime <- intersect(B, Bprime)
                                    symdif <- setdiff(BuBprime, BnBprime)
                                    res <- res + (cardBuBprime - 1) * prod(blvec[BnBprime]) * prod(-glvec[symdif])
                                }
                            }       
                        }
                        if (length(B) != 1) {
                            res <- res + 2 * (length(B) - 1) * prod(-glvec[B])
                        }
                    }   
                }
                
                return(res)
            }
            denom <- denom1(p, n, X, a, vecd, weight.choice)
            res <- 0.0
            Ip <- 1:p
            allBsets <- lapply(Ip, function(x) combn(p, x))
            for (indB1 in 1:p) {
                for (indB2 in 1:choose(p, indB1)) {
                    B <- allBsets[[indB1]][, indB2]
                    for (indBprime1 in 1:p) {
                        for (indBprime2 in 1:choose(p, indBprime1)) {
                            Bprime <- allBsets[[indBprime1]][, indBprime2]
                            term1 <- 0
                            for (j in 1:n) {
                                for (jprime in 1:n) {
                                    prod1 <- 1
                                    for (l in intersect(B, Bprime)) prod1 <- prod1 * betajjprimel(j, jprime, l, X, a, vecd, weight.choice)
                                    prod2 <- 1
                                    for (l in setdiff(B, Bprime)) prod2 <- prod2 * gammajl(j, l, X, a, vecd, weight.choice)
                                    prod3 <- 1
                                    for (l in setdiff(Bprime, B)) prod3 <- prod3 * gammajl(jprime, l, X, a, vecd, weight.choice)
                                    term1 <- term1 + prod1 * prod2 * prod3
                                }
                            }
                            term1 <- term1 / n ^ 2
                            term2 <- 0
                            for (j in 1:n) {
                                prod1 <- 1
                                for (l in intersect(B, Bprime)) {
                                    somme <- 0
                                    for (jprime in 1:n) {
                                        somme <- somme + betajjprimel(j, jprime, l, X, a, vecd, weight.choice)
                                    }
                                    somme <- somme / n
                                    prod1 <- prod1 * somme
                                }
                                prod2 <- 1
                                for (l in setdiff(B, Bprime)) prod2 <- prod2 * gammajl(j, l, X, a, vecd, weight.choice)
                                prod3 <- 1
                                for (l in setdiff(Bprime, B))  {
                                    somme <- 0
                                    for (jprime in 1:n) {
                                        somme <- somme + gammajl(jprime, l, X, a, vecd, weight.choice)
                                    }
                                    somme <- somme / n
                                    prod3 <- prod3 * somme
                                }
                                term2 <- term2 + prod1 * prod2 * prod3
                            }
                            term2 <- -2 * term2 / n
                            
                            prod1 <- 1
                            for (l in intersect(B, Bprime)) {
                                somme <- 0
                                for (j in 1:n) {
                                    for (jprime in 1:n) {
                                        somme <- somme + betajjprimel(j, jprime, l, X, a, vecd, weight.choice)
                                    }
                                }
                                somme <- somme / n ^ 2
                                prod1 <- prod1 * somme
                            }
                            prod2 <- 1
                            for (l in setdiff(B, Bprime))  {
                                somme <- 0
                                for (jprime in 1:n) {
                                    somme <- somme + gammajl(jprime, l, X, a, vecd, weight.choice)
                                }
                                somme <- somme / n
                                prod2 <- prod2 * somme
                            }
                            prod3 <- 1
                            for (l in setdiff(Bprime, B))  {
                                somme <- 0
                                for (jprime in 1:n) {
                                    somme <- somme + gammajl(jprime, l, X, a, vecd, weight.choice)
                                }
                                somme <- somme / n
                                prod3 <- prod3 * somme
                            }
                            term3 <- prod1 * prod2 * prod3
                            
                            res <- res + (-1) ^ (length(B) + length(Bprime)) * (term1 + term2 + term3)
                        }
                    }
                }
            }
        }
        if (version == 2) { # the one using the xi's
            xijjprimel <- function(j, jprime, l, X, a, vecd, weight.choice) {
                ind.begin.bloc.l <- if (l==1) 1 else sum(vecd[1:(l - 1)]) + 1
                prod <- 1
                if (weight.choice == 1) {
                    for (k in 1:vecd[l]) {
                        prod <- prod * exp(-(a * (X[j, ind.begin.bloc.l + k - 1] - X[jprime, ind.begin.bloc.l + k - 1])) ^ 2 / 2.0)
                    }
                }
                if (weight.choice == 2) {
                    for (k in 1:vecd[l]) {
                        prod <- prod / (1 + (a * (X[j, ind.begin.bloc.l + k - 1] - X[jprime, ind.begin.bloc.l + k - 1])) ^ 2)
                    }
                }
                if (weight.choice == 3) {
                    prod <- Inf
                }             
                if (weight.choice == 4) {
                    for (k in 1:vecd[l]) {
                        prod <- prod * (-2 * abs(X[j, ind.begin.bloc.l + k - 1] - X[jprime, ind.begin.bloc.l + k - 1]) + abs(X[j, ind.begin.bloc.l + k - 1] - X[jprime, ind.begin.bloc.l + k - 1] - 2 * a) + abs(X[j, ind.begin.bloc.l + k - 1] - X[jprime, ind.begin.bloc.l + k - 1] + 2 * a))/(4 * abs(a))
                    }
                }
                if (weight.choice == 5) {
                    for (k in 1:vecd[l]) {
                        tmp <- X[j, ind.begin.bloc.l + k - 1] - X[jprime, ind.begin.bloc.l + k - 1]
                        prod <- prod * (1 - (a * tmp) ^ 2) * exp(-(a * tmp) ^ 2 / 2)
                    }
                }
                return(prod)
            }

            denom2 <- function(p, n, X, a, vecd, weight.choice) {
                prod <- 1
                somme <- 0
                for (l in 1:p) {
                    xl <- 0
                    for (j in 1:n) {
                        for (jprime in 1:n) {
                            xl <- xl + xijjprimel(j, jprime, l, X, a, vecd, weight.choice)        
                        }
                    }
                    xl <- xl / n ^ 2
                    somme <- somme + 1 / xl
                    prod <- prod * xl
                }

                return(1 - (1 + somme - p) * prod)
            }

            denom <- denom2(p, n, X, a, vecd, weight.choice)
            
            term1 <- term2 <- 0
            term3 <- 1
            
            for (j in 1:n) {
                for (jprime in 1:n) {
                    prod <- 1
                    for (l in 1:p) {    
                        prod <- prod * xijjprimel(j, jprime, l, X, a, vecd, weight.choice)
                    }
                    term1 <- term1 + prod
                }
            }
            term1 <- term1 / n^2
            
            for (j in 1:n) {
                prod <- 1
                for (l in 1:p) {
                    somme <- 0
                    for (jprime in 1:n) {
                        somme <- somme + xijjprimel(j, jprime, l, X, a, vecd, weight.choice)
                    }
                    somme <- somme / n
                    prod <- prod * somme
                }
                term2 <- term2 + prod
            }
            term2 <- -2 * term2 / n
            
            for (l in 1:p) {
                somme <- 0
                for (j in 1:n) {
                    for (jprime in 1:n) {
                        somme <- somme + xijjprimel(j, jprime, l, X, a, vecd, weight.choice)        
                    }
                }
                somme <- somme / n ^ 2
                term3 <- term3 * somme
            }       
            res <- term1 + term2 + term3
        }
    }

    # Returns the value n * T_n(w) and the denominator of nT_n(w), namely H_n.
    stat <- n * res
    Hn <- denom
    if (norming) stat <- stat / Hn
    if (pval.comp) {
      pval.and.lambdas <- .mdcov.pval(stat, Hn, N = N, X = X, vecd = vecd, weight.choice = weight.choice, cubature = cubature, a = a, K = K, epsrel = epsrel, norming = norming, thresh.eigen = thresh.eigen)
    } else {
      pval.and.lambdas <- list(pval = NA, lambdas = NA)
    }
      
    return(list(mdcov = stat, Hn = Hn, pvalue = pval.and.lambdas$pval, lambdas = pval.and.lambdas$lambdas, a = a))
}




