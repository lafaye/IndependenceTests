#####################################
# Computation of the test statistic #
#####################################

# X : matrix of size n x q
# vecd: vector of length p giving the sizes q_1,...,q_p of each one of the p subvectors. 
# We have q=\sum_{l=1}^p q_l.
# a: real number (parameter in the w() weight function)
# weight.choice: weight.choice of the weight function as described in our paper.

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




