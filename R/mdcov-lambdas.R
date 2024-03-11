#######################################################
# Computation of eigenvalues of the integral operator #
#######################################################

.sinint <- function(y) {

    mask1 <- (0 <= y) & (y <= 4)
    mask2 <- (0 >= y) & (y >= -4)
    mask3 <- y > 4
    mask4 <- y < -4

    Si <- function(x) {x*(1. +
                x*x*(-4.54393409816329991e-2 +
                         x*x*(1.15457225751016682e-3 +
                                  x*x*(-1.41018536821330254e-5 +
                                           x*x*(9.43280809438713025e-8 +
                                                    x*x*(-3.53201978997168357e-10 +
                                                             x*x*(7.08240282274875911e-13 +
                                                                      x*x*(-6.05338212010422477e-16))))))))/ (1. + 
           x*x*(1.01162145739225565e-2 +
                    x*x*(4.99175116169755106e-5 + 
                             x*x*(1.55654986308745614e-7 +
                                      x*x*(3.28067571055789734e-10 +
                                               x*x*(4.5049097575386581e-13 + 
                                                        x*x*(3.21107051193712168e-16)))))))
                   }
    


    


   f <- function(x) {(1. + 
        (1 / x ^2)*(7.44437068161936700618e2 +
           (1 / x ^2)*(1.96396372895146869801e5 +
              (1 / x ^2)*(2.37750310125431834034e7 +
                 (1 / x ^2)*(1.43073403821274636888e9 +
                    (1 / x ^2)*(4.33736238870432522765e10 +
                       (1 / x ^2)*(6.40533830574022022911e11 +
                          (1 / x ^2)*(4.20968180571076940208e12 +
                             (1 / x ^2)*(1.00795182980368574617e13 +
                                (1 / x ^2)*(4.94816688199951963482e12 +
                                   (1 / x ^2)*(-4.94701168645415959931e11))))))))))) /
         (x*(1. +
              (1 / x ^2)*(7.46437068161927678031e2 +
                 (1 / x ^2)*(1.97865247031583951450e5 +
                    (1 / x ^2)*(2.41535670165126845144e7 +
                       (1 / x ^2)*(1.47478952192985464958e9 +
                          (1 / x ^2)*(4.58595115847765779830e10 +
                             (1 / x ^2)*(7.08501308149515401563e11 +
                                (1 / x ^2)*(5.06084464593475076774e12 +
                                   (1 / x ^2)*(1.43468549171581016479e13 +
                                      (1 / x ^2)*(1.11535493509914254097e13)))))))))))
                 }

    g <- function(x) {(1 / x ^2)*(1. +
          (1 / x ^2)*(8.1359520115168615e2 +
             (1 / x ^2)*(2.35239181626478200e5 +
                (1 / x ^2)*(3.12557570795778731e7 +
                   (1 / x ^2)*(2.06297595146763354e9 +
                      (1 / x ^2)*(6.83052205423625007e10 +
                         (1 / x ^2)*(1.09049528450362786e12 +
                            (1 / x ^2)*(7.57664583257834349e12 +
                               (1 / x ^2)*(1.81004487464664575e13 +
                                  (1 / x ^2)*(6.43291613143049485e12 +
                                     (1 / x ^2)*(-1.36517137670871689e12))))))))))) /
        (1. +
          (1 / x ^2)*(8.19595201151451564e2 +
             (1 / x ^2)*(2.40036752835578777e5 +
                (1 / x ^2)*(3.26026661647090822e7 +
                   (1 / x ^2)*(2.23355543278099360e9 +
                      (1 / x ^2)*(7.87465017341829930e10 +
                         (1 / x ^2)*(1.39866710696414565e12 +
                            (1 / x ^2)*(1.17164723371736605e13 +
                               (1 / x ^2)*(4.01839087307656620e13 +
                                  (1 / x ^2)*(3.99653257887490811e13))))))))))
                  }

    
    res <- rep(NA, length(y))
    res[mask1] <- Si(y[mask1])
    res[mask2] <- -Si(-y[mask2])
    x <- y[mask3]
    res[mask3] <- (pi / 2) - f(x) * cos(x) - g(x) * sin(x)
    x <- -y[mask4]
    res[mask4] <- -( (pi / 2) - f(x) * cos(x) - g(x) * sin(x) )
    return(res)
}

.Cnhatmat <- function(yMat, X, vecd) { # Matrix computed using the estimated covariance.
  # yMat: N x q
    N <- nrow(yMat)
    res <- rep(0.0, N * N)
    out <- .C("CnhatmatC", as.double(yMat), as.integer(nrow(yMat)), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.integer(length(vecd)), as.integer(vecd), res = as.complex(res), PACKAGE = "IndependenceTests")$res
    res <- matrix(out, nrow = N, ncol = N)
    return(res)
}


.mdcov.lambdas <- function(N = 200, X, vecd, weight.choice = 1, cubature = FALSE, a = 1, K = 100, thresh.eigen = 10 ^ -8) {


    
    q <- sum(vecd)
  # SOME WORK NEED TO BE DONE HERE TO BE SURE THAT the \aleavect{Y}_j (j=1,...,N) are generated with a density w(.)=\prod_{l=1}^p v(.) where v(.) 
  # is given by the 'weight.choice' argument (takes its value in {1,2,3,4,5}). See our draft paper on pages 15-17.
    if (!cubature) { # Monte-Carlo approach.
        if (weight.choice == 1 ) yMat <- matrix(rnorm(N * q, sd = a), nrow = N, ncol = q) # Check if the sd=a is allright!! Check also for the cubature approach!!
        if (weight.choice == 2 ) yMat <- matrix(urlaplace(N * q, location = 0, scale = a), nrow = N, ncol = q)
        if (weight.choice == 3 ) stop("Not done yet.") # TODO!!!
        if (weight.choice == 4 ) {
            tmp <- runif(N * q)
            res <- rep(NA, N * q)
            for (j in 1: (N * q)) {
                res[j] <- uniroot(function(x, u, a) {if (x == 0) (pi * abs(a) / 2  - a * sin(a * x)) / (pi * abs(a)) - u else (pi * abs(a) / 2 - (sin(a * x)) ^ 2 / x + a * .sinint(2 * a * x)) / (pi * abs(a)) - u}, interval = c(-10 ^ 12, 10 ^ 12), u = tmp[j], a = a, tol = .Machine$double.eps)$root
            }
            yMat <- matrix(res, nrow = N, ncol = q)
        }
        if (weight.choice == 5 ) {
            tmp <- runif(N * q)
            ord <- order(tmp)
            tmp <- sort(tmp)
            res <- rep(NA, N * q)
            CDF <- function(x, a) pnorm(x / a) - dnorm(x, sd = a) * x
            bound <- 1
            while (CDF(bound, a) < 1) bound <- bound + 1
            for (j in 1: (N * q)) {
              res[j] <- uniroot(function(x, u, a) {pnorm(x / a) - dnorm(x, sd = a) * x - u}, interval = c(-bound, bound), u = tmp[j], a = a, tol = .Machine$double.eps)$root
            }
            yMat <- matrix(res[ord], nrow = N, ncol = q)
        }
        weights <- rep(1 / N, N)   
    } else { # Cubature approach.
                                        # WORK HAS TO BE DONE HERE!!! Sometimes, the weights are not between 0 and 1 !!!
                                        #   tmp <- createIntegrationGrid("GQN", q, 4)  # Maximum value is 25 ... but then N can become very big and eigen() takes time!
#   yMat <- tmp$nodes
                                        #   weights <- tmp$weights

        if (weight.choice != 1) stop("Cubature formulas for these weight choices are not yet implemented (nor exist!).")
        
        cuba <- 1
        if (q == 2) {
            if (cuba == 1) N <- 44
            if (cuba == 2) N <- 44
            if (cuba == 3) N <- 172
        }
        if (q >= 3) N <- 2 ^ (q + 1) + 4 * q ^ 2
        tmp <- .C("cubature1", a = as.double(a), q = as.integer(q), N = as.integer(N), prec = as.double(10 ^ (-6)), MAXIT = as.integer(100), cuba = as.integer(cuba), respoids = as.double(rep(0.0, N)), resabs = as.double(rep(0.0, N * q)), PACKAGE = "IndependenceTests")
        yMat <- matrix(tmp$resabs, nrow = N, ncol = q) # Transposee de tmp$resabs??
        weights <- abs(tmp$respoids) # PROBLEME ICI!!!!!!! Des fois les poids sont < 0 (par exemple quand q = 10 et vecd=c(5,5) dans les exemples de notre papier)
        
  
#  N <- length(weights)
  
    }

#  Cmat <- Cmat(yMat,vecd)

#    tmp <- diag(sqrt(weights))
    # This step is quite time consuming!
#    Mat <- tmp %*% .Cnhatmat(yMat, X, vecd) %*%  tmp

# This instruction computes only the lower part of the below matrix 
# and stores the result into a vector (not a matrix) to be passed to AP below.
    Matlower <- .C("CnhatmatClower", as.double(yMat), as.integer(nrow(yMat)), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.integer(length(vecd)), as.integer(vecd), res = as.complex(rep(0.0, N * (N + 1) / 2)), sqrtweights = as.double(sqrt(weights)), PACKAGE = "IndependenceTests")$res
    
    if (K > N) K <- N

#    K <- N # WARNING: If I do not compute a sufficient number of these, then the p-value computed afterwards  
             # seems to be NOT correct! We will have to understand why!!

    # We compute the largest eigenvalue
#    eig.largest <- .Fortran("zhpevx", JOBZ = as.character("N"), RANGE = as.character("I"), UPLO = as.character("L"), 
#                             N = as.integer(N), AP = as.complex(Matlower), VL = as.double(0), 
#                             VU = as.double(0), IL = as.integer(N), IU = as.integer(N),
#                             ABSTOL = as.double(10 ^ (-10)), M = as.integer(0), W = as.double(rep(0, N)), 
#                             Z = as.complex(0), LDZ = as.integer(1), WORK = as.complex(rep(0, 2 * N)), 
#                             RWORK = as.double(rep(0, 7 * N)), IWORK = as.integer(rep(0, 5 * N)),
#                             IFAIL = as.integer(0), INFO = as.integer(0), PACKAGE = "IndependenceTests")$W[1]
  
  
    eig.largest <- .C("zhpevxC", JOBZ = as.character("N"), RANGE = as.character("I"), UPLO = as.character("L"), 
                             N = as.integer(N), AP = as.complex(Matlower), VL = as.double(0), 
                             VU = as.double(0), IL = as.integer(N), IU = as.integer(N),
                             ABSTOL = as.double(10 ^ (-10)), M = as.integer(0), W = as.double(rep(0, N)), 
                             Z = as.complex(0), LDZ = as.integer(1), WORK = as.complex(rep(0, 2 * N)), 
                             RWORK = as.double(rep(0, 7 * N)), IWORK = as.integer(rep(0, 5 * N)),
                             IFAIL = as.integer(0), INFO = as.integer(0), PACKAGE = "IndependenceTests")$W[1]

    # We compute all eigenvalues between 10 ^ -8 and the largest and keep only the K largest ones
#    eig.res <- sort(.Fortran("zhpevx", JOBZ = as.character("N"), RANGE = as.character("V"), UPLO = as.character("L"), 
#                             N = as.integer(N), AP = as.complex(Matlower), VL = as.double(thresh.eigen), 
#                             VU = as.double(eig.largest) + 1, IL = as.integer(0), IU = as.integer(0),
#                             ABSTOL = as.double(10 ^ (-10)), M = as.integer(0), W = as.double(rep(0, N)), 
#                             Z = as.complex(0), LDZ = as.integer(1), WORK = as.complex(rep(0, 2 * N)), 
#                             RWORK = as.double(rep(0, 7 * N)), IWORK = as.integer(rep(0, 5 * N)),
#                             IFAIL = as.integer(0), INFO = as.integer(0), PACKAGE = "IndependenceTests")$W, decreasing = TRUE)[1:K]

      eig.res <- sort(.C("zhpevxC", JOBZ = as.character("N"), RANGE = as.character("V"), UPLO = as.character("L"), 
                             N = as.integer(N), AP = as.complex(Matlower), VL = as.double(thresh.eigen), 
                             VU = as.double(eig.largest) + 1, IL = as.integer(0), IU = as.integer(0),
                             ABSTOL = as.double(10 ^ (-10)), M = as.integer(0), W = as.double(rep(0, N)), 
                             Z = as.complex(0), LDZ = as.integer(1), WORK = as.complex(rep(0, 2 * N)), 
                             RWORK = as.double(rep(0, 7 * N)), IWORK = as.integer(rep(0, 5 * N)),
                             IFAIL = as.integer(0), INFO = as.integer(0), PACKAGE = "IndependenceTests")$W, decreasing = TRUE)[1:K]

    # Because "zhpevx" outputs a W of length N, we disregard all the zeroes
    eig.res <- eig.res[eig.res != 0]
    
#  eig.res <- eigen(Mat,symmetric=TRUE,only.values=TRUE)$values
#  eig.res2 <- eigen(diag(sqrt(weights)) %*% Conj(Cmat) %*%  diag(sqrt(weights)),symmetric=TRUE)$values
  # There is a BUG in rARPACK. It does not seem to work with complex matrices.
  # This is bad because eigen() from base R computes ALL the eigen-values, 
  # which is quite time consuming.
  # require("rARPACK")
  # eig.res <- eigs_sym(Cmat %*%  diag(weights),k=trunc,which="LM",opts=list(retvec=TRUE))
#  lambdas <- eig.res$values
# If the proof in our paper is correct, we do not need anymore to compute the eigenvectors because 
# the p_k's = 1.
#  fks <- diag(1/sqrt(weights)) %*% eig.res$vectors
#  fks.minus <- diag(1/sqrt(weights)) %*% eig.res2$vectors
#  abspks <- apply(fks*fks.minus,MARGIN=2,FUN=function(x) abs(sum(weights*x)))
#abspks <- rep(1,N)
#return(list(half.lambdas=(lambdas/2)[1:trunc],abspks=abspks[1:trunc]))
    return(eig.res)
}
