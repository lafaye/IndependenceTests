######################################
# Computation of phihat and phitilde #
######################################

phiStandGauss <- function(vect) return(exp(-sum(vect ^ 2) / 2)) # NOT USED
  
phinhat <- function(vect, X, Cpp = TRUE) {
    if (Cpp) { # This one uses a C code.
        res <- .C("phinhatC", as.double(vect), as.double(X), as.integer(ncol(X)), as.integer(nrow(X)), res = as.complex(0.0), PACKAGE = "IndependenceTests")$res
    } else {
        res <- mean(exp(1i * X %*% as.matrix(vect)))
    }
  return(res)
}

phinlhat <- function(vect, l, X, vecd) { # NOT USED
  p <- length(vecd)
  if ((l > p) | (l < 1)) stop("l should be between 1 and p")
  ind.begin.bloc.l <- if (l == 1) 1 else sum(vecd[1:(l - 1)]) + 1
  res <- phinhat(as.matrix(vect[ind.begin.bloc.l:(ind.begin.bloc.l + vecd[l] - 1)]), X[, ind.begin.bloc.l:(ind.begin.bloc.l+vecd[l] - 1)])
  return(res)
}

phintilde <- function(vect, X, vecd) { # NOT USED
    p <- length(vecd)
    prod <- 1
    for (l in 1:p) {
      ind.begin.bloc.l <- if (l == 1) 1 else sum(vecd[1:(l - 1)]) + 1
      prod <- prod * phinhat(as.matrix(vect[ind.begin.bloc.l:(ind.begin.bloc.l + vecd[l] - 1)]), X[, ind.begin.bloc.l:(ind.begin.bloc.l + vecd[l] - 1)])
    }
    return(prod)
}

###########################################
# Computation of the covariance structure #
###########################################


Cst <- function(vecs, vect, vecd) { # True covariance structure.
  p <- length(vecd)
  prod1 <- prod2 <- 1
  somme <- 0
  for (l in 1:p) {
    ind.begin.bloc.l <- if (l == 1) 1 else sum(vecd[1:(l - 1)]) + 1
    vecsl <- as.matrix(vecs[ind.begin.bloc.l:(ind.begin.bloc.l + vecd[l] - 1)])
    vectl <- as.matrix(vect[ind.begin.bloc.l:(ind.begin.bloc.l + vecd[l] - 1)])
    tmp1 <- phiStandGauss(vecsl - vectl)
    tmp2 <- phiStandGauss(-vectl) * phiStandGauss(vecsl)
    prod1 <- prod1 * tmp1
    prod2 <- prod2 * tmp2
    somme <- somme + tmp1 / tmp2
  }
  res <- prod1 - prod2 * (1 - p + somme)
  return(res)  
}

Cmat <- function(yMat,vecd) { # Matrix computed using the true covariance.
  # yMat: N x q
    N <- nrow(yMat)
    res <- matrix(NA, nrow = N, ncol = N)
    for (i in 1:N) {
        for (j in 1:N) {
            res[i, j] <- Cst(as.vector(yMat[i, ]), as.vector(yMat[j, ]), vecd)
        }
    }
    return(res)
}

Cnhat <- function(vecs, vect, X, vecd, Cpp = TRUE) { # Estimated covariance structure.
    if (Cpp) { # This one uses a C code.
        res <- .C("CnhatC", as.double(vecs), as.double(vect), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.integer(length(vecd)), as.integer(vecd), res = as.complex(0.0), PACKAGE = "IndependenceTests")$res
    } else { # Kept for historical reasons.
        p <- length(vecd)
        prod1 <- prod2 <- 1
        somme <- 0
        for (l in 1:p) {
            ind.begin.bloc.l <- if (l == 1) 1 else sum(vecd[1:(l - 1)]) + 1
            vecsl <- as.matrix(vecs[ind.begin.bloc.l:(ind.begin.bloc.l + vecd[l] - 1)])
            vectl <- as.matrix(vect[ind.begin.bloc.l:(ind.begin.bloc.l + vecd[l] - 1)])
            Xl <- X[, ind.begin.bloc.l:(ind.begin.bloc.l + vecd[l] - 1)]
            tmp1 <- phinhat(vecsl - vectl, Xl)
            tmp2 <- phinhat(-vectl, Xl) * phinhat(vecsl, Xl)
            prod1 <- prod1 * tmp1
            prod2 <- prod2 * tmp2
            somme <- somme + tmp1 / tmp2
        }
        res <- prod1 - prod2 * (1 - p + somme)
    }
    return(res)
}


Cnhatmat <- function(yMat, X, vecd, Cpp = TRUE) { # Matrix computed using the estimated covariance.
  # yMat: N x q
    N <- nrow(yMat)
    if (Cpp) { # This one uses a C code.
        res <- rep(0.0, N * N)
        out <- .C("CnhatmatC", as.double(yMat), as.integer(nrow(yMat)), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.integer(length(vecd)), as.integer(vecd), res = as.complex(res), PACKAGE = "IndependenceTests")$res
        res <- matrix(out, nrow = N, ncol = N)
    } else {
        res <- matrix(NA, nrow = N, ncol = N)
        for (i in 1:N) {
            for (j in 1:N) {
                res[i, j] <- Cnhat(as.vector(yMat[i, ]), as.vector(yMat[j, ]), X, vecd, Cpp = FALSE)
            }
        }
    }
    return(res)
}
