\encoding{utf8}
\name{dependogram}
\alias{dependogram}
\title{Nonparametric tests of independence between random vectors}
\description{This function can be used for the following two problems:
  1) testing mutual independence between some 
  numerical random vectors, and 2) testing for serial independence of a multivariate
  stationary quantitative time series. The proposed test does not assume
  continuous marginals. It is valid for any probability distribution. It
is also invariant with respect to the affine general linear group of transformations on
the vectors. This test is based on a
characterization of mutual independence defined from probabilities
of half-spaces in a combinatorial formula of \enc{Möbius}{Mobius}. As such,
it is a natural generalization of tests of independence between
univariate random variables using the empirical distribution
function. Without the assumption that each vector is
one-dimensional with a continuous cumulative distribution
function, any test of independence can not be distribution free.
The critical values of the proposed test are thus computed with
the bootstrap which was shown to be consistent in this context.}
\usage{
dependogram(X, vecd.or.p, N = 10, B = 2000, alpha = 0.05, display =
            TRUE, graphics = TRUE, nbclus = 1)
}
\arguments{
  \item{X}{Data.frame or matrix with observations corresponding to rows and variables to columns.}
  \item{vecd.or.p}{For the mutual independence problem 1), a vector
    giving the sizes of each subvector. For the serial independence
    problem 2), an integer indicating the number of consecutive observations.}
  \item{N}{Integer. Number of points of the discretization to obtain
    directions on the sphere in order to evaluate the value of the test statistic.}
  \item{B}{Integer. Number of bootstrap samples. Note that \code{B} can
    be slightly modified if \code{nbclus} > 1}
  \item{alpha}{Double. Global significance level of the test.}
  \item{display}{Logical. TRUE to display values of the \eqn{A}-dependence statistics.}
\item{graphics}{Logical. TRUE to plot the dependogram.}
\item{nbclus}{Integer. Number of nodes in the cluster. Used only for parallel computations.}
}
\value{
 A list with the following components:

  In the mutual independence case:
  
  \item{norm.RnA}{Supremum norm (Kolmogorov). Test statistic is
  \eqn{\|R_{n,A}\|}{||Rna||} and is computed
  from the \enc{Möbius}{Mobius} independence half space processes \eqn{R_{n,A}}{RnA}.}
  \item{Rn}{Maximum value of \code{norm.RnA} over all subsets \eqn{A} of variables.}
  \item{rA}{Critical value of the bootstrap distribution of the test
  statistic \eqn{\|R_{n,A}\|}{||Rna||}.}
  \item{r}{Critical value of the bootstrap distribution of the test
  statistic \eqn{R_n}{Rn}.}
  \item{RnAsstar}{Matrix of size \eqn{(2 ^ p - p - 1)\times B}{(2 ^ p -
  p - 1) x B} which contains, for each of the \eqn{B} bootstrap samples, the statistics \code{norm.RnA} for all
  subsets \eqn{A} of variables.}

  In the serial case:

  \item{norm.SnA}{Supremum norm (Kolmogorov). Test statistic  is
  \eqn{\|S_{n,A}\|}{||Sna||} and is computed
  from the \enc{Möbius}{Mobius} independence half space processes \eqn{S_{n,A}}{SnA}.}
  \item{Sn}{Maximum value of \code{norm.SnA} over all subsets \eqn{A} of variables.}
  \item{sA}{Critical value of the bootstrap distribution of the test
  statistic \eqn{\|S_{n,A}\|}{||Sna||}.}
  \item{s}{Critical value of the bootstrap distribution of the test
  statistic \eqn{S_n}{Sn}.}
  \item{SnAsstar}{Matrix of size \eqn{(2 ^ {p - 1} - 1)\times B}{(2 ^ {p
  - 1} - 1) x B} which
  contains, for each of the \eqn{B} bootstrap samples, the statistics \code{norm.SnA} for all
  subsets \eqn{A} of variables.}
  
}
\author{Bilodeau M., Lafaye de Micheaux P.}

\references{Beran R., Bilodeau M., Lafaye de Micheaux
  P. (2007). Nonparametric tests of independence between random vectors, 
  Journal of Multivariate Analysis, 98, 1805-1824.
  }


\examples{
%% All these examples (but the last one) are taken from the paper
%% "Nonparametric tests of independence between random vectors"
%%  R. Beran, M. Bilodeau, P. Lafaye de Micheaux
%%  Journal of Multivariate Analysis 98 (2007) 1805-1824

# NOTE: In real applications, B should be set to at least 1000.

# Example 4.1: Test of mutual independence between four discrete Poisson
# variables. The pair (X1,X2) is independent of the pair (X3,X4), with
# each pair having a correlation of 3/4.
# NOTE: with B=1000, this one took 65s with nbclus=1 and 15s with nbclus=7 on my computer.
n <- 100
W1 <- rpois(n, 1)
W3 <- rpois(n, 1)
W4 <- rpois(n, 1)
W6 <- rpois(n, 1)
W2 <- rpois(n, 3)
W5 <- rpois(n, 3)
X1 <- W1 + W2
X2 <- W2 + W3
X3 <- W4 + W5
X4 <- W5 + W6
X <- cbind(X1, X2, X3, X4)
dependogram(X, vecd.or.p = c(1, 1, 1, 1), N = 10, B = 20, alpha = 0.05,
            display = TRUE, graphics = TRUE)

# Example 4.2: Test of mutual independence between three bivariate
#  vectors. The block-diagonal structure of the covariance matrix is
#  such that only the second and third subvectors are dependent.
# NOTE: with B=2000, this one took 3.8h with nbclus=1 on my computer.
n <- 50
mu <- rep(0,6)
Psi <- matrix(c(1, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0,
                0, 0, 1, 0,.4,.5,
                0, 0, 0, 1,.1,.2,
                0, 0,.4,.1, 1, 0,
                0, 0,.5,.2, 0, 1), nrow = 6, byrow = TRUE)
X <- mvrnorm(n, mu, Psi)
\donttest{
dependogram(X, vecd.or.p = c(2, 2, 2), N = 10, B = 20, alpha = 0.05,
            display = TRUE, graphics = TRUE)
}

# Example 4.3: Test of mutual independence between 4 dependent binary
# variables which are 2-independent (pairwise) and also 3-independent
# (any 3 of the 4 variables are mutually independent).
n <- 100
W <- sample(x = 1:8, size = n, TRUE)
X1 <- W \%in\% c(1, 2, 3, 5)
X2 <- W \%in\% c(1, 2, 4, 6)
X3 <- W \%in\% c(1, 3, 4, 7)
X4 <- W \%in\% c(2, 3, 4, 8)
X <- cbind(X1, X2, X3, X4)
dependogram(X, vecd.or.p = c(1, 1, 1, 1), N = 10, B = 20, alpha = 0.05,
            display = TRUE, graphics = TRUE)

# Example 4.4: Test of serial independence of binary sequences of zeros
# and ones. The sequence W is an i.i.d. sequence. The sequence Y is
# dependent at lag 3.
n <- 100 ; lag <- 3
W <- rbinom(n, 1, 0.8)
Y <- W[1:(n - lag)] * W[(1 + lag):n]
dependogram(W, vecd.or.p = 4, N = 10, B = 20, alpha = 0.05, display =
            TRUE, graphics = TRUE)
dependogram(Y, vecd.or.p = 4, N = 10, B = 20, alpha = 0.05, display =
            TRUE, graphics = TRUE)

# Example 4.5: Test of serial independence of sequences of directional
# data on the 2-dimensional sphere. The sequence W is an
# i.i.d. sequence. The sequence Y is dependent at lag 1.
# NOTE: with B=2000, this one took 7.9h with nbclus=1 on my computer.
n <- 75 ; lag <- 1
U <- matrix(rnorm(2 * n), nrow = n, ncol = 2)
W <- U[1:(n - lag),] + sqrt(2) * U[(1 + lag):n,]
Y <- W / apply(W, MARGIN = 1, FUN = function(x) {sqrt(x[1] ^ 2 + x[2] ^ 2)})
\donttest{
dependogram(Y, vecd.or.p = 3, N = 10, B = 20, alpha = 0.05, display =
            TRUE, graphics = TRUE)
}

# This one always gives the same value of the test statistic:
x <- rnorm(100)
dependogram(X = cbind(x, x), vecd.or.p = c(1, 1), N = 2, B = 2, alpha =
0.05, display = FALSE, graphics = FALSE, nbclus = 1)$Rn
# This is correct because this is equivalent to computing:
I <- 1:100
n <- 100
sqrt(n) * max(I/n - I ^ 2 / n ^ 2)

}
