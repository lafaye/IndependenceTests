\encoding{UTF-8}
\name{A.dep.tests}
\alias{A.dep.tests}
\title{Tests for mutual or serial independence between categorical variables}
\description{The tests are constructed from the \enc{Möbius}{Mobius} transformation
applied to the probability cells in a multi-way contingency table. The
Pearson chi-squared test of mutual independence is partitioned into
A-dependence statistics over all subsets A of variables. The goal of the
partition is to identify subsets of dependent variables when the mutual
independence hypothesis is rejected by the Pearson chi-squared test. The
methodology can be directly adapted to test for serial independence of d
successive observations of a
stationary categorical time series.

For categorical time series,
especially those of a nominal (non ordinal) nature, the user should be
aware that tests of serial independence obtained by methods suited to
quantitative sequences by quantification of the labels are not 
invariant to permutation of the labels contrary to the test described here.}
\usage{
A.dep.tests(Xmat, choice = 1, d = 0, m = d, freqname = "", type = "text")
}
\arguments{
  \item{Xmat}{Table, matrix or data-frame of the contingency table data,
    if \code{choice} = 1. Vector of the time series data, if
    \code{choice} = 2.}
  \item{choice}{Integer. 1 for mutual independence, 2 for serial independence.}
  \item{d}{Integer. Used only if \code{choice} = 2 for the number of successive obervations.}
  \item{m}{Integer. Maximum cardinality of subsets \eqn{A} for which an
    \eqn{A}-dependence statistic is required. This option is particularly
    useful for large values of \code{d}.}
  \item{freqname}{Character. Used only if \code{choice} = 1 and when \code{Xmat} is a
    matrix or a data-frame to identify the variable for the counts (frequencies).}
  \item{type}{"text" or "html"}
}
\value{Returns an object of class \code{list} containing the following
  components:
  \item{TA}{\eqn{A}-dependence statistics for each subset \eqn{A} of variables.}
  \item{fA}{degrees of freedom of the \eqn{A}-dependence statistics.}
  \item{pvalA}{\eqn{p}-values of the \eqn{A}-dependence statistics.}
  \item{X}{summary of the results.}
  \item{X2}{test statistic for mutual independence obtained by the
  sum of the \eqn{A}-dependence statistics, if \code{choice} = 1.}
  \item{Y2}{test statistic for serial independence obtained by the
  sum of the \eqn{A}-dependence statistics, if \code{choice} = 2.}
  \item{f}{number of degrees of freedom associated with the test
   statistic \code{X2} or \code{Y2}.}
  \item{pval}{the \eqn{p}-value associated with the test statistic \code{X2} or \code{Y2}.}
}
\author{Bilodeau M., Lafaye de Micheaux P.}

\references{Bilodeau M., Lafaye de Micheaux P. (2009). A-dependence
  statistics for mutual and serial independence of categorical
  variables, Journal of Statistical Planning and Inference, 139,
  2407-2419.

  Agresti A. (2002). Categorical data analysis, Wiley, p. 322
  
  Whisenant E.C., Rasheed B.K.A., Ostrer H., Bhatnagar  Y.M. (1991).
  Evolution and sequence analysis of a human Y-chromosomal
  DNA fragment, J. Mol. Evol., 33, 133-141.
}

\examples{
# Test of  mutual independence between 3 independent Bernoulli variables.

n <- 100
data <- data.frame(X1 = rbinom(n, 1, 0.3), X2 = rbinom(n, 1, 0.3) , X3 =
                   rbinom(n, 1, 0.3))
X <- table(data)
A.dep.tests(X)

# Test of mutual independence between 4 variables which are
# 2-independent and 3-independent, but are 4-dependent.

n <- 100
W <- sample(x = 1:8, size = n, TRUE)
X1 <- W \%in\% c(1, 2, 3, 5)
X2 <- W \%in\% c(1, 2, 4, 6)
X3 <- W \%in\% c(1, 3, 4, 7)
X4 <- W \%in\% c(2, 3, 4, 8)
data <- data.frame(X1, X2, X3, X4)
X <- table(data)
A.dep.tests(X)

# Test of serial independence of a nucleotide sequence of length
# 4156 described in Whisenant et al. (1991).

data(dna)
x2 <- dna[1]
for (i in 2:length(dna)) x2 <- paste(x2, dna[i], sep = "")
x <- unlist(strsplit(x2, ""))
x[x == "a" | x == "g"] <- "r"
x[x == "c" |  x== "t"] <- "y"
\donttest{
out <- A.dep.tests(x, choice = 2, d = 1501, m = 2)$TA[[1]]
plot(100:1500, out[100:1500], xlab = "lag j", ylab = "T(1,j+1)", pch = 19)
abline(h = qchisq(.995, df = 1))
}

# Analysis of a contingency table in Agresti (2002) p. 322 

data(highschool)
A.dep.tests(highschool, freqname = "count")

}
