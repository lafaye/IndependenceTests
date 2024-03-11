.mdcov.pval <- function(stat, Hn, N, X, vecd, weight.choice = 1, cubature = FALSE, a = 1, K = 100, epsrel = 10 ^ -6, norming, thresh.eigen = 10 ^ -8) {
    tmp <- .mdcov.lambdas(N, X, vecd, weight.choice = weight.choice, cubature, a, K = K, thresh.eigen = thresh.eigen)
#  pval <- davies(stat,lambda=sort(tmp$half.lambdas*c(1+tmp$abspks,1-tmp$abspks),decreasing=TRUE)[1:(2*trunc)])
    if (norming) {
        pval <- imhof(stat, lambda = tmp / Hn, epsrel = epsrel)$Qq
    } else {
        pval <- imhof(stat, lambda = tmp, epsrel = epsrel)$Qq
    }
    return(list(pval = pval, lambdas = tmp))
}
