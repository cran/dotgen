#' Methods for combining decorrelated summary statistics
#'
#' Decorrelates and combines per-variant genetic association test statistics, `Z`,
#'  given the correlation matrix among them, `C`.
#'
#' @details
#' 
#' These functions  first call [dot()]  to decorrelate the  genetic association
#' test  statistics and  then provide  various options  to combine  independent
#' statistics or corresponding P-values into the overall statistic and P-value.
#'
#' The two  rank truncated  tests (i.e.,  [dot_art()], [dot_rtp()])  require an
#' additional parameter `k` that specifes the number of smallest (decorrelated)
#' P-values to combine. By default, `k`  equals half of the number of variants.
#' The adaptive  rank truncation  method, [dot_arta()], determines  the optimal
#' truncation value between 1 and `k`.
#'
#' The truncated  product method,  [dot_tpm()], combines  P-values at  least as
#' small as `tau` (0.05 by default).  If  `tau` is equal to 1, then [dot_tpm()]
#' provides  the  same result  as  [dot_fisher()]  (i.e., Fisher's  method  for
#' combining  P-values). Similarly,  if `k`  is equal  to the  total number  of
#' tests, the results  of [dot_art()] and [dot_rtp()] will be  the same as that
#' of [dot_fisher()].
#'
#' Reference  (\strong{a})  below  details  how to  combine  decorrelated  test
#' statistics  or  P-values  via  [dot_art()],  [dot_rtp()]  and  [dot_arta()];
#' reference (\strong{b}) details [dot_tpm()] method.
#'
#' @param Z vector of association test statistics (i.e., Z-scores).
#' @param  C matrix of  correlation among the  test statistics, as  obtained by
#'     [cst()].
#' @param k combine `k` smallest (decorrelated) P-values.
#'
#' @param ... additional parameters
#'
#' @return a list of
#' \itemize{
#' \item{`X`:} {decorrelated  association statistics.}
#' \item{`H`:} {orthogonal transformation, such that `X = H %*% Z`.}
#' \item{`Y`:} {the overall combined statistic.}
#' \item{`P`:} {the P-value corresponding to \code{Y}.}
#' }
#'
#' @references
#' (a) \href{https://www.frontiersin.org/articles/10.3389/fgene.2019.01051/full}{
#' Vsevolozhskaya, O.   A., Hu, F., &  Zaykin, D.  V. (2019).   _Detecting weak
#' signals  by  combining  small  P-values  in  genetic  association  studies._
#' Frontiers in genetics, 10, 1051.}
#'
#' (b) \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.0042}{Zaykin,
#' D.    V.,   Zhivotovsky,    L.     A.,   Westfall,    P.    H.,   &    Weir,
#' B. S.  (2002). _Truncated  product method  for combining  P-values._ Genetic
#' Epidemiology, 22(2), 170-185.}
#'
#' @seealso [dot()]
#' 
#' @examples
#' ## get the test statistics and pre-calculated LD matrix
#' stt <- readRDS(system.file("extdata", 'art_zsc.rds', package="dotgen"))
#' sgm <- readRDS(system.file("extdata", 'art_ldm.rds', package="dotgen"))
#'
#' @name dot_sst
NULL

#' @describeIn dot_sst
#'
#' decorrelation followed by a Chi-square test.
#'
#' @examples
#'
#' ## decorrelated chi-square test
#' result <- dot_chisq(stt, sgm)
#' print(result$Y)  # 37.2854
#' print(result$P)  # 0.0003736988
#' @export
dot_chisq <- function(Z, C, ...)
{
    ret <- dot(Z, C, ...)        # decorrelate
    L <- ret$L                   # effective number of eigenvalues

    Y <- sum(ret$X^2)            # sum of squares
    P <- 1 - pchisq(Y, L)        # a single p-value
    c(list(P=P, Y=Y), ret)
}


#' @describeIn dot_sst
#'
#' decorrelated Fisher's combined P-value test.
#'
#' @examples
#'
#' ## decorrelated Fisher's combined P-value chi-square test
#' result <- dot_fisher(stt, sgm)
#' print(result$Y)  # 58.44147
#' print(result$P)  # 0.0002706851
#' @export
dot_fisher <- function(Z, C, ...)
{
    ret <- dot(Z, C, ...)               # decorrelate
    L <- ret$L                          # effective number of eigenvalues
    P <- 1 - pchisq(ret$X^2, df=1)      # decorrelated two-tail P-values

    Y <- -2 * sum(log(P))               # Fisher statistics
    P <- 1 - pchisq(Y, 2 * L)           # summarized P-value
    c(list(P=P, Y=Y), ret)
}

#' @describeIn dot_sst
#'
#' decorrelated Augmented Rank Truncated (ART) test. 
#'
#' @examples
#'
#' ## decorrelated augmented rank truncated (ART) test.
#' result <- dot_art(stt, sgm, k=6)
#' print(result$Y)  # 22.50976
#' print(result$P)  # 0.0006704994
#' @export
dot_art <- function(Z, C, k=NULL, ...)
{
    ret <- dot(Z, C, ...)               # decorrelate
    L <- ret$L                          # effective number of eigenvalues
    M <- ret$M                          # effective number of variants
    P <- 1 - pchisq(ret$X^2, df=1)      # decorrelated two-tail P-values
    P <- sort(P)                        # sort

    if(is.null(k))                      # k = M / 2 (default)
        k <- round(M * .5)              #
    k <- min(L, k)                      # k <= L

    if(any(P == 0))
        return(c(list(P=0, Y=1, k=k, ret)))

    ## summarized statistics: shape parameter
    S <- (k - 1) * (digamma(L + 1) - digamma(k))

    ## summarized statistics: summation
    Y <- 0
    Y <- Y + (k - 1) * log(P[k])
    Y <- Y - sum(log(P[1:(k - 1)]))
    Y <- Y + qgamma(1 - pbeta(P[k], k, L - k + 1), S)

    ## summarized statistics: P-value
    P <- 1 - pgamma(Y, k + S - 1)
    
    c(list(P=P, Y=Y), ret)
}


#' @describeIn dot_sst
#'
#' decorrelated Augmented Rank Truncated Adaptive (ARTA) test.
#'
#' @return
#' for Augmented Rank Truncated Adaptive (ARTA) test,
#' \itemize{
#' \item{k:} {the number of decorrelated P-values that were adaptively picked.}}
#'
#' @param w weight assigned to partial sums in ARTA implementation; default is 1.
#' 
#' @examples
#'
#' ## decorrelated Augmented Rank Truncated Adaptive (ARTA) test
#' result <- dot_arta(stt, sgm, k=6)
#' print(result$Y)  # -1.738662
#' print(result$k)  #  5 smallest P-values are retained
#' print(result$P)  #  0.003165 (varies)
#' @export
dot_arta <- function(Z, C, k=NULL, w=NULL, ...)
{
    ret <- dot(Z, C, ...)               # decorrelate
    L <- ret$L                          # effective number of eigenvalues
    M <- ret$M                          # effective number of variants
    P <- 1 - pchisq(ret$X^2, df=1)      # decorrelated two-tail P-values
    P <- sort(P)                        # sorted

    if(is.null(k))                      # k = L/2 (default)
        k <- round(M * 0.5)
    k <- min(k, L)                      # k <= L

    if(is.null(w))
        w <- rep(1, k)
    w <- w[1:k]                         # first k weights

    ## z_i = [(1 - p_i) / (1 - p_{i-1})]^(L - i + 1), i = 1 .. L, p_0 = 1
    z <- ((1 - P) / c(1, 1 - P[-M]))^(M:1)
    
    ## accumulate the top-ranking k statistics, apply weights
    p <- (1 - z)[1:k]
    q <- qnorm(p) * w

    if(max(q) == Inf)
        return(c(list(P=1, k=which(q==Inf)[1], Y=Inf), ret))
    
    sumQ <- cumsum(qnorm(p) * w)

    ## L(1) %*% diag(w^2) %*% U(1)
    pSg <- matrix(cumsum(w^2), k, k)
    pSg[lower.tri(pSg)] <- t(pSg)[lower.tri(pSg)]
    pCr <- cov2cor(pSg)

    sQ <- sumQ / sqrt(diag(pSg))

    Y <- max(sQ)
    P <- try(mvtnorm::pmvnorm(lower=rep(-Inf, k), upper=rep(max(sQ), k), sigma=pCr))
    if(inherits(P, 'try-error'))
    {
        print("bad")
        P <- 1
    }
    else
        P <- P[1]
    c(list(P=P, k=which.min(sQ), Y=Y), ret)
}


#' @describeIn dot_sst
#'
#' decorrelated Rank Truncated Product (RTP) test.
#'
#' @examples
#'
#' ## decorrelated Rank Truncated Product (RTP)
#' result <- dot_rtp(stt, sgm, k=6)
#' print(result$Y)  # 22.6757
#' print(result$P)  # 0.0007275518
#' @export
dot_rtp <- function(Z, C, k=NULL, ...)
{
    ret <- dot(Z, C, ...)               # decorrelate
    L <- ret$L                          # effective number of eigenvalues
    M <- ret$M                          # effective number of variants
    P <- 1 - pchisq(ret$X^2, df=1)      # decorrelated two-tail P-values
    P <- sort(P)                        # sorted
    if(is.null(k))                      # k = M / 2 (default)
        k <- round(M * .5)
    k <- min(L, k)                      # k <= L

    Y <- sum(-log(P[1:k]))
    P <- integrate(function(x, y, m, n)
    {
        1 - pgamma(log(qbeta(x, m + 1, n - m)) * m + y, m)
    },
    0, 1, Y, k, L)$va
    c(list(P=P, Y=Y), ret)
}


#' @describeIn dot_sst
#'
#' decorrelated Truncated Product Method (TPM) test.
#'
#' @examples
#'
#' ## decorrelated Truncated Product Method (TPM)
#' result <- dot_tpm(stt, sgm, tau=0.05)
#' print(result$Y)  #  1.510581e-08
#' print(result$k)  #  6 P-values <= tau
#' print(result$P)  #  0.0007954961
#'
#' @param tau combine (decorrelated) P-values no large than tau; default is 0.05.
#' @return
#' for Truncated Product Method (TPM),
#' \itemize{
#' \item{k:} {the number of decorrelated P-values \eqn{\le} \code{tau}.}}
#' @export
dot_tpm <- function(Z, C, tau=0.05, ...)
{
    if(is.null(tau))                    # k = L / 2 (default)
        tau <- 0.05

    ret <- dot(Z, C, ...)               # decorrelate
    L <- ret$L                          # effective number of eigenvalues
    P <- 1 - pchisq(ret$X^2, 1)         # decorrelated two-tail P-values
    tau <- min(tau, max(P))

    k <- sum(P <= tau)
    if(k == 0)
        return(c(list(P=1, Y=1, k=k), ret))
    if(any(P == 0))
        return(c(list(P=0, Y=1, k=k), ret))

    lw <- sum(log(P[P <= tau]))         # log(w), w = prod(p_i | p_i<=tau)
    
    ## sequence: k = 1 .. L
    K <- seq(L)

    ## log Pr(k) = log(L choose k) + log((1 - tau)^(L - k), k = 1..L
    lpk <- lchoose(L, K) + (L - K) * log(1 - tau)
    
    ## log Pr(W <= w | k)
    lpw <- double(L)
    klt <- K * log(tau)                 # k ln(tau), k = 1..L
    for(k in K[lw <= klt])              # case 1: ln(w) <= k ln(tau)
    {
        S <- seq(0, k - 1)
        ## sum_{s=0..k-1}(k log(tau) - log(w))^s / s!
        . <- S * log(klt[k] - lw) - lfactorial(S)
        . <- exp(.)
        lpw[k] <- lw + log(sum(.))
    }
    lpw[lw > klt] <- klt[lw > klt]      # case 2: ln(w)  > k ln(tau)
    
    ## sum Pr(k) Pr(W <= w|k), k = 1 .. L
    P <- sum(exp(lpk + lpw))
    Y <- exp(lw)
    
    c(list(P=P, Y=Y, k=k), ret)
}
