#' Decorrelation by Orthogonal Transformation (DOT)
#'
#' [dot()]  decorrelates   genetic  association  test  statistics   by  special
#' symmetric orthogonal transformation.
#'
#' @details
#' Genetic  association studies  typically provide  per-variant test  statistics
#' that can be  converted to asymptotically normal, signed  Z-scores. Once those
#' Z-scores are transformed to independent random variables, various methods can
#' be applied to combine them and obtain SNP-set overall association.
#' 
#' [dot()]  uses  per-variant  genetic  association  test  statistics  and  the
#' correlation among them to decorrelate Z-scores.
#'
#' To estimate the  correlation among genetic association  test statistics, use
#' [cst()].  If  P-values and  estimated effects  (i.e, beta  coefficients) are
#' given instead  of test statistics, [zsc()]  can be used to  recover the test
#' statistics (i.e., Z-scores).
#'
#' `tol.cor`: variants  with correlation too close  to 1 in absolute  value are
#' considered to be collinear  and only one of them will  be retained to ensure
#' that  the  LD  matrix  is   full-rank.   The  maximum  value  for  tolerable
#' correlation  is  1   -  `tol.cor`.  The  default  value   for  `tol.cor`  is
#' `sqrt(.Machine$double.eps)`.
#'
#' `tol.egv`: negative and close to  zero eigenvalues are truncated from matrix
#'  `D` in `H = EDE'`. The corresponding  columns of `E` are also deleted. Note
#'  the  the dimention  of the  square matrix  `H` does  not change  after this
#'  truncation. See DOT publication in the  reference below for more details on
#'  definitions  of `E`  and `D`  matrices.  The  default eigenvalue  tolerance
#'  value is `sqrt(.Machine$double.eps)`.
#' 
#' A number of methods are available for combining de-correlated P-values,
#' see [dot_sst] for details.
#'
#' @references
#' \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007819}{
#' Vsevolozhskaya, O. A., Shi, M., Hu, F., & Zaykin, D. V. (2020). DOT: Gene-set
#' analysis by combining decorrelated association statistics. PLOS Computational
#' Biology, 16(4), e1007819.}
#'
#' @param Z vector of association test statistics (i.e., Z-scores).
#' @param C correlation matrix among the association test statistics, as
#'     obtained by [cst()].
#' @param tol.cor tolerance threshold for the largest correlation absolute value.
#' @param tol.egv tolerance threshold for the smallest eigenvalue.
#' @param ... additional parameters.
#'
#' @return
#' a  list with
#' \itemize{
#' \item{`X`:} {association test statistics, de-correlated.}
#' \item{`H`:} {orthogonal transformation, such that `X = H %*% Z`.}
#' \item{`M`:} {effective number of variants after de-correlation.}
#' \item{`L`:} {effective number of eigenvalues after truncation.}
#' }
#' 
#' @seealso [cst()], [zsc()], [dot_sst]
#' @examples
#' ## get genotype and covariate matrices
#' gno <- readRDS(system.file("extdata", 'rs208294_gno.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#'
#' ## estimate the correlation among association test statistics
#' sgm <- cst(gno, cvr)
#'
#' ## get the result of genetic association analysis (P-values and effects)
#' res <- readRDS(system.file("extdata", 'rs208294_res.rds', package="dotgen"))
#'
#' ## recover Z-score statistics
#' stt <- with(res, zsc(P, BETA))
#'
#' ## decorrelate Z-scores by DOT
#' result <- dot(stt, sgm)
#' print(result$X)          # decorrelated statistics
#' print(result$H)          # orthogonal transformation
#'
#' ## sum of squares of decorrelated statistics is a chi-square
#' ssq <- sum(result$X^2)
#' pvl <- 1 - pchisq(ssq, df=result$L)
#'
#' print(ssq)            # sum of squares = 35.76306
#' print(pvl)            # chisq P-value =  0.001132132
#' @export
dot <- function(Z, C, tol.cor=NULL, tol.egv=NULL, ...)
{
    if(is.null(tol.cor))
        tol.cor <- sqrt(.Machine$double.eps)
    if(is.null(tol.egv))
        tol.egv <- sqrt(.Machine$double.eps)

    ## trim collinear variants
    m <- dvt(C, tol.cor)
    M <- sum(m)                      # effective number of variants
    C <- C[m, m]
    Z <- Z[m]

    ## get orthogonal transformation
    d <- nsp(C, eps=tol.egv, ...)
    H <- d$H                         # orthogonal transformation
    L <- d$L                         # effective number of eigenvalues
    
    X <- H %*% Z                     # decorrelated statistics

    list(H=H, X=X, L=L, M=M)
}

#' Calculate Z-scores from P-values and estimated effects
#'
#' [zsc()] recovers Z-scores from  P-values and corresponding effect directions
#' (or beta coefficients) reported by a genetic association analysis.
#'
#' @details
#'
#' For any  genetic variant,  its two-sided  P-value (\eqn{p})  and the  sign of
#' estimated effect (\eqn{\beta}) is used to recover the Z-score (\eqn{z}), that
#' is, \eqn{z = sign(\beta) \Phi^{-1}(p/2)}.
#'
#' @seealso [dot()]
#'
#' @param P vector of P-values.
#' @param BETA vector of effect directions or beta coefficients.
#' @return A vector of Z-scores.
#' @examples
#' ## result of per-variant analysis (P-values and estimated effects)
#' res <- readRDS(system.file("extdata", 'rs208294_res.rds', package="dotgen"))
#'
#' ## recover Z-score statistics
#' stt <- with(res, zsc(P, BETA))
#'
#' ## checking
#' stopifnot(all.equal(pnorm(abs(stt), lower.tail = FALSE) * 2, res$P))
#' 
#' @export
zsc <- function(P, BETA)
{
    qnorm(P / 2) * (sign(BETA))
}


#' Correlation among association test statistics
#'
#' Calculates the correlation among genetic association test statistics.
#'
#' @details
#' When no covariates are present in per-variant association analyses, that is,
#' `x==NULL`, correlation among test statistics  is the same as the correlation
#' among variants, `cor(g)`.
#'
#' With  covariates, correlation  among  test  statistics is  not  the same  as
#' `cor(g)`. In this case, [cst()] takes  the generalized inverse of the entire
#' correlation  matrix, `corr(cbind(g,  x))`, and  then inverts  back only  the
#' submtarix containing genotype variables, `g`.
#'
#' If Z-scores were calculated based on genotypes with some missing values, the
#' correlation among test statistics will be  reduced by the amount that can be
#' theoretically derived. It can be shown  that this reduced correlation can be
#' calculated by imputing  the missing values with the  averages of non-missing
#' values. Therefore, by default, [cst()]  fills missing values in each variant
#' with  the  average  of  non-missing  values  in  that  same  variant  (i.e.,
#' imputation  by  average, [imp_avg()]).  Other  imputation  methods are  also
#' available (see topic [imp] for other techniques that may improve power), but
#' note that  techniques other than the  imputation by average requires  one to
#' re-run  the  association  analyses  with  imputed  variants  to  ensure  the
#' correlation among new statistics (i.e.,  Z-scores) and the correlation among
#' imputed variants are identical. Otherwise, Type  I error may be inflated for
#' decorrelation-based methods.
#'
#' @param g matrix of  genotype, one row per  sample, one column  per variant,
#'     missing values allowed.
#' @param x matrix of covariates, one row per sample, no missing values allowed.
#'
#' @return Correlation matrix among association test statistics.
#'
#' @seealso [imp], [imp_avg()]
#' @examples
#' ## get genotype and covariate matrices
#' gno <- readRDS(system.file("extdata", 'rs208294_gno.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#'
#' ## correlation among association statistics, covariates involved
#' res <- cst(gno, cvr)
#' print(res[1:4, 1:4])
#'
#' ## genotype matrix with 2% randomly missing data
#' g02 <- readRDS(system.file("extdata", 'rs208294_g02.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#' res <- cst(g02, cvr)
#' print(res[1:4, 1:4])
#' 
#' @export
cst <- function(g, x=NULL)
{
    ## impute missed calls by naive average
    g <- imp_avg(g)
    
    if(is.null(x))
        r <- stats::cor(g)              # no covariate, use full cor
    else
    {
        r <- stats::cor(cbind(x, g))    # full cor
        i <- seq(ncol(x))               # index of covariates
        r <- stats::cov2cor(scp(r, i))  # cond cor
    }
    r
}
