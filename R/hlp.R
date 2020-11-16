#' Decolinearity by Variant Trimming 
#'
#' Retain only one of highly colinear variants.
#'
#' Decolinearity by variant trimming bring the correlation closer to full-rank,
#' Low rank LD correlation matrix causes inflated type I error for some summary
#' statistics.
#'
#' Tow or  more variants with  correlation higher  than \code{1 -  tol.cor} are
#' considered colienar, and  \code{dvt} marks all but one of  them for removal.
#' The    default   tolerence    of    high   correlation    (\code{tol.cor})is
#' \code{.Machine$double.eps}
#' 
#' @param C matrix of LD-correlation;
#' @param tol.cor tolerence of high correlation
#'
#' @return boolean mask of variants to retain.
#' @noRd
dvt <- function(C, tol.cor=NULL, ...)
{
    if(is.null(tol.cor))
        tol.cor <- sqrt(.Machine$double.eps)
    
    C <- abs(C)
    C[upper.tri(C, TRUE)] <- 0
    apply(C, 2, function(.) all(. < 1 - tol.cor))
}


#' Schur complement
#'
#' Given a full matrix \code{X} and a target block \code{C} within,
#' 
#' X = [A  B]
#'     [B' C],
#'
#' calculate the Schur complement \eqn{X/C = A - BC^{-1}B'}
#' 
#' @param X the whole matrix
#' @param C mask or index of the target block
#' @return the Schur complement of block C of matrix X
#'
#' @examples
#' ## get genotype and covariate matrices
#' gno <- readRDS(system.file("extdata", 'rs208294_gno.rds', package="dotgen"))
#' cvr <- readRDS(system.file("extdata", 'rs208294_cvr.rds', package="dotgen"))
#' 
#' ## overall correlation is a 56 x 56 positive definite matrix
#' X <- cor(cbind(gno, cvr))
#'
#' ## Schur complement for the covariate block
#' C <- 1:ncol(cvr) + ncol(gno)      # block of covariate
#' A <- 1:ncol(gno)                  # block of genotype
#' res1 <- solve(solve(X)[A, A])     # method 1
#' res2 <- scp(X, C)                 # method 2
#'
#' ## check equality
#' stopifnot(all.equal(res1, res2))  # must be TRUE
#' 
#' @noRd
scp <- function(X, C)
{
    ## C can be names, numbers, or Boolean masks, turn them into index
    i <- seq(ncol(X))
    names(i) <- colnames(X)
    i <- unname(i[C])

    ## A <- X[-C, -C]
    ## B <- X[-C, +C]
    ## B'<- X[+C, -C]
    ## C <- X[+C, +C]

    ## A - B C^{-1} B'
    X[-i, -i] - X[-i, +i] %*% solve(X[+i, +i], X[+i, -i])
}


#' Negative square root of a positive definite matrix
#'
#' Given a positive  semi-definite (PSD) matrix \code{X},  \code{nsp} gives the
#' negative square root of \code{X}, that is, \eqn{(R R')^{-1} = X}.
#'
#' For matrix \code{X}  that is not PSD, \code{nsp} trims  \code{X} and retains
#' \code{L} only positive eigenvalues.
#' 
#' @param X supposedly positive definite matrix
#' @return a list
#' \itemize{
#' \item{\code{W}: }{negative square root  of \code{X}}
#' \item{\code{L}: }{the effective number of positve eigenvalues}
#' }
#' @noRd
nsp <- function(X, L=NULL, eps=NULL, ...)
{
    if(is.null(eps))
        eps <- sqrt(.Machine$double.eps)

    ## eigen decomposition
    . <- eigen(X, TRUE)
    u <- .$vectors
    d <- .$values

    ## positive eigen values
    . <- d > d[1] * eps
    if(!all(.))
    {
        d <- d[  .]
        u <- u[, .]
    }
    L <- length(d)              # effective number of eigen
    
    ## square root
    d <- sqrt(1/d)
    H <- u %*% (d * t(u))       # U diag(d) U'
    H <- 0.5 * (H + t(H))
    
    list(H=H, L=L)
}
