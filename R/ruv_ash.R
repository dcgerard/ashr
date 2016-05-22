#########
## Wrapper for ash + RUV inflating variances by simple MLE estimate.
#########

#' Perform RUV to estimate confounders, estimate scale using simple
#' MLE, run ash.
#'
#' @param Y A matrix of numerics. These are the response variables
#'     where each column has its own variance. In a gene expression
#'     study, the rows are the individuals and the columns are the
#'     genes.
#' @param X A matrix of numerics. The covariates of interest.
#' @param k A non-negative integer.The number of unobserved
#'     confounders. If not specified and the R package sva is
#'     installed, then this function will estimate the number of
#'     hidden confounders using th methods of Buja and Eyuboglu
#'     (1992).
#' @param cov_of_interest A positive integer. The column number of
#'     covariate in X whose coefficients you want to apply ASH to.
#' @param ctl A vector of logicals of length \code{ncol(Y)}. If
#'     position i is \code{TRUE} then position i is considered a
#'     negative control.
#' @param ash_args A list of arguments to pass to ash. See
#'     \code{\link{ash.workhorse}} for details.
#' @param include_intercept A logical. If \code{TRUE}, then it will
#'     check \code{X} to see if it has an intercept term. If not, then
#'     it will add an intercept term. If \code{FALSE}, then \code{X}
#'     will be unchanged.
#'
#' @export
#'
#' @author David Gerard
#'
#' @references Andreas Buja and Nermin Eyuboglu. Remarks on parallel
#'     analysis. Multivariate behavioral research, 27(4):509–540,
#'     1992.
#'
ash_ruv <- function(Y, X, ctl, k = NULL, cov_of_interest = ncol(X), ash_args = list(),
                    include_intercept = TRUE) {


    if (!requireNamespace("ruv", quietly = TRUE)) {
        stop("R package ruv needs to be installed to run ash_ruv. To install, run in R:\n   install.packages(\"ruv\")")
    }

    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::assert_that(cov_of_interest >= 1 & cov_of_interest <= ncol(X))

    if (is.null(k)) {
        if (requireNamespace("sva", quietly = TRUE)) {
            message("Number of confounders not provided so being estimated with package sva.")
            k <- sva::num.sv(dat = t(Y), mod = X)
        } else {
            stop("If sva is not installed, then k needs to be provided. To install sva, run in R\n   source(\"https://bioconductor.org/biocLite.R\")\n   biocLite(\"sva\")")
        }
    }

    if (include_intercept) {
        X_scaled <- apply(X, 2, function(x) { x / sqrt(sum(x ^ 2)) })
        int_term <- rep(1, length = nrow(X)) / sqrt(nrow(X))

        any_int <- any(colSums((int_term - X_scaled) ^ 2) < 10 ^ -14)
        if (!any_int) {
            X <- cbind(X, rep(1, length = nrow(X)))
        }
    }

    ruvout <- ruv::RUV4(Y = Y, X = X[, cov_of_interest, drop = FALSE],
                        Z = X[, -cov_of_interest, drop = FALSE],
                        ctl = ctl, k = k)

    betahat <- ruvout$betahat
    sebetahat <- sqrt(ruvout$sigma2 * ruvout$multiplier)
    t_stat <- betahat / sebetahat
    assertthat::are_equal(t_stat, ruvout$t)

    multiplier <- mean(t_stat[ctl] ^ 2)

    sebetahat_scaled <- sebetahat * sqrt(multiplier)

    ash_args$betahat   <- betahat
    ash_args$sebetahat <- sebetahat_scaled
    ## ash_args$df        <- ruvout$df

    ash_out <- do.call(what = ash.workhorse, args = ash_args)
    ash_out$multiplier <- multiplier
    ash_out$ruv <- ruvout


    return(ash_out)
}
