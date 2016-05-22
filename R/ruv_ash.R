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

    ## Place desired covariate as last covariate
    X <- X[, c((1:ncol(X))[-cov_of_interest], cov_of_interest), drop = FALSE]
    cov_of_interest <- ncol(X)

    qr_x <- qr(X)
    ## multiply by sign so that it matches with beta_hat_ols
    Q <- qr.Q(qr_x, complete = TRUE) * sign(qr.R(qr_x)[cov_of_interest, cov_of_interest])
    Y_tilde <- crossprod(Q, Y)[cov_of_interest:nrow(Y), , drop = FALSE]  # discard first q-1 rows.

    ## Factor analysis using all but first row of Y_tilde
    pca_out <- pca_naive(Y = Y_tilde[2:nrow(Y_tilde), , drop = FALSE], r = k)
    alpha <- pca_out$Gamma
    sig_diag <- pca_out$Sigma

    ## absorb fnorm(X) into Y_tilde[1,], alpha, and sig_diag -------------------
    fnorm_x <- abs(qr.R(qr_x)[cov_of_interest, cov_of_interest])  ## since dealt with sign earlier
    betahat_ols <- t(Y_tilde[1, , drop = FALSE]) / fnorm_x ## this is betahat from OLS, called Y1 in CATE.
    alpha_scaled <- alpha / fnorm_x
    sig_diag_scaled <- sig_diag / (fnorm_x ^ 2) ## this is se of betahat ols if no confounders

    ## Use control genes to jointly estimate Z1 and variance scaling parameter.
    Yc <- betahat_ols[ctl, , drop = FALSE]
    if (k != 0) {
        alphac     <- alpha_scaled[ctl, , drop = FALSE]
        Sigmac_inv <- diag(1 / sig_diag_scaled[ctl])
        Z1 <- solve(t(alphac) %*% Sigmac_inv %*% alphac) %*% t(alphac) %*% Yc
        resid_mat <- Yc - alphac %*% Z1
        betahat <- betahat_ols - alpha_scaled %*% Z1
    } else {
        resid_mat <- Yc
        betahat <- betahat_ols
    }
    multiplier <- mean(resid_mat ^ 2 / sig_diag_scaled[ctl])

    ## run ASH
    sebetahat <- sqrt(sig_diag_scaled * multiplier)

    ash_args$betahat   <- betahat
    ash_args$sebetahat <- sebetahat
    ## ash_args$df        <- ruvout$df

    ash_out <- do.call(what = ash.workhorse, args = ash_args)
    ash_out$ruv <- list()
    ash_out$ruv$multiplier    <- multiplier
    ash_out$ruv$betahat_ols   <- betahat_ols
    ash_out$ruv$sebetahat_ols <- sqrt(sig_diag_scaled)
    ash_out$ruv$betahat       <- betahat
    ash_out$ruv$sebetahat     <- sebetahat
    ash_out$ruv$alphahat      <- alpha

    return(ash_out)
}



#' Basic PCA.
#'
#' Most if not all of code is from package \code{cate}. This is mostly
#' so people don't have to install sva and leapp if they want to use
#' it.
#'
#'
#' @param Y A matrix of numerics. The data.
#' @param r the rank.
pca_naive <- function (Y, r) {
    if(r == 0) {
        Gamma <- NULL
        Z <- NULL
        Sigma <- apply(Y, 2, function(x) mean(x ^ 2))
    } else {
        svd_Y <- svd(Y)
        Gamma <- svd_Y$v[, 1:r, drop = FALSE] %*% diag(svd_Y$d[1:r], r, r) / sqrt(nrow(Y))
        Z <- sqrt(nrow(Y)) * svd_Y$u[, 1:r, drop = FALSE]
        Sigma <- apply(Y - Z %*% t(Gamma), 2, function(x) mean(x ^ 2))
    }
    return(list(Gamma = Gamma, Z = Z, Sigma = Sigma))
}
