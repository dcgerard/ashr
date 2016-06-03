#' Use control genes to estimate hidden confounders and variance
#' inflation parameter, then run ASH.
#'
#' This function will perform a variant of Removing Unwanted Variation
#' 4-step (RUV4) (Gagnon-Bartsch et al, 2013) where the control genes
#' are used not only to estimate the hidden confounders, but to
#' estimate a variance inflation parameter. This variance inflation
#' step is akin to the "empirical null" approach of Efron
#' (2004). After this procedure, Adaptive SHrinkage (ASH) (Stephens,
#' 2016) is performed on the coefficient estimates and the inflated
#' standard errors.
#'
#' The model is \deqn{Y = XB + ZA + E,} where \eqn{Y} is a matrix of
#' responses (e.g. log-transformed gene expression levels), \eqn{X} is
#' a matrix of covariates, \eqn{B} is a matrix of coefficients,
#' \eqn{Z} is a matrix of unobserved confounders, \eqn{A} is a matrix
#' of unobserved coefficients of the unobserved confounders, and
#' \eqn{E} is the noise matrix where the elements are independent
#' Gaussian and each column shares a common variance. The rows of
#' \eqn{Y} are the observations (e.g. individuals) and the columns of
#' \eqn{Y} are the response variables (e.g. genes).
#'
#' This model is fit using a two-step approach proposed in
#' Gagnon-Bartsch et al (2013) and described in Wang et al (2015),
#' modified to include estimating a variance inflation
#' parameter. Rather than use OLS in the second step of this two-step
#' procedure, we estimate the coefficients using Adaptive SHrinkage
#' (ASH) (Stephens, 2016). In the current implementation, only the
#' coefficients of one covariate can be estimated using ASH. The rest
#' are regressed out using OLS.
#'
#' @param Y A matrix of numerics. These are the response variables
#'     where each column has its own variance. In a gene expression
#'     study, the rows are the individuals and the columns are the
#'     genes.
#' @param X A matrix of numerics. The covariates of interest.
#' @param k A non-negative integer.The number of unobserved
#'     confounders. If not specified and the R package sva is
#'     installed, then this function will estimate the number of
#'     hidden confounders using the methods of Buja and Eyuboglu
#'     (1992).
#' @param cov_of_interest A positive integer. The column number of
#'     covariate in X whose coefficients you want to apply ASH to.
#' @param ctl A vector of logicals of length \code{ncol(Y)}. If
#'     position i is \code{TRUE} then position i is considered a
#'     negative control. If \code{ctl = NULL} (the default) then ASH
#'     will be run on the OLS estimates and corresponding standard
#'     errors.
#' @param ash_args A list of arguments to pass to ash. See
#'     \code{\link{ash.workhorse}} for details.
#' @param include_intercept A logical. If \code{TRUE}, then it will
#'     check \code{X} to see if it has an intercept term. If not, then
#'     it will add an intercept term. If \code{FALSE}, then \code{X}
#'     will be unchanged.
#' @param gls A logical. Should we use generalized least squares
#'     (\code{TRUE}) or ordinary least squares (\code{FALSE}) for
#'     estimating the confounders? The OLS version is equivalent to
#'     using RUV to estimate the confounders.
#' @param likelihood Either \code{"normal"} or \code{"t"}. If
#'     \code{likelihood = "t"}, then the user may provide the degrees
#'     of freedom by including a \code{df} element in
#'     \code{ash_args}. If \code{ash_args$df} is \code{NULL} then the
#'     degrees of freedom will be the sample size minus the number of
#'     covariates minus \code{k}.
#' @param limmashrink A logical. Should we apply hierarchical
#'     shrinkage to the variances (\code{TRUE}) or not (\code{FALSE})?
#' @param fa_func A factor analysis function. The function must have
#'     as inputs a numeric matrix \code{Y} and a rank (numeric scalar)
#'     \code{r}. It must output a numeric matrix \code{alpha} and a
#'     numeric vector \code{sig_diag}. \code{alpha} is the estimate of
#'     the coefficients of the unobserved confounders, so it must be
#'     an \code{r} by \code{ncol(Y)} matrix. \code{sig_diag} is the
#'     estimate of the column-wise variances so it must be of length
#'     \code{ncol(Y)}. The default is the function \code{pca_naive}
#'     that just uses the first \code{r} singular vectors as the
#'     estimate of \code{alpha}. The estimated variances are just the
#'     column-wise mean square.
#' @param fa_args A list. Additional arguments you want to pass to
#'     fa_func.
#'
#'
#' @return Except for the list \code{ruv}, the values returned are the
#'     exact same as in \code{\link{ash.workhorse}}. See that function
#'     for more details. Elements in the \code{ruv} list are:
#'
#'     \code{multiplier} A numeric. The estimated variance inflation parameter.
#'
#'     \code{betahat_ols} A vector of numerics. The ordinary least
#'     squares estimates of the coefficients of the covariate of
#'     interest. This is when not including the estiamted confounding
#'     variables.
#'
#'     \code{sebetahat_ols} A vector of positive numerics. The
#'     pre-inflation standard errors of \code{ruv$betahat} (NOT
#'     \code{ruv$betahat_ols}).
#'
#'     \code{betahat} A vector of numerics. The ordinary least squares
#'     estimates of the coefficients of the covariate of interest WHEN
#'     YOU ALSO INCLUDE THE ESTIMATES OF THE UNOBSERVED CONFOUNDERS.
#'
#'     \code{sebetahat} A vector of positive numerics. This is equal
#'     to sqrt(ruv$sebethat_ols * ruv$multiplier). This is the
#'     post-inflation adjusted standard errors for \code{ruv$betahat}.
#'
#'     \code{tstats} A vector of numerics. The t-statistics for
#'     testing against the null hypothesis of the coefficient of the
#'     covariate of interest being zero.
#'
#'     \code{pvalues} A vector of numerics. The p-values of said test
#'     above.
#'
#'     \code{alphahat} A matrix of numerics. The estimates of the
#'     coefficients of the hidden confounders.
#'
#'     \code{input} A list of arguments sent to
#'     \code{\link{ash.workhorse}}.
#'
#'     \code{sigma2} A vector of positive numerics. The estimates of
#'     the variances.
#'
#'     \code{fnorm_x} A numeric. This is the diagonal element of
#'     \code{t(X) \%*\% X} that corresponds to the covariate of
#'     interest. Returned mostly for debugging reasons and may be
#'     removed in the future.
#'
#' @export
#'
#' @author David Gerard
#'
#' @references Gagnon-Bartsch, J., Laurent Jacob, and Terence
#'     P. Speed. "Removing unwanted variation from high dimensional
#'     data with negative controls."
#'     Berkeley: Department of Statistics. University of California
#'     (2013).
#'
#'     Andreas Buja and Nermin Eyuboglu. "Remarks on parallel
#'     analysis." Multivariate behavioral research, 27(4):509–540,
#'     1992.
#'
#'     Bradley Efron
#'     "Large-Scale Simultaneous Hypothesis Testing: The Choice of a Null
#'     Hypothesis",
#'     Journal of the American Statistical Association, 99:465,
#'     96-104, 2004.
#'
#'     Stephens, Matthew. "False Discovery Rates: A New Deal." bioRxiv
#'     (2016): 038216.
#'
#'     Wang, J., Zhao, Q., Hastie, T., & Owen, A. B
#'     "Confounder Adjustment in Multiple Hypotheses Testing."
#'     arXiv preprint arXiv:1508.04178 (2015).
#'
ash_ruv <- function(Y, X, ctl = NULL, k = NULL,
                    cov_of_interest = ncol(X), ash_args = list(),
                    include_intercept = TRUE, gls = TRUE,
                    likelihood = c("normal", "t"),
                    limmashrink = FALSE, fa_func = pca_naive,
                    fa_args = list()) {

    if (is.null(ctl)) {
        message("No control genes provided so just doing OLS then ASH.")
        k <- 0
        ctl <- rep(FALSE, length = ncol(Y))
    }

    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::are_equal(ncol(Y), length(ctl))
    assertthat::assert_that(is.logical(ctl))
    assertthat::assert_that(cov_of_interest >= 1 & cov_of_interest <= ncol(X))
    assertthat::assert_that(is.logical(gls))
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.list(ash_args))
    assertthat::assert_that(is.null(ash_args$betahat))
    assertthat::assert_that(is.null(ash_args$sebetahat))
    assertthat::assert_that(is.logical(limmashrink))
    assertthat::assert_that(is.list(fa_args))
    assertthat::assert_that(is.null(fa_args$Y))
    assertthat::assert_that(is.null(fa_args$r))
    assertthat::assert_that(is.function(fa_func))

    likelihood <- match.arg(likelihood)

    if (include_intercept) {
        X_scaled <- apply(X, 2, function(x) {
            x / sqrt(sum(x ^ 2))
        })
        int_term <- rep(1, length = nrow(X)) / sqrt(nrow(X))

        any_int <- any(colSums( (int_term - X_scaled) ^ 2) < 10 ^ -14)
        if (!any_int) {
            X <- cbind(X, rep(1, length = nrow(X)))
        }
    }

    if (is.null(k)) {
        if (requireNamespace("sva", quietly = TRUE)) {
            message("Number of confounders not provided so being estimated with package sva.")
            k <- sva::num.sv(dat = t(Y), mod = X)
        } else {
            stop("If sva is not installed, then k needs to be provided. To install sva, run in R\n   source(\"https://bioconductor.org/biocLite.R\")\n   biocLite(\"sva\")")
        }
    }

    assertthat::assert_that(k + ncol(X) < nrow(X))

    if (k >= sum(ctl) & k != 0) {
        stop("k is larger than the number of control genes so model not identified.\nReduce k or increase the number of control genes.\nYou can also try out succotashr. To install succotashr, run in R:\n    install.packages(\"devtools\")\n    devtools::install_github(\"dcgerard/succotashr\")")
    }

    ## Place desired covariate as last covariate
    X <- X[, c( (1:ncol(X))[-cov_of_interest], cov_of_interest), drop = FALSE]
    cov_of_interest <- ncol(X)

    qr_x <- qr(X)
    ## multiply by sign so that it matches with beta_hat_ols
    Q <- qr.Q(qr_x, complete = TRUE) *
        sign(qr.R(qr_x)[cov_of_interest, cov_of_interest])
    ## discard first q-1 rows.
    Y_tilde <- crossprod(Q, Y)[cov_of_interest:nrow(Y), , drop = FALSE]

    ## Factor analysis using all but first row of Y_tilde
    fa_args$Y <- Y_tilde[2:nrow(Y_tilde), , drop = FALSE]
    fa_args$r <- k
    fa_out    <- do.call(what = fa_func, args = fa_args)
    alpha     <- fa_out$alpha
    sig_diag  <- fa_out$sig_diag

    ## make sure the user didn't screw up the factor analysis.
    assertthat::assert_that(is.vector(sig_diag))
    assertthat::are_equal(length(sig_diag), ncol(Y))
    assertthat::assert_that(all(sig_diag > 0))
    if (k != 0) {
        assertthat::assert_that(is.matrix(alpha))
        assertthat::are_equal(nrow(alpha), k)
        assertthat::are_equal(ncol(alpha), ncol(Y))
    } else {
        assertthat::assert_that(is.null(alpha))
    }

    ## Shrink variances if desired.
    if (requireNamespace("limma", quietly = TRUE) & limmashrink) {
        limma_out <- limma::squeezeVar(var = sig_diag,
                                       df = nrow(X) - ncol(X) - k)
        sig_diag <- limma_out$var.post
    } else if (!requireNamespace("limma", quietly = TRUE) & limmashrink) {
        stop("limmashrink = TRUE but limma not installed. To install limma, run in R:\n    source(\"https://bioconductor.org/biocLite.R\")\n    biocLite(\"limma\")")
    }

    ## absorb fnorm(X) into Y_tilde[1,], alpha, and sig_diag -----------------
    ## since dealt with sign earlier
    fnorm_x <- abs(qr.R(qr_x)[cov_of_interest, cov_of_interest])
    ## this is betahat from OLS, called Y1 in CATE.
    betahat_ols <- t(Y_tilde[1, , drop = FALSE]) / fnorm_x
    alpha_scaled <- alpha / fnorm_x
    ## this is se of betahat_ols if no confounders
    sig_diag_scaled <- sig_diag / (fnorm_x ^ 2)

    ## Use control genes to jointly estimate Z1 and variance scaling parameter.
    Yc <- betahat_ols[ctl, , drop = FALSE]
    if (k != 0) {
        alphac       <- alpha_scaled[ctl, , drop = FALSE]
        sig_diag_inv <- 1 / sig_diag_scaled[ctl]
        if (gls) {
            Z1 <- crossprod(solve(crossprod(alphac, sig_diag_inv * alphac)),
                            crossprod(alphac, sig_diag_inv * Yc))
        } else {
            Z1 <- crossprod(solve(crossprod(alphac, alphac)),
                            crossprod(alphac, Yc))
        }
        resid_mat <- Yc - alphac %*% Z1
        betahat <- betahat_ols - alpha_scaled %*% Z1
    } else {
        resid_mat <- Yc
        betahat <- betahat_ols
    }

    ## similar to MLE to UMVUE adjustment, divide by degrees of freedom.
    if (sum(ctl) != 0) {
        multiplier <- mean(resid_mat ^ 2 / sig_diag_scaled[ctl]) *
            nrow(X) / (nrow(X) - k - ncol(X))
    } else {
        multiplier <- 1
    }

    ## run ASH
    sebetahat <- sqrt(sig_diag_scaled * multiplier)

    ash_args$betahat   <- c(betahat)
    ash_args$sebetahat <- sebetahat

    if (likelihood == "t" & is.null(ash_args$df)) {
        ash_args$df <- nrow(X) - k - ncol(X)
    } else if (likelihood == "normal" & !is.null(ash_args$df)) {
        message("likelihood = \"normal\" but ash_args$df not NULL. Ignoring ash_args$df.")
        ash_args$df <- NULL
    } ## else, ash_args$df = NULL gives normal likelihood

    ash_out <- do.call(what = ash.workhorse, args = ash_args)

    ## Output frequentist values.
    ash_out$ruv <- list()
    ash_out$ruv$multiplier    <- multiplier
    ash_out$ruv$betahat_ols   <- betahat_ols
    ash_out$ruv$sebetahat_ols <- sqrt(sig_diag_scaled)
    ash_out$ruv$betahat       <- betahat
    ash_out$ruv$sebetahat     <- sebetahat
    ash_out$ruv$tstats        <- betahat / sebetahat
    ash_out$ruv$pvalues       <- 2 * (stats::pt(q = -abs(ash_out$ruv$tstats),
                                                df = nrow(X) - k - ncol(X)))
    ash_out$ruv$alphahat      <- alpha
    ash_out$input             <- ash_args
    ash_out$ruv$sigma2        <- sig_diag
    ash_out$ruv$fnorm_x       <- fnorm_x


    return(ash_out)
}



#' Basic PCA.
#'
#' Most of this code is from the package \code{cate}. I corrected some
#' problems. Specifically, I allow \code{r = 0} and I included a few
#' needed \code{drop = FALSE} terms. I also divide by \code{nrow(Y) -
#' r} rather than by \code{nrow(Y)}.
#'
#'
#' @param Y A matrix of numerics. The data.
#' @param r the rank.
#'
#' @author David Gerard
pca_naive <- function (Y, r) {
    if (r == 0) {
        Gamma <- NULL
        Z <- NULL
        Sigma <- apply(Y, 2, function(x) mean(x ^ 2))
    } else {
        svd_Y <- svd(Y)
        Gamma <- svd_Y$v[, 1:r, drop = FALSE] %*% diag(svd_Y$d[1:r], r, r) /
            sqrt(nrow(Y))
        Z <- sqrt(nrow(Y)) * svd_Y$u[, 1:r, drop = FALSE]
        Sigma <- apply(Y - Z %*% t(Gamma), 2, function(x) sum(x ^ 2)) / (nrow(Y) - r)
    }
    return(list(alpha = Gamma, Z = Z, sig_diag = Sigma))
}


#' Perform RUV to estimate confounders, estimate scale using simple
#' MLE, run ash.
#'
#' This is the depricated version and may be removed at any moment.
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
#'
#' @author David Gerard
#'
#' @references Andreas Buja and Nermin Eyuboglu. Remarks on parallel
#'     analysis. Multivariate behavioral research, 27(4):509–540,
#'     1992.
#'
ash_ruv_old <- function(Y, X, ctl, k = NULL, cov_of_interest = ncol(X),
                        ash_args = list(), include_intercept = TRUE) {


    if (!requireNamespace("ruv", quietly = TRUE)) {
        stop("R package ruv needs to be installed to run ash_ruv_old. To install, run in R:\n   install.packages(\"ruv\")")
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
        X_scaled <- apply(X, 2, function(x) {
            x / sqrt(sum(x ^ 2))
        })
        int_term <- rep(1, length = nrow(X)) / sqrt(nrow(X))

        any_int <- any(colSums( (int_term - X_scaled) ^ 2) < 10 ^ -14)
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

    ash_out <- do.call(what = ash.workhorse, args = ash_args)
    ash_out$multiplier <- multiplier
    ash_out$ruv <- ruvout


    return(ash_out)
}
