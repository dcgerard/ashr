

#' Wrapper for SVA to ASH pipeline with MAD inflation.
#'
#' @inheritParams ash_ruv
#'
#' @export
ash_sva_mad <- function(Y, X, k = NULL,
                        cov_of_interest = ncol(X),
                        likelihood = c("normal", "t"), ash_args = list(),
                        include_intercept = TRUE,
                        limmashrink = FALSE) {

    if (!requireNamespace(package = "sva", quietly = TRUE)) {
        stop("sva needs to be installed to run this function")
    } else if (!requireNamespace(package = "limma", quietly = TRUE)) {
        stop("limma needs to be installed to run this function")
    }

    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::assert_that(cov_of_interest >= 1 & cov_of_interest <= ncol(X))
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.list(ash_args))
    assertthat::assert_that(is.null(ash_args$betahat))
    assertthat::assert_that(is.null(ash_args$sebetahat))

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
        message("Number of confounders not provided so being estimated with package sva.")
        k <- sva::num.sv(dat = t(Y), mod = X)
    }

    ## Place desired covariate as last covariate
    X <- X[, c( (1:ncol(X))[-cov_of_interest], cov_of_interest), drop = FALSE]
    cov_of_interest <- ncol(X)

    trash <- capture.output(svaout <- sva::sva(dat = t(as.matrix(Y)), mod = X, n.sv = k))
    Xsv   <- cbind(X, svaout$sv)

    limout <- limma::lmFit(object = t(Y), design = Xsv)
    ash_args$betahat <- limout$coefficients[, cov_of_interest]
    if (limmashrink) {
        svout <- limma::eBayes(limout)
        sebetahat <- sqrt(svout$s2.post) * svout$stdev.unscaled[, cov_of_interest]
    } else {
        sebetahat <- limout$stdev.unscaled[, cov_of_interest] * limout$sigma
    }

    tstats <- ash_args$betahat / sebetahat
    multiplier <- mad(tstats)

    ash_args$sebetahat <- sebetahat * multiplier

    if (likelihood == "t") {
        ash_args$df <- nrow(Xsv) - ncol(Xsv)
    }

    ash_out <- do.call(what = ash, args = ash_args)

    tfinal <- tstats / multiplier
    pfinal <- 2 * pt(-abs(tfinal), df = nrow(Xsv) - ncol(Xsv))

    ash_out$cal               <- list()
    ash_out$cal$limma_fit     <- limout
    ash_out$cal$input         <- ash_args
    ash_out$cal$multiplier    <- multiplier
    ash_out$cal$bethat_old    <- ash_args$betahat
    ash_out$cal$sebetahat_old <- sebetahat
    ash_out$cal$sebetahat_new <- ash_args$sebetahat
    ash_out$cal$sva           <- svaout
    ash_out$cal$tstats        <- tfinal
    ash_out$cal$pvalues       <- pfinal

    return(ash_out)
}


#' Wrapper for RUV4 to ASH pipeline with MAD inflation.
#'
#' @inheritParams ash_ruv
#' @param ruv_type Should we run RUV4 (\code{"ruv4"}) or RUV2
#'     (\code{"ruv2"})?
#'
#' @export
ash_ruv_mad <- function(Y, X, ctl, ruv_type = c("ruv4", "ruv2"),
                         k = NULL,
                         cov_of_interest = ncol(X),
                         likelihood = c("normal", "t"), ash_args = list(),
                         include_intercept = TRUE,
                         limmashrink = FALSE) {

    if (!requireNamespace(package = "ruv", quietly = TRUE)) {
        stop("ruv needs to be installed to run this function")
    } else if (!requireNamespace(package = "limma", quietly = TRUE)) {
        stop("limma needs to be installed to run this function")
    }

    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::assert_that(cov_of_interest >= 1 & cov_of_interest <= ncol(X))
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.list(ash_args))
    assertthat::assert_that(is.null(ash_args$betahat))
    assertthat::assert_that(is.null(ash_args$sebetahat))

    likelihood <- match.arg(likelihood)
    ruv_type   <- match.arg(ruv_type)


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

    ## Place desired covariate as last covariate
    X <- X[, c( (1:ncol(X))[-cov_of_interest], cov_of_interest), drop = FALSE]
    cov_of_interest <- ncol(X)

    if (is.null(k)) {
        if (requireNamespace("sva", quietly = TRUE)) {
            message("Number of confounders not provided so being estimated with package sva.")
            k <- sva::num.sv(dat = t(Y), mod = X)
        } else {
            stop("If sva is not installed, then k needs to be provided. To install sva, run in R\n   source(\"https://bioconductor.org/biocLite.R\")\n   biocLite(\"sva\")")
        }
    }

    assertthat::assert_that(k + ncol(X) < nrow(X))

    if (k >= sum(ctl)) {
        stop("k is larger than the number of control genes so model not identified.\nReduce k or increase the number of control genes.\nYou can also try out succotashr. To install succotashr, run in R:\n    install.packages(\"devtools\")\n    devtools::install_github(\"dcgerard/succotashr\")")
    }

    if(ruv_type == "ruv4") {
        ruvout <- ruv::RUV4(Y = Y, X = X[, cov_of_interest, drop = FALSE],
                            ctl = ctl, k = k, Z = X[, -cov_of_interest, drop = FALSE])
    } else if (ruv_type == "ruv2") {
        ruvout <- ruv::RUV2(Y = Y, X = X[, cov_of_interest, drop = FALSE],
                            ctl = ctl, k = k, Z = X[, -cov_of_interest, drop = FALSE])
    }

    sebetahat_old    <- sqrt(ruvout$sigma2 * ruvout$multiplier)
    ash_args$betahat <- ruvout$betahat

    multiplier <- mad(ruvout$t)

    ash_args$sebetahat <- sebetahat_old * multiplier

    ash_out <- do.call(what = ash, args = ash_args)

    tfinal <- ruvout$t / multiplier
    pfinal <- 2 * pt(-abs(tfinal), df = nrow(X) - ncol(X) - k)

    ash_out$cal               <- list()
    ash_out$cal$ruv_fit       <- ruvout
    ash_out$cal$input         <- ash_args
    ash_out$cal$multiplier    <- multiplier
    ash_out$cal$bethat_old    <- ash_args$betahat
    ash_out$cal$sebetahat_old <- sebetahat_old
    ash_out$cal$sebetahat_new <- ash_args$sebetahat
    ash_out$cal$tstats        <- tfinal
    ash_out$cal$pvalues       <- pfinal

    return(ash_out)
}

#' Voom-limma-ebayes-mad-ash pipeline.
#'
#' Yuck.
#'
#' @export
#'
#' @inheritParams ash_ruv
vlema <- function(Y, X, cov_of_interest = ncol(X),
                  likelihood = c("normal", "t"), ash_args = list(),
                  include_intercept = TRUE) {

    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("limma needs to be installed to run this function")
    }

    assertthat::assert_that(is.matrix(Y))
    assertthat::assert_that(is.matrix(X))
    assertthat::are_equal(nrow(Y), nrow(X))
    assertthat::assert_that(cov_of_interest >= 1 & cov_of_interest <= ncol(X))
    assertthat::assert_that(is.logical(include_intercept))
    assertthat::assert_that(is.list(ash_args))
    assertthat::assert_that(is.null(ash_args$betahat))
    assertthat::assert_that(is.null(ash_args$sebetahat))

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

    ## Place desired covariate as last covariate
    X <- X[, c( (1:ncol(X))[-cov_of_interest], cov_of_interest), drop = FALSE]
    cov_of_interest <- ncol(X)

    vout   <- limma::voom(counts = t(Y), design = X)
    limout <- limma::lmFit(vout)
    eout   <- limma::eBayes(limout)

    told <- eout$t[, cov_of_interest]

    multiplier <- mad(told)

    tstats  <- told / multiplier
    pvalues <- 2 * pt(-abs(tstats), df = ncol(X))

    sebetahat_old <- eout$stdev.unscaled[, cov_of_interest]  * sqrt(eout$s2.post)

    ash_args$betahat <- eout$coefficients[, cov_of_interest]

    assertthat::are_equal(args_val$betahat / sebetahat_old,
                          eout$t[, cov_of_interest])

    ash_args$sebetahat <- sebetahat_old * multiplier

    ash_out <- do.call(what = ash, args = ash_args)

    ash_out$cal               <- list()
    ash_out$cal$voom          <- vout
    ash_out$cal$limma         <- limout
    ash_out$cal$ebayes        <- eout
    ash_out$cal$input         <- ash_args
    ash_out$cal$multiplier    <- multiplier
    ash_out$cal$bethat_old    <- ash_args$betahat
    ash_out$cal$sebetahat_old <- sebetahat_old
    ash_out$cal$sebetahat_new <- ash_args$sebetahat
    ash_out$cal$tstats        <- tstats
    ash_out$cal$pvalues       <- pvalues

    return(ash_out)
}
