#' Find log-component density convolution when error distribution is a mixture.
#'
#' This function will find the log of the coefficients of the prior
#' mixing proportions in the ASH-likelihood when the model is
#' \eqn{betahat = g + errordist}, where \eqn{g} and \eqn{errordist}
#' are both mixtures of either normals or uniforms.
#'
#' There are four cases: (1) \code{g} and \code{errordist} are both of
#' class \code{normalmix}, (2) \code{g} is of class \code{normalmix}
#' and \code{errordist} is of class \code{unifmix}, (3) \code{g} is of
#' class \code{unifmix} and \code{errordist} is of class
#' \code{normalmix}, and (4) \code{g} and \code{errordist} are both of
#' class \code{unifmix}. All of these are supported. Though (2) and
#' (3) differ only in indexing.
#'
#' @param g A mixture density. Either of class \code{normalmix} or of
#'     class \code{unimix}. This is the prior.
#' @param errordist A list of objects of either class \code{normalmix}
#'     or \code{unimix}. The length of this list must be the length of
#'     \code{betahat}. \code{errordist[[i]]} is the \eqn{i}th error
#'     distribution of \code{betahat[i]}.
#' @param betahat A vector of numerics. The locations at which to
#'     evalutate the density.
#'
#' @return A matrix with row dimension \code{length(betahat)} and
#'     column dimension \code{length(errordist)}.
#'
#' @author David Gerard
log_compdens_conv_mix <- function(g, betahat, errordist) {
    class_e <- unique(sapply(errordist, class))
    class_g <- class(g)

    assertthat::are_equal(length(class_e), 1)
    assertthat::are_equal(length(class_g), 1)
    assertthat::assert_that(is.list(errordist))
    assertthat::are_equal(length(betahat), length(errordist))

    pisum <- sapply(errordist, FUN = function(x) { sum(x$pi) })
    assertthat::assert_that(all(pisum == 1))

    ## make sure error distribution doesn't have point mass on zero
    if (class_e == "normalmix") {
        assertthat::assert_that(!any(sapply(errordist, FUN = function(x) { any(x$sd == 0)})))
    } else if (class_e == "unimix") {
        assertthat::assert_that(!any(sapply(errordist, FUN = function(x) { any(x$a == 0 & x$b == 0)})))
    }

    if (class_g == "normalmix" & class_e == "normalmix") {
        matrix_llik <- log_compdens_conv_mix.normalnormal(g, betahat, errordist)
    } else if (class_g == "normalmix" & class_e == "unimix") {
        matrix_llik <- log_compdens_conv_mix.normaluni(g, betahat, errordist)
    } else if (class_g == "unimix" & class_e == "normalmix") {
        matrix_llik <- log_compdens_conv_mix.uninormal(g, betahat, errordist)
    } else if (class_g == "unimix" & class_e == "unimix") {
        matrix_llik <- log_compdens_conv_mix.uniuni(g, betahat, errordist)
    } else {
        stop("Error: either g or errordist is of an unsupported class.")
    }
    return(matrix_llik)
}

#' Normal-mixture - normal-mixture convolution.
#'
#' @inheritParams log_compdens_conv_mix
#'
#' @author David Gerard
log_compdens_conv_mix.normalnormal <- function(g, betahat, errordist) {
    class_e <- unique(sapply(errordist, class))
    class_g <- class(g)

    assertthat::are_equal(length(class_e), 1)
    assertthat::are_equal(length(class_g), 1)
    assertthat::are_equal(class_e, "normalmix")
    assertthat::are_equal(class_g, "normalmix")
    assertthat::are_equal(length(betahat), length(errordist))

    n <- length(betahat)
    p <- length(g$mean)
    matrix_llik <- matrix(NA, nrow = p, ncol = n)
    for (index in 1:n) {
        matrix_llik[, index] <- nnconv_dense(g = g,
                                             errordist = errordist[[index]],
                                             x = betahat[index])
    }
    return(matrix_llik)
}


#' Normal-normal convolution log-density evaluated at a single point.
#'
#' @param g A \code{normalmix} object. This is the prior mixture.
#' @param errordist A \code{normalmix} object. This is the error mixture.
#' @param x A single numeric. This is the data.
#'
#' @author David Gerard
#'
#' @seealso \code{\link{log_compdens_conv_mix.normalnormal}}
nnconv_dense <- function(g, errordist, x) {
    assertthat::are_equal(class(g), "normalmix")
    assertthat::are_equal(class(errordist), "normalmix")
    assertthat::are_equal(length(x), 1)

    new_means <- outer(errordist$mean, g$mean, FUN = "+")
    new_sd    <- sqrt(outer(errordist$sd ^ 2, g$sd ^ 2, FUN = "+"))

    log_vals <- stats::dnorm(x = x, mean = new_means, sd = new_sd, log = TRUE)

    if (length(dim(log_vals)) == 2) {
        max_log_vals <- apply(log_vals, 2, max)
        like_vals <- log(colSums(exp(log_vals - max_log_vals) * errordist$pi)) + max_log_vals
    } else if (length(dim(log_vals)) == 0 ) {
        ## working with a scalar
        like_vals <- log_vals
    } else {
        stop("Contact David Gerard. This is a scenario I wasn't expecting")
    }

    return(like_vals)
}


#' Normal-mixture - uniform-mixture convolution.
#'
#' @inheritParams log_compdens_conv_mix
#'
#' @author David Gerard
log_compdens_conv_mix.normaluni <- function(g, betahat, errordist) {
    class_e <- unique(sapply(errordist, class))
    class_g <- class(g)

    assertthat::are_equal(length(class_e), 1)
    assertthat::are_equal(length(class_g), 1)
    assertthat::are_equal(class_e, "unimix")
    assertthat::are_equal(class_g, "normalmix")
    assertthat::are_equal(length(betahat), length(errordist))

    n <- length(betahat)
    p <- length(g$mean)
    matrix_llik <- matrix(NA, nrow = p, ncol = n)
    for (index in 1:n) {
        matrix_llik[, index] <- nuconv_dense(g = g,
                                             errordist = errordist[[index]],
                                             x = betahat[index])
    }
    return(matrix_llik)
}


#' Normal-uniform convolution log-density evaluated at a single point.
#'
#' @param g A \code{normalmix} object. This is the prior mixtrue.
#' @param errordist A \code{unimix} object. This is the error mixture.
#' @param x A single numeric. This is the data.
#'
#' @author David Gerard
#'
#' @seealso \code{\link{log_compdens_conv_mix.normaluni}}
nuconv_dense <- function(g, errordist, x) {
    assertthat::are_equal(class(g), "normalmix")
    assertthat::are_equal(class(errordist), "unimix")
    assertthat::are_equal(length(x), 1)


    left_means  <- outer(errordist$b, g$mean - x, FUN = "+")
    right_means <- outer(errordist$a, g$mean - x, FUN = "+")

    sd_mat <- matrix(rep(g$sd, length(errordist$b)), ncol = length(g$sd), byrow = TRUE)

    lp1 <- stats::pnorm(q = left_means, sd = sd_mat, log.p = TRUE)
    lp2 <- stats::pnorm(q = right_means, sd = sd_mat, log.p = TRUE)

    max_log <- apply(cbind(apply(lp1, 2, max), apply(lp2, 2, max)), 1, max)

    lp1_sub <- t(lp1) - max_log
    lp2_sub <- t(lp2) - max_log

    diff_vec <- 1 / (errordist$b - errordist$a)

    like_vals <- log(colSums((diff_vec * errordist$pi) * t(exp(lp1_sub) - exp(lp2_sub)))) +
        max_log

    ## deal with place where point mass at zero
    which_zero <- g$sd == 0
    assertthat::assert_that(sum(which_zero) <= 1)
    like_vals[which_zero] <- log(sum(errordist$pi * (x >= errordist$a & x <= errordist$b) /
                                     (errordist$b - errordist$a)))

    return(like_vals)
}


#' Uniform-mixture - normal-mixture convolution.
#'
#' @inheritParams log_compdens_conv_mix
#'
#' @author David Gerard
log_compdens_conv_mix.uninormal <- function(g, betahat, errordist) {
    class_e <- unique(sapply(errordist, class))
    class_g <- class(g)

    assertthat::are_equal(length(class_e), 1)
    assertthat::are_equal(length(class_g), 1)
    assertthat::are_equal(class_e, "normalmix")
    assertthat::are_equal(class_g, "unimix")
    assertthat::are_equal(length(betahat), length(errordist))

    n <- length(betahat)
    p <- length(g$pi)
    matrix_llik <- matrix(NA, nrow = p, ncol = n)
    for (index in 1:n) {
        matrix_llik[, index] <- unconv_dense(g = g,
                                             errordist = errordist[[index]],
                                             x = betahat[index])
    }
    return(matrix_llik)
}


#' Uniform-normal convolution log-density evaluated at a single point.
#'
#' @param g A \code{unimix} object. This is the prior mixtrue.
#' @param errordist A \code{normalmix} object. This is the error mixture.
#' @param x A single numeric. This is the data.
#'
#' @author David Gerard
#'
#' @seealso \code{\link{log_compdens_conv_mix.uninormal}}
unconv_dense <- function(g, errordist, x) {
    assertthat::are_equal(class(g), "unimix")
    assertthat::are_equal(class(errordist), "normalmix")
    assertthat::are_equal(length(x), 1)

    left_means  <- outer(errordist$mean, g$b - x, FUN = "+")
    right_means <- outer(errordist$mean, g$a - x, FUN = "+")

    sd_mat <- matrix(rep(errordist$sd, length(g$b)), ncol = length(g$pi), byrow = FALSE)

    lp1 <- stats::pnorm(q = left_means, sd = sd_mat, log.p = TRUE)
    lp2 <- stats::pnorm(q = right_means, sd = sd_mat, log.p = TRUE)

    max_log <- apply(cbind(apply(lp1, 2, max), apply(lp2, 2, max)), 1, max)

    lp1_sub <- t(lp1) - max_log
    lp2_sub <- t(lp2) - max_log

    diff_vec <- 1 / (g$b - g$a)

    like_vals <- log(colSums((errordist$pi) * t(exp(lp1_sub) - exp(lp2_sub)))) +
        max_log + log(diff_vec)

    ## deal with values where prior is point mass on zero.
    which_zero <- g$a == 0 & g$b == 0
    assertthat::assert_that(sum(which_zero) <= 1)
    zeropart <- stats::dnorm(x = x, mean = errordist$mean, sd = errordist$sd, log = TRUE)
    max_zeropart <- max(zeropart)
    like_vals[which_zero] <- log(sum(errordist$pi * exp(zeropart - max_zeropart))) + max_zeropart

    return(like_vals)
}



#' Uniform-mixture - uniform-mixture convolution.
#'
#' @inheritParams log_compdens_conv_mix
#'
#' @author David Gerard
log_compdens_conv_mix.uniuni <- function(g, betahat, errordist) {
    class_e <- unique(sapply(errordist, class))
    class_g <- class(g)

    assertthat::are_equal(length(class_e), 1)
    assertthat::are_equal(length(class_g), 1)
    assertthat::are_equal(class_e, "unimix")
    assertthat::are_equal(class_g, "unimix")
    assertthat::are_equal(length(betahat), length(errordist))

    n <- length(betahat)
    p <- length(g$pi)
    matrix_llik <- matrix(NA, nrow = p, ncol = n)
    for (index in 1:n) {
        matrix_llik[, index] <- uuconv_dense(g = g,
                                             errordist = errordist[[index]],
                                             x = betahat[index])
    }
    return(matrix_llik)
}


#' Uniform-uniform convolution log-density evaluated at a single point.
#'
#' @param g A \code{unimix} object. This is the prior mixtrue.
#' @param errordist A \code{unimix} object. This is the error mixture.
#' @param x A single numeric. This is the data.
#'
#' @author David Gerard
#'
#' @seealso \code{\link{log_compdens_conv_mix.uniuni}}
uuconv_dense <- function(g, errordist, x) {
    assertthat::are_equal(class(g), "unimix")
    assertthat::are_equal(class(errordist), "unimix")
    assertthat::are_equal(length(x), 1)

    firstmat <- outer(rep(1, length = length(errordist$a)), g$b)
    secondmat <- outer(x - errordist$a, rep(1, length = length(g$b)))
    minval <- matrix(apply(cbind(c(firstmat), c(secondmat)),
                           1, min), nrow = length(errordist$a))

    firstmat <- outer(rep(1, length = length(errordist$b)), g$a)
    secondmat <- outer(x - errordist$b, rep(1, length = length(g$a)))
    maxval <- matrix(apply(cbind(c(firstmat), c(secondmat)),
                           1, max), nrow = length(errordist$a))

    numerator <- outer(errordist$b - errordist$a, g$b - g$a, FUN = "*")

    like_dens_mat <- (minval - maxval) / numerator
    like_dens_mat[like_dens_mat < 0] <- 0

    like_vals <- log(colSums(errordist$pi * like_dens_mat))

    ## deal with place where point mass at zero
    which_zero <- g$a == 0 & g$b == 0
    assertthat::assert_that(sum(which_zero) <= 1)
    like_vals[which_zero] <- log(sum(errordist$pi * (x >= errordist$a & x <= errordist$b) /
                                     (errordist$b - errordist$a)))
    return(like_vals)
}
