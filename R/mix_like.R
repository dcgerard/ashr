#' Find log-component density convolution when error distribution is a
#' mixture.
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




#' Generate the posterior distribution when the error distribution is
#' a mixture.
#'
#' @param g The prior distribution. Either of class \code{"normalmix"}
#'     or \code{"unimix"}.
#' @param betahat The observations, a vector of numerics.
#' @param errordist A list of mixture distributions.
#'
#' @return A list of arrays containing the parameters for the posterior.
#' The first dimension is the observations, the second is the mixture
#' components for the prior, and the third is the mixture components of
#' the error.
#'
#'     \code{weights} The mixing proportions
#'
#'     \code{means} The location parameter if the posterior is either
#'     a mixture of normals or truncated normals.
#'
#'     \code{variances} The scale parameter if the posterior is either
#'     a mixture of normals or truncated normals.
#'
#'     \code{lower} The lower bound of support if the posterior is
#'     either a mixture of truncated normals of uniforms.
#'
#'     \code{upper} The upper bound of support if the posterior is
#'     either a mixture of truncated normals of uniforms.
#'
#' @author David Gerard
post_mix_dist <- function(g, betahat, errordist) {
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
        post_dist <- post_mix_dist.normalnormal(g, betahat, errordist)
    } else if (class_g == "normalmix" & class_e == "unimix") {
        post_dist <- post_mix_dist.normaluni(g, betahat, errordist)
    } else if (class_g == "unimix" & class_e == "normalmix") {
        post_dist <- post_mix_dist.uninormal(g, betahat, errordist)
    } else if (class_g == "unimix" & class_e == "unimix") {
        post_dist <- post_mix_dist.uniuni(g, betahat, errordist)
    } else {
        stop("Error: either g or errordist is of an unsupported class.")
    }
    return(post_dist)
}

#' Normal-mixture - normal-mixture posterior calculation.
#'
#'
#' @inheritParams post_mix_dist
#'
#' @author David Gerard
post_mix_dist.normalnormal <- function(g, betahat, errordist) {
    n  <- length(betahat)
    K  <- length(g$pi)
    Lj <- length(errordist[[1]]$pi)

    assertthat::are_equal(n, length(errordist))

    carray   <- array(NA, dim = c(n, K, Lj))
    weights_array  <- array(NA, dim = c(n, K, Lj))
    means_array    <- array(NA, dim = c(n, K, Lj))
    variances_array <- array(NA, dim = c(n, K, Lj))

    ## compute carray
    for (seindex in 1:n) {
        cerr  <- errordist[[seindex]]
        cmean <- outer(g$mean, cerr$mean, FUN = "+")
        csd   <- sqrt(outer(g$sd ^ 2, cerr$sd ^2, FUN = "+"))
        carray[seindex, , ] <- stats::dnorm(x = betahat[seindex], mean = cmean, sd = csd)
    }


    which_pointmass <- g$sd == 0

    ## compute weights_array, means_array, and variances_array
    for (seindex in 1:n) {
        cerr  <- errordist[[seindex]]
        uncalibrated_weights_array <- outer(g$pi, cerr$pi, FUN = "*") * carray[seindex, , ]
        weights_array[seindex, , ] <- uncalibrated_weights_array / sum(uncalibrated_weights_array)

        current_var  <- 1 / outer(1 / g$sd ^ 2, 1 / cerr$sd ^ 2, FUN = "+")
        current_mean <- outer(g$mean / g$sd ^ 2, (betahat[seindex] - cerr$mean) /
                                                 cerr$sd ^ 2, FUN = "+") * current_var
        current_mean[which_pointmass, ] <- g$mean[which_pointmass]

        means_array[seindex, , ]     <- current_mean
        variances_array[seindex, , ] <- current_var
    }

    post_list <- list(weights = weights_array,
                      means = means_array,
                      variances = variances_array)
    class(post_list) <- "normalmix_array"
    return(post_list)
}

#' Normal-mixture - uniform-mixture posterior calculation.
#'
#'
#' @inheritParams post_mix_dist
#'
#' @author David Gerard
post_mix_dist.normaluni <- function(g, betahat, errordist) {
    n  <- length(betahat)
    K  <- length(g$pi)
    Lj <- length(errordist[[1]]$pi)

    assertthat::are_equal(class(g), "normalmix")
    assertthat::are_equal(unique(sapply(errordist, class)), "unimix")
    assertthat::are_equal(n, length(errordist))

    carray          <- array(NA, dim = c(n, K, Lj))
    weights_array   <- array(NA, dim = c(n, K, Lj))
    means_array     <- array(NA, dim = c(n, K, Lj))
    variances_array <- array(NA, dim = c(n, K, Lj))
    lower_array     <- array(NA, dim = c(n, K, Lj))
    upper_array     <- array(NA, dim = c(n, K, Lj))

    which_pointmass <- g$sd == 0

    ## compute carray
    for (seindex in 1:n) {
        cerr <- errordist[[seindex]]
        left_vals <- outer(g$mean, cerr$b - betahat[seindex], FUN = "+") / g$sd
        right_vals <- outer(g$mean, cerr$a - betahat[seindex], FUN = "+") / g$sd

        current_c <- t(t(stats::pnorm(q = left_vals) -
                         stats::pnorm(q = right_vals)) / (cerr$b - cerr$a))

        current_c[which_pointmass, ] <-
            (betahat[seindex] >= cerr$a & betahat[seindex] <= cerr$b) / (cerr$b - cerr$a)

        carray[seindex, , ] <- current_c
    }

    ## means, variances, weights, lower and upper
    for (seindex in 1:n) {
        cerr  <- errordist[[seindex]]
        uncalibrated_weights_array <- outer(g$pi, cerr$pi, FUN = "*") * carray[seindex, , ]
        weights_array[seindex, , ] <- uncalibrated_weights_array / sum(uncalibrated_weights_array)

        current_mean <- matrix(rep(g$mean, Lj), nrow = K)
        current_var  <- matrix(rep(g$sd ^ 2, Lj), nrow = K)
        current_lower <- matrix(rep(betahat[seindex] - cerr$b, K), nrow = K, byrow = TRUE)
        current_upper <- matrix(rep(betahat[seindex] - cerr$a, K), nrow = K, byrow = TRUE)

        means_array[seindex, , ]     <- current_mean
        variances_array[seindex, , ] <- current_var
        lower_array[seindex, , ]     <- current_lower
        upper_array[seindex, , ]     <- current_upper
    }

    post_list <- list(weights = weights_array,
                      means = means_array,
                      variances = variances_array,
                      lower = lower_array,
                      upper = upper_array)
    class(post_list) <- "truncnormalmix_array"
    return(post_list)
}


#' Uniform-mixture - normal-mixture posterior calculation.
#'
#'
#' @inheritParams post_mix_dist
#'
#' @author David Gerard
post_mix_dist.uninormal <- function(g, betahat, errordist) {
    n  <- length(betahat)
    K  <- length(g$pi)
    Lj <- length(errordist[[1]]$pi)

    assertthat::are_equal(class(g), "unimix")
    assertthat::are_equal(unique(sapply(errordist, class)), "normalmix")
    assertthat::are_equal(n, length(errordist))

    carray          <- array(NA, dim = c(n, K, Lj))
    weights_array   <- array(NA, dim = c(n, K, Lj))
    means_array     <- array(NA, dim = c(n, K, Lj))
    variances_array <- array(NA, dim = c(n, K, Lj))
    lower_array     <- array(NA, dim = c(n, K, Lj))
    upper_array     <- array(NA, dim = c(n, K, Lj))

    which_pointmass <- g$a == 0 & g$b == 0

    ## compute carray
    for (seindex in 1:n) {
        cerr <- errordist[[seindex]]
        left_vals <- t(outer(cerr$mean - betahat[seindex], g$b, FUN = "+") / cerr$sd)
        right_vals <- t(outer(cerr$mean - betahat[seindex], g$a, FUN = "+") / cerr$sd)

        current_c <- (stats::pnorm(left_vals) - stats::pnorm(right_vals)) / (g$b - g$a)

        current_c[which_pointmass, ] <- stats::dnorm(x = betahat[seindex],
                                                     mean = cerr$mean,
                                                     sd = cerr$sd)
        carray[seindex, , ] <- current_c
    }

    ## means, variances, weights, lower and upper
    for (seindex in 1:n) {
        cerr  <- errordist[[seindex]]
        uncalibrated_weights_array <- outer(g$pi, cerr$pi, FUN = "*") * carray[seindex, , ]
        weights_array[seindex, , ] <- uncalibrated_weights_array / sum(uncalibrated_weights_array)

        current_mean  <- matrix(rep(betahat[seindex] - cerr$mean, K), nrow = K, byrow = TRUE)
        current_var   <- matrix(rep(cerr$sd ^ 2, K), nrow = K, byrow = TRUE)
        current_lower <- matrix(rep(g$a, Lj), nrow = K, byrow = FALSE)
        current_upper <- matrix(rep(g$b, Lj), nrow = K, byrow = FALSE)

        means_array[seindex, , ]     <- current_mean
        variances_array[seindex, , ] <- current_var
        lower_array[seindex, , ]     <- current_lower
        upper_array[seindex, , ]     <- current_upper
    }

    post_list <- list(weights = weights_array,
                      means = means_array,
                      variances = variances_array,
                      lower = lower_array,
                      upper = upper_array)
    class(post_list) <- "truncnormalmix_array"
    return(post_list)
}


#' Uniform-mixture - uniform-mixture posterior calculation.
#'
#'
#' @inheritParams post_mix_dist
#'
#' @author David Gerard
post_mix_dist.uniuni <- function(g, betahat, errordist) {
    n  <- length(betahat)
    K  <- length(g$pi)
    Lj <- length(errordist[[1]]$pi)

    assertthat::are_equal(class(g), "unimix")
    assertthat::are_equal(unique(sapply(errordist, class)), "unimix")
    assertthat::are_equal(n, length(errordist))

    carray          <- array(NA, dim = c(n, K, Lj))
    weights_array   <- array(NA, dim = c(n, K, Lj))
    lower_array     <- array(NA, dim = c(n, K, Lj))
    upper_array     <- array(NA, dim = c(n, K, Lj))

    which_pointmass <- g$a == 0 & g$b == 0

    ## compute carray
    for (seindex in 1:n) {
        cerr <- errordist[[seindex]]

        topleft <- matrix(apply(cbind(c(matrix(rep(cerr$b, K), nrow = K, byrow = TRUE)),
                               c(matrix(rep(betahat[seindex] - g$a, Lj), nrow = K))),
                               1, min), nrow = K)

        topright <- matrix(apply(cbind(c(matrix(rep(cerr$a, K), nrow = K, byrow = TRUE)),
                               c(matrix(rep(betahat[seindex] - g$b, Lj), nrow = K))),
                               1, max), nrow = K)

        denom <- outer(g$b - g$a, cerr$b - cerr$a, FUN = "*")
        current_c <- (topleft - topright) / denom

        current_c[which_pointmass, ] <-
            (betahat[seindex] >= cerr$a & betahat[seindex] <= cerr$b) / (cerr$b - cerr$a)

        current_c[current_c < 0] <- 0

        carray[seindex, , ] <- current_c
    }


    ## weights, lower, and upper
    for (seindex in 1:n) {
        cerr  <- errordist[[seindex]]
        uncalibrated_weights_array <- outer(g$pi, cerr$pi, FUN = "*") * carray[seindex, , ]
        weights_array[seindex, , ] <- uncalibrated_weights_array / sum(uncalibrated_weights_array)

        lower_current <- matrix(apply(cbind(c(matrix(rep(betahat[seindex] - cerr$b, K),
                                                     nrow = K, byrow = TRUE)),
                                            c(matrix(rep(g$a, Lj), nrow = K, byrow = FALSE))),
                                      1, max), nrow = K)

        upper_current <- matrix(apply(cbind(c(matrix(rep(betahat[seindex] - cerr$a, K),
                                                     nrow = K, byrow = TRUE)),
                                            c(matrix(rep(g$b, Lj), nrow = K, byrow = FALSE))),
                                      1, min), nrow = K)

        lower_array[seindex, , ] <- lower_current
        upper_array[seindex, , ] <- upper_current


        lower_array[seindex, which_pointmass, ] <- 0
        upper_array[seindex, which_pointmass, ] <- 0
    }
    assertthat::assert_that(all(upper_array >= lower_array))

    post_list <- list(weights = weights_array,
                      lower = lower_array,
                      upper = upper_array)

    class(post_list) <- "unimix_array"
    return(post_list)
}


#' Compute mean from a mix_array.
#'
#' @param mixdist Of class \code{unimix_array},
#'     \code{normalmix_array}, or \code{truncnormalmix_array}. These
#'     are the classes of the output from \code{\link{post_mix_dist}}.
#'
#' @author David Gerard
mix_mean_array <- function(mixdist) {
      UseMethod("mix_mean_array")
}
mix_mean_array.default <- function(mixdist) {
    stop(paste("Error: Invalid class", class(mixdist), "for argument in", match.call()))
}
mix_mean_array.normalmix_array <- function(mixdist) {
    return(apply(mixdist$means * mixdist$weights, 1, sum))
}
mix_mean_array.truncnormalmix_array <- function(mixdist) {
    which_pointmass <- mixdist$variances == 0 | (mixdist$lower == 0 & mixdist$upper == 0)

    sdarray <- sqrt(mixdist$variances)
    alpha <- (mixdist$lower - mixdist$means) / sdarray
    beta  <- (mixdist$upper - mixdist$means) / sdarray

    ## numerically unstable way
    ## Z <- stats::pnorm(beta) - stats::pnorm(alpha)
    ## num <- stats::dnorm(alpha) - stats::dnorm(beta)
    ## postmeans <- mixdist$mean + num / Z * sdarray
    ## postmeans[which_pointmass] <- 0
    ## outvec1 <- apply(postmeans * mixdist$weights, 1, sum)
    ## badmeans <- postmeans == Inf

    ## numerically stable way
    which_switch <- alpha > 0
    Z                  <- array(NA, dim = dim(alpha))
    num                <- array(NA, dim = dim(alpha))
    Z[!which_switch]   <- stats::pnorm(beta[!which_switch]) - stats::pnorm(alpha[!which_switch])
    num[!which_switch] <- stats::dnorm(alpha[!which_switch]) - stats::dnorm(beta[!which_switch])
    Z[which_switch]    <- stats::pnorm(-alpha[which_switch]) - stats::pnorm(-beta[which_switch])
    num[which_switch]  <- stats::dnorm(-alpha[which_switch]) - stats::dnorm(-beta[which_switch])
    postmeans <- mixdist$mean + num / Z * sdarray
    postmeans[which_pointmass] <- 0
    outvec <- apply(postmeans * mixdist$weights, 1, sum)

    ## truncnorm method
    ## postmeans <- truncnorm::etruncnorm(a = mixdist$lower, b = mixdist$upper,
    ##                                    mean = mixdist$means, sd = sqrt(mixdist$variances))
    ## postmeans <- array(postmeans, dim = dim(mixdist$means))

    return(outvec)
}
mix_mean_array.unimix_array <- function(mixdist) {
    return(apply((mixdist$upper + mixdist$lower) / 2 * mixdist$weights, 1, sum))
}

#' Compute the probability of zero.
#'
#' @inheritParams mix_mean_array
#'
#' @author David Gerard
mix_probzero_array <- function(mixdist) {
    UseMethod("mix_probzero_array")
}
mix_probzero_array.default <- function(mixdist) {
    stop(paste("Error: Invalid class", class(mixdist), "for argument in", match.call()))
}
mix_probzero_array.normalmix_array <- function(mixdist) {
    which_pointmass <- mixdist$variances == 0
    return(apply(mixdist$weights * which_pointmass, 1, sum))
}
mix_probzero_array.truncnormalmix_array <- function(mixdist) {
    which_pointmass <- mixdist$variances == 0 | (mixdist$lower == 0 & mixdist$upper == 0)
    return(apply(mixdist$weights * which_pointmass, 1, sum))
}
mix_probzero_array.unimix_array <- function(mixdist) {
    which_pointmass <- mixdist$lower == 0 & mixdist$upper == 0
    return(apply(mixdist$weights * which_pointmass, 1, sum))
}



#' Compute SD from a mix_array.
#'
#' @param mixdist Of class \code{unimix_array},
#'     \code{normalmix_array}, or \code{truncnormalmix_array}. These
#'     are the classes of the output from \code{\link{post_mix_dist}}.
#'
#'
#' @return The posterior standard deviations of the location
#'     parameters.
#'
#' @author David Gerard
mix_sd_array <- function(mixdist) {
      UseMethod("mix_sd_array")
}
mix_sd_array.normalmix_array <- function(mixdist) {
    second_moment <- apply(mixdist$weights * (mixdist$means ^ 2 + mixdist$variances), 1, sum)
    first_moment2 <- apply(mixdist$weights * mixdist$means, 1, sum) ^ 2
    return(sqrt(second_moment - first_moment2))
}
mix_sd_array.truncnormalmix_array <- function(mixdist) {
    which_pointmass <- mixdist$variances == 0 | (mixdist$lower == 0 & mixdist$upper == 0)

    sdarray <- sqrt(mixdist$variances)
    alpha <- (mixdist$lower - mixdist$means) / sdarray
    beta  <- (mixdist$upper - mixdist$means) / sdarray
    which_switch <- alpha > 0
    Z                  <- array(NA, dim = dim(alpha))
    num                <- array(NA, dim = dim(alpha))
    Z[!which_switch]   <- stats::pnorm(beta[!which_switch]) - stats::pnorm(alpha[!which_switch])
    num[!which_switch] <- stats::dnorm(alpha[!which_switch]) - stats::dnorm(beta[!which_switch])
    Z[which_switch]    <- stats::pnorm(-alpha[which_switch]) - stats::pnorm(-beta[which_switch])
    num[which_switch]  <- stats::dnorm(-alpha[which_switch]) - stats::dnorm(-beta[which_switch])
    postmeans <- mixdist$mean + num / Z * sdarray
    postmeans[which_pointmass] <- 0

    ## truncnorm method
    ## postmeans <- truncnorm::etruncnorm(a = mixdist$lower, b = mixdist$upper,
    ##                                    mean = mixdist$means, sd = sqrt(mixdist$variances))
    ## postmeans <- array(postmeans, dim = dim(mixdist$means))

    postvars <- (1 + (alpha * stats::dnorm(alpha) - beta * stats::dnorm(beta)) / Z -
                 ((stats::dnorm(alpha) - stats::dnorm(beta)) / Z) ^ 2) * mixdist$variances
    postvars[which_pointmass] <- 0

    ## truncnorm method
    ## postvars <- truncnorm::vtruncnorm(a = mixdist$lower, b = mixdist$upper,
    ##                                   mean = mixdist$means, sd = sqrt(mixdist$variances))
    ## postvars <- array(postvars, dim = dim(mixdist$means))

    second_moment <- apply(mixdist$weights * (postmeans ^ 2 + postvars), 1, sum)
    first_moment2 <- apply(mixdist$weights * postmeans, 1, sum) ^ 2
    return(sqrt(second_moment - first_moment2))
}
mix_sd_array.unimix_array <- function(mixdist) {
    which_pointmass <- mixdist$lower == 0 & mixdist$upper == 0
    postmeans <- (mixdist$upper + mixdist$lower) / 2
    postvars  <- (mixdist$upper - mixdist$lower) ^ 2 / 12
    postmeans[which_pointmass] <- 0
    postvars[which_pointmass] <- 0
    second_moment <- apply(mixdist$weights * (postmeans ^ 2 + postvars), 1, sum)
    first_moment2 <- apply(mixdist$weights * postmeans, 1, sum) ^ 2
    return(sqrt(second_moment - first_moment2))
}


#' Compute cdf from a mix_array.
#'
#' @param mixdist Of class \code{unimix_array},
#'     \code{normalmix_array}, or \code{truncnormalmix_array}. These
#'     are the classes of the output from \code{\link{post_mix_dist}}.
#' @param q The quantile.
#'
#' @return The probability of being less than or equal to q.
#'
#' @author David Gerard
mix_cdf_array <- function(mixdist, q) {
      UseMethod("mix_cdf_array")
}
mix_cdf_array.normalmix_array <- function(mixdist, q) {
    pless <- stats::pnorm(q = q, mean = mixdist$means, sd = sqrt(mixdist$variances))
    return(apply(pless * mixdist$weights, 1, sum))
}
mix_cdf_array.truncnormalmix_array <- function(mixdist, q) {
    which_pointmass <- mixdist$variances == 0 | (mixdist$lower == 0 & mixdist$upper == 0)

    sdarray <- sqrt(mixdist$variances)
    alpha <- (mixdist$lower - mixdist$means) / sdarray
    beta  <- (mixdist$upper - mixdist$means) / sdarray
    xi    <- (q - mixdist$mean) / sdarray
    which_switch <- alpha > 0
    Z                  <- array(NA, dim = dim(alpha))
    num                <- array(NA, dim = dim(alpha))
    Z[!which_switch]   <- stats::pnorm(beta[!which_switch]) - stats::pnorm(alpha[!which_switch])
    num[!which_switch] <- stats::pnorm(xi[!which_switch]) - stats::pnorm(alpha[!which_switch])
    Z[which_switch]    <- stats::pnorm(-alpha[which_switch]) - stats::pnorm(-beta[which_switch])
    num[which_switch]  <- stats::pnorm(-alpha[which_switch]) - stats::pnorm(-xi[which_switch])
    pless <- num / Z
    pless[which_pointmass] <- (q >= 0) * 1

    ## truncnorm package way
    ## pless <- truncnorm::ptruncnorm(q = q, a = mixdist$lower, b = mixdist$upper,
    ##                                mean = mixdist$means, sd = sqrt(mixdist$variances))
    ## pless <- array(pless, dim = dim(mixdist$means))
    ## pless[which_pointmass] <- (q >= 0) * 1

    return(apply(pless * mixdist$weights, 1, sum))
}
mix_cdf_array.unimix_array <- function(mixdist, q) {
    which_pointmass <- mixdist$lower == 0 & mixdist$upper == 0
    pless <- stats::punif(q = q, min = mixdist$lower, max = mixdist$upper)
    pless[which_pointmass] <- (q >= 0) * 1
    return(apply(pless * mixdist$weights, 1, sum))
}


#' Calculate log-likelihood.
#'
#' @author David Gerard
#'
#' @inheritParams log_compdens_conv_mix
calc_loglik_array <- function(g, betahat, errordist) {
    matrix_llik <- log_compdens_conv_mix(g, betahat, errordist)
    maxmat <- apply(matrix_llik, 2, max)
    matrix_llik_norm <- t(t(matrix_llik) - maxmat)
    llike <- sum(log(colSums(exp(matrix_llik_norm) * g$pi)) + maxmat)
    return(llike)
}

#' Calculate the log-likelihood when all null.
#'
#' @author David Gerard
#'
#' @inheritParams log_compdens_conv_mix
calc_nulllik_array <- function(betahat, errordist) {
    g <- ashr::normalmix(pi = 1, mean = 0, sd = 0)
    matrix_llik <- log_compdens_conv_mix(g, betahat, errordist)
    return(sum(matrix_llik))
}



#' Generate an object of class \code{normalmix} that is an
#' approximation to a t-distribution.
#'
#' This is based on the representation of the t as an inverse-gamma
#' scaled-mixture of Gaussians.
#'
#' @param mu The mean of the t-distribution.
#' @param sig The scale parameter of the t-distribution.
#' @param df The degrees of freedom of the t-distribution.
#' @param gridsize The number of mixture components to use. The larger
#'     the more accurate the approximation, but the higher the
#'     computational load --- especially if you intend to use this for
#'     \code{\link{ash.workhorse}}.
#'
#' @export
#'
#' @author David Gerard
t_to_mix <- function(mu, sig, df, gridsize = 20) {
    pgrid <- seq(1 / 10000, 1 - 1 / 10000, length = gridsize + 1)

    shape_param <- df / 2
    rate_param  <- df * sig ^ 2 / 2

    mean_grid      <- rep(mu, length = gridsize)
    temp_grid      <- stats::qgamma(p = pgrid, shape = shape_param, rate = rate_param)
    precision_grid <- (temp_grid[2:length(temp_grid)] +
                       temp_grid[1:(length(temp_grid) - 1)]) / 2

    weight_grid <- rep(1 / gridsize, length = gridsize)

    tapprox <- normalmix(pi = weight_grid, mean = mean_grid, sd = 1 / sqrt(precision_grid))

    return(tapprox)
}



#' Generate an object of class \code{normalmix} that is an
#' approximation to a Laplace distribution.
#'
#' This is based on the representation of the Laplace distribution as
#' an exponential scale mixtures of normals.
#'
#' @param mu The location of the Laplace distribution.
#' @param sig The scale parameter of the Laplace distribution.
#' @param gridsize The number of mixture components to use. The larger
#'     the more accurate the approximation, but the higher the
#'     computational load --- especially if you intend to use this for
#'     \code{\link{ash.workhorse}}.
#'
#' @export
#'
#' @author David Gerard
laplace_to_mix <- function(mu, sig, gridsize = 20) {
    if (!requireNamespace("VGAM", quietly = TRUE)) {
        stop("VGAM must be installed to use laplace_to_mix")
    }
    pgrid <- seq(1 / 10000, 1 - 1 / 10000, length = gridsize + 1)

    mean_grid <- rep(mu, length = gridsize)
    temp_grid <- VGAM::qrayleigh(p = pgrid, scale = sig ^ 2)
    var_grid  <- (temp_grid[2:length(temp_grid)] +
                  temp_grid[1:(length(temp_grid) - 1)]) / 2

    weight_grid <- rep(1 / gridsize, length = gridsize)

    lapprox <- normalmix(pi = weight_grid, mean = mean_grid, sd = sqrt(var_grid))

    return(lapprox)
}

#' Random draw from a mixture of normals.
#'
#' @param n the number of samples to draw.
#' @param mixdense An object of class \code{normalmix}.
#'
#' @author David Gerard
#'
#' @export
rnormalmix <- function(n, mixdense) {
    assertthat::are_equal(class(mixdense), "normalmix")
    k <- length(mixdense$pi)
    which_norm <- sample(1:k, size = n, prob = mixdense$pi, replace = TRUE)
    samp <- stats::rnorm(n = n, mean = mixdense$mean[which_norm], sd = mixdense$sd[which_norm])
    return(samp)
}

#' CDF function for a mixture of normals.
#'
#' @param q The quantile.
#' @param mixdense An object of class \code{normalmix}.
#'
#' @author David Gerard
#'
#' @export
pnormalmix <- function(q, mixdense) {
    assertthat::are_equal(class(mixdense), "normalmix")
    pout <- sum(mixdense$pi * stats::pnorm(q = q, mean = mixdense$mean, sd = mixdense$sd))
    return(pout)
}

#' Density function for a mixture of normals.
#'
#' @param x The location to calculate the density.
#' @param mixdense An object of class \code{normalmix}.
#'
#' @author David Gerard
#'
#' @export
dnormalmix <- function(x, mixdense) {
    assertthat::are_equal(class(mixdense), "normalmix")
    dout <- sum(mixdense$pi * stats::dnorm(x = x, mean = mixdense$mean, sd = mixdense$sd))
    return(dout)
}

#' Quantile function for a mixture of normals.
#'
#' @param p The probability at which to calculate the quantile.
#' @param mixdense An object of class \code{normalmix}.
#'
#' @author David Gerard
#'
#' @export
qnormalmix <- function(p, mixdense) {
    psub <- function(q, p, mixdense) {
        return(pnormalmix(q = q, mixdense = mixdense) - p)
    }
    rootout <- stats::uniroot(f = psub, interval = c(-100, 100), p = p, mixdense = mixdense)
    return(rootout)
}

#' Random draw from a mixture of uniforms.
#'
#' @param n the number of samples to draw.
#' @param mixdense An object of class \code{unimix}.
#'
#' @author David Gerard
#'
#' @export
runimix <- function(n, mixdense) {
    assertthat::are_equal(class(mixdense), "unimix")
    k <- length(mixdense$pi)
    which_uni <- sample(1:k, size = n, prob = mixdense$pi, replace = TRUE)
    samp <- stats::runif(n = n, min = mixdense$a[which_uni], max = mixdense$b[which_uni])
    return(samp)
}
