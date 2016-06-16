library(ashr)
context("Postcalc for Mixlike")

test_that("post_mix_dist.normalnormal", {
    gpi <- c(0.5, 0.3, 0.2)
    gmean <- c(0, 0, 1)
    gsd <- c(0, 1, 2)
    g <- normalmix(pi = gpi, mean = gmean, sd = gsd)

    epi <- c(0.3, 0.7)
    emean <- c(0, 0)
    esd <- c(2, 3)
    n <- 11
    errordist <- list()
    for (index in 1:n) {
        errordist[[index]] <- normalmix(pi = epi, mean = emean, sd = esd)
    }

    betahat <- rnorm(n)

    postout <- post_mix_dist.normalnormal(g = g, betahat = betahat, errordist = errordist)
    expect_true(all(abs(apply(postout$weights, 1, sum) - 1) < 10 ^ -14))
    meanout      <- mix_mean_array(postout)
    probzero_out <- mix_probzero_array(postout)
    expect_true(all(probzero_out >= 0 & probzero_out <= 1))
    sdout        <- mix_sd_array(postout)
    pless        <- mix_cdf_array(postout, 0)
    expect_true(all(pless >= 0 & pless <= 1))
}
)

test_that("post_mix_dist.normaluni", {
    set.seed(76)
    gpi <- c(0.5, 0.3, 0.2)
    gmean <- c(0, 0, 1)
    gsd <- c(0, 1, 2)
    g <- normalmix(pi = gpi, mean = gmean, sd = gsd)

    epi    <- c(0.3, 0.7)
    elower <- c(-2, -1)
    eupper <- c(2, 1)
    n <- 11
    errordist <- list()
    for (index in 1:n) {
        errordist[[index]] <- unimix(pi = epi, a = elower, b = eupper)
    }

    betahat <- rnorm(n)

    postout <- post_mix_dist.normaluni(g = g, betahat = betahat, errordist = errordist)

    expect_true(all(abs(apply(postout$weights, 1, sum) - 1) < 10 ^ -14))
    meanout      <- mix_mean_array(postout)

    ## truncnorm way to calculate means
    tpostmeans <- truncnorm::etruncnorm(a = postout$lower, b = postout$upper,
                                        mean = postout$means, sd = sqrt(postout$variances))
    tpostmeans <- array(tpostmeans, dim = dim(postout$means))
    tpostmeans[, 1, ] <- 0
    tmeanout <- apply(tpostmeans * postout$weights, 1, sum)
    expect_equal(meanout, tmeanout)

    probzero_out <- mix_probzero_array(postout)
    expect_true(all(probzero_out >= 0 & probzero_out <= 1))
    sdout        <- mix_sd_array(postout)

    ## truncnorm method
    tpostvars <- truncnorm::vtruncnorm(a = postout$lower, b = postout$upper,
                                       mean = postout$means, sd = sqrt(postout$variances))
    tpostvars <- array(tpostvars, dim = dim(postout$means))
    tpostvars[, 1, ] <- 0
    second_moment <- apply(postout$weights * (tpostmeans ^ 2 + tpostvars), 1, sum)
    first_moment2 <- apply(postout$weights * tpostmeans, 1, sum) ^ 2
    tsdout <- sqrt(second_moment - first_moment2)
    expect_equal(sdout, tsdout)


    pless        <- mix_cdf_array(postout, 0)
    expect_true(all(pless >= 0 & pless <= 1))
}
)


test_that("post_mix_dist.uninormal", {
    set.seed(98)
    gpi <- c(0.5, 0.3, 0.2)
    ga  <- c(-5, -2, 0)
    gb <- c(7, 3, 0)
    g <- unimix(pi = gpi, a = ga, b = gb)

    epi <- c(0.3, 0.7)
    emean <- c(0, 0)
    esd <- c(2, 3)
    n <- 11
    errordist <- list()
    for (index in 1:n) {
        errordist[[index]] <- normalmix(pi = epi, mean = emean, sd = esd)
    }

    betahat <- rnorm(n)

    postout <- post_mix_dist.uninormal(g = g, betahat = betahat, errordist = errordist)

    expect_true(all(abs(apply(postout$weights, 1, sum) - 1) < 10 ^ -14))

    meanout      <- mix_mean_array(postout)

    ## truncnorm way to calculate means
    tpostmeans <- truncnorm::etruncnorm(a = postout$lower, b = postout$upper,
                                        mean = postout$means, sd = sqrt(postout$variances))
    tpostmeans <- array(tpostmeans, dim = dim(postout$means))
    tpostmeans[, 3, ] <- 0
    tmeanout <- apply(tpostmeans * postout$weights, 1, sum)
    expect_equal(meanout, tmeanout)

    probzero_out <- mix_probzero_array(postout)
    expect_true(all(probzero_out >= 0 & probzero_out <= 1))

    sdout        <- mix_sd_array(postout)

    ## truncnorm method
    tpostvars <- truncnorm::vtruncnorm(a = postout$lower, b = postout$upper,
                                       mean = postout$means, sd = sqrt(postout$variances))
    tpostvars <- array(tpostvars, dim = dim(postout$means))
    tpostvars[, 3, ] <- 0
    second_moment <- apply(postout$weights * (tpostmeans ^ 2 + tpostvars), 1, sum)
    first_moment2 <- apply(postout$weights * tpostmeans, 1, sum) ^ 2
    tsdout <- sqrt(second_moment - first_moment2)
    expect_equal(sdout, tsdout)

    pless        <- mix_cdf_array(postout, 0)
    expect_true(all(pless >= 0 & pless <= 1))
}
)

test_that("post_mix_dist.uniuni", {
    gpi <- c(0.5, 0.3, 0.2)
    ga  <- c(-5, -2, 0)
    gb <- c(7, 3, 0)
    g <- unimix(pi = gpi, a = ga, b = gb)


    epi    <- c(0.3, 0.7)
    elower <- c(-2, -1)
    eupper <- c(2, 1)
    n <- 11
    errordist <- list()
    for (index in 1:n) {
        errordist[[index]] <- unimix(pi = epi, a = elower, b = eupper)
    }

    betahat <- rnorm(n)

    postout <- post_mix_dist.uniuni(g = g, betahat = betahat, errordist = errordist)

    expect_true(all(abs(apply(postout$weights, 1, sum) - 1) < 10 ^ -14))
    meanout      <- mix_mean_array(postout)
    probzero_out <- mix_probzero_array(postout)
    expect_true(all(probzero_out >= 0 & probzero_out <= 1))
    sdout        <- mix_sd_array(postout)
    pless        <- mix_cdf_array(postout, 0)
    expect_true(all(pless >= 0 & pless <= 1))
}
)


test_that("giving normal error dist and mixnorm prior will result in same as regular ashr", {
    set.seed(199)
    p <- 11
    beta <- rnorm(p)
    sebetahat <- sqrt(rchisq(p, df = 5) / 5)
    y <- beta + rnorm(p, sd = sebetahat)

    errordist <- list()
    for (index in 1:p) {
        errordist[[index]] <- normalmix(pi = 1, mean = 0, sd = sebetahat[index])
    }

    ashout <- ash.workhorse(betahat = y, sebetahat = sebetahat, mixcompdist = "normal",
                            outputlevel = 4)

    eout   <- ash.workhorse(betahat = y, errordist = errordist, mixcompdist = "normal")

    postout <- post_mix_dist(g = eout$fitted.g,
                             betahat = y,
                             errordist = errordist)
    postmean <- mix_mean_array(mixdist = postout)

    expect_true(max(c(abs(ashout$flash.data$comp_postprob - t(postout$weights[, , 1])))) <
                10 ^ -14)
    expect_true(max(c(abs(ashout$flash.data$comp_postmean - t(postout$means[, , 1])))) <
                10 ^ -14)
    expect_true(max(c(abs(ashout$flash.data$comp_postmean2 -
                          t(postout$variances[, , 1] + postout$means[, , 1] ^ 2)))) <
                10 ^ -14)
    expect_equal(eout$fitted.g, ashout$fitted.g)
    expect_equal(eout$PosteriorMean, ashout$PosteriorMean)
    expect_equal(eout$PosteriorSD, ashout$PosteriorSD)
    expect_equal(eout$PosteriorSD, ashout$PosteriorSD)
    expect_equal(postmean, eout$PosteriorMean)
    expect_equal(eout$lfdr, ashout$lfdr)
    expect_equal(eout$lfsr, ashout$lfsr)
    expect_equal(eout$loglik, ashout$loglik)
    expect_equal(eout$logLR, ashout$logLR)
}
)

test_that("giving normal error dist and unimix prior will result in same as regular ashr", {
    set.seed(810)
    p <- 11
    beta <- rnorm(p)
    sebetahat <- sqrt(rchisq(p, df = 5) / 5)
    y <- beta + rnorm(p, sd = sebetahat)

    errordist <- list()
    for (index in 1:p) {
        errordist[[index]] <- normalmix(pi = 1, mean = 0, sd = sebetahat[index])
    }

    ashout <- ash.workhorse(betahat = y, sebetahat = sebetahat, mixcompdist = "uniform",
                            outputlevel = 4)

    eout   <- ash.workhorse(betahat = y, errordist = errordist, mixcompdist = "uniform")

    postout <- post_mix_dist(g = eout$fitted.g,
                             betahat = y,
                             errordist = errordist)
    postmean <- mix_mean_array(mixdist = postout)
    postsd   <- mix_sd_array(mixdist = postout)
    expect_true(max(c(abs(ashout$flash.data$comp_postprob - t(postout$weights[, , 1])))) <
                10 ^ -14)

    expect_equal(eout$fitted.g, ashout$fitted.g)
    expect_equal(eout$PosteriorMean, ashout$PosteriorMean)
    expect_equal(postmean, eout$PosteriorMean)
    expect_equal(eout$PosteriorSD, ashout$PosteriorSD)
    expect_equal(eout$PosteriorSD, ashout$PosteriorSD)
    expect_equal(eout$lfdr, ashout$lfdr)
    expect_equal(eout$lfsr, ashout$lfsr)
    expect_equal(eout$loglik, ashout$loglik)
    expect_equal(eout$logLR, ashout$logLR)
}
)
