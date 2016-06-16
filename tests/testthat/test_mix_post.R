library(ashr)
context("Postcalc for Mixlike")

test_that("post_mix_dist.normalnormal", {
    gpi <- c(0.5, 0.3, 0.2)
    gmean <- c(0, 0, 1)
    gsd <- c(0, 1, 2)
    g <- ashr::normalmix(pi = gpi, mean = gmean, sd = gsd)

    epi <- c(0.3, 0.7)
    emean <- c(0, 0)
    esd <- c(2, 3)
    n <- 11
    errordist <- list()
    for (index in 1:n) {
        errordist[[index]] <- ashr::normalmix(pi = epi, mean = emean, sd = esd)
    }

    betahat <- rnorm(n)

    postout <- post_mix_dist.normalnormal(g = g, betahat = betahat, errordist = errordist)
    expect_true(all(abs(apply(postout$weights, 1, sum) - 1) < 10 ^ -14))
    meanout      <- mix_mean_array(postout)
    probzero_out <- mix_probzero_array(postout)
    sdout        <- mix_sd_array(postout)
    pless        <- mix_cdf_array(postout, 0)
    expect_true(all(pless >= 0 & pless <= 1))
}
)

test_that("post_mix_dist.normaluni", {
    gpi <- c(0.5, 0.3, 0.2)
    gmean <- c(0, 0, 1)
    gsd <- c(0, 1, 2)
    g <- ashr::normalmix(pi = gpi, mean = gmean, sd = gsd)

    epi    <- c(0.3, 0.7)
    elower <- c(-2, -1)
    eupper <- c(2, 1)
    n <- 11
    errordist <- list()
    for (index in 1:n) {
        errordist[[index]] <- ashr::unimix(pi = epi, a = elower, b = eupper)
    }

    betahat <- rnorm(n)

    postout <- post_mix_dist.normaluni(g = g, betahat = betahat, errordist = errordist)

    expect_true(all(abs(apply(postout$weights, 1, sum) - 1) < 10 ^ -14))
    meanout      <- mix_mean_array(postout)
    probzero_out <- mix_probzero_array(postout)
    sdout        <- mix_sd_array(postout)
    pless        <- mix_cdf_array(postout, 0)
    expect_true(all(pless >= 0 & pless <= 1))
}
)


test_that("post_mix_dist.uninormal", {
    gpi <- c(0.5, 0.3, 0.2)
    ga  <- c(-5, -2, 0)
    gb <- c(7, 3, 0)
    g <- ashr::unimix(pi = gpi, a = ga, b = gb)

    epi <- c(0.3, 0.7)
    emean <- c(0, 0)
    esd <- c(2, 3)
    n <- 11
    errordist <- list()
    for (index in 1:n) {
        errordist[[index]] <- ashr::normalmix(pi = epi, mean = emean, sd = esd)
    }

    betahat <- rnorm(n)

    postout <- post_mix_dist.uninormal(g = g, betahat = betahat, errordist = errordist)

    expect_true(all(abs(apply(postout$weights, 1, sum) - 1) < 10 ^ -14))

    meanout      <- mix_mean_array(postout)
    probzero_out <- mix_probzero_array(postout)
    sdout        <- mix_sd_array(postout)
    pless        <- mix_cdf_array(postout, 0)
    expect_true(all(pless >= 0 & pless <= 1))
}
)

test_that("post_mix_dist.uniuni", {
    gpi <- c(0.5, 0.3, 0.2)
    ga  <- c(-5, -2, 0)
    gb <- c(7, 3, 0)
    g <- ashr::unimix(pi = gpi, a = ga, b = gb)


    epi    <- c(0.3, 0.7)
    elower <- c(-2, -1)
    eupper <- c(2, 1)
    n <- 11
    errordist <- list()
    for (index in 1:n) {
        errordist[[index]] <- ashr::unimix(pi = epi, a = elower, b = eupper)
    }

    betahat <- rnorm(n)

    postout <- post_mix_dist.uniuni(g = g, betahat = betahat, errordist = errordist)

    expect_true(all(abs(apply(postout$weights, 1, sum) - 1) < 10 ^ -14))
    meanout      <- mix_mean_array(postout)
    probzero_out <- mix_probzero_array(postout)
    sdout        <- mix_sd_array(postout)
    pless        <- mix_cdf_array(postout, 0)
    expect_true(all(pless >= 0 & pless <= 1))
}
)
