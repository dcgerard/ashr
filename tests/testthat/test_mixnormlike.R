library(ashr)
context("\nMixNormLikelihood")

test_that("estimate_mixprop correctly handles mixtures of normals likelihood", {
    p <- 10
    beta <- rnorm(p)
    sebetahat <- sqrt(rchisq(p, df = 5))
    betahat <- rnorm(beta, sebetahat)

    ## assumed likelihood
    pi_vals  <- c(0.5, 0.3, 0.2)
    mean_seq <- c(0, 0, 0)
    sd_seq   <- c(0, 1 , 2)

    errordist <- list()
    for (index in 1:p) {
        errordist[[index]] <- ashr::normalmix(pi = pi_vals, mean = mean_seq, sd = sd_seq)
    }
}
)

test_that("log_compdens_conv_mix.normalnormal function works", {
    set.seed(699)
    p <- 11
    betahat <- rnorm(p)

    ## assumed likelihood
    pi_vals  <- c(0.5, 0.3, 0.2)
    mean_seq <- c(0, 0, 0)
    sd_seq   <- c(0.5, 1 , 2)
    errordist <- list()
    for (index in 1:p) {
        errordist[[index]] <- ashr::normalmix(pi = pi_vals, mean = mean_seq, sd = sd_seq)
    }
    pi_vals  <- rep(0.2, length = 5)
    mean_seq <- rep(0, length = 5)
    sd_seq   <- seq(0, 4, length = 5)
    g <- ashr::normalmix(pi = pi_vals, mean = mean_seq, sd = sd_seq)
    lout <- log_compdens_conv_mix.normalnormal(g = g, betahat = betahat, errordist = errordist)

    expect_equal(length(betahat), ncol(lout))
    expect_equal(length(g$pi), nrow(lout))

    ## test same as log_compdens_conv.normalmix
    errordist <- list()
    pi_vals <- 1
    mean_seq <- rep(0, length = p)
    sd_seq <- rchisq(p, 10) / 10
    for (index in 1:p) {
        errordist[[index]] <- ashr::normalmix(pi = pi_vals, mean = mean_seq[index],
                                              sd = sd_seq[index])
    }

    lout1 <- log_compdens_conv_mix.normalnormal(g = g, betahat = betahat,
                                                errordist = errordist)
    lout2 <-log_compdens_conv.normalmix(m = g, x = betahat, s = sd_seq, v = NULL)
    expect_equal(lout1, lout2)

    ## test point mass at zero where I know results
    normal1 <- ashr::normalmix(pi = 1, mean = 0, sd = 0)
    normal2 <- ashr::normalmix(pi = 1, mean = 0, sd = 1)
    mylike <- nnconv_dense(g = normal1, errordist = normal2, x = 0)
    knownlike <- log(dnorm(0, normal2$mean, normal2$sd))
    expect_equal(mylike, knownlike)

    ## test point mass at zero where I know results
    normal1 <- ashr::normalmix(pi = 1, mean = 0, sd = 0)
    normal2 <- ashr::normalmix(pi = c(0.5, 0.5), mean = c(0, 0), sd = c(1, 2))
    mylike <- nnconv_dense(g = normal1, errordist = normal2, x = 0)
    knownlike <- log(sum(dnorm(0, normal2$mean, normal2$sd) * normal2$pi))
    expect_equal(mylike, knownlike)
}
)


test_that("log_compdens_conv_mix.normaluni works", {
    set.seed(543)
    p <- 11
    betahat <- rnorm(p)

    ## assumed likelihood
    pi_vals  <- c(0.5, 0.3, 0.2)
    a_seq <- c(-4, -2, -1)
    b_seq <- c(4, 2, 1)
    errordist <- list()
    for (index in 1:p) {
        errordist[[index]] <- ashr::unimix(pi = pi_vals, a = a_seq, b = b_seq)
    }
    pi_vals  <- rep(0.1, length = 5)
    mean_seq <- rep(0, length = 5)
    sd_seq   <- seq(0, 4, length = 5)
    g <- ashr::normalmix(pi = pi_vals, mean = mean_seq, sd = sd_seq)
    lout <- log_compdens_conv_mix.normaluni(g = g, betahat = betahat, errordist = errordist)

    expect_equal(ncol(lout), p)
    expect_equal(nrow(lout), length(g$pi))

    ## simple example that I can check by hand
    uni_part  <- ashr::unimix(pi = 1, a = -1, b = 1)
    norm_part <- ashr::normalmix(pi = 1, mean = 0, sd = 1)
    x <- 0
    my_val <- nuconv_dense(g = norm_part, errordist = uni_part, x = x)

    like_val <- (pnorm(q = uni_part$b, mean = x - norm_part$mean, sd = norm_part$sd) -
                 pnorm(q = uni_part$a, mean = x - norm_part$mean, sd = norm_part$sd)) /
        (uni_part$b - uni_part$a)

    expect_equal(like_val, exp(my_val))

    ## Point mass at 0 example where I know results
    uni_part  <- ashr::unimix(pi = 1, a = -1, b = 1)
    norm_part <- ashr::normalmix(pi = 1, mean = 0, sd = 0)
    x <- 0
    my_val <- nuconv_dense(errordist = uni_part, g = norm_part, x = x)
    like_val <- log(1 / (uni_part$b - uni_part$a))
    expect_equal(like_val, my_val)
}
)


test_that("unconv_dense works", {
    set.seed(33)
    p <- 11
    betahat <- rnorm(p)

    ## assumed likelihood
    pi_vals  <- c(0.5, 0.3, 0.2)
    mean_seq <- c(0, 0, 0)
    sd_seq   <- c(1, 2, 3)
    errordist <- list()
    for (index in 1:p) {
        errordist[[index]] <- ashr::normalmix(pi = pi_vals, mean = mean_seq,
                                              sd = sd_seq)
    }
    pi_vals <- rep(0.1, length = 5)
    a_seq   <- seq(-5, 0, length = 5)
    b_seq   <- seq(5, 0, length = 5)
    g <- ashr::unimix(pi = pi_vals, a = a_seq, b = b_seq)
    lout <- log_compdens_conv_mix.uninormal(g = g, betahat = betahat, errordist = errordist)

    expect_equal(ncol(lout), p)
    expect_equal(nrow(lout), length(g$pi))


    ## Example where I know results
    uni_part  <- ashr::unimix(pi = 1, a = -1, b = 1)
    norm_part <- ashr::normalmix(pi = 1, mean = 0, sd = 1)
    x <- 0
    my_val <- unconv_dense(g = uni_part, errordist = norm_part, x = x)

    like_val <- (pnorm(q = uni_part$b, mean = x - norm_part$mean, sd = norm_part$sd) -
                 pnorm(q = uni_part$a, mean = x - norm_part$mean, sd = norm_part$sd)) /
        (uni_part$b - uni_part$a)

    expect_equal(like_val, exp(my_val))

    ## Point mass at 0 example where I know results
    uni_part  <- ashr::unimix(pi = 1, a = 0, b = 0)
    norm_part <- ashr::normalmix(pi = 1, mean = 0, sd = 1)
    x <- 0
    my_val <- unconv_dense(g = uni_part, errordist = norm_part, x = x)
    like_val <- dnorm(x, norm_part$mean, norm_part$sd)
    expect_equal(like_val, exp(my_val))

}
)


test_that("uuconv_dense works", {
    set.seed(889)
    k1 <- 5
    k2 <- 7
    a_seq1 <- seq(-5, 0, length = k1)
    b_seq1 <- seq(5, 0, length = k1)
    a_seq2 <- seq(-4, -1, length = k2)
    b_seq2 <- seq(3, 6, length = k2)
    pi_vals1 <- abs(rnorm(k1))
    pi_vals1 <- pi_vals1 / sum(pi_vals1)
    pi_vals2 <- abs(rnorm(k2))
    pi_vals2 <- pi_vals2 / sum(pi_vals2)
    g <- ashr::unimix(pi = pi_vals1, a = a_seq1, b = b_seq1)

    p <- 11
    betahat <- rnorm(p)
    errordist <- list()
    for (index in 1:p) {
        errordist[[index]] <- ashr::unimix(pi = pi_vals2, a = a_seq2, b = b_seq2)
    }

    uuconv_dense(g = g, errordist = errordist[[1]], x = 0)

    lout <- log_compdens_conv_mix.uniuni(g = g, betahat = betahat, errordist = errordist)

    expect_equal(nrow(lout), length(g$pi))
    expect_equal(ncol(lout), p)

    ## some examples where I know the results
    g <- ashr::unimix(pi = 1, a = 0, b = 1)
    errordist <- ashr::unimix(pi = 1, a = 0, b = 1)
    expect_equal(uuconv_dense(g = g, errordist = errordist, x = 0.7),
                 log(0.7))
    expect_equal(uuconv_dense(g = g, errordist = errordist, x = 1.4),
                 log(0.6))


    ## Point mass at 0 example where I know results
    g  <- ashr::unimix(pi = 1, a = 0, b = 0)
    errordist <- ashr::unimix(pi = 1, a = -1, b = 1)
    x <- 0
    my_val <- uuconv_dense(g = g, errordist = errordist, x = x)
    like_val <- log(1 / (errordist$b - errordist$a))
    expect_equal(like_val, my_val)
}
)



test_that("t_to_mix is approximately a t", {
    mu <- 0
    sig <- 1
    df <- 4

    g <- t_to_mix(mean = mu, sd = sig, df = df, gridsize = 20)

    samp <- rnormalmix(n = 10000, mixdense = g)

    pseq <- seq(0.01, 0.99, length = 1000)
    empq <- quantile(samp, probs = pseq)
    theoq <- qt(p = pseq, df = df)

    plot(theoq, empq)
}
)
