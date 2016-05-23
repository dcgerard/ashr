library(ashr)
context("\nRUVASH")

test_that("ash_ruv add interecept", {
    set.seed(68)
    n <- 10
    p <- 20
    k <- 3
    cov_of_interest <- k
    X <- matrix(stats::rnorm(n * k), nrow = n)
    beta <- matrix(stats::rnorm(k * p), nrow = k)
    beta[, 1:round(p/2)] <- 0
    ctl <- beta[cov_of_interest, ] == 0
    E <- matrix(stats::rnorm(n * p), nrow = n)
    Y <- X %*% beta + E

    num_sv <- 2

    ashout1   <- ash_ruv(Y = Y, X = X, ctl = ctl, k = num_sv, include_intercept = TRUE)
    X <- cbind(X, rep(1, length = nrow(X)))
    ashout2 <- ash_ruv(Y = Y, X = X, ctl = ctl, k = num_sv, cov_of_interest = k,
                       include_intercept = TRUE)
    ashout3 <- ash_ruv(Y = Y, X = X, ctl = ctl, k = num_sv, cov_of_interest = k,
                       include_intercept = FALSE)

    expect_equal(ashout1$fitted.g$pi, ashout2$fitted.g$pi)
    expect_equal(ashout2$fitted.g$pi, ashout3$fitted.g$pi)
    expect_equal(ashout1$fitted.g$pi, ashout3$fitted.g$pi)
}
)


test_that("ash_ruv with k=0 same as ols + ash", {
    set.seed(50)
    n <- 10
    p <- 20
    k <- 3
    cov_of_interest <- k
    X <- matrix(stats::rnorm(n * k), nrow = n)
    beta <- matrix(stats::rnorm(k * p), nrow = k)
    beta[, 1:round(p/2)] <- 0
    ctl <- beta[cov_of_interest, ] == 0
    E <- matrix(stats::rnorm(n * p), nrow = n)
    Y <- X %*% beta + E
    num_sv <- 0

    ruvash_out <- ash_ruv(Y = Y, X = X, ctl = ctl, k = num_sv, cov_of_interest <- k,
                          include_intercept = FALSE)


    xtx_inv <- solve(t(X) %*% X)
    betahat_ols <- xtx_inv %*% t(X) %*% Y
    resid_mat <- X %*% betahat_ols - Y
    col_mse <- colSums(resid_mat ^ 2) / (n - k)
    sebetahat <- sqrt(diag(xtx_inv)[cov_of_interest] * col_mse)
    mult <- sqrt(mean((betahat_ols[cov_of_interest, ctl] / sebetahat[ctl]) ^ 2))

    ashout <- ash.workhorse(betahat = betahat_ols[cov_of_interest, ],
                            sebetahat = sebetahat * mult)

    expect_equal(ruvash_out$ruv$multiplier, mult^ 2)
    expect_equal(ruvash_out$fitted.g$pi, ashout$fitted.g$pi)
}
)


test_that("ash_ruv and ash_ruv_old give same results when using ols", {
    set.seed(68)
    n <- 10
    p <- 20
    k <- 3
    cov_of_interest <- k
    X <- matrix(stats::rnorm(n * k), nrow = n)
    beta <- matrix(stats::rnorm(k * p), nrow = k)
    beta[, 1:round(p/2)] <- 0
    ctl <- beta[cov_of_interest, ] == 0
    E <- matrix(stats::rnorm(n * p), nrow = n)
    Y <- X %*% beta + E

    num_sv <- 2

    ruvash_out <- ash_ruv(Y = Y, X = X, ctl = ctl, k = num_sv,
                            include_intercept = FALSE, gls = FALSE)
    ruvold_out <- ash_ruv_old(Y = Y, X = X, ctl = ctl, k = num_sv,
                              include_intercept = FALSE)

    expect_equal(c(as.matrix(ruvash_out$ruv$betahat)), c(ruvold_out$ruv$betahat))
    expect_equal(ruvash_out$fitted.g$pi, ruvold_out$fitted.g$pi)
}
)


test_that("t-likelihood in ash_ruv works", {
    set.seed(68)
    n <- 10
    p <- 20
    k <- 3
    cov_of_interest <- k
    X <- matrix(stats::rnorm(n * k), nrow = n)
    beta <- matrix(stats::rnorm(k * p), nrow = k)
    beta[, 1:round(p/2)] <- 0
    ctl <- beta[cov_of_interest, ] == 0
    E <- matrix(stats::rnorm(n * p), nrow = n)
    Y <- X %*% beta + E
    num_sv <- 2


    ashout <- ash_ruv(Y = Y, X = X, ctl = ctl, k = num_sv,
                      include_intercept = TRUE, likelihood = "t")
    ashoutnorm <- ash_ruv(Y = Y, X = X, ctl = ctl, k = num_sv,
                          include_intercept = TRUE, likelihood = "normal")
    expect_true(all(ashout$fitted.g$pi != ashoutnorm$fitted.g$pi))
}
)
