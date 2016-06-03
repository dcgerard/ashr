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

    ruvash_out <- ash_ruv(Y = Y, X = X, ctl = ctl, k = num_sv, cov_of_interest = k,
                          include_intercept = FALSE)


    xtx_inv <- solve(t(X) %*% X)
    betahat_ols <- xtx_inv %*% t(X) %*% Y
    resid_mat <- X %*% betahat_ols - Y
    col_mse <- colSums(resid_mat ^ 2) / (n - k)
    sebetahat <- sqrt(diag(xtx_inv)[cov_of_interest] * col_mse)
    mult <- sqrt(mean((betahat_ols[cov_of_interest, ctl] / sebetahat[ctl]) ^ 2) *
                 n / (n - k))

    ashout <- ash.workhorse(betahat = betahat_ols[cov_of_interest, ],
                            sebetahat = sebetahat * mult)

    expect_equal(ruvash_out$ruv$multiplier, mult^ 2)
    expect_equal(ruvash_out$fitted.g$pi, ashout$fitted.g$pi)
}
)


test_that("ash_ruv and ash_ruv_old give same results when using ols", {
    set.seed(68)
    n <- 11
    p <- 19
    k <- 3
    cov_of_interest <- k
    X <- matrix(stats::rnorm(n * k), nrow = n)
    beta <- matrix(stats::rnorm(k * p), nrow = k)
    beta[, 1:round(p/2)] <- 0
    ctl <- beta[cov_of_interest, ] == 0
    E <- matrix(stats::rnorm(n * p), nrow = n)
    Y <- X %*% beta + E

    num_sv <- 2

    if (requireNamespace("ruv", quietly = TRUE)) {
        ruv4_out <- ruv::RUV4(Y = Y, X = X[, cov_of_interest, drop = FALSE],
                              Z = X[, -cov_of_interest, drop = FALSE],
                              ctl = ctl, k = num_sv)

        ruvash_out <- ash_ruv(Y = Y, X = X, ctl = ctl, k = num_sv,
                              include_intercept = FALSE, gls = FALSE)
        ruvold_out <- ash_ruv_old(Y = Y, X = X, ctl = ctl, k = num_sv,
                                  include_intercept = FALSE)


        ## ruv4_out$multiplier is term from design INCLUDING
        ## unobserved covariates.  1 / ruvash_out$ruv$fnorm_x ^ 2 is
        ## the term from the design NOT including unobserved
        ## covariates.

        ## make sure I understand RUV4 output
        names(ruv4_out)
        XZW <- cbind(ruv4_out$X, ruv4_out$Z, ruv4_out$W)
        expect_equal(solve(t(XZW) %*% XZW)[1, 1], c(ruv4_out$multiplier))
        betahat_xzw <- solve(t(XZW) %*% XZW) %*% t(XZW) %*% Y
        sig2_xzw <- colSums((Y - XZW %*% betahat_xzw) ^ 2) / (nrow(XZW) - ncol(XZW))
        expect_equal(solve(t(XZW) %*% XZW)[1, 1], c(ruv4_out$multiplier))
        expect_equal(betahat_xzw[1, ], c(ruv4_out$betahat))
        expect_equal(sig2_xzw, ruv4_out$sigma2)

        expect_equal(solve(t(X) %*% X)[cov_of_interest, cov_of_interest],
                     1 / ruvash_out$ruv$fnorm_x ^ 2)

        expect_equal(c(ruv4_out$betahat), c(ruvash_out$ruv$betahat))
        expect_equal(ruvash_out$ruv$sigma2, ruv4_out$sigma2)

        expect_equal(c(as.matrix(ruvash_out$ruv$betahat)), c(ruvold_out$ruv$betahat))
        expect_equal(ruvold_out$ruv$sigma2, ruv4_out$sigma2)
    }
}
)


test_that("tregress_em increases likelihood", {
    set.seed(871)
    p  <- 21
    k  <- 5
    nu <- 5

    alpha <- matrix(stats::rnorm(p * k), nrow = p)
    Z     <- matrix(stats::rnorm(k), ncol = 1)
    sig_diag <- stats::rchisq(p, df = 2)
    E <- matrix(stats::rt(p, df = nu), ncol = 1) * sqrt(sig_diag)

    Y <- alpha %*% Z + E

    lambda_init <- 1
    Z_init <- rep(0, length = ncol(alpha))
    zlambda <- c(Z_init, lambda_init)

    itermax <- 20
    llike_vec <- rep(NA, length = itermax)
    llike_vec[1] <- tregress_obj(zlambda = zlambda, Y = Y, alpha = alpha,
                                 sig_diag = sig_diag, nu = nu)

    for (index in 2:itermax) {
        zlambda <- tregress_fix(zlambda = zlambda, Y = Y, alpha = alpha,
                                sig_diag = sig_diag, nu = nu)
        llike_vec[index] <- tregress_obj(zlambda = zlambda, Y = Y, alpha = alpha,
                                         sig_diag = sig_diag, nu = nu)
    }
    expect_true(all(llike_vec[2:length(llike_vec)] >= llike_vec[1:(length(llike_vec) - 1)]))

    oout <- stats::optim(par = c(Z_init, lambda_init), fn = tregress_obj, Y = Y,
                         alpha = alpha, sig_diag = sig_diag, nu = nu,
                         control = list(fnscale = -1, maxit = 5000))

    expect_equal(llike_vec[index], oout$value, tol = 10 ^ -5)

    tregress_obj(zlambda = oout$par, Y = Y, alpha = alpha, sig_diag = sig_diag, nu = nu)
    tregress_obj(zlambda = zlambda, Y = Y, alpha = alpha, sig_diag = sig_diag, nu = nu)

    expect_equal(zlambda, oout$par, tol = 10 ^ -3)

}
)

test_that("when given no control genes, same as OLS + ASH", {
    set.seed(68)
    n <- 10
    p <- 20
    k <- 3
    cov_of_interest <- k
    X <- matrix(stats::rnorm(n * k), nrow = n)
    beta <- matrix(stats::rnorm(k * p), nrow = k)
    beta[, 1:round(p/2)] <- 0
    ctl <- NULL
    E <- matrix(stats::rnorm(n * p), nrow = n)
    Y <- X %*% beta + E
    num_sv <- 2

    ruvash_out <- ash_ruv(Y = Y, X = X, ctl = ctl, k = num_sv,
                      include_intercept = FALSE, likelihood = "normal",
                      gls = FALSE)

    xtxi <- solve(t(X) %*% X)
    betahat_ols <-  xtxi %*% t(X) %*% Y
    sigma2 <- colSums((Y - X %*% betahat_ols) ^ 2) / (n - k)
    sebetahat_ols <- sqrt(sigma2 * xtxi[cov_of_interest, cov_of_interest])
    ash_out <- ashr::ash(betahat = betahat_ols[cov_of_interest, ],
                         sebetahat = sebetahat_ols)

    expect_equal(c(ruvash_out$ruv$betahat_ols), betahat_ols[cov_of_interest, ])
    expect_equal(c(ruvash_out$ruv$betahat), betahat_ols[cov_of_interest, ])
    expect_equal(ruvash_out$ruv$sebetahat, ruvash_out$ruv$sebetahat_ols)
    expect_equal(ruvash_out$ruv$sebetahat, sebetahat_ols)
    expect_equal(ash_out$fitted.g, ruvash_out$fitted.g)
}
)


test_that("ruvash with t-likelihood converges to ruvash with normal likelihood", {
    set.seed(31)
    nseq <- c(11, 19, 59, 103)

    p <- 23
    k <- 3
    num_sv <- 5

    ## see if t results converge to normal results
    tmult_vec <- rep(NA, length = length(nseq))
    nmult_vec <- rep(NA, length = length(nseq))

    tz_mat <- matrix(NA, nrow = num_sv, ncol = length(nseq))
    nz_mat <- matrix(NA, nrow = num_sv, ncol = length(nseq))

    for (index in 1:length(nseq)) {
        n <- nseq[index]
        cov_of_interest <- k
        X <- matrix(stats::rnorm(n * k), nrow = n)
        beta <- matrix(stats::rnorm(k * p), nrow = k)
        beta[, 1:round(p/2)] <- 0
        ctl <- beta[cov_of_interest, ] == 0
        E <- matrix(stats::rnorm(n * p), nrow = n)
        Y <- X %*% beta + E


        ruvash_t <- ash_ruv(Y = Y, X = X, ctl = ctl, k = num_sv, cov_of_interest = k,
                            include_intercept = FALSE, likelihood = "t")
        ruvash_n <- ash_ruv(Y = Y, X = X, ctl = ctl, k = num_sv, cov_of_interest = k,
                            include_intercept = FALSE, likelihood = "normal")
        tmult_vec[index] <- ruvash_t$ruv$multiplier
        nmult_vec[index] <- ruvash_n$ruv$multiplier
        tz_mat[, index]  <- ruvash_t$ruv$Z1
        nz_mat[, index]  <- ruvash_n$ruv$Z1
    }

    err <- colSums((tz_mat - nz_mat) ^ 2)
    expect_true(all(err[1:(length(err) - 1)] > err[2:length(err)]))
    merr <- (tmult_vec - nmult_vec) ^ 2
    expect_true(all(merr[1:(length(merr) - 1)] > merr[2:length(merr)]))
}
)
