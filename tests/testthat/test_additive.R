library(ashr)
context("Additive Inflation")

test_that("additive_update_zl1 and additive_update_l2 increase likelihood", {
    set.seed(81)
    p <- 19
    k <- 7

    alpha    <- matrix(stats::rnorm(p * k), nrow = p)
    Z        <- matrix(stats::rnorm(k), ncol = 1)
    sig_diag <- stats::rchisq(p, df = 5) / 5
    E        <- matrix(stats::rnorm(p, sd = sqrt(2 * sig_diag + 1)), ncol = 1)
    Y <- alpha %*% Z + E

    Zinit <- matrix(0, nrow = k, ncol = 1)
    lambda1_init <- 1
    lambda2_init <- 0

    numiter <- 20
    llike_vec <- rep(NA, numiter)

    llike_vec[1] <- additive_obj(Z = Zinit, lambda1 = lambda1_init,
                                 lambda2 = lambda2_init, Y = Y,
                                 alpha = alpha, sig_diag = sig_diag)

    Zcurrent <- Zinit
    lambda1_current <- lambda1_init
    lambda2_current <- lambda2_init

    for (index in 2:numiter) {
        aout <- additive_update_zl1(Z = Zcurrent,
                                    lambda1 = lambda1_current,
                                    lambda2 = lambda2_current,
                                    Y = Y, alpha = alpha,
                                    sig_diag = sig_diag)
        Zcurrent <- aout$Znew
        lambda1_current <- aout$lambda1_new
        lambda2_current <- additive_update_l2(Z = Zcurrent,
                                              lambda1 = lambda1_current,
                                              lambda2 = lambda2_current,
                                              Y = Y, alpha = alpha,
                                              sig_diag = sig_diag)
        llike_vec[index] <- additive_obj(Z = Zcurrent, lambda1 = lambda1_current,
                                         lambda2 = lambda2_current, Y = Y,
                                         alpha = alpha, sig_diag = sig_diag)
    }

    expect_true(all(llike_vec[1:(numiter - 1)] <= llike_vec[2:numiter]))


    llike_vec2 <- rep(NA, length = numiter)
    zlambda <- c(Zinit, lambda1_init, lambda2_init)
    llike_vec2[1] <- additive_obj_wrapper(zlambda = zlambda,
                                              Y = Y, alpha = alpha,
                                              sig_diag = sig_diag)
    for (index in 2:numiter) {
        zlambda <- additive_fix_wrapper(zlambda = zlambda, Y = Y,
                                        alpha = alpha,
                                        sig_diag = sig_diag)
        llike_vec2[index] <- additive_obj_wrapper(zlambda = zlambda,
                                                  Y = Y,
                                                  alpha = alpha,
                                                  sig_diag = sig_diag)
    }

    expect_equal(llike_vec, llike_vec2)

}
)

test_that("ash_ruv works with additive inflation", {
    set.seed(99)
    n <- 11
    p <- 31
    k <- 2
    q <- 3

    cov_of_interest <- k
    ctl <- rep(FALSE, length = p)
    ctl[1:10] <- TRUE

    alpha <- matrix(stats::rnorm(q * p), nrow = q)
    Z     <- matrix(stats::rnorm(n * q), nrow = n)
    beta  <- matrix(stats::rnorm(k * p), nrow = k)
    beta[cov_of_interest, 1:20] <- 0
    X     <- matrix(stats::rnorm(n * k), nrow = n)
    sig_diag <- stats::rchisq(p, df = 5) / 5
    E     <- matrix(stats::rnorm(n * p), nrow = n) %*% diag(sqrt(sig_diag))

    Y <- X %*% beta + Z %*% alpha + E

    ruvash_out <- ash_ruv(Y = Y, X = X, ctl = ctl, k = q, cov_of_interest = cov_of_interest,
                          additive_inflate = TRUE)
    expect_error(ash_ruv(Y = Y, X = X, ctl = ctl, k = q, cov_of_interest = cov_of_interest,
                          additive_inflate = TRUE, likelihood = "t"))
}
)
