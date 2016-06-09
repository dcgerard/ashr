library(ashr)
context("\nMAD adjustments")

test_that("ash_sva_mat works", {
    set.seed(91)
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

    svaout  <- ash_sva_mad(Y = Y, X = X, k = num_sv)
    ruvout <- ash_ruv_mad(Y = Y, X = X, k = num_sv, ctl = ctl)

    poi_out <- matrix(rpois(prod(dim(Y)), lambda = exp(c(X %*% beta))), nrow = nrow(Y))
    vlema_out <- vlema(Y = poi_out, X = X)
}
)
