test_that("generate INAR series", {
    N <- 500
    a1 <- 0.5
    lam <- 2
    set.seed(1234)
    x1 <- genINAR(N, a = a1, par = lam, arrival = "poisson")
    set.seed(1234)
    x2 <- genINAR(N, a = a1, par = lam, arrival = "poisson")
    expect_equal(x1, x2)
})
