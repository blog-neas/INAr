test_that("estimate identical values", {
    set.seed(1234)
    s1 <- genINAR(500, a = 0.5, par = 2, arrival = "poisson")$X
    mod1 <- INARfit(s1,1,"poisson")
    mod2 <- INARfit(s1,1,"poisson")
    expect_equal(mod1, mod2)
})
