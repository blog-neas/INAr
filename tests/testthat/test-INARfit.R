test_that("estimate identical values", {
    set.seed(1245)
    s1 <- genINAR(500, a = 0.5, par = 2, arrival = "poisson")$X
    mod1 <- INAR(s1, order = 1, arrival = "poisson", method = "CLS")
    mod2 <- INAR(s1, order = 1, arrival = "poisson", method = "CLS")
    expect_equal(mod1, mod2)

    set.seed(1245)
    s2 <- genINAR(500, a = c(0.3,0.2), par = c(2,0.66), arrival = "negbin")$X
    mod3 <- INAR(s2, order = 2, arrival = "negbin", method = "CLS")
    mod4 <- INAR(s2, order = 2, arrival = "negbin", method = "CLS")
    expect_equal(mod3, mod4)
})
