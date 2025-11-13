# SOSPESO PER ADESSO IN quanto sto generando nuoe function di stima
# test_that("estimate identical values", {
#     set.seed(1245)
#     s1 <- genINAR(500, a = 0.5, par = 2, inn = "poi")$X
#     mod1 <- INAR(s1, p = 1, inn = "poi", method = "CLS")
#     mod2 <- INAR(s1, p = 1, inn = "poi", method = "CLS")
#     expect_equal(mod1, mod2)
#
#     set.seed(1245)
#     s2 <- genINAR(500, a = c(0.3,0.2), par = c(2,0.66), inn = "negbin")$X
#     mod3 <- INAR(s2, p = 2, inn = "negbin", method = "CLS")
#     mod4 <- INAR(s2, p = 2, inn = "negbin", method = "CLS")
#     expect_equal(mod3, mod4)
#
#     set.seed(1245)
#     s3 <- genINAR(500, a = 0.25, par = 2, inn = "poi")$X
#     mod5 <- INAR(s3, p = 1, inn = "poi", method = "YW")
#     mod6 <- INAR(s3, p = 1, inn = "poi", method = "YW")
#     expect_equal(mod5, mod6)
# })
