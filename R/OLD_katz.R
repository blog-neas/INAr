#
# OLD SCRIPT FOR KATZ DISTRIBUTION
# This code will be deleted in the next versions
#

# dkatz <- function(x, a, b) {
#     stopifnot(x >= 0, x%%1 == 0)
#     stopifnot(a > 0, b < 1)
#
#     pmf <- ifelse(x < 0,0,choose(a/b + x - 1, x) * (1-b)^(a/b) * (b)^x)
#     return(pmf)
# }
#
# pkatz <- function(x, a, b){
#     stopifnot(x >= 0, x%%1 == 0)
#     stopifnot(a > 0, b < 1)
#
#     cdf <- sapply(x,function(xi)
#         tail(cumsum(dkatz(0:xi, a, b)), 1)
#         )
#     return(cdf)
# }
#
# rkatz <- function(n, a, b){
#     stopifnot(n > 0, n%%1 == 0)
#     stopifnot(a > 0, b < 1)
#
#     u <- runif(n)
#     samples <- sapply(u, function(ui) {
#         x <- 0
#         cdf_value <- dkatz(x, a, b)
#         while (cdf_value < ui) {
#             x <- x + 1
#             cdf_value <- cdf_value + dkatz(x, a, b)
#         }
#         return(x)
#     })
#     return(samples)
# }
#
# qkatz <- function(p, a, b){
#     stopifnot(all(p >= 0 & p <= 1))
#     stopifnot(a > 0, b < 1)
#
#     quantiles <- sapply(p, function(pi) {
#         x <- 0
#         cdf_value <- dkatz(x, a, b)
#         while (cdf_value < pi) {
#             x <- x + 1
#             cdf_value <- cdf_value + dkatz(x, a, b)
#         }
#         return(x)
#     })
#     return(quantiles)
# }
