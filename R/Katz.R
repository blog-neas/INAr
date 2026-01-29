#' probability mass function for the Katz distribution
#'
#' This function computes the probability mass function for the Katz distribution.
#'
#' @param x A vector of non-negative integers representing the number of occurrences.
#' @param a A positive parameter of the Katz distribution.
#' @param b A positive parameter of the Katz distribution.
#' @return A vector of probabilities corresponding to the input values.
#' @examples
#' katz_pmf(0:10, a = 2, b = 3)
#' @export
dkatz <- function(x, a, b) {
    stopifnot(x >= 0, x%%1 == 0)
    stopifnot(a > 0, b < 1)
    # if(x < 0) return(0)

    pmf <- ifelse(x < 0,0,choose(a/b + x - 1, x) * (1-b)^(a/b) * (b)^x)
    return(pmf)
}

#' cumulative distribution function for the Katz distribution
#'
#' This function computes the cumulative distribution function for the Katz distribution.
#'
#' @param x A vector of non-negative integers representing the number of occurrences.
#' @param a A positive parameter of the Katz distribution.
#' @param b A positive parameter of the Katz distribution.
#' @return A vector of cumulative probabilities corresponding to the input values.
#' @examples
#' katz_cdf(0:10, a = 2, b = 3)
#' @export
pkatz <- function(x, a, b){
    stopifnot(x >= 0, x%%1 == 0)
    stopifnot(a > 0, b < 1)

    cdf <- sapply(x,function(xi)
        tail(cumsum(dkatz(0:xi, a, b)), 1)
        )
    return(cdf)
}

#' random generation for the Katz distribution
#'
#' This function generates random samples from the Katz distribution.
#'
#' @param n The number of random samples to generate.
#' @param a A positive parameter of the Katz distribution.
#' @param b A positive parameter of the Katz distribution.
#' @return A vector of random samples from the Katz distribution.
#' @examples
#' katz_rnd(10, a = 2, b = 3)
#' @export
rkatz <- function(n, a, b){
    stopifnot(n > 0, n%%1 == 0)
    stopifnot(a > 0, b < 1)

    u <- runif(n)
    samples <- sapply(u, function(ui) {
        x <- 0
        cdf_value <- dkatz(x, a, b)
        while (cdf_value < ui) {
            x <- x + 1
            cdf_value <- cdf_value + dkatz(x, a, b)
        }
        return(x)
    })
    return(samples)
}
# rkatz(1000, a = 2, b = 0.3) |> table() |> barplot()

#' Quantile function for the Katz distribution
#'
#' This function computes the quantile function for the Katz distribution.
#'
#' @param p A vector of probabilities (between 0 and 1).
#' @param a A positive parameter of the Katz distribution.
#' @param b A positive parameter of the Katz distribution.
#'
#' @return A vector of quantiles corresponding to the input probabilities.
#' @examples
#' katz_quantile(c(0.1, 0.5, 0.9), a = 2, b = 3)
#' @export
qkatz <- function(p, a, b){
    stopifnot(all(p >= 0 & p <= 1))
    stopifnot(a > 0, b < 1)

    quantiles <- sapply(p, function(pi) {
        x <- 0
        cdf_value <- dkatz(x, a, b)
        while (cdf_value < pi) {
            x <- x + 1
            cdf_value <- cdf_value + dkatz(x, a, b)
        }
        return(x)
    })
    return(quantiles)
}
# qkatz(c(0.1, 0.5, 0.9), a = 2, b = 0.3)
# pkatz(c(0, 1, 3, 6), a = 2, b = 0.3)
