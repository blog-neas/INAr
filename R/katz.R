#' Katz probability mass function
#'
#' @param x Vector of quantiles.
#' @param a Katz parameter, a >= 0.
#' @param b Katz parameter, b < 1.
#' @param log Logical; if TRUE, probabilities are returned on log scale.
#' @param tol Numerical tolerance.
#'
#' @details
#' The Katz family is defined by Sun and McCabe as:
#'
#'$$
#' p_{j+1} / p_j = (a + b j) / (j + 1), j = 0, 1, 2, ...
#'$$
#'
#' with a > 0 and b < 1.
#'
#' Special cases:
#' - b < 0: Binomial, with m = -a / b and p = b / (b - 1);
#' - b = 0: Poisson, with lambda = a;
#' - 0 < b < 1: Negative Binomial, with r = a / b and p = b.
#'
#' @export
dkatz <- function(x, a, b, log = FALSE, tol = 1e-10) {
    info <- check_katz_par(a, b, tol = tol)

    x0 <- x
    x <- floor(x)

    out_log <- rep(-Inf, length(x))

    valid <- is.finite(x0) & x0 >= 0 & x0 == x

    if (!any(valid)) {
        return(if (log) out_log else exp(out_log))
    }

    xv <- x[valid]

    if (info$type == "degenerate") {
        out_log[valid] <- ifelse(xv == 0, 0, -Inf)
    }

    if (info$type == "poisson") {
        out_log[valid] <- dpois(
            x = xv,
            lambda = a,
            log = TRUE
        )
    }

    if (info$type == "negative-binomial") {
        r <- a / b
        p <- b

        out_log[valid] <- dnbinom(
            x = xv,
            size = r,
            prob = 1 - p,
            log = TRUE
        )
    }

    if (info$type == "binomial") {
        m <- info$N
        p <- b / (b - 1)

        out_log[valid] <- dbinom(
            x = xv,
            size = m,
            prob = p,
            log = TRUE
        )
    }

    if (log) out_log else exp(out_log)
}

#' Katz cumulative distribution function
#'
#' @param q Vector of quantiles.
#' @param a Katz parameter.
#' @param b Katz parameter.
#' @param lower.tail Logical; if TRUE, returns P(X <= q).
#' @param log.p Logical; if TRUE, probabilities are returned on log scale.
#' @param tol Numerical tolerance.
#'
#' @export
pkatz <- function(q,
                  a,
                  b,
                  lower.tail = TRUE,
                  log.p = FALSE,
                  tol = 1e-10) {
    info <- check_katz_par(a, b, tol = tol)

    q <- floor(q)

    if (info$type == "degenerate") {
        ans <- ifelse(q < 0, 0, 1)

        if (!lower.tail) ans <- 1 - ans
        if (log.p) ans <- log(ans)

        return(ans)
    }

    if (info$type == "poisson") {
        return(ppois(
            q = q,
            lambda = a,
            lower.tail = lower.tail,
            log.p = log.p
        ))
    }

    if (info$type == "negative-binomial") {
        r <- a / b
        p <- b

        return(pnbinom(
            q = q,
            size = r,
            prob = 1 - p,
            lower.tail = lower.tail,
            log.p = log.p
        ))
    }

    if (info$type == "binomial") {
        m <- info$N
        p <- b / (b - 1)

        return(pbinom(
            q = q,
            size = m,
            prob = p,
            lower.tail = lower.tail,
            log.p = log.p
        ))
    }
}

#' Katz quantile function
#'
#' @param p Vector of probabilities.
#' @param a Katz parameter.
#' @param b Katz parameter.
#' @param lower.tail Logical; if TRUE, probabilities are P(X <= x).
#' @param log.p Logical; if TRUE, probabilities are supplied on log scale.
#' @param tol Numerical tolerance.
#'
#' @export
qkatz <- function(p,
                  a,
                  b,
                  lower.tail = TRUE,
                  log.p = FALSE,
                  tol = 1e-10) {
    info <- check_katz_par(a, b, tol = tol)

    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p

    if (any(p < 0 | p > 1, na.rm = TRUE)) {
        stop("p must be in [0, 1].", call. = FALSE)
    }

    if (info$type == "degenerate") {
        return(rep(0L, length(p)))
    }

    if (info$type == "poisson") {
        return(qpois(
            p = p,
            lambda = a,
            lower.tail = TRUE,
            log.p = FALSE
        ))
    }

    if (info$type == "negative-binomial") {
        r <- a / b
        pp <- b

        return(qnbinom(
            p = p,
            size = r,
            prob = 1 - pp,
            lower.tail = TRUE,
            log.p = FALSE
        ))
    }

    if (info$type == "binomial") {
        m <- info$N
        pp <- b / (b - 1)

        return(qbinom(
            p = p,
            size = m,
            prob = pp,
            lower.tail = TRUE,
            log.p = FALSE
        ))
    }
}

#' Random generation from the Katz distribution
#'
#' @param n Number of observations.
#' @param a Katz parameter.
#' @param b Katz parameter.
#' @param tol Numerical tolerance.
#'
#' @export
rkatz <- function(n, a, b, tol = 1e-10) {
    if (length(n) != 1 || n < 0 || n != floor(n)) {
        stop("n must be a non-negative integer.", call. = FALSE)
    }

    info <- check_katz_par(a, b, tol = tol)

    if (info$type == "degenerate") {
        return(rep(0L, n))
    }

    if (info$type == "poisson") {
        return(rpois(
            n = n,
            lambda = a
        ))
    }

    if (info$type == "negative-binomial") {
        r <- a / b
        p <- b

        return(rnbinom(
            n = n,
            size = r,
            prob = 1 - p
        ))
    }

    if (info$type == "binomial") {
        m <- info$N
        p <- b / (b - 1)

        return(rbinom(
            n = n,
            size = m,
            prob = p
        ))
    }
}

#' Check Katz parameters according to Sun and McCabe
#'
#' The Katz family is defined by
#'
#' p_{j+1} / p_j = (a + b j) / (j + 1), j = 0, 1, 2, ...
#'
#' with a > 0 and b < 1.
#'
#' Special cases:
#' - b < 0: Binomial-type, finite support, underdispersion;
#' - b = 0: Poisson;
#' - 0 < b < 1: Negative-Binomial-type, overdispersion.
#'
#' @noRd
check_katz_par <- function(a, b, tol = 1e-10) {
    if (!is.finite(a) || !is.finite(b)) {
        stop("Katz parameters must be finite.", call. = FALSE)
    }

    if (a < 0) {
        stop("Katz parameter 'a' must be non-negative.", call. = FALSE)
    }

    if (b >= 1) {
        stop("Katz parameter 'b' must be smaller than 1.", call. = FALSE)
    }

    if (a == 0) {
        if (b != 0) {
            stop(
                "Degenerate Katz distribution with a = 0 is only allowed with b = 0.",
                call. = FALSE
            )
        }

        return(list(
            type = "degenerate",
            support = 0L,
            N = 0L
        ))
    }

    if (abs(b) <= tol) {
        return(list(
            type = "poisson",
            support = NULL,
            N = Inf
        ))
    }

    if (b > 0 && b < 1) {
        return(list(
            type = "negative-binomial",
            support = NULL,
            N = Inf
        ))
    }

    if (b < 0) {
        N_raw <- -a / b
        N <- round(N_raw)

        if (N < 1 || abs(N_raw - N) > tol) {
            stop(
                "For b < 0, the Katz distribution has finite support and requires ",
                "m = -a / b to be a positive integer.",
                call. = FALSE
            )
        }

        return(list(
            type = "binomial",
            support = 0:N,
            N = N
        ))
    }
}

# # Esempi
# set.seed(1)
# b <- 0.4
# a <- 2
#
# x_over <- rkatz(10000, a = a, b = b)
#
# mean(x_over);var(x_over);var(x_over) / mean(x_over)
#
# X <- rkatz(10000, a = a, b = b)
# dens <- dkatz(0:max(X), a = a, b = b)
# plots <- data.frame(x = 0:max(X), dens = dens)
# plot(table(X) / length(X), type = "h")
# points(x = 0:max(X), y =table(X) / length(X), pch = 16)
# lines(plots, col = "red", lty = 2)
# points(plots, col = "red", pch = 16)
#
# #
# # RQEUI
# set.seed(1)
#
# b <- 0
# a <- 3
#
# x_poi <- rkatz(10000, a = a, b = b)
#
# mean(x_poi);var(x_poi);var(x_poi) / mean(x_poi)
#
# X <- rkatz(10000, a = a, b = b)
# dens <- dkatz(0:max(X), a = a, b = b)
# plots <- data.frame(x = 0:max(X), dens = dens)
# plot(table(X) / length(X), type = "h")
# points(x = 0:max(X), y =table(X) / length(X), pch = 16)
# lines(plots, col = "red", lty = 2)
# points(plots, col = "red", pch = 16)
#
# #
# # RUNDER
# N <- 10
# p <- 0.3
#
# b <- -p / (1 - p)
# a <- p * (N + 1) / (1 - p)
#
# x_under <- rkatz(10000, a = a, b = b)
#
# mean(x_under);var(x_under);var(x_under) / mean(x_under)
#
# X <- rkatz(10000, a = a, b = b)
# dens <- dkatz(0:max(X), a = a, b = b)
# plots <- data.frame(x = 0:max(X), dens = dens)
# plot(table(X) / length(X), type = "h")
# points(x = 0:max(X), y =table(X) / length(X), pch = 16)
# lines(plots, col = "red", lty = 2)
# points(plots, col = "red", pch = 16)

