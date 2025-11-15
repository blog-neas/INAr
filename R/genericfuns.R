#' Print the results of INAR tests.
#'
#' This is a modified version of the print function for class `htest` from the package `stats`.
#' @rdname INARtest
#' @method print INARtest
#'
#' @param x, an `INARtest` object.
#' @param digits, number of digits.
#' @param prefix, prefix for the output.
#' @param ..., additional options.
#' @export
print.INARtest <- function (x, digits = getOption("digits"), prefix = "\t", ...){
    cat("\n")
    cat(strwrap(x$method, prefix = prefix), sep = "\n")
    cat("\n")
    cat("data:  ", x$data.name, "\n", sep = "")
    out <- character()
    out2 <- character()
    if (!is.null(x$statistic))
        out <- c(out, paste(names(x$statistic), "=", format(x$statistic,
                                                            digits = max(1L, digits - 2L))))
    if (!is.null(x$parameter))
        out <- c(out, paste(names(x$parameter), "=", format(x$parameter,
                                                            digits = max(1L, digits - 2L))))
    if (!is.null(x$p.value)) {
        fp <- format.pval(x$p.value, digits = max(1L, digits -
                                                      3L))
        out <- c(out, paste("p-value", if (startsWith(fp, "<")) fp else paste("=",
                                                                              fp)))
    }
    if (!is.null(x$statistic.boot))
        out2 <- c(out2, paste(names(x$statistic.boot), "=", format(x$statistic.boot,
                                                            digits = max(1L, digits - 2L))))
    if (!is.null(x$parameter.boot))
        out2 <- c(out2, paste(names(x$parameter.boot), "=", format(x$parameter.boot,
                                                            digits = max(1L, digits - 2L))))
    if (!is.null(x$p.value.boot)) {
        fp.boot <- format.pval(x$p.value.boot, digits = max(1L, digits -
                                                      3L))
        out2 <- c(out2, paste("p-value bootstrap", if (startsWith(fp.boot, "<")) fp.boot else paste("=",
                                                                              fp.boot)))
    }

    cat(strwrap(paste(out, collapse = ", ")),strwrap(paste(out2, collapse = ", ")), sep = "\n")

    if (!is.null(x$alternative)) {
        cat("alternative hypothesis: ")
        if (!is.null(x$null.value)) {
            if (length(x$null.value) == 1L) {
                alt.char <- switch(x$alternative, two.sided = "not equal to",
                                   less = "less than", greater = "greater than")
                cat("true ", names(x$null.value), " is ", alt.char,
                    " ", x$null.value, "\n", sep = "")
            }
            else {
                cat(x$alternative, "\nnull values:\n", sep = "")
                print(x$null.value, digits = digits, ...)
            }
        }
        else cat(x$alternative, "\n", sep = "")
    }
    if (!is.null(x$conf.int)) {
        cat(format(100 * attr(x$conf.int, "conf.level")), " percent confidence interval:\n",
            " ", paste(format(x$conf.int[1:2], digits = digits),
                       collapse = " "), "\n", sep = "")
    }
    if (!is.null(x$estimate)) {
        cat("sample estimates:\n")
        print(x$estimate, digits = digits, ...)
    }
    cat("\n")
    invisible(x)
}


#' Summarizing INAR(p) Models
#'
#' @rdname INAR
#' @method summary INAR
#'
#' @param object, an `INAR` object
#' @param ..., additional options
#' @export
summary.INAR <- function (object, ...){
    #
    #
    # TO DO:
    # - add the mean ( = intercept) into the parameter estimates

    # PRESA DA SUMMARY.LM
    #
    # z <- object
    # p <- z$rank
    # rdf <- z$df.residual
    # if (p == 0) {
    #     r <- z$residuals
    #     n <- length(r)
    #     w <- z$weights
    #     if (is.null(w)) {
    #         rss <- sum(r^2)
    #     }
    #     else {
    #         rss <- sum(w * r^2)
    #         r <- sqrt(w) * r
    #     }
    #     resvar <- rss/rdf
    #     ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    #     class(ans) <- "summary.lm"
    #     ans$aliased <- is.na(coef(object))
    #     ans$residuals <- r
    #     ans$df <- c(0L, n, length(ans$aliased))
    #     ans$coefficients <- matrix(NA_real_, 0L, 4L, dimnames = list(NULL,
    #                                                                  c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))
    #     ans$sigma <- sqrt(resvar)
    #     ans$r.squared <- ans$adj.r.squared <- 0
    #     ans$cov.unscaled <- matrix(NA_real_, 0L, 0L)
    #     if (correlation)
    #         ans$correlation <- ans$cov.unscaled
    #     return(ans)
    # }
    # if (is.null(z$terms))
    #     stop("invalid 'lm' object:  no 'terms' component")
    # if (!inherits(object, "lm"))
    #     warning("calling summary.lm(<fake-lm-object>) ...")
    # Qr <- qr.lm(object)
    # n <- NROW(Qr$qr)
    # if (is.na(z$df.residual) || n - p != z$df.residual)
    #     warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    # r <- z$residuals
    # f <- z$fitted.values
    # w <- z$weights
    # if (is.null(w)) {
    #     mss <- if (attr(z$terms, "intercept"))
    #         sum((f - mean(f))^2)
    #     else sum(f^2)
    #     rss <- sum(r^2)
    # }
    # else {
    #     mss <- if (attr(z$terms, "intercept")) {
    #         m <- sum(w * f/sum(w))
    #         sum(w * (f - m)^2)
    #     }
    #     else sum(w * f^2)
    #     rss <- sum(w * r^2)
    #     r <- sqrt(w) * r
    # }
    # resvar <- rss/rdf
    # if (is.finite(resvar) && resvar < (mean(f)^2 + var(c(f))) *
    #     1e-30)
    #     warning("essentially perfect fit: summary may be unreliable")
    # p1 <- 1L:p
    # R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    # se <- sqrt(diag(R) * resvar)
    # est <- z$coefficients[Qr$pivot[p1]]
    # tval <- est/se
    # ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    # ans$residuals <- r
    # ans$coefficients <- cbind(Estimate = est, `Std. Error` = se,
    #                           `t value` = tval, `Pr(>|t|)` = 2 * pt(abs(tval), rdf,
    #                                                                 lower.tail = FALSE))
    # ans$aliased <- is.na(z$coefficients)
    # ans$sigma <- sqrt(resvar)
    # ans$df <- c(p, rdf, NCOL(Qr$qr))
    # if (p != attr(z$terms, "intercept")) {
    #     df.int <- if (attr(z$terms, "intercept"))
    #         1L
    #     else 0L
    #     ans$r.squared <- mss/(mss + rss)
    #     ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n -
    #                                                          df.int)/rdf)
    #     ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
    #                         numdf = p - df.int, dendf = rdf)
    # }
    # else ans$r.squared <- ans$adj.r.squared <- 0
    # ans$cov.unscaled <- R
    # dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,
    #                                                            1)]
    # if (correlation) {
    #     ans$correlation <- (R * resvar)/outer(se, se)
    #     dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    #     ans$symbolic.cor <- symbolic.cor
    # }
    # if (!is.null(z$na.action))
    #     ans$na.action <- z$na.action
    # class(ans) <- "summary.lm"
    # ans
    warning("summary.INAR is not yet implemented. Returning the original object.")
    invisible(object)
}

#' Plotting INAR(p) Models
#'
#' summary method for class `INAR`.
#' @rdname INAR
#' @method plot INAR
#'
#' @param object, an `INAR` object
plot.INAR <- function (object){
    warning("plot.INAR is not yet implemented. Returning the original object.")
    invisible(object)
}
# methods(summary)

#' Get INAR(p) fitted values
#'
#' summary method for class `INAR`.
#' @rdname INAR
#' @method fitted INAR
#'
#' @param object, an `INAR` object
#' @param ..., additional options
#' @export
fitted.INAR <- function(object, ...){
    warning("fitted.INAR is not yet implemented. Returning the original object.")
    invisible(object)
}



# Set methods (S4 style) ...................
# setMethod("print", "INAR", print.INAR)
# setMethod("summary", "INAR", summary.INAR)

