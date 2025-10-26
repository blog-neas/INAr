#' Fitting INAR(p) Models
#'
#' @rdname INAR
#' @method print INAR
#'
#' @param x, an `INAR` object.
#' @param digits, number of digits.
#' @param se, whether to show standard error or not.
#' @param ..., additional options
#' @export
print.INAR <- function (x, digits = max(3L, getOption("digits") - 3L), se = TRUE, ...) {
    #
    #
    # TO DO:
    # - add the mean ( = intercept) into the parameter estimates

    cat(x$call0,"\nCall:", deparse(x$call, width.cutoff = 75L),"", sep = "\n")
    if (length(x$coef)) {
        cat("Thinning Parameters:\n")
        alphas <- round(x$coef$alphas, digits = digits)
        if (se && NROW(x$var.coef$alphas)) {
            ses <- rep.int(0, length(alphas))
            ses[x$mask$alphas] <- round(sqrt(diag(x$var.coef$alphas)), digits = digits)
            alphas <- matrix(alphas, 1L, dimnames = list(NULL, names(alphas)))
            alphas <- rbind(alphas, s.e. = ses)
        }
        print.default(alphas, print.gap = 2)

        cat("\nInnovation Parameters,",toupper(x$inn),"Distribution:\n")
        pars <- round(x$coef$pars, digits = digits)
        if (se && NROW(x$var.coef$pars)) {
            ses <- rep.int(0, length(pars))
            ses[x$mask$pars] <- round(sqrt(diag(x$var.coef$pars)), digits = digits)
            pars <- matrix(pars, 1L, dimnames = list(NULL, names(pars)))
            pars <- rbind(pars, s.e. = ses)
        }
        print.default(pars, print.gap = 2)
    }

    # QUESTA PARTE ANDRA' MESSA IN SUMMARY:
    cat("\n - sigma^2 estimated as ", format(x$momentsINN["sigma2"], digits = digits),
        ":  log likelihood = ", format(round(x$loglik, 2L)),
        ",  aic = ", format(round(x$aic, 2L)),",  bic = ", format(round(x$bic, 2L)),"\n",sep = "")
    if (length(x$SMCtest)) {
            cat(" - SMC Significance test = ",x$SMCtest$stat,
                "\n      p-value = ",x$SMCtest$pval,", bootstrap p-value = ",x$SMCtest$pvalboot,"\n",sep = "")
    }
    invisible(x)
}


#' Summarizing INAR(p) Models
#'
#' @rdname INAR
#' @method summary INAR
#'
#' @param object, an `INAR` object
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
    object
}

#' Plotting INAR(p) Models
#'
#' summary method for class `INAR`.
#' @rdname INAR
#' @method plot INAR
#'
#' @param object, an `INAR` object
plot.INAR <- function (object){
    object
    #
    #
    #
    #
}
# methods(summary)

#' Get INAR(p) fitted values
#'
#' summary method for class `INAR`.
#' @rdname INAR
#' @method fitted INAR
#'
#' @param object, an `INAR` object
#' @export
fitted.INAR <- function(object, ...){
    # fitted <- genINAR(object$n, a = object$coef$alphas, par = object$coef$pars, inn = object$inn, burnout = 500)$X
    # return(fitted)
    object
}



# Set methods (S4 style) ...................
# setMethod("print", "INAR", print.INAR)
# setMethod("summary", "INAR", summary.INAR)

