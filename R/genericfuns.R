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


#' Summary of INAR forecast
#'
#' summary method for class `INARforecast`.
#' @rdname INARforecast
#' @method summary INARforecast
#'
#' @param object an `INARforecast` object
#' @param ..., additional options
#'
#' @return the original object, invisibly.
#' @export
summary.INARforecast <- function(object, ...){

    cat("\nINAR forecast\n")
    cat("Call:\n")
    print(object$call)
    cat("Horizon:", object$n.ahead, "step ahead(s) \n")
    warning("summary.INARforecast is not yet implemented. Returning the original object.")
    invisible(object)

    #
    # if (x$method == "bootstrap") {
    #     cat("Bootstrap replications:", x$B, "\n")
    #     cat("Prediction interval level:", x$level, "\n\n")
    #
    #     tab <- data.frame(
    #         h = seq_len(x$n.ahead),
    #         mean = x$mean,
    #         median = x$median,
    #         lower = x$lower,
    #         upper = x$upper
    #     )
    #
    # } else {
    #
    #     cat("\n")
    #
    #     tab <- data.frame(
    #         h = seq_len(x$n.ahead),
    #         mean = x$mean
    #     )
    # }
    #
    # print(round(tab, digits), row.names = FALSE)
    #
    # invisible(x)
}

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
    stopifnot(inherits(object, "INAR"))

    fitted <- INARfitted_cpp(object$data, object$mINN, object$alphas)
    return(fitted)
}


#' Forecast method for INAR(p) models
#' predict method for class `INAR`.
#'
#' @rdname INAR
#' @method predict INAR
#' @importFrom stats median quantile
#'
#' @param object an object of class "INAR".
#' @param n.ahead forecast horizon.
#' @param type either "mean" or "bootstrap".
#' @param B number of bootstrap trajectories.
#' @param level prediction interval level.
#' @param seed optional random seed.
#' @param ... further arguments.
#'
#' @return an object of class "INARforecast".
#' @export
predict.INAR <- function(object,
                         n.ahead = 1,
                         type = c("mean", "bootstrap"),
                         B = 999,
                         level = 0.95,
                         seed = NULL,
                         ...) {

    type <- match.arg(type)

    stopifnot(B > 1, n.ahead >= 1)
    stopifnot(inherits(object, "INAR"))

    if(type == "mean"){
        forecast <- INARforecast_cpp(object$data , object$mINN, object$alphas, h = n.ahead, B = 0)$forecast
        lower <- NA
        upper <- NA
        forecastmedian <- NA
    }else{
        if(!is.null(seed)) set.seed(seed)
        Bdata <- INARforecast_cpp(object$data , object$mINN, object$alphas, h = n.ahead, B = B)
        forecast <- Bdata$forecast
        paths <- Bdata$paths
        forecastmedian <- apply(paths, 1, median)
        lower <- apply(paths, 1, quantile, probs = (1 - level) / 2)
        upper <- apply(paths, 1, quantile, probs = 1 - ( 1 - level) / 2)
    }

    OUT <- list(
        forecast = forecast,
        forecastmedian = forecastmedian,
        lower = lower,
        upper = upper,
        call = object$call,
        n.ahead = n.ahead,
        type = type,
        B = B,
        level = level
    )
    class(OUT) <- "INARforecast"
    return(OUT)
}


#' Print method for INAR forecast
#' Print method for class `INARforecast`.
#'
#' @rdname INARforecast
#' @method print INARforecast
#' @export
#' @importFrom utils head
#'
#' @param x an object of class "INARforecast".
#' @param ... further arguments.
#'
#' @return the original object, invisibly.
#' @export
print.INARforecast <- function(x, ...) {

    cat("\nINAR forecast\n")
    cat("Horizon:", x$n.ahead, "step ahead(s)\n")
    cat("Method:", x$type, "\n\n")

    if (x$type == "mean") {

        result <- data.frame(
            step = seq_len(x$n.ahead),
            forecast = x$forecast
        )

    } else {

        result <- data.frame(
            step = seq_len(x$n.ahead),
            forecast = x$forecast,
            median = x$forecastmedian,
            lower = x$lower,
            upper = x$upper
        )
    }

    print(result, row.names = FALSE)

    invisible(x)
}


#' Plotting INAR(p) Models
#'
#' summary method for class `INAR`.
#' @rdname INAR
#' @method plot INAR
#'
#' @importFrom ggplot2 ggplot geom_line labs
#'
#' @param x, an `INAR` object
#' @param ..., additional options
#' @return the original object, invisibly.
#' @export
plot.INAR <- function (x, ...){

    data <- x$data
    df <- data.frame(time = length(data), value = data) # , resid = data$residuals, stdresid = data$stdresiduals)

    gg1 <- ggplot(mapping = aes(x=df$time, y=df$value)) +
        geom_line(color = "darkorange", linewidth = 1) +
        labs(x = "Time",
             y = NULL)

    gg1

    # gg2 <- ggacf(df$value, lag.max = 20)
    # gg3 <- ggpacf(df$value, lag.max = 20)
    #
    # gg1 / (gg2 + gg3) # richiede l'import di patchwork

    # gg4 <- ggplot(df, aes(x = value)) +
    #     geom_histogram(aes(y = ..density..), bins = max(df$value), fill = "steelblue", color = "black", alpha = 0.7) +
    #     geom_density(color = "darkorange", size = 1) +
    #     theme_minimal() +
    #     labs(x = "Value",
    #          y = "Density")
    #
    # gg5 <- ggplot(df, aes(x = time, y = resid)) +
    #     geom_line(color = "steelblue", linewidth = 1) +
    #     theme_minimal() +
    #     labs(x = "Time",
    #          y = "Residuals")

    invisible(x)
}



#' Plotting forecast INAR(p) Models
#'
#' plot method for class `INARforecast`.
#' @rdname INARforecast
#' @method plot INARforecast
#'
#' @param x, an `INARforecast` object
#' @param ..., additional options
#' @return the original object, invisibly.
#' @export
plot.INARforecast <- function (x, ...){
    warning("plot.INARforecast is not yet implemented. Returning the original object.")
    invisible(x)
}



#' Plot bootstrap distribution of an INAR test
#'
#' @method plot INARtest
#'
#' @importFrom ggplot2 ggplot stat_ecdf geom_density stat_function labs theme_bw
#' @importFrom stats ks.test pnorm
#'
#' @param x an object of class `INARtest`.
#' @param which a character string specifying the type of plot to produce. Options are "density" for a density plot or "ecdf" for an empirical cumulative distribution function (ECDF) plot.
#' @param ... additional arguments.
#'
#' @return A `ggplot` object.
#' @export
plot.INARtest <- function(x, which = "density", ...) {

    if (!inherits(x, "INARtest")) {
        stop("'x' must be an object of class 'INARtest'.", call. = FALSE)
    }

    if (x$test =="smc") {
        # inserire parametrizzazioni specifici per SMC
        # se serve, altrimenti eliminare questi if

        if (is.null(x$bootvec)) {
            stop(
                "No bootstrap replications were found in the object.",
                "Run SMCtest(..., saveboot = TRUE).",
                call. = FALSE
            )
        }

    } else {
        stop(
            "This plot method is currently implemented only for SMC tests.",
            call. = FALSE
        )
    }


    # controllo per presenza di replicazioni bootstrap finite
    boot <- as.numeric(x$bootvec)
    boot <- boot[is.finite(boot)]

    if (length(boot) < 2L) {
        stop(
            "At least two finite bootstrap replications are required.",
            call. = FALSE
        )
    }

    # normal_mean <- unname(x$normal.reference["mean"])
    # normal_sd <- unname(x$normal.reference["sd"])
    # reference vc va cambiata solo per test di dispersion, per il resto è una Normale Standard
    ks_result <- ks.test(boot,"pnorm",mean = 0,sd = 1,exact = FALSE)

    ks_statistic <- unname(ks_result$statistic)
    ks_pvalue <- ks_result$p.value

    plot_subtitle <- sprintf(
        "Kolmogorov-Smirnov test: D = %.4f, p-value = %s",
        ks_statistic,
        format.pval(ks_pvalue,digits = 4,eps = 0.0001)
    )

    plot_data <- data.frame(
        statistic = boot
    )


    if(which == "density") {
        gg <- ggplot(plot_data,aes(x = statistic)) +
            # density of a standard normal distribution:
            stat_function(fun = dnorm,args = list(mean = 0,sd = 1),
                          linewidth = 0.9,colour = "darkred",linetype = "dashed") +
            geom_density(linewidth = 0.9,colour = "navy") +
            labs(
                title = "Bootstrap distribution",
                subtitle = plot_subtitle,
                x = paste(toupper(x$test),"Bootstrap statistics"),
                y = NULL
            )
    }else if(which == "ecdf") {
        gg <- ggplot(plot_data,aes(x = statistic)) +
            stat_ecdf(geom = "step",linewidth = 0.9,colour = "navy") +
            stat_function(fun = pnorm,args = list(mean = 0,sd = 1),
                linewidth = 0.9,colour = "darkred",linetype = "dashed") +
            labs(
                title = "Bootstrap distribution",
                subtitle = plot_subtitle,
                x = paste(toupper(x$test),"Bootstrap statistics"),
                y = NULL
                )
    }
    return(gg)
}


# Set methods (S4 style) ...................
# setMethod("print", "INAR", print.INAR)
# setMethod("summary", "INAR", summary.INAR)

