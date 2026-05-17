

#' ACF plot for INAR data
#'
#' @description
#' Plot the autocorrelation function (ACF) for INAR data. This is a modified version of the `pacf` function from the package `stats`, adapted to handle INAR data and to return a ggplot object for visualization.
#'
#' @param x, observed series
#' @param lag.max, maximum lag to display
#' @param plot, if TRUE, returns a ggplot object; if FALSE, returns
#' the ACF object
#' @param ..., additional options
#'
#' @return A ggplot object if `plot = TRUE`; otherwise an object of class `"acf"`.
#'
#' @importFrom stats acf qnorm
#' @importFrom ggplot2 ggplot aes geom_hline geom_col
#'
#' @export
ggacf <- function(x, lag.max = NULL, plot = TRUE, ...)
{
    object <- acf(x, lag.max = lag.max, plot = FALSE, ...)

    ci <-  0.95
    clim <- qnorm((1 + ci)/2)/sqrt(object$n.used)

    df <- data.frame("lags" = object$lag[-1], "acf" = object$acf[-1])

    gg <- ggplot(mapping = aes(x = df$lags, y = df$acf)) +
        geom_hline(yintercept = 0) +
        geom_col(stat = "identity", fill = "steelblue", width = 0.15) +
        geom_hline(yintercept = c(-clim,clim), linetype = "dashed", color = "gray30")

    if(plot){
        return(gg)
    }else{
        return(object)
    }
}

#' PACF plot for INAR data
#'
#' @description
#' Plot the partial autocorrelation function (PACF) for INAR data. This is a modified version of the `pacf` function from the package `stats`, adapted to handle INAR data and to return a ggplot object for visualization.
#'#'
#' @param x, observed series
#' @param lag.max, maximum lag to display
#' @param plot, if TRUE,  returns a ggplot object; if FALSE, returns
#' the ACF object
#' @param ..., additional options
#'
#' @return A ggplot object if `plot = TRUE`; otherwise an object of class `"acf"`.
#'
#' @importFrom stats pacf qnorm
#' @importFrom ggplot2 ggplot aes geom_hline geom_col
#'
#' @export
ggpacf <- function(x, lag.max = NULL, plot = TRUE, ...)
{
    object <- pacf(x, lag.max = lag.max, plot = FALSE, ...)

    ci <-  0.95
    clim <- qnorm((1 + ci)/2)/sqrt(object$n.used)
    df <- data.frame("lags" = object$lag, "pacf" = object$acf)

    gg <- ggplot(mapping = aes(x = df$lags, y = df$pacf)) +
        geom_hline(yintercept = 0) +
        geom_col(stat = "identity", fill = "steelblue", width = 0.15) +
        geom_hline(yintercept = c(-clim,clim), linetype = "dashed", color = "gray30")

    if(plot){
        return(gg)
    }else{
        return(object)
    }
}

