
#' Matrix of the first p lagged observations
#'
#' Internal function
#'
#' @param x, observed series
#' @param p, number of lags
#' @param cut, if TRUE saves the matrix with a smaller number of rows
#' @noRd
lagmat <- function(x, p, cut = TRUE){
    n <- length(x)
    stopifnot(n > p)
    Xmat <- matrix(0,nrow = n, ncol = p + 1)
    Xmat[,1] <- 1
    for(j in 1:p){
        # Xmat[,j+1] <- c(rep(0,j),x[1:(n - j)])
        Xmat[,j+1] <- c(rep(NA, j), x)[1:n]

    }
    if(cut) Xmat <- Xmat[(p+1):n,]

    return(Xmat)
}
