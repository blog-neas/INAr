# Experimental functions for MINAR(p) models, with poisson innovations.

#' Generate MINAR(p) models with different innovations
#'
#' Experimental function
#'
#' @param n integer, number of observations, length of the series.
#' @param a vector, thinning parameters. The lenght of this vector defines the number of lags `p` of the INAR(p) process.
#' @param par vector, parameters related with the model, see the details section.
#' @param inn character, innovation distribution. Default value is `"poi"`, alternative values are `"negbin"` for Negative Binomial, `"genpoi"` for Generalized Poisson and `"katz"` for the Katz family.
#' @param burnout integer, number of starting observations to discard. Set to 500 by default.
#'
#' @examples
#'
#' @export
genMINAR <- function(n, m, A_list, lambdas, inn="poi", burnout=500){
    stopifnot(inn %in% c("poi"))
    lags <- length(A_list)
    inn <- tolower(inn)
    s <- n + burnout

    resid_ <- matrix(NA,s,m)

    # Check for stationarity
    A_sum <- Reduce("+", A_list)
    rho <- max(Mod(eigen(A_sum)$values))
    A_mat <- do.call(cbind, A_list)
    if(rho >= 1){
        stop(paste("The resulting process is NOT stationary: specrtal radius =", round(rho, 3)))
    }

    X <- matrix(0, nrow = s, ncol = m)

    # Poisson marginal
    if(inn == "poi"){
        # stopifnot(length(lambdas) == m)

        resid_ <- matrix(0, nrow = s, ncol = m)
        for(j in 1:m){
            resid_[,j] <- rpois(s, lambdas[j])
        }
    }
    else{
        stop("please specify one of the available distributions", call. = FALSE)
    }

    X <- MINARp_gen_cpp(resid_, A_mat)

    X_sim <- X[(burnout+1):s,]

    return(X_sim)
}

# example
# n <- 1000
# m <- 3
# A1 <- matrix(c(0.3,0.1,0.2, 0.2,0.4,0.1, 0.1,0.2,0.3), nrow = m, byrow = TRUE)
# A2 <- matrix(c(0.1,0.05,0.1, 0.05,0.1,0.05, 0.1,0.05,0.1), nrow = m, byrow = TRUE)
# A_list <- list(A1, A2)
# lambdas <- c(2.5, 5, 1)
# set.seed(1234)
# sim_MINAR <- genMINAR(n, m, A_list, lambdas, inn = "poi", burnout = 500)
# sim_MINAR |> cor()
#
# x <- rpois(10000,2.5)
# y <- rpois(10000,5)
# z <- rpois(10000,1)
# RES <- cbind(x,y,z)
# A1 <- matrix(c(0.5,0,0,0,0.4,0,0,0,0.2),nrow=3,ncol=3)
# A2 <- matrix(c(0.2,0,0,0,0.3,0,0,0,0.4),nrow=3,ncol=3)
# A <- cbind(A1,A2)
# set.seed(1234)
# MINARp_gen_cpp(RES,A) |> cor()
# sim_MINAR <- genMINAR(10000, 3, list(A1,A2), c(2.5,5,1), inn = "poi", burnout = 500)
# sim_MINAR |> cor()



#' Yule-Walker for MINAR(p) alphas parameter estimation
#'
#' Experimental function
#'
#' @param X, matrix of m observed series
#' @param p, number of lags
#' @param inn, distribution of the innovation process
#' @param ..., additional parameters
#'
#' @importFrom stats acf
#' @importFrom stats var
#' @importFrom RcppML nnls
#'
#' @examples
#'
#' @export
estimMYW <- function(X, p, inn = "poi", ...) {
    stopifnot(inn %in% c("poi"))

    X <- as.matrix(X)
    n <- nrow(X)
    m <- ncol(X)

    # Media campionaria
    mu_hat <- colMeans(X)

    # Funzione per la covarianza campionaria a lag h
    Gamma_hat <- function(h) {
        X1 <- X[(h + 1):n, , drop = FALSE]
        X2 <- X[1:(n - h), , drop = FALSE]

        X1c <- sweep(X1, 2, mu_hat)
        X2c <- sweep(X2, 2, mu_hat)

        t(X1c) %*% X2c / (n - h)
    }

    # Stima delle matrici Gamma(0), ..., Gamma(p)
    Gamma_list <- lapply(0:p, Gamma_hat)

    # Costruzione di Gamma = [Gamma(1) ... Gamma(p)]
    Gamma_mat <- do.call(cbind, Gamma_list[2:(p + 1)])

    # Costruzione della matrice Toeplitz a blocchi
    G_block <- matrix(0, nrow = m * p, ncol = m * p)

    for (i in 1:p) {
        for (j in 1:p) {
            lag_ij <- abs(i - j)
            G_block[((i - 1) * m + 1):(i * m),
                    ((j - 1) * m + 1):(j * m)] <- Gamma_list[[lag_ij + 1]]
        }
    }

    # Stima Yule-Walker
    A_hat <- Gamma_mat %*% solve(G_block)

    # Separazione delle matrici A_1, ..., A_p
    A_list <- vector("list", p)
    names(A_list) <- paste0("A", 1:p)

    for (k in 1:p) {
        A_list[[k]] <- A_hat[, ((k - 1) * m + 1):(k * m)]
    }

    # Stima lambda (innovazioni Poisson)
    lambda_hat <- mu_hat
    for (k in 1:p) {
        lambda_hat <- lambda_hat - A_list[[k]] %*% mu_hat
    }

    OUT <- list("alphas" = A_list,
                "lambdas" = as.vector(lambda_hat),
                "meanX" =  mu_hat, "varX" = apply(X,2,var)
    )
    return(OUT)
}


# example on the use of estimMYW
# set.seed(123)
# X <- matrix(rpois(10000000*3, lambda = 5), ncol = 3)
# estimMYW(X, p = 1, inn = "poi")
