#' Innovations map
#'
#' @description
#' Internal data frame used to map innovations among back-end functions
#'
#' The data frame contains the following columns:
#' - `inn`: frontend string used to specify the innovations
#' - `inn_names`: full name of the innovations
#' - `inn_num`: numeric code used in the C++ back-end functions
#' - `smc`: logical code to check if the innovation is available for the SMC test
#' - `hmc`: logical code to check if the innovation is available for the HMC test
#'
#' @keywords internal
#' @noRd
info_inn <- data.frame(
    inn = c("poi", "negbin", "genpoi", "katz","none"),
    inn_name = c("Poisson", "Negative Binomial", "Generalized Poisson", "Katz","none"),
    inn_num = c(1,2,3,4,NA),
    smc = c(TRUE, TRUE, FALSE, FALSE,NA),
    hmc = c(FALSE, FALSE, FALSE, FALSE,NA)
)
# usethis::use_data(info_inn, internal = TRUE)


#' Generate test output list
#'
#' @description
#' Internal function to generate a list containing the information about the test
#'
#' @param LIST, list containing the parameters of the test. Its arguments are:
#' - `test`: type of test, either "smc" (Sun-McCabe) or "hmc" (Harrison-McCabe)
#' - `inn`: type of innovations, either "poi" (Poisson), "negbin" (Negative Binomial), "genpoi" (Generalized Poisson) or "katz" (Katz)
#' - `method`: if NA the test is not bootstrap, otherwise it is either "par" or "semipar"
#' - `B`: number of bootstrap replications
#'
#' @keywords internal
#' @noRd
get_info <- function(LIST){
    id_inn <- which(info_inn$inn == LIST$inn)
    innovation <- info_inn$inn_name[id_inn]
    method_name <- if(LIST$B > 0 & !is.na(LIST$method) & LIST$method %in% c("par","semipar")){ifelse(LIST$method == "par", "Parametric", "Semi-parametric")}else{""}

    if(LIST$test == "smc" | LIST$test == "hmc"){
        OUT <- list(
            method = paste(method_name,
                           ifelse(LIST$test == "smc", "Sun-McCabe", "Harrison-McCabe"),
                           ifelse(LIST$B > 0, "bootstrap", ""),
                           "test with", innovation, "innovations.",
                           ifelse(LIST$B > 0, paste(LIST$B, "bootstrap replications."), "")),
            inn_name = innovation,
            inn_num = info_inn$inn_num[id_inn],
            check = ifelse(LIST$test == "smc", info_inn$smc[id_inn], info_inn$hmc[id_inn])
        )
    }else if(LIST$test == "zi_pv"){
        OUT <- list(
            method = paste("Test for zero-inflation in INAR(1) models - PV version."), # with", innovation, "innovations."),
            inn_name = innovation,
            inn_num = info_inn$inn_num[id_inn],
            check = TRUE
        )
    }else if(LIST$test == "zi_vdb"){
        OUT <- list(
            method = paste("Test for zero-inflation in INAR(1) models - VDB version."), # with", innovation, "innovations."),
            inn_name = innovation,
            inn_num = info_inn$inn_num[id_inn],
            check = TRUE
        )
    }else if(LIST$test == "di"){
        OUT <- list(
            method = paste("Test for dispersion in INAR(1) models."), #  with", innovation, "innovations."),
            inn_name = innovation,
            inn_num = info_inn$inn_num[id_inn],
            check = TRUE
        )
    }else if(LIST$test == "zidi_pv"){
        OUT <- list(
            method = paste("Chi-squared test for joint zero-inflation and dispersion in INAR(1) models - PV version."), # with", innovation, "innovations."),
            inn_name = innovation,
            inn_num = info_inn$inn_num[id_inn],
            check = TRUE
        )
    }else if(LIST$test == "zidi_vdb"){
        OUT <- list(
            method = paste("Chi-squared test for joint zero-inflation and dispersion in INAR(1) models - VDB version."), # with", innovation, "innovations."),
            inn_name = innovation,
            inn_num = info_inn$inn_num[id_inn],
            check = TRUE
        )
    }

    if(!is.null(LIST$parameter)) OUT$parameter <- LIST$parameter
    if(is.null(LIST$statistic)){OUT$statistic <- c(T = NA)}else{OUT$statistic <- LIST$statistic}
    if(is.null(LIST$p.value)){OUT$p.value = NA}else{OUT$p.value <- LIST$p.value}
    if(is.null(LIST$data.name)){OUT$data.name = NA}else{OUT$data.name <- LIST$data.name}

    class(OUT) <- "htest"
    return(OUT)
}


#' Matrix of the first p lagged observations
#'
#' Internal function
#'
#' @param x, observed series
#' @param p, number of lags
#' @param cut, if TRUE saves the matrix with a smaller number of rows
#'
#' @keywords internal
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

# #' Partial values necessary to compute the variance of zero inflation estimator
# #'
# #' Internal function
# #'
# #' @param mu, parameter to compute the value
# #' @param alpha, parameter to compute the value
# #'
# #' @keywords internal
# #' @noRd
# loopseries <- function(mu, alpha) {
# 	go <- TRUE
# 	sumloop <- 0
# 	sumfin <- 0
# 	j <- 0
#
# 	while (go) {
# 		j <- j + 1
# 		sumfin <- sumfin + ((mu^j) / gamma(j + 1)) * (alpha^j - alpha^(j*2))
#
# 		if (abs(sumfin - sumloop) < 1e-9) go <- FALSE
#
# 		sumloop <- sumfin
# 	}
#
# 	return(sumfin)
# }


#' Parameter transformation and back-transformation
#'
#' @description
#' Internal function to transform and back-transform parameters
#'
#' @importFrom stats qlogis
#' @keywords internal
#' @noRd
par_transform <- function(a, par, inn){
    a_t <- qlogis(a)
    if(inn == "poi"){
        stopifnot(length(par) == 1)
        par_t <- log(par)
    }
    theta <- c(a_t,par_t)
    return(theta)
}

#' Parameter transformation and back-transformation
#'
#' @description
#' Internal function to transform and back-transform parameters
#'
#' @importFrom stats plogis
#' @keywords internal
#' @noRd
par_back <- function(theta, inn){
    n <- length(theta)
    if(inn == "poi"){
        m <- 1 # nr of innovation parameters
        a <- plogis(theta[1:(n-m)])
        attr(a,"names") <- paste0("a",1:(n-m))
        par <- exp(theta[(n-m+1):n])
    }
    return(list("alphas"=a,"par"=par))
}
