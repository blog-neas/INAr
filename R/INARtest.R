#' Perform Sun-McCabe INAR tests
#'
#' @param X vector, thinning parameters. The length of this vector defines the number of lags `p` of the INAR(p) process.
#' @param inn character, the innovation distribution to be used in the test. The options are "poi", "negbin" and "genpoi".
#' @param B integer, the number of bootstrap samples to be generated. The default value is 0, corresponding to no bootstrap samples.
#' @param method character, the type of test to be performed, the alternatives are "par" or "semipar".
#' @return A number.
#'
#' @details
#' The function performs the Sun-McCabe (SMC) tests for the INAR(1) model.
#' It is possible to run an exact (B = 0) or a bootstrap (B > 0) test.
#' SMC tests are based on Poisson, negative binomial and generalized Poisson
#' distributions, in all the considered cases, both parametric and semiparametric
#' methods are available.
#' It is possible to perform the exact tests using the original data (B=0) or
#' using bootstrap samples (B > 0).
#' @examples
#' # ....... examples .....
#' # SMCtest(X = rpois(100,2), type = "semiparametric", B = 0)
#'
#' @export
SMCtest <- function(X, inn = "poi", B = 0, method = NA){
    if(!all(X == as.integer(X))) X <- as.integer(X)
    if(!is.integer(B)) B <- as.integer(B)

    stopifnot(inn %in% info_inn$inn, B >= 0)
    stopifnot(is.na(method) | method %in% c("par","semipar"))
    if(B > 0 & is.na(method)) stop("Specify the method: 'par' or 'semipar'")
    stopifnot(all(X == as.integer(X)), !all(X == X[1]), is.integer(B))
    if(inn == "negbin" & var(X) <= mean(X)){
        stop("Underdispersion detected: Negative Binomial test not applicable (var <= mean).\n Try switching to generalized Poisson distribution.");
    }

    OUT <- get_info(list(
        test = "smc",
        inn = inn,
        method = method,
        B = B
    ))

    if(B > 0 & method == "par" & !OUT$check) stop(paste("The selected innovation distribution (",inn,") is not available for the parametric SMC bootstrap test",sep=""))
    OUT$data.name = deparse1(substitute(X)) # deparse(args(X))

    smc_est <- SMC_Cpp(X, OUT$inn_num)

    if(B > 0){
        # perform bootstrap test
        if(method == "par"){
            smc_boot <- SMC_parBOOT_Cpp(X, B, OUT$inn_num)
        }else if(method == "semipar"){
            smc_boot <- SMC_semiparBOOT_Cpp(X, B, OUT$inn_num)
        }else{
            stop("Specify a correct test type: 'par' or 'semipar'")
        }

        # add Boot results to OUT
        OUT$statistic <- c(T = mean(smc_boot))
        OUT$p.value <- mean(abs(smc_boot) > abs(smc_est[1]), na.rm = TRUE)
    }else{
        OUT$statistic <- c(T = smc_est[1])
        OUT$p.value <- smc_est[2]
    }
    return(OUT)
}

# example usage:
# set.seed(1926)
# x <- rpois(100000,2)
# SMCtest(x, inn = "katz", B = 0)
# SMCtest(x, inn = "negbin", B = 0, method = "semipar")


#' Perform the over(under)-Dispersion Index test for the INAR(1) (Weiss et al, 2019)
#'
#' @param X vector, a series of counts.
#' @return ........
#'
#' @details
#' The function performs the zero inflation test for the INAR(1) model.
#' @importFrom stats pnorm
#' @examples
#' # ....... examples .....
#' # DItest(X = rpois(100,2))
#'
#' @export
DItest <- function(X){ # , inn = "poi"
    if(!all(X == as.integer(X))) X <- as.integer(X)

    n <- length(X)
    est <- estimYW(X, p = 1)
    alpha <- abs(est$alphas)

    # sotto H0 ho sqrt(T)*(Id - 1)/s --> N(0,1)
    stderr <- sqrt(2 * (1 + alpha^2)/(1 - alpha^2))
    DIstat <- sqrt(n)*(est$varX / est$meanX - 1)/stderr

    pval <- 2 * (1 - pnorm(abs(DIstat)))

    OUT <- get_info(list(
        data.name = deparse1(substitute(X)),
        test = "di",
        statistic = c(T = DIstat),
        parameter = c("a1" = alpha),
        p.value = pval,
        method = NA,
        B = 0
    ))

    return(OUT)
}

# DItest(rpois(1000,2))
# DItest(genINAR(1000,a = 0.5, par = 2,arrival = "poisson")$X)
# DItest(genINAR(1000,a = 0.5, par = c(2,0.7),arrival = "negbin")$X)


#' Perform the zero inflation test for the INAR(1) (Weiss et al, 2019)
#'
#' @param X vector, a series.
#' @param type character, one of the possible zero inflation tests, ... ("pv") or ... ("vdb").
#' @return ........
#'
#' @details
#' The function performs the zero inflation test for the INAR(1) model.
#' @examples
#' # ....... examples .....
#' # ZItest(X = rpois(100,2))
#'
#' @export
ZItest <- function(X,type = "pv"){ # , inn = "poi"
    stopifnot(type %in% c("pv","vdb"))
    if(!all(X == as.integer(X))) X <- as.integer(X)
    # stopifnot(inn %in% info_inn$inn)

    n <- length(X)
    est <- estimYW(X, p = 1)
    alpha <- abs(est$alphas)
    p0 <- sum(X == 0) / n

    if(type == "pv"){
        ser1 <- series_varpv(est$meanX, alpha)
        se_pv <- sqrt((exp(est$meanX) - 1)/(est$meanX^2) - (1/est$meanX)*((1 + alpha)/(1 - alpha)) + (2/est$meanX^2)*ser1)
        ZIpv <- sqrt(n)*((1 + (log(p0)/est$meanX))/se_pv)
        stat <- ZIpv
        pval <- 2*(1 - pnorm(abs(ZIpv)))
    }else if(type == "vdb"){
        ser1 <- series_varpv(est$meanX, alpha)
        se_vdb <- sqrt(exp(est$meanX) - 1 - est$meanX*((1+alpha)/(1-alpha)) + 2*ser1)
        ZIvdb <- sqrt(n)*((p0*exp(est$meanX) - 1)/se_vdb)
        stat <- ZIvdb
        pval <- 2*(1 - pnorm(abs(ZIvdb)))
    }

    OUT <- get_info(list(
        data.name = deparse1(substitute(X)),
        test = paste("zi",type,sep="_"),
        statistic = c(T = stat),
        p.value = pval,
        parameter = c("a1" = alpha),
        inn = "none",
        method = NA,
        B = 0
    ))

    return(OUT)
}

# ZItest(rpois(1000,2),"pv")
# ZItest(rpois(1000,2),"vdb")
#
# ZItest(genINAR(1000,a = 0.5, par = 2,arrival = "poisson")$X,"pv")
# ZItest(genINAR(1000,a = 0.5, par = 2,arrival = "poisson")$X,"vdb")
# ZItest(genINAR(1000,a = 0.5, par = c(2,0.7),arrival = "negbin")$X,"pv")
# ZItest(genINAR(1000,a = 0.5, par = c(2,0.7),arrival = "negbin")$X,"vdb")


# #' Perform Harrison-McCabe INAR tests
# #'
# #' @param X vector, thinning parameters. The length of this vector defines the number of
# #' lags `p` of the INAR(p) process.
# #' @param inn character, the innovation distribution to be used in the test.
# #' The options are "poi", "negbin" and "genpoi".
# #' @param B integer, the number of bootstrap samples to be generated.
# #' The default value is 0, corresponding to no bootstrap samples.
# #' @param method character, the type of test to be performed,
# #' the alternatives are "par" or "semipar".
# #' @return A number.
# #'
# #' @details
# #' The function performs the Harrison-McCabe (HMC) tests for the INAR
# #' model. It is possible to run an exact (B = 0) or a bootstrap (B > 0) test.
# #' HMC tests are based on Poisson, negative binomial and generalized Poisson
# #' distributions, in all the considered cases, both parametric and semiparametric
# #' methods are available.
# #' It is possible to perform the exact tests using the original data (B=0)
# #' or using bootstrap samples (B > 0).
# #' @examples
# #' # ....... examples .....
# #' # HMCtest(X = rpois(100,2), type = "semiparametric", B = 0)
# #'
# #' @export
# #'
# HMCtest <- function(X, inn = "poi", B = 0, method = NA){
#     if(!all(X == as.integer(X))) X <- as.integer(X)
#     if(!is.integer(B)) B <- as.integer(B)
#
#     stopifnot(inn %in% info_inn$inn, B >= 0)
#     stopifnot(is.na(method) | method %in% c("par","semipar"))
#     if(B > 0 & is.na(method)) stop("Specify the method: 'par' or 'semipar'")
#     stopifnot(all(X == as.integer(X)), !all(X == X[1]), is.integer(B))
#     if(inn == "negbin" & var(X) <= mean(X)){
#         stop("Underdispersion detected: Negative Binomial test not applicable (var <= mean).\n Try switching to generalized Poisson distribution.");
#     }
#
#     OUT <- get_info(list(
#         test = "hmc",
#         inn = inn,
#         method = method,
#         B = B
#     ))
#
#     if(B > 0 & method == "par" & !OUT$check) stop(paste("The selected innovation distribution (",inn,") is not available for the parametric HMC bootstrap test",sep=""))
#     OUT$data.name = deparse1(substitute(X)) # deparse(args(X))
#
#     hmc_est <- HMC_Cpp(X, OUT$inn_num)
#
#     if(B > 0){
#         # perform bootstrap test
#         if(method == "par"){
#             hmc_boot <- HMC_parBOOT_Cpp(X, B, OUT$inn_num)
#         }else if(method == "semipar"){
#             hmc_boot <- HMC_semiparBOOT_Cpp(X, B, OUT$inn_num)
#         }else{
#             stop("Specify a correct test type: 'par' or 'semipar'")
#         }
#
#         # add Boot results to OUT
#         OUT$statistic <- c(T = mean(hmc_boot))
#         OUT$p.value <- mean(abs(hmc_boot) > abs(hmc_est[1]), na.rm = TRUE)
#     }else{
#         OUT$statistic <- c(T = hmc_est[1])
#         OUT$p.value <- hmc_est[2]
#     }
#     return(OUT)
# }
