# require(gamlss.dist)
# require(RNGforGPD)
# require(skellam)
# require(good)

# Rcpp::sourceCpp("src/INAR1_gen.cpp")

#' Generate INAR(1) models with different innovations
#'
#' @param n integer, number of observations, length of the series.
#' @param par vector, parameters related with the model, see the details section.
#' @param arrival character, arrival distribution. Default value is `"poisson"`, alternative values are `"negbin"` for Negative Binomial, `"genpoi"` for Generalized Poisson and `"katz"` for the Katz family.
#' @param burnout integer, number of starting observations to discard, set to 500 by default.
#' @details
#' The function generates observations drawn from an INAR(1) model
#' \deqn{X_t = \alpha {*} X_{t-1} + \varepsilon_t}
#' where \eqn{*} is the binomial thinning operator and \eqn{\alpha} is the thinning parameter \insertCite{Steutel1979,al1987first}{INAr}.
#' The first value of the `par` vector is the thinning operator \eqn{\alpha}, while the remaining depend mainly on the type of arrival distribution. The arrivals already implemented are listed below.
#'
#' **Poisson**. Innovations are drawn from a Poisson distribution \eqn{Poi(\lambda)}
#' \deqn{P(\varepsilon_t = k)= \dfrac{\lambda^k e^{-\lambda}}{k!},}
#' The positive real number \eqn{\lambda} is equal to the expected value and also to the variance.The parameter is then equal to:
#' * `par` = c(\eqn{\alpha},\eqn{\lambda})
#'
#' **Negative Binomial**. Innovations are drawn from a Poisson distribution \eqn{Poi(\lambda)}
#' * `par` = ...
#'
#' **Generalized Poisson**. Innovations are \eqn{ \varepsilon_t \sim Poi(\lambda)}
#' * `par` = ...
#'
#' **Katz**. Innovations are \eqn{ \varepsilon_t \sim Poi(\lambda)}
#' * `par` = ...
#'
#' By the way, gatto gatto miao miao [genINAR()].
#' @references
#'   \insertAllCited{}
#' @return A number.
#' @examples
#' # Generate 500 an INAR(1) model with thinning \eqn{\alpha}=0.5 and Poisson(2) innovations
#' set.seed(1234)
#' pars <- c(0.5,2)
#' genINAR(500, pars, arrival = "poisson")
#'
#' @export
genINAR <- function(n,a,par,arrival="poisson",burnout=500){
    # n <- 500
    # a <- c(0.2,0.9,0.8)
    # par <- 2
    # arrival="poisson"
    # burnout=500
    stopifnot(is.vector(a))
    stopifnot(all(a >= 0))
    lags <- length(a)
    arrival <- tolower(arrival)

    s <- n + burnout
    mat <- matrix(0, nrow = lags, ncol = lags)
    mat[row(mat) - 1 == col(mat)] <- 1
    mat[1,] <- a
    ev <- eigen(mat)$values
    stopifnot(all(abs(ev) < 1))

    # -9999 = valori NA che saranno sovrascritti!
    sim_X <- rep(-99,s)

    # Poisson marginal
    if(arrival=="poisson"){
        stopifnot(length(par)==1)
        l_ <- unname(par)
        resid_ <- rpois(s,l_)

    }
    else if(arrival=="overdisp_poisson"){
        stopifnot(length(par)==2)

        l_ <- unname(par[1])
        d_ <- unname(par[2])
        stopifnot(d_ > 1)

        resid_ <- rnbinom(s, size=(l_/(d_-1)), mu=l_)
    }
    else if(arrival=="zip"){
    stopifnot(length(par)==2)

    l_ <- unname(par[1])
    sig_ <- unname(par[2])

    resid_ <- gamlss.dist::rZIP(s, mu = l_, sigma = sig_)
  }
  else if(arrival=="bimodal_poisson"){
    stopifnot(length(par)==3)

    lam1_ <- unname(par[1])
    lam2_ <- unname(par[2])
    mixp_ <- unname(par[3])

    u_rand <- runif(s)

    selettore <- u_rand < mixp_
    resid_ <- rep(NA,s)
    resid_[selettore] <- rpois(sum(selettore),lam1_)
    resid_[!selettore] <- rpois(sum(!selettore),lam2_)
  }
  else if(arrival=="negbin"){
    stopifnot(length(par)==2)
    g_ <- unname(par[1]) # size, gamma
    b_ <- unname(par[2]) # transf. prob, beta

    p.compl_ <- b_/(1+b_)
    # dato che rnbinom prende le prob. invertite uso il complemento ad 1
    # della vera formula, che sarebbe: p_ <- 1/(1+b_)

    resid_ <- rnbinom(s,g_,p.compl_)
  }
  else if(arrival=="discunif"){
    stopifnot(length(par)==2)

    minu_ <- unname(par[1])
    maxu_ <- unname(par[2])
    stopifnot(minu_ < maxu_)

    resid_ <- round(runif(s,minu_,maxu_),0)
  }
  else if(arrival=="binomial"){
    stopifnot(length(par)==2)

    enne_ <- unname(par[1]) # size
    p_ <- unname(par[2]) # prob

    resid_ <- rbinom(s,enne_,p_)
  }
  else if(arrival=="genpoi"){
    stopifnot(length(par)==2)

    l_ <- unname(par[1]) # lambda, come in sunmccabe
    k_ <- unname(par[2]) # kappa, come in sunmccabe

    # methods: five methods including Inversion, Branching, Normal-Approximation, Build-Up, and Chop-Down.
    # All five methods come from Demirtas (2017).
    # When lambda equals to 0, it is the ordinary Poisson distribution, so there is no need to specify the method.
    # "Branching" only works when lambda is positive. When theta is less than 10, the "Normal-Approximation" may not be reliable.

    # attenzione al metodo utilizzato!
    resid_ <- RNGforGPD::GenUniGpois(theta = l_, lambda = k_,s, method = "Branching",...)$data
  }
  else if(arrival=="katz"){
    #
    # come fare???
    # resid_ <-
    #
  }
  else if(arrival=="truncnorm"){
    mu_ <- unname(par[1])
    sig_ <- unname(par[2])

    vals <- round(rnorm(s,mu_,sig_),0)
    resid_ <- ifelse(vals < 0,0,vals)
  }
  else if(arrival=="truncskel"){
    lam1_ <- unname(par[1])
    lam2_ <- unname(par[2])

    vals <- skellam::rskellam(s, lam1_, lam2_)
    resid_ <- ifelse(vals < 0,0,vals)
  }
  else if(arrival=="good"){
    z_ <- unname(par[1])
    s_ <- unname(par[2])

    resid_ <- good::rgood(s, z_, s_)
  }
  else{
    stop("please specify one of the available distributions")
  }

  # generazione del campione
  a_ <- unname(a) # alpha

  ##### VERSIONE C++ #####
  # sim_X <- INAR1_cpp(resid_,a_) # per INAR1
  sim_X <- INARp_cpp(resid_,a_)

  dati_sim <- data.frame(X=sim_X[(burnout+1):s],res=resid_[(burnout+1):s])
  return(dati_sim)
}

# zzz <- rbinom(1000,40,0.7)
# asd <- function(A,B){
#   OUT <- A
#   for(i in 2:length(A)){
#     OUT[i] <- rbinom(1,OUT[i-1],B) + A[i]
#   }
#   return(OUT)
# }
#
# microbenchmark::microbenchmark(INAR1_cpp(zzz,0.7),asd(zzz,0.7),times = 100L)


# veloce esempio --------------------------------------------------------
# N <- 500
# y <- genINAR(N,0.9,2,arrival="poisson")$X
# plot(y)
# N <- 500
# par <- c("a"=0.7,"lambda"=2)
# set.seed(1926)
# sim <- genINAR(N,par,arrival="poisson")$X
# plot(sim,type="b", main = "Poisson")
#
# #
# par <- c("a"=0.7,"lambda"=1.5,"disp"=2)
# set.seed(1926)
# sim <- genINAR(N,par,arrival="overdisp_poisson")$X
# plot(sim,type="b", main = "Overdispersed Poisson")
#
# #
# par <- c("a"=0.7,"mu"=4,"sigma"=0.33)
# set.seed(1926)
# sim <- genINAR(N,par,arrival="zip")$X
# plot(sim,type="b", main = "ZIP")
#
# #
# par <- c("a"=0.7,"lam1"=2,"lam2"=10,"mixprob"=0.66)
# set.seed(1926)
# sim <- genINAR(N,par,arrival="bimodal_poisson")$X
# plot(sim,type="b", main = "Bimodal Poisson")
#
# #
# par <- c("a"=0.7,"g"=2,"b"=2/3)
# set.seed(1926)
# sim <- genINAR(N,par,arrival="negbin")$X
# plot(sim,type="b", main = "Negative Binomial")
#
# #
# par <- c("a"=0.7,"min"=0,"max"=9)
# set.seed(1926)
# sim <- genINAR(N,par,arrival="discunif")$X
# plot(sim,type="b", main = "Discrete Uniform")
#
# #
# par <- c("a"=0,"n"=5,"p"=0.4)
# set.seed(1926)
# sim <- genINAR(N,par,arrival="binomial")$X
# plot(sim,type="b", main = "Binomial")
#
# #
# par <- c("a"=0.7,"lambda"=2,"kappa"=0.9)
# set.seed(1926)
# sim <- genINAR(N,par,arrival="genpoi", details=FALSE)$X
# plot(sim,type="b", main = "Generalized Poisson")
#
# # sim <- genINAR("katz") # non ancora implememntato
#
# #
# par <- c("a"=0.7,"mu"=0,"sigma"=6)
# set.seed(1926)
# sim <- genINAR(N,par,arrival="truncnorm")$X
# plot(sim,type="b", main = "Truncated Normal")
#
# #
# par <- c("a"=0.7,"lam1"=2,"lam2"=10)
# set.seed(1926)
# sim <- genINAR(N,par,arrival="trunkskel")$X
# plot(sim,type="b", main = "Truncated Skellam")


# zzz <- skellam::rskellam(1000, 2,8)
#
# hist(zzz)
