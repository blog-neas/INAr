#' Generate INAR(p) models with different innovations
#'
#' @param n integer, number of observations, length of the series.
#' @param a vector, thinning parameters. The lenght of this vector defines the number of lags `p` of the INAR(p) process.
#' @param par vector, parameters related with the model, see the details section.
#' @param inn character, innovation distribution. Default value is `"poi"`, alternative values are `"negbin"` for Negative Binomial, `"genpoi"` for Generalized Poisson and `"katz"` for the Katz family.
#' @param burnout integer, number of starting observations to discard. Set to 500 by default.
# #' @param ... Additional arguments passed to the functions generating the random numbers.
#' @details
#' The function generates `n` observations drawn from an INAR(p) model
#' \deqn{X_t = \alpha_1 {*} X_{t-1} + \ldots + \alpha_p {*} X_{t-p} + \varepsilon_t}
#' where \eqn{*} is the binomial thinning operator and \eqn{\alpha_p} are the thinning parameters \insertCite{Steutel1979,al1987first}{INAr}.
#' The values of the thinning parameters \eqn{\alpha_i} are the elements of the `a` vector, for \eqn{i=1, \ldots, p}, while in `par` vector are specified the parameters related with the innovation distribution `\varepsilon`. The available processes for the innovations are listed below.
#'
#' **Poisson**. Innovations are drawn from a Poisson distribution \eqn{Poi(\lambda)}
#' \deqn{P(\varepsilon_t = k)= \dfrac{\lambda^k e^{-\lambda}}{k!},}
#' for \eqn{k \in \mathbb{N}}. The positive real number \eqn{\lambda} is equal to the expected value and also to the variance. The parameter is then equal to:
#' * `par` = c(\eqn{\lambda})
#'
#' **Negative Binomial**. Innovations are drawn from a Negative Binomial distribution \eqn{Poi(\lambda)}
#' \deqn{P(\varepsilon_t = k) = \dfrac{\Gamma(k+\gamma)}{\Gamma(\gamma)k!}p^k (1-p)^\gamma ,}
#' for \eqn{k \in \mathbb{N}}. The positive real number \eqn{\gamma} is the size and \eqn{p \in (0,1)} is the probability. The parameter is then equal to:
#' * `par` = c(\eqn{\gamma,p})
#'
#' **Generalized Poisson**. Innovations are \eqn{\varepsilon_t \sim Poi(\lambda)}
#' * `par` = ...
#'
#' **Katz**. Innovations are \eqn{\varepsilon_t \sim Poi(\lambda)}
#' * `par` = ...
#'
#' By the way, gatto gatto miao miao [genINAR()].
#' @references
#'   \insertAllCited{}
#' @return A number.
#' @importFrom stats rpois
#' @importFrom stats rgeom
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats rbinom
#' @importFrom stats rnbinom
#' @importFrom gamlss.dist rZIP
#' @importFrom skellam rskellam
#' @importFrom HMMpa rgenpois
#' @importFrom VGAM ryules
#' @importFrom tolerance rpoislind
#' @importFrom RNGforGPD GenUniGpois
#' @examples
#' # Generate 500 observations from an INAR(1) model with
#' # thinning parameter \eqn{\alpha}=0.5 and Poisson(2) innovations
#' set.seed(1234)
#' lam <- 2
#' genINAR(500, a = 0.5, par = lam, inn = "poi")
#'
#' @export
genINAR <- function(n, a, par, inn="poi", burnout=500){
    stopifnot(is.vector(a) ,all(a >= 0), sum(a) < 1)
    lags <- length(a)
    inn <- tolower(inn)
    s <- n + burnout
    resid_ <- rep(NA,s)

    # Check for stationarity
    mat <- matrix(0, nrow = lags, ncol = lags)
    mat[row(mat) - 1 == col(mat)] <- 1
    mat[1,] <- a
    ev <- eigen(mat)$values
    stopifnot(all(abs(ev) < 1))

    # Poisson marginal
    if(inn == "poi"){
        stopifnot(length(par) == 1)
        l_ <- unname(par)
        resid_ <- rpois(s, l_)
    }
    else if(inn == "overdisp_poisson"){
        stopifnot(length(par) == 2)

        l_ <- unname(par[1])
        d_ <- unname(par[2])
        stopifnot(d_ > 1)

        resid_ <- rnbinom(s, size=(l_/(d_-1)), mu=l_)
    }
    else if(inn == "zip"){
        stopifnot(length(par)==2)

        l_ <- unname(par[1])
        sig_ <- unname(par[2])

        # gamlss.dist::rZIP
        resid_ <- rZIP(s, mu = l_, sigma = sig_)
    }
    else if(inn == "bimodal_poisson"){
        stopifnot(length(par)==3)

        lam1_ <- unname(par[1])
        lam2_ <- unname(par[2])
        mixp_ <- unname(par[3])

        u_rand <- runif(s)

        sel <- u_rand < mixp_
        resid_[sel] <- rpois(sum(sel),lam1_)
        resid_[!sel] <- rpois(sum(!sel),lam2_)
    }
    else if(inn == "negbin"){
        stopifnot(length(par) == 2)

        # # VECCHIA PARAMETRIZZAZIONE, CON GAMMA E BETA
        # g_ <- unname(par[1]) # size, gamma
        # b_ <- unname(par[2]) # transf. prob, beta
        #
        # p.compl_ <- b_/(1+b_)
        # # dato che rnbinom prende le prob. invertite uso il complemento ad 1
        # # della vera formula, che sarebbe: p_ <- 1/(1+b_)
        # resid_ <- rnbinom(s,g_,p.compl_

        g_ <- unname(par[1]) # size, gamma
        p_ <- unname(par[2]) # prob successo

        p.compl_ <- 1-p_
        resid_ <- rnbinom(s,g_,p.compl_)

    }
    else if(inn == "geom"){
        stopifnot(length(par) == 2)

        pi_ <- unname(par[1])
        stopifnot(pi_ > 0, pi_ < 1)

        resid_ <- rgeom(s, pi_)
    }
    # else if(inn == "disc_unif"){
    #     stopifnot(length(par) == 1)
    #
    #     max_ <- unname(par[1])
    #
    #     resid_ <- sample(0:max_,s,replace = TRUE)
    # }
    else if(inn == "binomial"){
    stopifnot(length(par) == 2)

    enne_ <- unname(par[1])
    p_ <- unname(par[2])

    resid_ <- rbinom(s,enne_,p_)
    }
    else if(inn == "genpoi"){
        stopifnot(length(par) == 2)

        # lambda and kappa as in Sun and McCabe
        l_ <- unname(par[1])
        k_ <- unname(par[2])


        stopifnot(l_ > 0)
        stopifnot(k_ > -1, k_ < 1)

        # if(k_ < 0){
        #     m <- l_ + 4:1000*k_
        #     m <- ifelse(m <= 0, NA, m)
        #     m <- which.min(m) + 3 # perche' parto da 4
        #     lowb <- max(-1,-l_/m)
        #     stopifnot(k_ > lowb)
        # }

        # methods: five methods including Inversion, Branching, Normal-Approximation, Build-Up, and Chop-Down.
        # All five methods come from Demirtas (2017).
        # When lambda equals to 0, it is the ordinary Poisson distribution, so there is no need to specify the method.
        # "Branching" only works when lambda is positive. When theta is less than 10, the "Normal-Approximation" may not be reliable.

        # attenzione al metodo utilizzato!
        # HMMpa::rgenpois
        # RNGforGPD::GenUniGpois
        # resid_ <- GenUniGpois(theta = l_, lambda = k_,s, method = "Branching",...)$data
        # RNGforGPD::GenUniGpois
        resid_ <- GenUniGpois(l_, k_, s, method = ifelse(l_ > 10, "Normal-Approximation", "Inversion"), details = FALSE)$data
    }
    else if(inn=="katz"){
        stopifnot(length(par) == 2)
        a_ <- unname(par[1])
        b_ <- unname(par[2])

        resid_ <- rkatz(s, a_, b_)
  }
  else if(inn == "truncnorm"){
        mu_ <- unname(par[1])
        sig_ <- unname(par[2])

        vals <- round(rnorm(s, mu_, sig_),0)
        resid_ <- ifelse(vals < 0, 0, vals)
  }
  else if(inn == "truncskel"){
    lam1_ <- unname(par[1])
    lam2_ <- unname(par[2])

    # skellam::rskellam
    vals <- rskellam(s, lam1_, lam2_)
    resid_ <- ifelse(vals < 0, 0, vals)
  }
    # good package is not available anymore
    # else if(inn=="good"){
    #     z_ <- unname(par[1])
    #     s_ <- unname(par[2])
    #
    #     # good::rgood
    #     resid_ <- rgood(s, z_, s_)
    # }
    else if(inn=="yule"){
        rho_ <- unname(par[1])
        stopifnot(rho_ > 0)

        # VGAM::ryules
        resid_ <- ryules(s,rho_)
    }
    else if(inn=="zeta"){
        esse_ <- unname(par[1])
        stopifnot(esse_ > 0)

        # VGAM::rzeta
        resid_ <- VGAM::rzeta(s, esse_)
    }
    else if(inn=="poislind"){
        theta_ <- unname(par[1])
        stopifnot(theta_ > 0)

        # tolerance::rpoislind
        resid_ <- tolerance::rpoislind(s,theta_)
    }
    else if(inn=="mix_bin"){
        stopifnot(length(par)==5)

        enne1_ <- unname(par[1]) # size
        pb1_ <- unname(par[2]) # prob
        enne2_ <- unname(par[3]) # size, gamma
        pb2_ <- unname(par[4]) # prob successo
        mixp_ <- unname(par[5])

        selettore <- runif(s) < mixp_

        resid_[selettore] <-  rbinom(sum(selettore),enne1_,pb1_)
        resid_[!selettore] <- rbinom(sum(!selettore),enne2_,pb2_)
    }
    else if(inn=="mix_bin_negbin"){
        stopifnot(length(par)==5)

        enne_ <- unname(par[1]) # size
        pb_ <- unname(par[2]) # prob
        g_ <- unname(par[3]) # size, gamma
        pnb_ <- unname(par[4]) # prob successo
        mixp_ <- unname(par[5])

        selettore <- runif(s) < mixp_
        p.compl_ <- 1-pnb_

        resid_[selettore] <-  rbinom(sum(selettore),enne_,pb_)
        resid_[!selettore] <- rnbinom(sum(!selettore),g_,p.compl_)
    }
    else{
    stop("please specify one of the available distributions", call. = FALSE)
  }

  # generazione del campione
  a_ <- unname(a) # alpha

  ##### VERSIONE C++ #####
  sim_X <- INARp_cpp(resid_,a_)

  dati_sim <- data.frame(X=sim_X[(burnout+1):s],res=resid_[(burnout+1):s])
  return(dati_sim)
}

# veloce esempio --------------------------------------------------------
# N <- 500
# y <- genINAR(N,0.9,2,inn="poi")$X
# plot(y)
# N <- 500
# par <- c("a"=0.7,"lambda"=2)
# set.seed(1926)
# sim <- genINAR(N,par,inn="poi")$X
# plot(sim,type="b", main = "Poisson")
#
# #
# par <- c("a"=0.7,"lambda"=1.5,"disp"=2)
# set.seed(1926)
# sim <- genINAR(N,par,inn="overdisp_poisson")$X
# plot(sim,type="b", main = "Overdispersed Poisson")
#
# #
# par <- c("a"=0.7,"mu"=4,"sigma"=0.33)
# set.seed(1926)
# sim <- genINAR(N,par,inn="zip")$X
# plot(sim,type="b", main = "ZIP")
#
# #
# par <- c("a"=0.7,"lam1"=2,"lam2"=10,"mixprob"=0.66)
# set.seed(1926)
# sim <- genINAR(N,par,inn="bimodal_poisson")$X
# plot(sim,type="b", main = "Bimodal Poisson")
#
# #
# par <- c("a"=0.7,"g"=2,"b"=2/3)
# set.seed(1926)
# sim <- genINAR(N,par,inn="negbin")$X
# plot(sim,type="b", main = "Negative Binomial")
#
# #
# par <- c("a"=0.7,"min"=0,"max"=9)
# set.seed(1926)
# sim <- genINAR(N,par,inn="discunif")$X
# plot(sim,type="b", main = "Discrete Uniform")
#
# #
# par <- c("a"=0,"n"=5,"p"=0.4)
# set.seed(1926)
# sim <- genINAR(N,par,inn="binomial")$X
# plot(sim,type="b", main = "Binomial")
#
# #
# par <- c("a"=0.7,"lambda"=2,"kappa"=0.9)
# set.seed(1926)
# sim <- genINAR(N,par,inn="genpoi", details=FALSE)$X
# plot(sim,type="b", main = "Generalized Poisson")
#
# # sim <- genINAR("katz") # non ancora implememntato
#
# #
# par <- c("a"=0.7,"mu"=0,"sigma"=6)
# set.seed(1926)
# sim <- genINAR(N,par,inn="truncnorm")$X
# plot(sim,type="b", main = "Truncated Normal")
#
# #
# par <- c("a"=0.7,"lam1"=2,"lam2"=10)
# set.seed(1926)
# sim <- genINAR(N,par,inn="trunkskel")$X
# plot(sim,type="b", main = "Truncated Skellam")


# zzz <- skellam::rskellam(1000, 2,8)
#
# hist(zzz)
