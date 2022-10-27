#
#
# Sun McCabe test and Bootstrapped Sun McCabe test
#

#' Sun McCabe score statistic to test for dependence in an integer autoregressive process
#'
#' @param x data vector
#' @param order p, the order of the INAR(p) process
#' @param arrival distribution of the innovation process
#'
#' @returns
#' A list with class "htest" containing the following components:
#' * statistic the value of the Score test statistic.
#' * p.value the p-value for the test.
#' * null.value the specified hypothesized value of the thinning operator under the null hypothesis.
#' * alternative a character string describing the alternative hypothesis.
#' * method a character string indicating what type of Score test was performed.
#' * data.name a character string giving the assumption about the dgp of the data.
#' @details
#' This is a one-sided test, where the null hypothesis is \eqn{\alpha = 0} and under the alternative hypothesis the true thinning operator is greater than 0. The test statistics is computed according to the arrival distribution, see \insertCite{sun2013score;textual}{INAr} for additional details.
#' @references
#'   \insertAllCited{}
#' @export
SMC.test <- function(x, order=1, arrival = "poisson") {

    if(!(is.integer(x))) { stop("Error: Data should be integers"); }
    if(min(x) < 0) { stop("Error: Data should be non-negative"); }
    if(order!=1) { stop("Only INAR(1) process available at the moment"); }
    disp <- var(x)/mean(x)


    #Set description of test
    method      <- "Sun McCabe score test for INAR dependence.";

    #Set null and alternative hypotheses
    null.value  <- 0;
    attr(null.value, "names") <- "thinning operator";
    alternative <- "greater";

    #Calculate test statistics
    met <- 0
    if(tolower(arrival)=="poisson"){
        met <- 1
        #Set description of data
        data.name   <- paste0("INAR counts from Poisson INAR(",order,") dgp.")
    } else if(tolower(arrival)=="negbin"){
        met <- 2
        data.name   <- paste0("INAR counts from Negative Binomial INAR(",order,") dgp.")
        stopifnot('Data must be overdispersed when arrival = "negbin"' = disp > 1)
    }

    if(met==0) { stop('Wrong arrivals specified. Available options: "poisson", "negbin".'); }

    # Call C++ routine
    tval <- sunMC_Cpp(x, method = met)

    # questo non ci serve, qui volendo potremmo restituire gli alpha_hat
    # estimate    <- sum(M)/sum(N);
    # attr(estimate, "names") <- "proportion parameter";
    statistic   <- tval[1];
    attr(statistic, "names") <- "S";

    #Calculate p-value
    p.value     <- tval[2]
    attr(p.value, "names") <- NULL;

    #Create htest object
    TEST        <- list(method = method, data.name = data.name,
                        null.value = null.value, alternative = alternative, # estimate = estimate
                        statistic = statistic, p.value = p.value);
    class(TEST) <- "htest";
    TEST;
}

#
# SMC.test(1:100, arrival = "poisson")

#' Sun McCabe score statistic to test for dependence in an integer autoregressive process
#'
#' @param x data vector.
#' @param order p, the order of the INAR(p) process.
#' @param arrival distribution of the innovation process.
#' @param procedure the bootstrap procedure to use, either "semiparametric" or "parmetric" .
#' @param B number of Bootstrap replications.
#'
#' @returns
#' A list with class "htest" containing the following components:
#' * statistic the value of the bootstrap Score test statistic.
#' * p.value the p-value for the test.
#' * null.value the specified hypothesized value of the thinning operator under the null hypothesis.
#' * alternative a character string describing the alternative hypothesis.
#' * method a character string indicating what type of bootstrap Score test was performed.
#' * data.name a character string giving the assumption about the dgp of the data.
#' * estimate the value of the Score test statistic.
#' @details
#' This is a one-sided test, where the null hypothesis is \eqn{\alpha = 0} and under the alternative hypothesis the true thinning operator is greater than 0. The observed Score test statistic is computed according to the arrival distribution, see \insertCite{sun2013score;textual}{INAr} for additional details, while the bootstrap is performed according to \insertCite{palazzo2022semiparametric;textual}{INAr}.
#' @references
#'   \insertAllCited{}
#' @export
SMCboot.test <- function(x, order=1, arrival = "poisson", procedure, B = 499) {

    if(!(is.integer(x))) { stop("Error: Data should be integers"); }
    if(min(x) < 0) { stop("Error: Data should be non-negative"); }
    if(order!=1) { stop("Only INAR(1) process available at the moment"); }
    disp <- var(x)/mean(x)

    #Set description of test
    method      <- paste("Bootstrap Sun McCabe score test for INAR dependence,",tolower(procedure),"procedure. B =",B,"replications.");

    #Set null and alternative hypotheses
    null.value  <- 0;
    attr(null.value, "names") <- "thinning operator";
    alternative <- "greater";

    #Calculate test statistics
    met <- 0
    if(tolower(arrival)=="poisson"){
        met <- 1
        #Set description of data
        data.name   <- paste0("INAR counts from Poisson INAR(",order,") dgp.")
    } else if(tolower(arrival)=="negbin"){
        met <- 2
        data.name   <- paste0("INAR counts from Negative Binomial INAR(",order,") dgp.")
        stopifnot('Data must be overdispersed when arrival = "negbin"' = disp > 1)
    }

    if(met==0) { stop('Wrong arrivals specified. Available options: "poisson", "negbin".'); }

    # Call C++ routine
    tval <- sunMC_Cpp(x, method = met)
    Smc <- tval[1]
    estimate    <- Smc
    attr(estimate, "names") <- "obs. Score statistic";

    # BOOT Call C++
    if(tolower(procedure) == "parametric"){
        SmcB <- sunMC_parBOOT_Cpp(x, B, met)
    }
    else if(tolower(procedure) == "semiparametric"){
        SmcB <- sunMC_semiparBOOT_Cpp(x, B, met)
    } else{ stop('Wrong procedure chosen. Available options: "parametric", "semiparametric".') }

    statistic   <- mean(SmcB);
    attr(statistic, "names") <- "S_boot";

    #Calculate p-value
    p.value     <- mean(abs(SmcB) > abs(Smc))
    attr(p.value, "names") <- "bootstrap"

    #Create htest object
    TEST        <- list(method = method, data.name = data.name,
                        null.value = null.value, alternative = alternative,
                        estimate = estimate,
                        statistic = statistic, p.value = p.value)
    class(TEST) <- "htest"
    TEST
}
#
# SMCboot.test(1:100, arrival = "poisson", procedure = "semiparametric")

# https://stats.stackexchange.com/questions/441651/how-do-you-program-a-custom-hypothesis-test-in-r
#
# data:  COUNTS successes from TRIALS trials
# z = 2.5988, p-value = 0.009355
# alternative hypothesis: true dispersion parameter is greater than 0
# sample estimates:
#     proportion parameter
# 0.4359756
#
# TRIALS <- c(30, 32, 40, 28, 29, 35, 30, 34, 31, 39);
# COUNTS <- c( 9, 10, 22, 15,  8, 19, 16, 19, 15, 10);
#
# #Apply Tarone's test to the example data
# TEST <- Tarone.test(TRIALS, COUNTS);
# TEST;
#
# Tarone.test <- function(N, M) {
#
#     #Check validity of inputs
#     if(!(all(N == as.integer(N)))) { stop("Error: Number of trials should be integers"); }
#     if(min(N) < 1) { stop("Error: Number of trials should be positive"); }
#     if(!(all(M == as.integer(M)))) { stop("Error: Count values should be integers"); }
#     if(min(M) < 0) { stop("Error: Count values cannot be negative"); }
#     if(any(M > N)) { stop("Error: Observed count value exceeds number of trials"); }
#
#     #Set description of test and data
#     method      <- "Tarone's Z test";
#     data.name   <- paste0(deparse(substitute(M)), " successes from ",
#                           deparse(substitute(N)), " trials");
#
#     #Set null and alternative hypotheses
#     null.value  <- 0;
#     attr(null.value, "names") <- "dispersion parameter";
#     alternative <- "greater";
#
#     #Calculate test statistics
#     estimate    <- sum(M)/sum(N);
#     attr(estimate, "names") <- "proportion parameter";
#     S           <- ifelse(estimate == 1, sum(N),
#                           sum((M - N*estimate)^2/(estimate*(1 - estimate))));
#     statistic   <- (S - sum(N))/sqrt(2*sum(N*(N-1)));
#     attr(statistic, "names") <- "z";
#
#     #Calculate p-value
#     p.value     <- 2*pnorm(-abs(statistic), 0, 1);
#     attr(p.value, "names") <- NULL;
#
#     #Create htest object
#     TEST        <- list(method = method, data.name = data.name,
#                         null.value = null.value, alternative = alternative,
#                         estimate = estimate, statistic = statistic, p.value = p.value);
#     class(TEST) <- "htest";
#     TEST; }
