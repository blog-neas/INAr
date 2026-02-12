#' @name cryptosporidiosis
#'
#' @title Counts of cryptosporidiosis infect
#'
#' @description This data contain a time series expressing the weekly counts of new infections in Germany (period 2002–2008).
#'
#' @docType data
#'
#' @usage data(cryptosporidiosis)
#'
#' @format An object of class `"tibble"`;
#'
#' @keywords datasets
#'
#' @references Robert-Koch-Institut (2016)
#' @references Example 5.1.7, Weiß, C. H. (2018). An introduction to discrete-valued time series. John Wiley & Sons.
#'
#' @source SurvStat@@RKI 2.0 repository, <https://survstat.rki.de/>
#'
#' @examples
#' data(cryptosporidiosis)
#' \donttest{
#' plot(cryptosporidiosis$X,type="b")}
"cryptosporidiosis"
