#' @name downloads
#' 
#' @title Download counts
#'
#' @description This data contain a time series expressing the daily number of downloads of a TEX editor for the period from June 2006 to February 2007 (T = 267)
#'
#' @docType data
#'
#' @usage data(downloads)
#'
#' @format An object of class `"tibble"`;
#'
#' @keywords datasets
#'
#' @references Weiß, C. H. (2018). An introduction to discrete-valued time series. John Wiley & Sons.
#' @references Weiß, C. H. (2008). Thinning operations for modeling time series of counts—a survey. AStA Advances in Statistical Analysis, 92(3), 319-341.
#'
#' @source Wiley repository, <www.wiley.com/go/weiss/discrete-valuedtimeseries>
#'
#' @examples
#' data(downloads)
#' \donttest{
#' plot(downloads$X,type="b")}
"downloads"
