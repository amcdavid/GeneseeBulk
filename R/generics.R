#' Extract results from a differential expression fit
#'
#' @param object an object containing differential expression results
#' @param ... optional parameters
#'
#' @return `tibble`
#' @export
gresults = function(object, ...){
  UseMethod("gresults")
}

#' @export
#' @describeIn gresults names of coefficients in [gresults()]
gresults_names = function(object){
  UseMethod("gresults_names")
}
