#' Testing and analysis of the CIMseq and method.
#'
#' Description
#'
#' \tabular{ll}{ Package: \tab CIMseq\cr Type: \tab Package\cr
#' Version: \tab 1.0\cr Date: \tab 2016-02-28\cr License: \tab GPL-3\cr }
#'
#' @name CIMseq.testing-package
#' @aliases CIMseq.testing-package CIMseq.testing
#' @docType package
#' @author Author: Jason T. Serviss
#' @keywords package
#'
#' @import methods
#' @import ggplot2
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
