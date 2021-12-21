#' SummarizedExperiment with Nebula results
#'
#' A subclass that provides no additional slots,
#' but exists to implement polymorphic behavior of `tidy`, `results` and `resultsNames`
#' @seealso [run_nebula()]
#' @export
#' @importClassesFrom DESeq2 DESeqDataSet
setClass('NebulaSE', contains = 'DESeqDataSet')

