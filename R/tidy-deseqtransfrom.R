unrowname = function(x) {
  rownames(x) = NULL
  x
}

biobroom_finish = function (x) {
  x = unrowname(x)
  opt = getOption("biobroom.return", default = "tbl_df")
  if (opt == "tbl_df") {
    tibble::as_tibble(x)
  }
  else if (opt == "tbl_dt") {
    tibble::as_tibble(x)
  }
  else if (opt == "data.table") {
    data.table::as.data.table(tibble::as_tibble(x))
  }
  else if (opt == "data.frame") {
    as.data.frame(x)
  }
  else {
    stop(paste("Invalid biobroom.return format", opt))
  }
}



#' Tidying methods for DESeq2 DESeqTransform objects
#'
#' This reshapes a DESeq2 transform object (output from `rlog()` and
#' `varianceStabilizingTransform()` functions) into a tidy format.
#'
#' @param x DESeqTransform object
#' @param colData whether colData should be included in the tidied output
#' for those in the DESeqTransform object.
#' @param ... extra arguments (not used)
#'
#' @details \code{colDat=TRUE} adds covariates from colData to the data frame.
#' @note Unclear what the license for this code might be.  It's copied directly from https://gist.github.com/tavareshugo/3973461a7daf8a43e65e3566d5deed14
#'
#' @return the result is a data frame with the columns
#'   \item{gene}{gene ID}
#'   \item{sample}{sample ID}
#'   \item{log2Count}{transformed normalized counts in this gene in this sample on a log2 scale}
#'
#' If \code{colData = TRUE}, it also merges this with the columns present
#' in \code{colData(x)}.
#'
#' @name DESeq2_tidiers
#' @importClassesFrom  DESeq2 DESeqTransform
#'
#' @examples
#'
#' # From DESeq2 documentation
#'
#' if (require("DESeq2")) {
#'     dds = makeExampleDESeqDataSet(betaSD = 1)
#'
#'     # vst transform
#'     vst_norm = varianceStabilizingTransformation(dds)
#'
#'     # rlog transform
#'     rlog_norm = rlog(dds)
#'
#'     tidy(rlog_norm)
#'     tidy(vst_norm)
#'
#'     # With design included
#'     tidy(rlog_norm, colData=TRUE)
#'     tidy(vst_norm, colData=TRUE)
#' }
#'
#' @method tidy DESeqTransform
#' @importFrom dplyr "%>%"
#' @export
#' @author Hugo Tavares https://gist.github.com/tavareshugo/3973461a7daf8a43e65e3566d5deed14
tidy.DESeqTransform = function(x, colData = FALSE, ...){
  ellipsis::check_dots_empty()
  if(!requireNamespace(biobroom)){
    stop("please install `biobroom` package")
  }

  # Tidy the normalised expression data within the object
  expressions = tibble::as_tibble(assay(x), rownames = "gene")

  ret = expressions %>%
    tidyr::gather(sample, log2Count, -gene) %>%
    dplyr::mutate(sample = as.character(sample))

  if (colData) {
    cdat = data.frame(SummarizedExperiment::colData(x), stringsAsFactors=FALSE)
    rownames(cdat) = colnames(x)
    ret = cbind(
      ret[, c('gene', 'sample')],
      cdat[ret$sample,,drop=FALSE],
      log2Count = ret$log2Count) %>%
      unrowname()
  }
  biobroom_finish(ret)
}
