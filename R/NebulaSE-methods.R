#' @name coerce
#' @title coerce
#' @aliases coerce,DESeqDataSet,NebulaSE-method
#' @rdname NebulaSE-class
setAs('DESeqDataSet', 'NebulaSE', function(from) {
  value = new('NebulaSE')
  for (what in slotNames(value)) {
    slot(value, what) = slot(from, what)
  }
  value
})

.nebula_terms = function(object){
  stopifnot(inherits(object, 'NebulaSE'))
  S4Vectors::metadata(object)$resultsNames
}

#' Tidying methods for ad-hoc `g_nebula` class
#'
#' @param x output from [run_nebula()]
#' @param ... not used
#'
#' @return `tibble` with columns matching [biobroom::tidy.DESeqResults()]
#' @export
#' @importFrom generics tidy
#' @importFrom dplyr rename group_by mutate ungroup
#' @method tidy NebulaSE
#' @autoglobal
#' @describeIn NebulaSE-class return tidy results
tidy.NebulaSE = function(x, ...){
  ellipsis::check_dots_empty()
  rd_info = mcols(mcols(x))
  betas = mcols(x)[rd_info$type %in% c('beta', 'convergence')] %>%
    as.data.frame()
  betas[['idx']] = seq_len(nrow(betas))
  fit_summary =  betas %>% tidyr::pivot_longer(c(-idx, -convergence),
                                             names_to = 'name')
  terms = rd_info[which(rd_info$type == 'beta'),c('name', 'term', 'quantity')] %>% as.data.frame()
  fit_summary = left_join(fit_summary, terms, by = 'name')
  out = fit_summary %>%
    tidyr::pivot_wider(c(idx, term, convergence), names_from = quantity) %>%
    rename(estimate = logFC, stderror = se, p.value = p) %>%
    group_by(term) %>%
    mutate(p.adjusted = p.adjust(p.value, method = 'fdr')) %>%
    ungroup()
  add_gene_name(out, x)
}

#' @importFrom dplyr left_join
#' @autoglobal
add_gene_name = function(results, obj){
  gene_tab = tibble(gene = rownames(obj), idx = seq_len(nrow(obj)))
  left_join(results, gene_tab, by = 'idx') %>% dplyr::select(-idx)
}

#' @export
#' @describeIn NebulaSE-class S4 tidy
setMethod('tidy', signature = c(x = 'NebulaSE'), tidy.NebulaSE)


#' @export
#' @describeIn gresults-g_nebula return the names of the model coefficients.
#' @method gresults_names NebulaSE
gresults_names.NebulaSE = .nebula_terms

#' Return DESeq2-style `results`
#'
#' One and only one of `contrast` or `name` must be specified.
#' @param object output from [run_nebula()]
#' @param contrast `numeric`, giving linear combination of coefficients
#' @param name a name of a term.
#' @param ... unused
#' @return `tibble`
#' @importFrom stats as.formula model.matrix na.omit p.adjust pchisq prcomp rnorm sd terms update
#' @importFrom methods as new slot "slot<-" slotNames
#' @importFrom dplyr filter
#' @rdname gresults-g_nebula
#' @autoglobal
#' @export
gresults.NebulaSE = function(object, contrast, name, ...) {
  ellipsis::check_dots_empty()
  if (!missing(contrast)) {
    if (!is.numeric(contrast) || length(contrast) != length(gresults_names(object))) {
      stop("`contrast` must be numeric with same length `gresults_names(object)`.")
    }
    result = do_nebula_contrasts(object, contrast)
  } else if (!missing(name)) {
    result = tidy(object) %>% filter(.data[['term']] == name)
  } else {
    stop("Must specify only one of `contrast` or `name`.")
  }
  result = result %>% dplyr::rename(log2FoldChange = estimate,
                         pvalue = p.value,
                         padj = p.adjusted) %>%
    as.data.frame()
  rownames(result) = result[['gene']]
  return(result)
}

#' @autoglobal
do_nebula_contrasts = function(object, contrast){
  # genes x p
  betas = as.matrix(object$summary[which(object$terms$quantity == 'logFC')])
  # p x 1
  dim(contrast) = c(ncol(betas), 1)
  # genes x 1
  beta_rot = betas %*% contrast
  if(is.null(cov_df <- object$covariance)) stop("Set `covariance = TRUE` when running nebula.")
  sigma2_rot = rep(NA_real_, nrow(cov_df))
  for(i in seq_along(sigma2_rot)) {
    this_sigma = new("dspMatrix", uplo = "L", x = as.numeric(cov_df[i,]),
                     Dim = c(ncol(betas), ncol(betas)), Dimnames = list(NULL,NULL),
                     factors = list())
    sigma2_rot[i] = t(contrast) %*% this_sigma %*% contrast
  }
  return(tibble(gene_id = seq_along(sigma2_rot),
                gene = object$summary[['gene']],
                term = 'specified contrast',
                convergence = object$convergence,
                estimate = as.numeric(beta_rot),
                stderror = sqrt(sigma2_rot),
                p.value = 1 - pchisq(estimate^2/sigma2_rot, df = 1),
                p.adjusted = p.adjust(p.value, method = 'fdr')
  ))
}

