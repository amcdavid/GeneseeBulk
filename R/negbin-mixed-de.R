#' Differential Expression with Negative Binomial Mixed Models
#'
#' @param dge_ [DESeq2::DESeqDataSet()] object
#' @param assay_ which assay in `dge_` to use
#' @param formula_ [formula()] with a simple random intercept
#' @param subset_ an R expression, as a `character`, evaluated in the context of `colData(dge_)` to subset the columns (samples)
#' @param method `character`
#' @param ... passed to [lme4::glmer.nb()] or [nebula::nebula()]
#'
#' @return an object with a `tidy`, `results` and `resultNames` method
#' @importFrom stringr str_c fixed
#' @importFrom SummarizedExperiment assay colData
#' @importFrom tibble tibble
#'
#' @export
run_glmm_nb_design = function(dge_, formula_, assay_ = 'counts', subset_, method = c('nebula', 'glmer.nb'), ...){
  method = match.arg(method, c('nebula', 'glmer.nb'))
  dge_ = subset_dge(subset_, dge_)
  if(is.null(DESeq2::sizeFactors(dge_))){
    message("Estimating size factors.")
    dge_ = DESeq2::estimateSizeFactors(dge_)
  }
  n_gene = nrow(assay(dge_, assay_))
  mf = as.data.frame(colData(dge_))

  if(method == 'glmer.nb'){
    if(!requireNamespace('lme4') || !requireNamespace('broom.mixed')) stop("Please install lme4 and broom.mixed.")
    formula_ = update(y ~ . + offset(log(sizeFactor)), formula_)
    safe_mm = purrr::safely(run_one_glmm, otherwise = tibble(converged = FALSE))
    purrr::map_dfr(seq_len(n_gene), safe_mm,
                   mf = mf, formula_ = formula_, exprs = assay(dge_, assay_),
                   method = method, ...)
  } else if(method == 'nebula'){
    if(!requireNamespace('nebula')) stop("Please install nebula.")
    run_nebula(dge_, formula_, ...)
  } else{
    stop("Bad `method`")
  }

}

run_one_glmm = function(i, mf, formula_, exprs, method, ...){
  y = exprs[i,,drop = TRUE]
  environment(formula_) = environment()
  if(method == 'glmer.nb'){
    fit = lme4::glmer.nb(formula_, data = mf, ...)
  }
  df = tidy(fit, effect = 'fixed')
  df
}


#' Do a set of formulaes have a random effect?
#'
#' @param formulae vector that can be coerced to formula (e.g. `character()`)
#' @return `logical`
#' @export
#'
#' @examples
#' has_re(list(~ 1 + (1 |x), "~ treatment + (1|group)"))
has_re = function(formulae){
  purrr::map_lgl(formulae, function(x){
    f = as.formula(x)
    tt = try(extract_re_var(f), silent = TRUE)
    if(inherits(tt, 'try-error')) return(FALSE)
    return(TRUE)
  })

}



extract_re_var = function(formula_){
  termNames = labels(terms(formula_))
  termNames = labels(terms(formula_))
  if(any(bad <- stringr::str_detect(termNames, '[^\\w|: ]'))) {
    stop("Only alphanumeric variable names are supported in formulae currently.  Illegal terms: ",
         paste(termNames[bad], collapse = ', '), '.')
  }
  hasRE = stringr::str_detect(termNames, stringr::fixed('|'))

  ## collapse all variables into something that can be used for model.frame
  if(sum(hasRE) != 1) stop("Must specify one (and only one) random effect.")
  re_id = stringr::str_match(termNames[hasRE], '^1\\s?\\|[\\s]*([\\w]+)\\s*$')[,2]
  if(is.na(re_id)) stop("Couldn't parse ", termNames[hasRE], '.  Only random intercepts and alphanumeric variable names are supported.')
  ## save portion of formula that contained random effects
  fe_form =  paste(termNames[!hasRE], collapse='+')
  if(stringr::str_trim(fe_form)=='') fe_form = '1'

  ## re_id: random effect variable name
  ## fe_form: the actual formula specifying the fixed effects
  ## All are character vectors of length 1.
  list(re_id = re_id, fe_form = str_c("~ ", fe_form))
}

drop_na = function(x) x[!is.na(x)]

#' Run a negative binomial mixed model regression using Nebula
#'
#' @param dge_ DESeq2 object
#' @param formula_ formula
#' @param ... passed to [nebula::nebula()]
#'
#' @return NebulaSE-class
#' @autoglobal
#' @export
run_nebula = function(dge_, formula_, ...){
  parse_form = extract_re_var(formula_)
  dge_ = dge_[,order(colData(dge_)[[parse_form$re_id]])]
  id = colData(dge_)[[parse_form$re_id]]
  mm = model.matrix(as.formula(parse_form$fe_form), data = colData(dge_))
  counts = assay(dge_, 'counts')
  offset = DESeq2::sizeFactors(dge_)
  fit = Genesee::call_intercalate_left(nebula::nebula, ..., count = counts, id = id, pred = mm, offset = offset, extra = list(covariance = TRUE))
  # somehow lose rowData when we coerce -- not sure why
  dge_ = as(dge_, 'NebulaSE')
  stopifnot(nrow(fit$summary) == nrow(dge_))
  # may want to index with these
  fit$summary = fit$summary %>% dplyr::select(-gene, -gene_id)
  results = S4Vectors::DataFrame(fit$summary,
                      'names<-'(fit$overdispersion, str_c('overdisp_', names(fit$overdispersion))),
                      convergence = unname(fit$convergence),
                      fit$covariance)
  snames = names(results)
  terms = stringr::str_match(snames, '^(logFC|se|p)_(.+?)$')
  colnames(terms) = c('name', 'quantity', 'term')
  SummarizedExperiment::mcols(results) = cbind(S4Vectors::DataFrame(terms),
                         type = rep(c('beta', 'overdispersion', 'convergence', 'covariance'),
                                    times = c(ncol(fit$summary), ncol(fit$overdispersion), 1, ncol(fit$covariance))))
  SummarizedExperiment::rowData(dge_) = cbind(SummarizedExperiment::rowData(dge_), results)
  S4Vectors::metadata(dge_) = c(S4Vectors::metadata(dge_),
                                   list(resultsNames = drop_na(unique(terms[,'term']))))
  dge_
}
