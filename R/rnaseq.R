#' Run PCA and return a plot
#'
#' @param object `SummarizedExperiment`
#' @param ntop use this many variable genes
#' @param print_plot `logical`
#' @param aes_point a call to `aes` evaluated in the context of the colData for the `geom_point`, eg `aes(color = treatment)`.
#' @param aes_label a call to `aes` evaluated in the context of the colData for the `geom_textrepel`, eg `aes(label = sample_id)`.
#' @return list with elements `pca_dat`, containing `colData` plus PC1, PC2, `pca` output of `prcomp` and `plot`.
#' @import ggplot2
#' @autoglobal
#' @export
se_pca_plot = function(object, ntop = 500, aes_point, aes_label, print_plot = TRUE){
  rv <- MatrixGenerics::rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  pca_dat <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], colData(object))
  attr(pca_dat, "percentVar") <- percentVar[1:2]

  if(!missing(aes_point)){
    p = ggplot(pca_dat, aes(x = PC1, y = PC2)) + geom_point(aes_point)  + theme_minimal()
  } else{
    p = ggplot(pca_dat, aes(x = PC1, y = PC2)) + geom_point()  + theme_minimal()
  }
  if(!missing(aes_label)){
    p = p + ggrepel::geom_text_repel(aes_label, size = 2)
  }
  if(print_plot) print(p)

  invisible(list(pca_dat = pca_dat, pca = pca, plot = p))
}

#' @importFrom dplyr group_by mutate min_rank ungroup
#' @autoglobal
plotLoadings = function(pca, ntop = 10){
  pca_df = tibble::as_tibble(pca$rotation[,1:2], rownames = 'gene') %>% tidyr::gather(,,-gene)
  pca_df = pca_df %>% group_by(key, -sign(value)) %>% mutate(rank = min_rank(-abs(value))) %>% filter(rank < ntop)
  pca_df = pca_df %>% ungroup() %>% mutate(genef = factor(gene), genef = forcats::fct_reorder(genef, -value))
  ggplot(pca_df, aes(x = genef, y = value, fill = factor(`-sign(value)`))) + facet_wrap(~key, scales = 'free_x') + geom_col() + theme_minimal() + theme(axis.text.x = element_text(angle = 90), legend.position = 'none') + xlab('Gene')
}

#' @autoglobal
se_sparsepca_plot = function(data,
                             SPC_args = list(sumabsv = 4),
                             ggord_args = list(vec_ext = 1, ext = 1.2),
                             grp_in = 'treatment', color_by = NULL, label_by = NULL) {
  datx = t(assay(data))
  sdatx = scale(datx, scale = FALSE)
  ss = Genesee::call_intercalate_right(PMA::SPC, x = sdatx, K=2, orth = TRUE, extra = SPC_args)
  v = ss$v
  rownames(v) = make.names(colnames(sdatx))
  vs = v[rowSums(abs(v)) > 0, ]
  vs = vs + rnorm(length(vs), sd = sd(vs) / 4)
  uu  = as.data.frame(scale(ss$u))
  vs = scale(vs, center = FALSE)
  rownames(uu) = rownames(sdatx)
  ord1 = Genesee::call_intercalate_right(ggord,
    obs = uu,
    vecs = as.data.frame(vs),
    txt = 2.5,
    size = 1,
    arrow = NULL,
    grp_in = SummarizedExperiment::colData(data)[[grp_in]],
    extra = ggord_args)

  cu = cbind(ss$u, as.data.frame(colData(data)))
  ord2 = ggplot(cu, aes(x = `1`, y = `2`, color = !!color_by,label = !!label_by)) +
    geom_text() + theme_minimal()
  SummarizedExperiment::rowData(data) = cbind(SummarizedExperiment::rowData(data), as.data.frame(ss$v))
  invisible(list(ord1 = ord1, ord2 =ord2, data = data))
}

clamp = function(x, modulus = 5) {
  x[x < -modulus] = -modulus
  x[x > modulus] = modulus
  x
}

#' @autoglobal
modified_volc = function(dd, name, heatmap_genes = NULL, heatmap_top_group = 'group', heatmap_max_gene = Inf, sample_subset_idx = TRUE, title = '', subtitle = "", ...){
  if(!missing(name) && !any(name == gresults_names(dd)))
    stop("bad coefficient ", name, ". Options are ", paste0(gresults_names(dd), collapse = ', '), '.')
  res = gresults(dd, name = name, ...) %>% as.data.frame() %>% tibble::rownames_to_column('SYMBOL') %>%
    mutate(sign = sign(log2FoldChange)) %>% group_by(sign) %>%
    mutate(p_rank = rank(pvalue), label = ifelse(p_rank < 50, SYMBOL, ''))
  res[is.na(res$padj),'padj'] = 1
  trans = DESeq2::vst(dd)
  if(is.null(heatmap_genes)) heatmap_genes = filter(res, padj < .1, p_rank< heatmap_max_gene)$SYMBOL
  sub = trans[heatmap_genes,sample_subset_idx]
  if(nrow(sub)>0){
  SummarizedExperiment::assay(sub, 'zscore') = t(scale(t(assay(sub))))
  if(!('set' %in% names(SummarizedExperiment::rowData(sub)))){
    SummarizedExperiment::rowData(sub)$set = NA
  }
  catmat2 = SummarizedExperiment::rowData(sub)[, 'set', drop = FALSE]
  h = ComplexHeatmap::Heatmap(assay(sub, 'zscore'),
              top_annotation =  ComplexHeatmap::HeatmapAnnotation(df = as.data.frame(colData(sub)[, heatmap_top_group]), which = 'column'),
              name = 'Z-scored\nNormalized Exp', row_names_gp = grid::gpar(fontsize = 4), clustering_distance_rows = 'spearman', clustering_distance_columns = 'spearman', column_names_gp = grid::gpar(fontsize = 4))# + HeatmapAnnotation(catmat2, which = 'row')
  print(h)
  } else{
    message("No significant comparisons in ", title, "(", subtitle, ").")
  }
  ggplot(res, aes(x = log2FoldChange, y  = -clamp(log10(pvalue), 20), color = padj < .1)) + geom_point() + ggrepel::geom_text_repel(data = filter(res, p_rank < 25), mapping = aes(label = label), size = 2)  + geom_vline(lty =2, xintercept = 0) + theme_minimal() + scale_color_discrete('FDR', labels = c('TRUE' = '<10%', 'FALSE' = '>=10%')) + ggtitle(title, subtitle = subtitle)
}

deseq_within = function(dge_, stratify_by){
  strata = unique(colData(dge_)[[stratify_by]])
  res = purrr::map(strata, function(.x){
    dges = dge_[,colData(dge_)[[stratify_by]] == .x]
    DESeq2::DESeq(dges)
  })
  names(res) = strata
  res
}


#' @importFrom SummarizedExperiment colData rowData "colData<-" "rowData<-" mcols
subset_dge = function(subset_, dge_) {
  if(!missing(subset_) && !is.null(subset_) && !is.na(subset_)){
    cd = as.data.frame(SummarizedExperiment::colData(dge_))
    cd[['__idx__']] = seq_len(nrow(cd))
    subset_ = parse(text = subset_, n = 1)[[1]]
    cd = dplyr::filter(cd, !!subset_)
    dge_ = dge_[,cd[['__idx__']]]
    cdd = SummarizedExperiment::colData(dge_)
    for(j in seq_len(ncol(cdd))){
      if(is.factor(cdd[,j])) SummarizedExperiment::colData(dge_)[,j] = droplevels(cdd[,j])
    }
  }
  dge_
}

#' Run DESeq2 with a specified formula, on a specified subset of data
#' @inheritParams run_glmm_nb_design
#' @return [DESeq2::DESeqDataSet()]
#' @export
run_deseq_design = function(dge_, formula_, subset_){
  dge_ = subset_dge(subset_, dge_)
  DESeq2::design(dge_) = formula_
  DESeq2::DESeq(dge_)
}

gresults.DESeqDataSet = function(object, ...) {
  DESeq2::results(object,...)
}

gresults_names.DESeqDataSet = function(object){
  DESeq2::resultsNames(object)
}
