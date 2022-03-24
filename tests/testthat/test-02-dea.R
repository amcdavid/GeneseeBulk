local({
  Genesee::create_exampleproject(skeleton_args = list(authors = 'you and me', project_type = 'RNA', investigator = 'alligator', project_title = 'bulk', navigate_rawdata = FALSE))
  cache_path = file.path(system.file('..', 'tests', 'testthat', 'known_cache', package = 'GeneseeBulk', mustWork = TRUE))

  test_that("Can run DEA with fixed effects",{
    rmarkdown::render('02_qc_dea.Rmd',
                      params = list(input_root = file.path(cache_path, '01load_se.rds'),
                                    se_subset = NULL,
                                    design_csv = 'extradata/fe_design.csv',
                                    terms_annotation_csv = 'extradata/fe_design_terms.csv',
                                    sample_id = 'sample_id',
                                    plot_covariates = c('group', 'IFN_beta'),
                                    organism = 'human',
                                    gsea_level = 0.05,
                                    compareCluster_args = list(pvalueCutoff = 1.0, qvalueCutoff = 1.0)
                      ), output_dir = 'reports')
  })

})

