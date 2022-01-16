local({
  Genesee::create_exampleproject(skeleton_args = list(authors = 'you and me', project_type = 'RNA', investigator = 'alligator', project_title = 'bulk', navigate_rawdata = FALSE))
  cache_path = file.path(system.file('..', 'tests', 'testthat', 'known_cache', package = 'GeneseeBulk', mustWork = TRUE))

  test_that("Can create SummarizedExperiment",{
    rmarkdown::render('01_load.Rmd', params = list( input_tsv = "rawdata/deSeq2_subcounts.txt",
                                                    input_project = NULL,
                                                    sample_sheet_csv= "sample_sheet.csv",
                                                    output_root= "refined/01load"),  output_dir = 'reports', quiet = TRUE)
    expect_equal(nrow(dge), nrow(counts))
    saveRDS(dge, file.path(cache_path, '01load_se.rds'))
  })

})

