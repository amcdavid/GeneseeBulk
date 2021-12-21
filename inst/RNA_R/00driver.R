main_directory='/path/to/Projects/'
investigator='investigator'
project_title='{{Genesee.title}}'
authors = '{{Genesee.authors}}'
type = 'RNA'

devtools::load_all("path/to/GeneseeSC")
## or
# library(GeneseeSC)

#### set up directories and create sample file
geneseesc_skeleton(main_directory, investigator, project_title, authors = authors, project_type = type)

rmarkdown::render('01_load.Rmd', params = list( input_tsv = "deSeq_counts.txt",
                                                input_project = NULL,
                                                sample_sheet_csv= "all_samples_file.csv",
                                                output_root= "refined/01load"),  output_dir = 'reports')

rmarkdown::render('02_qc_dea.Rmd',
                  params = list(input_root = 'refined/01load_se.rds',
                                se_subset = NULL,
                                design_csv = 'extradata/deseq_design.csv',
                                sample_id = 'sample_id',
                                plot_covariates = c(NULL)
                  ), output_dir = 'reports')
