main_directory='/path/to/Projects/'
investigator='investigator'
project_title='{{Genesee.title}}'
authors = '{{Genesee.authors}}'
type = 'RNA'

devtools::load_all("path/to/GeneseeBulk")
devtools::load_all("path/to/Genesee")

## or
# library(GeneseeBulk)
# library(Genesee)

#### set up directories and create sample file
geneseesc_skeleton(main_directory, investigator, project_title, authors = authors, project_type = type)

rmarkdown::render('01_load.Rmd', params = list( input_tsv = "deSeq_counts.txt",
                                                input_project = NULL,
                                                sample_sheet_csv= "all_samples_file.csv",
                                                output_root= "refined/01load"),  output_dir = 'reports')

rmarkdown::render('02_qc_dea.Rmd',
                  params = list(input_root = 'refined/01load_se.rds',
                                se_subset = NULL,
                                design_csv = 'extradata/fe_design.csv',
                                sample_id = 'sample_id',
                                plot_covariates = c("group", "tp", "mouse_group"),
                                dea_method = 'deseq2',
                                cache = TRUE
                                output_root = "refined/02_qc_dea"
                  ), output_dir = 'reports',
                  output_file = '02_qc_dea.html'
)
