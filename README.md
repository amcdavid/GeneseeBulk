
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GeneseeBulk

<!-- badges: start -->
<!-- badges: end -->

The goal of GeneseeBulk is to facilitate analysis (differential
expression) of Bulk RNASeq.

## Installation

You can install the development version of GeneseeBulk from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("amcdavid/GeneseeBulk")
```

## Example

The Genesee\* family of packages is organized around project templates
that copy over various parametrized markdown scripts.  
As such, there’s some code required to set overall parameters for the
project.

``` r
library(GeneseeBulk)
dir.create(tmp <- tempfile())
proj_directory = Genesee::genesee_skeleton(genesee_root = tmp,
                          investigator = 'Plato', 
                          project_title = 'Crito',
                          project_type = 'RNA',
                          init_git = FALSE)
#> Creating /var/folders/33/vpk59z855zzbrczltjnc52qh0000gn/T//RtmpR1YZSG/file7f0a2c93dfb4/Plato
#> Creating /var/folders/33/vpk59z855zzbrczltjnc52qh0000gn/T//RtmpR1YZSG/file7f0a2c93dfb4/Plato/Crito
#> Creating /var/folders/33/vpk59z855zzbrczltjnc52qh0000gn/T//RtmpR1YZSG/file7f0a2c93dfb4/Plato/Crito/rawdata
#> Creating /var/folders/33/vpk59z855zzbrczltjnc52qh0000gn/T//RtmpR1YZSG/file7f0a2c93dfb4/Plato/Crito/refined
#> Now symlink or copy raw data to /var/folders/33/vpk59z855zzbrczltjnc52qh0000gn/T//RtmpR1YZSG/file7f0a2c93dfb4/Plato/Crito/rawdata and run make_sample_sheet.
```

After that, you would edit `00driver.R` to customize how the markdown
are run.

``` r
project_title='Crito'
authors = 'default'
type = 'RNA'

# devtools::load_all("path/to/GeneseeBulk")
# devtools::load_all("path/to/Genesee")

## or
library(GeneseeBulk)
library(Genesee)

#### set up directories and create sample file
geneseesc_skeleton(main_directory, investigator, project_title, authors = authors, project_type = type)

rmarkdown::render('01_load.Rmd', params = list( input_tsv = "deSeq_counts.txt",
                                                input_project = NULL,
                                                sample_sheet_csv= "all_samples_file.csv",
                                                output_root= "refined/01load"),  output_dir = 'reports')
```

### Templates

At the moment, there are 2 template available:

``` r
list.files(proj_directory, pattern = '*.Rmd')
#> [1] "01_load.Rmd"   "02_qc_dea.Rmd"
```
