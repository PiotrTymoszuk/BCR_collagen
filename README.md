# BCR_collagen
Modeling of biochemical relapse-free survival in PCa with expression of collagen-related genes

## Figures, Tables, and Supplementary Material

The transcriptomic part of the main manuscript text is available as a [MS Word file](https://github.com/PiotrTymoszuk/BCR_collagen/blob/main/paper/manuscript_figures_tables.docx). 
Supplementary Material for the transcriptomic part of the manuscript is available as a [MS Word file](https://github.com/PiotrTymoszuk/BCR_collagen/blob/main/paper/manuscript_supplement.docx).

PDF files with main paper Figures are available [here](https://github.com/PiotrTymoszuk/BCR_collagen/tree/main/paper/figures).
PDF files with Supplementary Figures are available [here](https://github.com/PiotrTymoszuk/BCR_collagen/tree/main/paper/supplementary%20figures).

## Terms of use

Please cite the repository and the peer-reviewed publication, when available. The analysis pipeline is [GPL-3 licensed](https://github.com/PiotrTymoszuk/BCR_collagen/blob/main/LICENSE). 
The raw data files will be made upon request to the senior study authors, [Dr. Isabel Heidegger](mailto:isabel-maria.heidegger@i-med.ac.at) and [Dr. Maria Frantzi](mailto:frantzi@mosaiques.de).

## Usage

The following development packages are required to run the pipeline:

```r
## handling of tabular data
devtools::install_github('PiotrTymoszuk/trafo')

## script sourcing
devtools::install_github('PiotrTymoszuk/soucer') 

## exploratory data analysis and statistical hypothesis testing
devtools::install_github('PiotrTymoszuk/ExDA') 
devtools::install_github('PiotrTymoszuk/fastTest')
devtools::install_github('PiotrTymoszuk/microViz')
devtools::install_github('PiotrTymoszuk/clustTools')

## network analyses
devtools::install_github('PiotrTymoszuk/graphExtra')

## management of figures and tables in Rmd documents
devtools::install_github('PiotrTymoszuk/figur')

## handling and analysis of Cox proportional hazard models
devtools::install_github('PiotrTymoszuk/coxExtensions')

```

Source 'exec.R' to launch the entire pipeline:

```r

source('exec.R')

```

## Contact

The repository maintainer is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com). Data requests should be addressed to [Dr. Isabel Heidegger](mailto:isabel-maria.heidegger@i-med.ac.at) and [Dr. Maria Frantzi](mailto:frantzi@mosaiques.de)

<br>
