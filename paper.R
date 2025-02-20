# Figures, Tables, and Supplementary Material for the paper 

# tools -------

  library(tidyverse)
  library(exda)
  library(figur)
  library(rmarkdown)
  library(knitr)
  library(bookdown)
  library(flextable)
  library(writexl)
  library(soucer)
  library(stringi)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(pathview)
  library(coxExtensions)
  library(rlang)
  library(survival)
  
  explore <- exda::explore
  select <- dplyr::select
  reduce <- purrr::reduce
  set_rownames <- trafo::set_rownames
  extract <- clustTools::extract
  
  insert_head()
  
  c('./tools/globals.R', 
    './tools/functions.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# Figures, Tables, and Supplementary Material -------
  
  insert_msg('Figures, Tables, Supplements')
  
  c('./paper scripts/figures.R', 
    './paper scripts/supplementary_figures.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# END ------
  
  insert_tail()