# This script imports prostate cancer data from GEO 
# (GSE16560, GSE54460, GSE70768 and GSE70769, gse220095)
# and the TCGA and DKFZ RNA seq data obtained from cBioportal
# Expression values log2(x) (Microarray) or log2(x + 1) (RNAseq) transformed.
#
# 

# toolbox ----

  library(tidyverse)
  library(readxl)
  library(rlang)
  library(trafo)
  library(stringi)
  
  library(htGLMNET)
  library(fastTest)
  library(exda)
  library(microViz)

  library(GEOquery)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  
  library(furrr)
  library(soucer)

  insert_head()

  explore <- exda::explore
  select <- dplyr::select
  reduce <- purrr::reduce
  set_rownames <- trafo::set_rownames

  c('./tools/globals.R', 
    './tools/functions.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# executing data import scripts, from scratch if not done before, saving the raw data ----
  
  insert_msg('Reading the expression and clinical data')

  for(i in globals$analysis_studies[globals$analysis_studies != 'geo_pool']) {
    
    access_cache(cache_path = paste0('./data/', i, '.RData'), 
                 script_path = paste0('./import scripts/', 
                                      globals$study_labels[i], '.R'),
                 message = paste('Loading cleared data for:', i))
    
  }

# Treatment of gene expression with COMBAT and the pooled training cohort -------
  
  insert_msg('COMBAT and cohort pooling')

  access_cache(cache_path = './data/combat.RData',
               script_path = './import scripts/combat.R', 
               message = 'Loading cached COMBAT results')
  
  ## appending the data sets with the ComBat estimates
  
  for(i in globals$analysis_studies[globals$analysis_studies != 'geo_pool']) {
    
    assign_expr <- paste0(i, '$combat <- combat$expression$', i) %>% 
      parse_expr
    
    eval(assign_expr)
    
  }
  
  ## the pooled data set
  
  geo_pool <- list()
  
  geo_pool$studies <- c('gse54460', 'gse70768', 'gse70769', 'gse220095')
  
  pool_expr <- geo_pool$studies %>% 
    map(~paste0(.x, ' = ', .x)) %>% 
    paste(collapse = ', ') %>% 
    paste0('list(', ., ')') %>% 
    parse_expr
  
  geo_pool$clinic <- pool_expr %>% 
    eval %>% 
    map(~.x$clinic) %>% 
    map2(., names(.), ~mutate(.x, study = .y)) %>% 
    reduce(full_rbind) %>% 
    mutate(study = factor(study, geo_pool$studies))
  
  geo_pool$combat <- combat$expression[geo_pool$studies] %>% 
    reduce(rbind)
  
  rm(combat)

# END -----
  
  rm(i, assign_expr, pool_expr)
  
  insert_tail()
