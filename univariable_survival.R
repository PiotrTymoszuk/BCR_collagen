# Uni-variable analysis of BCR-free survival as a function of dichtomized 
# expression of the collagen-related genes in the prostate cancer tissue. 

# tools ------

  library(tidyverse)
  library(trafo)
  library(rlang)
  library(stringi)

  library(survival)
  library(survminer)
  library(kmOptimizer)

  library(furrr)
  library(soucer)
  
  library(ggtext)
  
  explore <- exda::explore
  select <- dplyr::select
  reduce <- purrr::reduce
  set_rownames <- trafo::set_rownames
  
  c('./tools/globals.R', 
    './tools/functions.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# analysis globals ------
  
  insert_msg('Analysis globals')
  
  c('./univariable survival scripts/globals.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# modeling and visualization scripts --------
  
  insert_msg('Modeling and visualization scripts')
  
  access_cache(cache_path = './cache/uni_cut.RData', 
               script_path = './univariable survival scripts/univariable_rfs.R', 
               message = 'Loading cached univairable modeling results')
  
  ## visualization scripts
  
  c('./univariable survival scripts/plots.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# END --------
  
  insert_tail()