# Co-expression networks analyses for the collagen-related genes

# tools --------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)

  library(igraph)
  library(graphExtra)
  
  library(ggrepel)
  library(ggtext)
  library(RColorBrewer)
  
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
  
# analysis scripts --------
  
  insert_msg('Analysis scripts')
  
  c('./network scripts/networks.R') %>% 
    source_all(message = TRUE)

# END ------
  
  insert_tail()