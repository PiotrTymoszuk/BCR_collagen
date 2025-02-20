# Exploratory analysis of expression of collagen-related genes shared by all 
# transcriptomic datasets
#
# 1) Characteristic of the cohorts and samples
#
# 2) Distribution of the collagen-related gene variables 
# (Shapiro-Wilk test) and their information content 
# (element frequency, variance and Gini coefficient)
#
# 3) Co-expression analysis of the genes in the cancer tissue. This analysis 
# step involves calculation of Euclidean distances between the genes and 
# multi-dimensional scaling as well as by a PCA
#
# 4) Co-expression network analysis for the collagen-related genes 

# tools ------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)
  library(readxl)
  
  library(exda)
  library(fastTest)
  library(microViz)
  library(clustTools)
  
  library(igraph)
  library(graphExtra)

  library(survival)
  library(survminer)
  
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
  
# Exploration analysis scripts --------
  
  insert_msg('Exploration scripts')
  
  ## general distribution stats
  ## coexpression analysis and co-expression networks
  
  c('./exploration scripts/cohorts.R', 
    './exploration scripts/distribution.R', 
    './exploration scripts/coexpression.R', 
    './exploration scripts/networks.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
  ## characteristic of the pooled GEO cohort
  
  c('./exploration scripts/pooled_geo.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# END -------
  
  insert_tail()