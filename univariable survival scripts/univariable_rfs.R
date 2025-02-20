# Univariable analysis of overall survival in dichotomous gene strata. 
# significant association with survival is assumed for: 
# pFDR < 0.05 and Cox's HR > 0 (unfavorable marker), and pFDR < 0.05 
# and Cox's HR < 0 (favorable marker).

  insert_head()
  
# container -------
  
  uni_cut <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis data frames --------
  
  insert_msg('Analysis data frames')
  
  ## variables and ComBat-adjusted log2 expression values
  ## no Z-scores: after ComBat, the distribution of expression of the genes
  ## of interest is similar
  
  uni_cut$variables <- globals$genes

  uni_cut$data <- uni_globals$data

  ## minimal number of observations in strata: 25% of the cohort

  uni_cut$n_min <- floor(0.25 * map_dbl(uni_cut$data, nrow))

# Cutoff finding --------
  
  insert_msg('Cutoff finding')
  
  ## working in a safely mode, since some genes do not yield
  ## any sound cutoffs, e.g. due to constant expression
  
  uni_cut$cutoff_obj <- 
    list(x = uni_cut$data, 
         y = uni_cut$n_min) %>% 
    pmap(function(x, y) uni_cut$variables %>% 
           future_map(~safely(find_cutoff)(data = x, 
                                           time = 'rfs_months', 
                                           event = 'relapse', 
                                           variable = .x, 
                                           min_n = y, 
                                           .parallel = FALSE), 
                      .options = furrr_options(seed = TRUE)) %>% 
           set_names(uni_cut$variables))
  
  uni_cut$cutoff_obj <- uni_cut$cutoff_obj %>% 
    map(map, ~.x$result) %>% 
    map(compact)
  
# Testing summary -------
  
  insert_msg('Testing summary and significant genes')
  
  ## in case there are more cutoffs, the first is selected
  ## log-rank test summary: this won't be shown in the paper
  
  uni_cut$test <- uni_cut$cutoff_obj %>% 
    map(map, summary) %>% 
    map(map, ~.x[1, ]) %>% 
    map(compress, names_to = 'gene_symbol') %>% 
    map(re_adjust)

# Univariable Cox modeling, Cox versus low expressors -----
  
  insert_msg('Cox PH modeling, high versus low expressors')
  
  uni_cut$modeling <- uni_cut$cutoff_obj %>% 
    map(future_map_dfr, 
        hr_from_cut, 
        .options = furrr_options(seed = TRUE))
  
  ## formatting the Cox PH modeling results: 
  ## FDR correction and identification of significant 
  ## favorable and unfavorable markers
  
  uni_cut$modeling <- uni_cut$modeling %>% 
    map(re_adjust) %>% 
    map(mutate, 
        marker = ifelse(p_adjusted >= 0.05, 'ns', 
                        ifelse(hr < 1, 'favorable', 'unfavorable')), 
        marker = factor(marker, 
                        c('unfavorable', 'favorable', 'ns')))

# Significant effects --------
  
  insert_msg('Significant effects')
  
  ## significant effects in single cohorts
  
  uni_cut$significant <- uni_cut$modeling %>% 
    map(filter, 
        marker %in% c('unfavorable', 'favorable')) %>%
    map(blast, marker) %>% 
    transpose %>% 
    map(map, ~.x$gene_symbol)

  ## common significant genes: shared by all cohorts
  
  uni_cut$common_significant <- uni_cut$significant %>% 
    map(reduce, intersect)

# Caching the results ------
  
  insert_msg('Caching the results')
  
  uni_cut$variables <- NULL
  uni_cut$data <- NULL
  
  uni_cut <- compact(uni_cut)
  
  save(uni_cut, file = './cache/uni_cut.RData')
  
# END -------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()