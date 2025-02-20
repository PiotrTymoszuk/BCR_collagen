# Analysis globals for univariable analysis of BCR-free survival. 
# Those are data frames with ComBat-adjusted log2 estimates of expression 
# levels

  insert_head()

# container ------
  
  uni_globals <- list()
  
# BCR-free survival and ComBat expression -------
  
  insert_msg('Survival and gene expression')
  
  ## cohorts of interest
  
  uni_globals$cohorts <- c('geo_pool', 'tcga', 'dkfz')
  
  ## survival: plan survival in months, adding one month 
  ## to BCR-free survival to avoid zero survival times
  
  uni_globals$survival <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$clinic) %>% 
    map(function(x) if('relapse' %in% names(x)) x else NULL) %>% 
    compact %>% 
    map(select, sample_id, relapse, rfs_months) %>% 
    map(mutate, 
        rfs_months = rfs_months + 1)

  ## expression: the pooled GEO cohort, TCGA and DKFZ, 
  ## adding second-order terms (used by regularized expression), 
  ## the expression is left as log2 expression, not Z-scores
  
  uni_globals$variables <- globals$genes
  
  uni_globals$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$combat) %>% 
    map(select, sample_id, all_of(uni_globals$variables))

  ## merging with the survival data
  
  uni_globals$data <- 
    map2(uni_globals$survival[uni_globals$cohorts], 
         uni_globals$expression[uni_globals$cohorts], 
         inner_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x)))
  
# Total numbers of observations and deaths --------
  
  insert_msg('Total number of observations and deaths')
  
  uni_globals$n_numbers$n_total <- uni_globals$data %>% 
    map_dbl(nrow)
  
  uni_globals$n_numbers$n_events <- uni_globals$data %>% 
    map(filter, relapse == 1) %>% 
    map_dbl(nrow)
  
# END -------
  
  insert_tail()