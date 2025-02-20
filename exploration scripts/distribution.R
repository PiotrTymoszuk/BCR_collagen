# Distribution tests for the collagen variables (ComBat-adjusted log2 expression):  
# Normality check by Shapiro-Wilk test, distribution statistics
# (means, variances, Gini indexes, and fractions of unique values)

  insert_head()
  
# container --------
  
  distr <- list()
  
# expression data frames -------
  
  insert_msg('Expression data')
  
  ## variables of interest
  
  distr$variables <- globals$genes
  
  ## ComBat-adjusted log2 expression estimates
  
  distr$data <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$combat) %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, all_of(distr$variables))
  
# Normality by Shapiro-Wilk test --------
  
  insert_msg('Normality by Shapiro-Wilk test')
  
  ## violations assumed for W < 0.9 and p < 0.05
  
  distr$normality <- distr$data %>% 
    map(f_shapiro_test) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'variable') %>%
    map(mutate, 
        violated = ifelse(w < 0.9 & p_value < 0.05,
                          'yes', 'no')) %>% 
    map(as_tibble)
  
  ## violations and number of violations: 
  ## potentially problematic variables for regression models are 
  ## MMP13, COL2A1, COL9A1, COL11A2 - normality is violated in 
  ## four cohorts and more

  distr$norm_violations <- distr$normality %>% 
    map(filter, violated == 'yes') %>% 
    map(~.x$variable)
  
  distr$number_violations <- distr$norm_violations %>% 
    count_features
  
# Variance and information content ------
  
  insert_msg('Variance and Gini coefficients')
  
  ## at a glance, there are no uninformative variables with 
  ## near-zero Gini indexes or many non-unique values
  
  distr$stats <- distr$data %>% 
    map(map_dfc, minimum_shift) %>% 
    map(distr_stats) %>% 
    map(arrange, gini_coef)
    
# END ------
  
  distr <- distr[c("normality", "stats")]
  
  insert_tail()