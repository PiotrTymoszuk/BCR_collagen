# Development of the collagen score in the pooled GEO training cohort
# The procedure: RIDGE with the initial explanatory variable set of 55
# collagen pathway genes (1st and 2nd order), 
# response: relapse-free survival, family: Cox.
# Pre-processing of the explanatory variables: normalization Z/score
# Lambda finding in 200-repeat 10-fold CV
# The Collagen Score is defined as linear predictor score of the training
# model. 
# The model explanatory factors are the first- and second-order log2 expression 
# values of the collagen-related genes.

  insert_head()
  
# container ------
  
  ridge_surv <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals -------
  
  insert_msg('Globals')

  ## analysis tables: obtained from globals
  ## inclusion of the second order terms
  
  ridge_surv$variables <- 
    c(surv_globals$variables, 
      surv_globals$sq_variables)
  
  ridge_surv$data <- surv_globals$data

  ridge_surv$data <- ridge_surv$data %>% 
    map(~filter(.x, complete.cases(.x))) 

  ## survival objects
  
  ridge_surv$y <- ridge_surv$data %>% 
    map(~Surv(.x$scaled_rfs_months, .x$relapse))

  ## matrices of normalized explanatory variables
  
  ridge_surv$x <- ridge_surv$data %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(~.x[ridge_surv$variables]) %>% 
    map(as.matrix)
  
  ## CV folds 

  ridge_surv$folds <- surv_globals$folds
  
# Tuning of the lambda parameter ------
  
  insert_msg('Lambda tuning')
  
  ridge_surv$lambda_tune <- 
    tune_glmnet(x = ridge_surv$x$geo_pool, 
                y = ridge_surv$y$geo_pool, 
                fold_ids = ridge_surv$folds, 
                type.measure = 'default', 
                family = 'cox', 
                alpha = 0, 
                standardize = FALSE) ## Z-scores provided!

# Calculating the linear predictor scores for the training and test cohorts -------
  
  insert_msg('Calculating the collagen scores')
  
  ## predictions
  
  ridge_surv$score_tbl <- ridge_surv$x %>% 
    map(predict, object = ridge_surv$lambda_tune) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(set_names, c('sample_id', 'collagen_score')) %>% 
    map(as_tibble)
  
  ## appending with the survival information
  
  ridge_surv$score_tbl <- 
    map2(ridge_surv$score_tbl, 
         map(ridge_surv$data, 
             ~.x[c('sample_id', 'scaled_rfs_months', 'relapse')]), 
         left_join, by = 'sample_id')

# Building univariable Cox models ------
  
  insert_msg('Uni-variable Cox models')
  
  # working with metaprogramming, to get the entire
  # data sets kept in place with the models
  
  ridge_surv$models <- ridge_surv$score_tbl %>%
    map(~call2('coxph', 
               formula = Surv(scaled_rfs_months, relapse) ~ collagen_score, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., 
         ridge_surv$score_tbl, 
         as_coxex)

# Characteristic of the Cox model: assumptions, fit stats and inference -----
  
  insert_msg('Characteristic of collagen score in the training cohort')

  ## assumptions: mild violation of the PHZ assumption for the DKFZ cohort
  
  ridge_surv$assumptions <- ridge_surv$models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  ridge_surv$stats <- ridge_surv$models %>% 
    map(summary, type = 'fit')
  
  ## inference: estimates expressed as hazard ratios
  
  ridge_surv$inference <- ridge_surv$models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  ridge_surv[c("assumptions", "stats", "inference")] <- 
    ridge_surv[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'geo_pool', 'training', 'test'))
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')

  ridge_surv$calibration <- ridge_surv$models %>% 
    future_map(calibrate.coxex, 
               n = 3, 
               labels = c('low', 'int', 'high'), 
               .options = furrr_options(seed = TRUE))
  
  ridge_surv$global_cal <- ridge_surv$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo_pool', 'training', 'test'))

  ridge_surv$brier_scores <- ridge_surv$models %>% 
    map(surv_brier)
    
# Differences in survival between the collagen score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  ## median survival in the tertiles and Peto-Peto test
  
  ridge_surv$tertile_stats <- ridge_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  ridge_surv$tertile_test <- ridge_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo_pool', 'training', 'test'), 
           p_value = pval) %>% 
    re_adjust

# training model estimates ------
  
  insert_msg('Estimates of the training model')
  
  ridge_surv$coefs <- ridge_surv$lambda_tune$model %>% 
    coef %>% 
    as.matrix %>% 
    as.data.frame %>% 
    rownames_to_column('variable') %>% 
    set_names(c('variable', 'coef')) %>% 
    mutate(exp_coef = exp(coef)) %>% 
    filter(coef != 0) %>% 
    as_tibble
  
# Caching the results -------
  
  insert_msg('Caching the results')

  ridge_surv$x <- NULL
  ridge_surv$y <- NULL
  ridge_surv$data <- NULL
  ridge_surv$n_rep <- NULL
  ridge_surv$folds <- NULL
  
  ridge_surv <- compact(ridge_surv)
  
  save(ridge_surv, file = './cache/ridge_surv.RData')

# END -----
  
  plan('sequential')
  
  insert_tail()