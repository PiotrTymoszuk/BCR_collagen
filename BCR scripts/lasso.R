# Development of the collagen score in the pooled GEO cohort
# The procedure: LASSO with the initial explanatory variable set of 55
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
  
  lasso_surv <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals -------
  
  insert_msg('Globals')

  ## analysis tables: obtained from globals
  ## inclusion of the second order terms
  
  lasso_surv$variables <- 
    c(surv_globals$variables, 
      surv_globals$sq_variables)
  
  lasso_surv$data <- surv_globals$data

  lasso_surv$data <- lasso_surv$data %>% 
    map(~filter(.x, complete.cases(.x))) 

  ## survival objects
  
  lasso_surv$y <- lasso_surv$data %>% 
    map(~Surv(.x$scaled_rfs_months, .x$relapse))

  ## matrices of normalized explanatory variables
  
  lasso_surv$x <- lasso_surv$data %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(~.x[lasso_surv$variables]) %>% 
    map(as.matrix)
  
  ## CV folds 

  lasso_surv$folds <- surv_globals$folds
  
# Tuning of the lambda parameter ------
  
  insert_msg('Lambda tuning')
  
  lasso_surv$lambda_tune <- 
    tune_glmnet(x = lasso_surv$x$geo_pool, 
                y = lasso_surv$y$geo_pool, 
                fold_ids = lasso_surv$folds, 
                type.measure = 'default', 
                family = 'cox', 
                alpha = 1, 
                standardize = FALSE) ## Z-scores provided!

# Calculating the linear predictor scores for the training and test cohorts -------
  
  insert_msg('Calculating the collagen scores')
  
  ## predictions
  
  lasso_surv$score_tbl <- lasso_surv$x %>% 
    map(predict, object = lasso_surv$lambda_tune) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(set_names, c('sample_id', 'collagen_score')) %>% 
    map(as_tibble)
  
  ## appending with the survival information
  
  lasso_surv$score_tbl <- 
    map2(lasso_surv$score_tbl, 
         map(lasso_surv$data, 
             ~.x[c('sample_id', 'scaled_rfs_months', 'relapse')]), 
         left_join, by = 'sample_id')

# Building univariable Cox models ------
  
  insert_msg('Uni-variable Cox models')
  
  # working with metaprogramming, to get the entire
  # data sets kept in place with the models
  
  lasso_surv$models <- lasso_surv$score_tbl %>%
    map(~call2('coxph', 
               formula = Surv(scaled_rfs_months, relapse) ~ collagen_score, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., 
         lasso_surv$score_tbl, 
         as_coxex)

# Characteristic of the Cox model: assumptions, fit stats and inference -----
  
  insert_msg('Characteristic of collagen score in the training cohort')

  ## assumptions: mild violation of the PHZ assumption for the DKFZ cohort
  
  lasso_surv$assumptions <- lasso_surv$models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  lasso_surv$stats <- lasso_surv$models %>% 
    map(summary, type = 'fit')
  
  ## inference: estimates expressed as hazard ratios
  
  lasso_surv$inference <- lasso_surv$models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  lasso_surv[c("assumptions", "stats", "inference")] <- 
    lasso_surv[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'geo_pool', 'training', 'test'))
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')

  lasso_surv$calibration <- lasso_surv$models %>% 
    future_map(calibrate.coxex, 
               n = 3, 
               labels = c('low', 'int', 'high'), 
               .options = furrr_options(seed = TRUE))
  
  lasso_surv$global_cal <- lasso_surv$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo_pool', 'training', 'test'))

  lasso_surv$brier_scores <- lasso_surv$models %>% 
    map(surv_brier)
    
# Differences in survival between the collagen score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  ## median survival in the tertiles and Peto-Peto test
  
  lasso_surv$tertile_stats <- lasso_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  lasso_surv$tertile_test <- lasso_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo_pool', 'training', 'test'), 
           p_value = pval) %>% 
    re_adjust

# training model estimates ------
  
  insert_msg('Estimates of the training model')
  
  lasso_surv$coefs <- lasso_surv$lambda_tune$model %>% 
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

  lasso_surv$x <- NULL
  lasso_surv$y <- NULL
  lasso_surv$data <- NULL
  lasso_surv$n_rep <- NULL
  lasso_surv$folds <- NULL
  
  lasso_surv <- compact(lasso_surv)
  
  save(lasso_surv, file = './cache/lasso_surv.RData')

# END -----
  
  plan('sequential')
  
  insert_tail()