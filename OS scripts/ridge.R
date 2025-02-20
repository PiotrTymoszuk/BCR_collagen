# Development of the collagen score in the training portion
# of the GSE16560 cohort
# The procedure: RIDGE with the initial explanatory variable set of 55
# collagen pathway genes (1st and 2nd order), 
# response: overall survival, family: Cox.
# Pre-processing of the explanatory variables: normalization Z/score
# Lambda finding in 200-repeat 10-fold CV
# The Collagen Score is defined as linear predictor score of the training
# model. 
# The model explanatory factors are the first- and second-order log2 expression 
# values of the collagen-related genes.
# Evaluation in the training collective and the GSE16560 cohort. 

  insert_head()
  
# container ------
  
  ridge_os <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals -------
  
  insert_msg('Globals')

  ## analysis tables: obtained from globals
  ## inclusion of the second order terms
  
  ridge_os$variables <- 
    c(os_globals$variables, 
      os_globals$sq_variables)
  
  ridge_os$data <- os_globals$data

  ## survival objects
  
  ridge_os$y <- ridge_os$data %>% 
    map(~Surv(.x$os_months, .x$death))

  ## matrices of normalized explanatory variables
  
  ridge_os$x <- ridge_os$data %>%  
    map(column_to_rownames, 'sample_id') %>% 
    map(~.x[ridge_os$variables]) %>% 
    map(as.matrix)
  
  ## CV folds 

  ridge_os$folds <- os_globals$folds
  
# Tuning of the lambda parameter ------
  
  insert_msg('Lambda tuning')
  
  ridge_os$lambda_tune <- 
    tune_glmnet(x = ridge_os$x$training, 
                y = ridge_os$y$training, 
                fold_ids = ridge_os$folds, 
                type.measure = 'default', 
                family = 'cox', 
                alpha = 0, 
                standardize = FALSE) ## Z-scores provided!

# Calculating the linear predictor scores for the training and test cohorts -------
  
  insert_msg('Calculating the collagen scores')
  
  ## predictions
  
  ridge_os$score_tbl <- ridge_os$x %>% 
    map(predict, object = ridge_os$lambda_tune) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(set_names, c('sample_id', 'collagen_score')) %>% 
    map(as_tibble)
  
  ## appending with the survival information
  
  ridge_os$score_tbl <- 
    map2(ridge_os$score_tbl, 
         map(ridge_os$data, 
             ~.x[c('sample_id', 'os_months', 'death')]), 
         left_join, by = 'sample_id')

# Building univariable Cox models ------
  
  insert_msg('Uni-variable Cox models')
  
  # working with metaprogramming, to get the entire
  # data sets kept in place with the models
  
  ridge_os$models <- ridge_os$score_tbl %>%
    map(~call2('coxph', 
               formula = Surv(os_months, death) ~ collagen_score, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., 
         ridge_os$score_tbl, 
         as_coxex)

# Characteristic of the Cox model: assumptions, fit stats and inference -----
  
  insert_msg('Characteristic of collagen score in the training cohort')

  ## assumptions: met
  
  ridge_os$assumptions <- ridge_os$models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  ridge_os$stats <- ridge_os$models %>% 
    map(summary, type = 'fit')
  
  ## inference: estimates expressed as hazard ratios
  
  ridge_os$inference <- ridge_os$models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  ridge_os[c("assumptions", "stats", "inference")] <- 
    ridge_os[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'tcga', 'training', 'test'))
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')

  ridge_os$calibration <- ridge_os$models %>% 
    future_map(calibrate.coxex, 
               n = 3, 
               labels = c('low', 'int', 'high'), 
               .options = furrr_options(seed = TRUE))
  
  ridge_os$global_cal <- ridge_os$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'tcga', 'training', 'test'))

  ridge_os$brier_scores <- ridge_os$models %>% 
    map(surv_brier)
    
# Differences in survival between the collagen score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  ## median survival in the tertiles and Peto-Peto test
  
  ridge_os$tertile_stats <- ridge_os$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  ridge_os$tertile_test <- ridge_os$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'tcga', 'training', 'test'), 
           p_value = pval) %>% 
    re_adjust

# training model estimates ------
  
  insert_msg('Estimates of the training model')
  
  ridge_os$coefs <- ridge_os$lambda_tune$model %>% 
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

  ridge_os$x <- NULL
  ridge_os$y <- NULL
  ridge_os$data <- NULL
  ridge_os$n_rep <- NULL
  ridge_os$folds <- NULL
  
  ridge_os <- compact(ridge_os)
  
  save(ridge_os, file = './cache/ridge_os.RData')

# END -----
  
  plan('sequential')
  
  insert_tail()