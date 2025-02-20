# Development of the collagen score in the training portion
# of the GSE16560 cohort
# The procedure: Elastic Net with the initial explanatory variable set of 55
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
  
  lasso_os <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals -------
  
  insert_msg('Globals')

  ## analysis tables: obtained from globals
  ## inclusion of the second order terms
  
  lasso_os$variables <- 
    c(os_globals$variables, 
      os_globals$sq_variables)
  
  lasso_os$data <- os_globals$data

  ## survival objects
  
  lasso_os$y <- lasso_os$data %>% 
    map(~Surv(.x$os_months, .x$death))

  ## matrices of normalized explanatory variables
  
  lasso_os$x <- lasso_os$data %>%  
    map(column_to_rownames, 'sample_id') %>% 
    map(~.x[lasso_os$variables]) %>% 
    map(as.matrix)
  
  ## CV folds 

  lasso_os$folds <- os_globals$folds
  
# Tuning of the lambda parameter ------
  
  insert_msg('Lambda tuning')
  
  lasso_os$lambda_tune <- 
    tune_glmnet(x = lasso_os$x$training, 
                y = lasso_os$y$training, 
                fold_ids = lasso_os$folds, 
                type.measure = 'default', 
                family = 'cox', 
                alpha = 1, 
                standardize = FALSE) ## Z-scores provided!

# Calculating the linear predictor scores for the training and test cohorts -------
  
  insert_msg('Calculating the collagen scores')
  
  ## predictions
  
  lasso_os$score_tbl <- lasso_os$x %>% 
    map(predict, object = lasso_os$lambda_tune) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(set_names, c('sample_id', 'collagen_score')) %>% 
    map(as_tibble)
  
  ## appending with the survival information
  
  lasso_os$score_tbl <- 
    map2(lasso_os$score_tbl, 
         map(lasso_os$data, 
             ~.x[c('sample_id', 'os_months', 'death')]), 
         left_join, by = 'sample_id')

# Building univariable Cox models ------
  
  insert_msg('Uni-variable Cox models')
  
  # working with metaprogramming, to get the entire
  # data sets kept in place with the models
  
  lasso_os$models <- lasso_os$score_tbl %>%
    map(~call2('coxph', 
               formula = Surv(os_months, death) ~ collagen_score, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., 
         lasso_os$score_tbl, 
         as_coxex)

# Characteristic of the Cox model: assumptions, fit stats and inference -----
  
  insert_msg('Characteristic of collagen score in the training cohort')

  ## assumptions: met
  
  lasso_os$assumptions <- lasso_os$models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  lasso_os$stats <- lasso_os$models %>% 
    map(summary, type = 'fit')
  
  ## inference: estimates expressed as hazard ratios
  
  lasso_os$inference <- lasso_os$models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  lasso_os[c("assumptions", "stats", "inference")] <- 
    lasso_os[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'tcga', 'training', 'test'))
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')

  lasso_os$calibration <- lasso_os$models %>% 
    future_map(calibrate.coxex, 
               n = 3, 
               labels = c('low', 'int', 'high'), 
               .options = furrr_options(seed = TRUE))
  
  lasso_os$global_cal <- lasso_os$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'tcga', 'training', 'test'))

  lasso_os$brier_scores <- lasso_os$models %>% 
    map(surv_brier)
    
# Differences in survival between the collagen score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  ## median survival in the tertiles and Peto-Peto test
  
  lasso_os$tertile_stats <- lasso_os$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  lasso_os$tertile_test <- lasso_os$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'tcga', 'training', 'test'), 
           p_value = pval) %>% 
    re_adjust

# training model estimates ------
  
  insert_msg('Estimates of the training model')
  
  lasso_os$coefs <- lasso_os$lambda_tune$model %>% 
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

  lasso_os$x <- NULL
  lasso_os$y <- NULL
  lasso_os$data <- NULL
  lasso_os$n_rep <- NULL
  lasso_os$folds <- NULL
  
  lasso_os <- compact(lasso_os)
  
  save(lasso_os, file = './cache/lasso_os.RData')

# END -----
  
  plan('sequential')
  
  insert_tail()