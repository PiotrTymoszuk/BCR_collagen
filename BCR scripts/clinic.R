# Survival gradient boosted machines with clinical predictors of 
# BCR risk. 
# Explanatory variables: Gleason (ISUP), log(PSA) at diagnosis and 
# pT stage (3+ and 2 versus 1).
# Response: BCR-free survival, with the time variable subjected to min/max 
# scaling. 
# Tuning andf training in the pooled GEO cohort. 
# The tuning criterion was minimization of the model deviance in 10-fold 
# cross-validation. 
# The trained model generates 
# predictor scores ('clinic_score') for the training collective and test cohorts 
# (TCGA and DKFZ). 
# Evaluation: evaluation of performance of univariable Cox models of BCR-free 
# survival as a function of the clinical predictor score. 
# Variable importance: hazard ratios for the explanatory factors
#
# Performance of the 'clinical' model is compared with performance of the GBM 
# model that combines the clinical predictors with expression levels of the 
# collagen-related genes

  insert_head()
  
# container -------
  
  surv_clin <- list()

# analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## analysis data: from the globals
  
  surv_clin$data <- 
    map2(map(surv_globals$data, 
             ~.x[, c('sample_id', 'scaled_rfs_months', 'relapse')]), 
         surv_globals$clinic, 
         inner_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(column_to_rownames, 'sample_id')
  
  ## tuning grid for the GBM model: partially optimized per hand
  
  surv_clin$tune_grid <- surv_globals$gbm_grid
  
# CV tuning in the GEO data set -------
  
  insert_msg('CV tuning in the GEO data set')
  
  set.seed(12345)
  
  surv_clin$tuning <- gbm_tune(data = surv_clin$data$geo_pool, 
                                time_variable = 'scaled_rfs_months', 
                                event_variable = 'relapse', 
                                n_folds = 10, 
                                tune_grid = surv_clin$tune_grid, 
                                distribution = 'coxph')
  
# Training of the GBM model in the GEO cohort ------
  
  insert_msg('Training the GEO model')
  
  set.seed(12345)
  
  surv_clin$gbm_model <- 
    gbm(formula = Surv(scaled_rfs_months, relapse) ~ ., 
        data = surv_clin$data$geo_pool, 
        distribution = 'coxph', 
        n.trees = surv_clin$tuning$best_tune$n.trees[1],
        shrinkage = surv_clin$tuning$best_tune$shrinkage[1],
        interaction.depth = surv_clin$tuning$best_tune$interaction.depth[1],
        n.minobsinnode = surv_clin$tuning$best_tune$n.minobsinnode[1],
        cv.folds = 10) 
  
  ## the optimal iteration number for making predictions
  ## as specified by the minimal CV deviance
  
  surv_clin$best_iter <- gbm.perf(surv_clin$gbm_model, method = 'cv')
  
# Predictions and score tables --------
  
  insert_msg('Predictions and score tables')
  
  surv_clin$score_tbl <- surv_clin$data %>% 
    map(predict, 
        object = surv_clin$gbm_model, 
        n.trees = surv_clin$best_iter) %>% 
    map2(surv_clin$data, ., 
         ~mutate(.x[c('scaled_rfs_months', 'relapse')], 
                 clinic_score = .y)) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(as_tibble)
  
# Univariable Cox models -------
  
  insert_msg('Univariable Cox models')
  
  surv_clin$cox_models <- surv_clin$score_tbl %>%
    map(~call2('coxph', 
               formula = Surv(scaled_rfs_months, relapse) ~ clinic_score, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., 
         surv_clin$score_tbl, 
         as_coxex)
  
# Characteristic of the Cox models: assumptions, fit stats and inference -----
  
  insert_msg('Characteristic of collagen scores')
  
  ## assumptions: met
  
  surv_clin$assumptions <- surv_clin$cox_models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  surv_clin$stats <- surv_clin$cox_models %>% 
    map(summary, type = 'fit')
  
  ## inference
  
  surv_clin$inference <- surv_clin$cox_models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  surv_clin[c("assumptions", "stats", "inference")] <- 
    surv_clin[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'geo_pool', 'training', 'test'))
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')
  
  surv_clin$calibration <- surv_clin$cox_models %>% 
    map(calibrate, 
        n = 3, 
        labels = c('low', 'int', 'high'))
  
  surv_clin$global_cal <- surv_clin$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo_pool', 'training', 'test'))
  
  surv_clin$brier_scores <- surv_clin$cox_models %>% 
    map(surv_brier)
  
# Differences in survival between the score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  surv_clin$tertile_stats <- surv_clin$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  surv_clin$tertile_test <- surv_clin$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo_pool', 'training', 'test')) %>% 
    re_adjust(p_variable = 'pval')
  
# Variable importance ------
  
  insert_msg('Variable importance')
  
  surv_clin$importance <- surv_clin$gbm_model %>% 
    summary(n.trees = surv_clin$best_iter) %>% 
    select(var, rel.inf) %>% 
    set_names(c('variable', 'rel_influence')) %>% 
    as_tibble

# Caching the results -------
  
  insert_msg('Caching the results')

  surv_clin$data <- NULL
  
  surv_clin <- compact(surv_clin)
  
  save(surv_clin, file = './cache/surv_clin.RData')
  
# END ------

  insert_tail()