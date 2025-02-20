# functions used specifically for support vector machines and 
# gradient boosted machine modeling of survival.
# 

# tools --------

  library(tidyverse)
  library(trafo)
  library(rlang)

  library(survival)
  library(survminer)
  library(survivalsvm)
  library(coxExtensions)
  library(gbm)

# Predictor scores and performance stats of a SVM survival model ------

  svm_score <- function(prediction, 
                        data, 
                        time_variable = 'rfs_months', 
                        event_variable = 'relapse', 
                        n_cuts = 3, 
                        labels = c('low', 'int', 'high')) {
    
    ## computes the SVM score (e.g. regression predictors or ranks)
    
    score_tbl <- 
      tibble(sample_id = rownames(data), 
             !!time_variable := data[[time_variable]], 
             !!event_variable := data[[event_variable]], 
             svm_score = prediction$predicted[1, ])
    
    cuts <- (1:(n_cuts - 1))/n_cuts
    
    quants <- quantile(score_tbl$svm_score, cuts, na.rm = TRUE)
    
    score_tbl %>% 
      mutate(score_cuts = cut(svm_score, 
                              c(-Inf, quants, Inf), 
                              labels))
    
    
  }

# Out-of-fold prediction stats for SVM and GBM ------

  svm_oof <- function(data, 
                      time_variable = 'rfs_months', 
                      event_variable = 'relapse', 
                      folds, ...) {
    
    ## given a data set, a formula as well as a list of assignments 
    ## to the training subsets ('folds' argument), the function constructs 
    ## a SVM model in the training portions of the data set and computes 
    ## the predictor scores for the test portion. 
    ## For such out-of-fold predictor scores, C-indexes are calculated
    
    ## fold assignment list and data split, model formulas -------
    
    data_splits <- folds %>% 
      map(~list(train = data[.x, ], 
                test = data[-.x, ])) %>% 
      transpose
    
    mod_formula <- paste0('Surv(', time_variable, ', ', 
                          event_variable, ') ~.') %>% 
      as.formula

    
    cox_formula <- paste0('Surv(', time_variable, ', ', 
                          event_variable, ') ~ svm_score') %>% 
      as.formula
    
    ## construction of SVM models and predictions --------
    
    models <- data_splits$train %>% 
      future_map(function(x) survivalsvm(formula = mod_formula, 
                                         data = x, ...), 
                 .options = furrr_options(seed = TRUE))
    
    
    predictions <- 
      list(object = models, 
           newdata = data_splits$test) %>% 
      pmap(predict)
    
    score_tbl <- 
      map2(predictions, 
           data_splits$test, 
           svm_score, 
           time_variable = time_variable, 
           event_variable = event_variable)
    
    ## computation of the concordance indexes and Brier scores --------
    
    oof_stats <- score_tbl %>% 
      map(concordance, 
          object = cox_formula) %>% 
      map_dbl(~.x$concordance) %>% 
      compress(names_to = 'fold_id', 
               values_to = 'c_index')
    
    ## output: average stats ------
    
    if(stri_detect(oof_stats$fold_id[[1]], fixed = '.')) {
      
      oof_stats <- oof_stats %>% 
        mutate(repetition = stri_split_fixed(fold_id, 
                                             pattern = '.', 
                                             simplify = TRUE)[, 1], 
               fold = stri_split_fixed(fold_id, 
                                       pattern = '.', 
                                       simplify = TRUE)[, 2])
      
      summary_stats <- oof_stats %>% 
        blast(repetition, .skip = TRUE) %>% 
        map(select, -fold_id, -fold) %>% 
        map(colMeans) %>% 
        reduce(rbind) %>% 
        colMeans
      
    } else {
      
      oof_stats <- oof_stats %>% 
        mutate(repetition = NA, 
               fold = fold_id)
      
      summary_stats <- oof_stats %>% 
        select(-fold, -fold_id, -repetition) %>% 
        colMeans
      
    }
    
    summary_stats <- summary_stats %>% 
      as.matrix %>% 
      t %>% 
      as.data.frame %>% 
      as_tibble
    
    list(oof = oof_stats, 
         summary = summary_stats)
    
  }
  
  gbm_oof <- function(data, 
                      time_variable = 'rfs_months',
                      event_variable = 'relapse', 
                      n_folds = 10, ...) {
    
    ## OOF C-index for a GBM model
    
    mod_formula <- paste0('Surv(', time_variable, 
                          ', ', event_variable, ') ~ .') %>% 
      as.formula
    
    gbm_model <- gbm(formula = mod_formula, 
                     data = data, 
                     cv.folds = n_folds, ...)
    
    min(gbm_model$cv.error)
    
  }

# tuning of survival SVM and GBM models ---------
  
  svm_tune <- function(data, 
                       time_variable = 'rfs_months', 
                       event_variable = 'relapse', 
                       folds, 
                       tune_grid, ...) {
    
    ## cross-validation tuning of survival SVM models

    tune_cond <- paste0('cond_', 1:nrow(tune_grid))
    
    tune_stats <- tune_grid %>% 
      pmap(svm_oof, 
           data = data, 
           time_variable = time_variable, 
           event_variable = event_variable, 
           folds = folds, ...) %>% 
      set_names(tune_cond) %>% 
      transpose
    
    ## summary stats and the best tune
    
    tune_grid <- as_tibble(tune_grid) %>% 
      mutate(tune_id = tune_cond)
    
    summary_stats <- tune_stats$summary %>% 
      compress(names_to = 'tune_id')
    
    summary_stats <- inner_join(tune_grid, summary_stats, by = 'tune_id')

    best_tune <- summary_stats %>% 
      filter(.data[['c_index']] == max(.data[['c_index']]))
    
    best_tune <- best_tune[1, ]
    
    best_folds <- tune_stats[['oof']][[best_tune$tune_id]]
    
    ## output 
    
    list(summary = summary_stats, 
         best_tune = best_tune, 
         best_oof = best_folds)
    
  }
  
  gbm_tune <- function(data, 
                       time_variable = 'rfs_months',
                       event_variable = 'relapse', 
                       n_folds = 10, 
                       tune_grid, ...) {
    
    ## CV tuning of GBM cox models.
    
    ## construction of models, extraction of the errors -------
    
    tune_cond <- paste0('cond_', 1:nrow(tune_grid))

    cv_errors <- tune_grid %>% 
      pmap_dbl(gbm_oof, 
               data = data, 
               time_variable = time_variable, 
               event_variable = event_variable, 
               n_folds = n_folds, ...)
    
    cv_errors <- set_names(cv_errors, tune_cond) %>% 
      compress(names_to = 'condition', 
               values_to = 'cv_deviance')
    
    tune_grid <- tune_grid %>% 
      mutate(condition = tune_cond)
    
    summary <- left_join(tune_grid, cv_errors, by = 'condition')
    
    best_tune <- summary %>% 
      filter(cv_deviance == min(cv_deviance, na.rm = TRUE))
    
    best_tune <- best_tune[1, ]
    
    list(summary = summary, 
         best_tune = best_tune)
    
  }
  
# Permutation importance of survival SVM --------
  
  svm_importance <- function(data, 
                             time_variable = 'rfs_months', 
                             event_variable = 'relapse', 
                             folds, ...) {
    
    ## permutation importance of explanatory variables. 
    ## The permutation measure is OOF delta in C-index between the full model
    ## and each of the model with random re-shuffling of an explanatory variable
    
    ## re-shuffled data sets ------
    
    data_lst <- list()
    
    data_lst$full <- data
    
    expl_variables <- 
      names(data)[!names(data) %in% c(time_variable, event_variable)]
    
    for(i in expl_variables) {
      
      data_lst[[i]] <- data %>% 
        mutate(!!i := sample(.data[[i]], size = nrow(data), replace = TRUE))
      
    }
    
    ## construction of the models, and computation of CV stats ----------
    
    cv_stats <- data_lst %>% 
      map(svm_oof, 
          time_variable = time_variable, 
          event_variable = event_variable, 
          folds = folds, ...) %>% 
      map(~.x$summary)
    
    cv_stat_tbl <- cv_stats %>% 
      compress(names_to = 'variable') %>% 
      mutate(delta_c_index = cv_stats[['full']]$c_index - c_index)
      
    cv_stat_tbl
    
  }
  
# END ------
