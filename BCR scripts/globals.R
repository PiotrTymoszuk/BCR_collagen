# COMBAT expression estimates used for multi-parameter modeling of survival
# in th pooled GEO, TCGA, and DKFZ cohorts. 
#
# BCR-free survival (+1 months to avoid zero survival times) is min/max scaled.
# Clinical information (Gleason/ISUP, pT stage) is extracted 
# as well and will be used for construction of a clinical model. 
# Unfortunately, information on pre-operative PSA in the TCGA cohort is 
# provided only for a small minority of patients (n = 185 with 19 BCR cases), 
# which excludes it as a predictor of BCR-free survival. 
# Instead we include ComBat-adjusted log2 expression of PSA (KLK3 gene) in the 
# clinical predictor set (references: DOI 10.3325/cmj.2020.61.450 and 
# DO 10.1038/s41391-020-00283-3)

  insert_head()
  
# container ------
  
  surv_globals <- list()
  
# input data: COMBAT estimates and survival times -----
  
  insert_msg('COMBAT expression estimates and survival times')
  
  ## cohorts of interest
  
  surv_globals$cohorts <- c('geo_pool', 'tcga', 'dkfz')
  
  ## survival: min/max scaled survival times
  
  surv_globals$survival <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$clinic) %>% 
    map(function(x) if('relapse' %in% names(x)) x else NULL) %>% 
    compact %>% 
    map(select, sample_id, relapse, rfs_months) %>% 
    map(mutate, 
        rfs_months = rfs_months + 1, 
        scaled_rfs_months = minMax(rfs_months) + 0.01) %>% 
    map(select, -rfs_months)
  
  ## genes of interest: first and second order terms
  
  surv_globals$variables <- globals$genes
  
  surv_globals$sq_variables <- 
    paste0(surv_globals$variables, '_sq')
  
  ## expression: the pooled GEO cohort, TCGA and DKFZ, 
  ## adding second-order terms (used by regularized expression), 
  ## calculating Z-scores
  
  surv_globals$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$combat) %>% 
    map(select, sample_id, all_of(surv_globals$variables))
  
  for(i in surv_globals$variables) {
    
    surv_globals$expression <- surv_globals$expression %>% 
      map(mutate, 
          !!paste0(i, '_sq') := .data[[i]]^2)
    
  }
  
  for(i in names(surv_globals$expression)) {
    
    surv_globals$expression[[i]][c(surv_globals$variables, surv_globals$sq_variables)] <- 
      surv_globals$expression[[i]][c(surv_globals$variables, surv_globals$sq_variables)] %>% 
      map_dfc(zScores)
    
  }

  ## merging with the survival data

  surv_globals$data <- 
    map2(surv_globals$survival[surv_globals$cohorts], 
         surv_globals$expression[surv_globals$cohorts], 
         inner_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x)))
  
# N numbers -------
  
  insert_msg('N numbers')
  
  ## n numbers: total and events
  
  surv_globals$n_numbers <- surv_globals$data %>% 
    map(~tibble(n_total = nrow(.x), 
                n_events = sum(.x$relapse))) %>% 
    compress(names_to = 'cohort')
  
  ## ready-to-use plot captions
  
  surv_globals$n_tags <- surv_globals$n_numbers %>% 
    mutate(plot_cap = paste0('total: n = ', n_total, 
                             ', events: n = ', n_events)) %>% 
    .$plot_cap %>% 
    set_names(surv_globals$n_numbers$cohort)

# CV folds used for tuning of the GLMNET models -----
  
  insert_msg('CV folds')
  
  set.seed(1234)
  
  surv_globals$n_rep <- 200
  
  surv_globals$folds <- 1:surv_globals$n_rep %>% 
    map(function(x) createFolds(y = factor(surv_globals$data$geo_pool$relapse), 
                                k = 10, 
                                list = FALSE, 
                                returnTrain = TRUE)) %>% 
    set_names(paste0('rep_', 1:surv_globals$n_rep))
  
# Clinical data for the expression data sets ------
  
  insert_msg('Clinical data for the expression data sets')
  
  surv_globals$clinic <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor')
  
  surv_globals$clinic <- 
    surv_globals$clinic[surv_globals$cohorts] %>%
    map(select, sample_id, gleason_simple, pt_stage) %>% 
    map(~filter(.x, complete.cases(.x)))
  
  ## some minimal wrangling: 
  ## pooling of pT stages

  surv_globals$clinic <- surv_globals$clinic %>% 
    map(mutate, 
        #psa_diagnosis = log(psa_diagnosis + 1), 
        #psa_diagnosis = zScores(psa_diagnosis), 
        pt_stage = car::recode(pt_stage, 
                              "'T3' = 'T3+'; 'T4' = 'T4+'"), 
        pt_stage = factor(pt_stage, c('T1', 'T2', 'T3+')))
  
  ## merging with the KLK3 expression levels

  surv_globals$psa <- list(geo_pool = geo_pool, 
                           tcga = tcga, 
                           dkfz = dkfz) %>% 
    map(~.x$combat[c('sample_id', 'KLK3')]) %>% 
    map(mutate, KLK3 = zScores(KLK3))
  
  surv_globals$clinic <- 
    map2(surv_globals$clinic, 
         surv_globals$psa, 
         left_join, by = 'sample_id')
  
# Tuning grid for the GBM, RF, and SVM models --------
  
  insert_msg('Tuning grid for the ML models')
  
  ## SVM 
  
  surv_globals$svm_grid <- 
    expand.grid(gamma.mu = c(0.0001, 0.0025, 0.005, 0.01, 
                             0.015, 0.02, 0.025, 0.05, 
                             seq(0.1, 1, by = 0.1)), 
                type = c('vanbelle1'), 
                diff.meth = c('makediff3'), 
                kernel = c('add_kernel'), 
                stringsAsFactors = FALSE)
  
  ## Random Forest
  
  surv_globals$rf_grid <- 
    expand.grid(mtry = seq(2, 27, by = 2), 
                splitrule = c('logrank', 'bs.gradient'), 
                nsplit = c(1, 2, 5), 
                nodesize = c(5, 10, 15, 20), 
                stringsAsFactors = FALSE)
  
  ## GBM
  
  surv_globals$gbm_grid <- 
    expand.grid(n.trees = c(500, 1000),
                shrinkage = c(0.01, 0.02, 0.03, 0.04, 0.05),
                interaction.depth = c(2, 3, 4),
                n.minobsinnode = c(2, 5))
  
# END -----

  surv_globals$survival <- NULL
  surv_globals$expression <- NULL
  
  surv_globals <- compact(surv_globals)
  
  insert_tail()