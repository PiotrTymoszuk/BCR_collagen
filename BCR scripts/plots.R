# Plots for the results of multi-parameter survival modeling

  insert_head()
  
# container ------
  
  surv_plots <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# Plotting globals --------
  
  insert_msg('Plotting globals')
  
  ## names of models with the expression only
  
  surv_plots$expression_models <- 
    c("ridge", "elnet", "lasso", "svm", "rf", "gbm")
  
# Tuning process -------
  
  insert_msg('Plots of the tuning process')
  
  ## for the regularized Cox models
  
  surv_plots$tuning[c('ridge', 'elnet', 'lasso')] <- 
    list(ridge_surv, elnet_surv, lasso_surv) %>% 
    map(~.x$lambda_tune) %>% 
    map(plot)
    
  surv_plots$tuning[c('ridge', 'elnet', 'lasso')] <- 
    list(x = surv_plots$tuning[c('ridge', 'elnet', 'lasso')], 
         y = globals$algo_labels[c("ridge", "elnet", "lasso")]) %>% 
    pmap(function(x, y) x + 
           globals$common_theme + 
           labs(title = y, 
                y = 'deviance, repeated CV'))

  ## for SVM: 
  
  surv_plots$tuning$svm <- svm_surv$tuning$summary %>% 
    mutate(best = ifelse(c_index == max(c_index), 
                         'yes', 'no')) %>% 
    ggplot(aes(x = gamma.mu, 
               y = c_index, 
               fill = best)) + 
    geom_point(shape = 21, 
               size = 2) + 
    scale_x_continuous(trans = 'log', 
                       labels = function(x) signif(x, 2)) + 
    scale_fill_manual(values = c(no = 'steelblue', 
                                 yes = 'coral3'), 
                      name = 'best tune') + 
    globals$common_theme + 
    labs(title = globals$algo_labels["svm"], 
         subtitle = paste('Best tune: \u03B3 =', 
                          signif(svm_surv$tuning$best_tune$gamma.mu, 3)), 
         x = expression(gamma),
         y = 'C-index, CV')
  
  ## for RF: 
  
  surv_plots$tuning$rf <- rf_surv$tuning$summary %>% 
    mutate(best = ifelse(c_index == max(c_index), 
                         'yes', 'no')) %>% 
    ggplot(aes(x = mtry, 
               y = c_index, 
               color = splitrule, 
               fill = best)) + 
    geom_line(aes(group = splitrule)) + 
    geom_point(shape = 21, 
               size = 1, 
               color = 'black') + 
    scale_fill_manual(values = c(no = 'steelblue', 
                                yes = 'coral3'), 
                     name = 'best tune') + 
    facet_grid(nsplit ~ nodesize) + 
    globals$common_theme + 
    labs(title = globals$algo_labels["rf"], 
         subtitle = paste0('Best tune: ', 
                           'mtry = ', signif(rf_surv$tuning$best_tune$mtry[1], 2), 
                           ', splitrule: ', rf_surv$tuning$best_tune$splitrule[1], 
                           ', n splits = ', rf_surv$tuning$best_tune$nsplit[1], 
                           ', node size = ', rf_surv$tuning$best_tune$nodesize[1]), 
         x = 'mtry', 
         y = 'C-index, CV')
  
  ## for the GBM models
  
  surv_plots$tuning[c('gbm', 'gbm_clinic', 'gbm_combi')] <- 
    list(gbm_surv, surv_clin, surv_combi) %>% 
    map(~.x$tuning$summary) %>% 
    map(mutate, 
        best = ifelse(cv_deviance == min(cv_deviance), 
                      'yes', 'no'))
    
  surv_plots$tuning[c('gbm', 'gbm_clinic', 'gbm_combi')] <- 
    list(x = surv_plots$tuning[c('gbm', 'gbm_clinic', 'gbm_combi')], 
         y = globals$algo_labels[c("gbm", "gbm_clinic", "gbm_combi")], 
         z = list(gbm_surv, surv_clin, surv_combi) %>% 
           map(~.x$tuning$best_tune)) %>% 
    pmap(function(x, y, z) x %>% 
           ggplot(aes(x = shrinkage, 
                      y = cv_deviance, 
                      color = factor(n.trees), 
                      fill = best)) + 
           geom_line(aes(group = factor(n.trees))) + 
           geom_point(shape = 21, 
                      size = 1, 
                      color = 'black') + 
           scale_fill_manual(values = c(no = 'steelblue', 
                                        yes = 'coral3'), 
                             name = 'best tune') + 
           facet_grid(interaction.depth ~ n.minobsinnode) + 
           globals$common_theme + 
           labs(title = y, 
                subtitle = paste0('Best tune: ', 
                                  'N trees = ', z$n.trees[1], 
                                  ', shrinkage = ', z$shrinkage[1], 
                                  ', interaction depth = ', z$interaction.depth[1], 
                                  ', minimal node size = ', z$n.minobsinnode[1]), 
                x = 'N trees', 
                y = 'Deviance, CV'))

# Performance stats for the algorithms: C-index and IBS, expression models ------
  
  insert_msg('Performance of the algorithms: C-index and IBS, expression models')
  
  ## comparison of the cohorts
  
  surv_plots$stat_plots <- 
    list(stats = surv_summary$stats[surv_plots$expression_models], 
         plot_subtitle = globals$algo_labels[surv_plots$expression_models] %>% 
           paste('algorithm')) %>% 
    pmap(plot_surv_stats)
  
  ## comparison of the algorithms within single cohorts, expression-only models
  
  surv_plots$algorithm_stat_plots <- 
    surv_summary$stats[surv_plots$expression_models] %>% 
    compress(names_to = 'algorithm') %>% 
    mutate(cohort = factor(cohort, surv_summary$stats[[1]]$cohort)) %>% 
    blast(cohort) %>% 
    list(stats = ., 
         plot_title = globals$study_labels[names(.)]) %>% 
    pmap(plot_surv_stats, 
         palette = set_names(globals$algo_colors, 
                             unname(globals$algo_labels)), 
         labels = globals$algo_labels, 
         label_variable = 'algorithm', 
         color_variable = 'algorithm', 
         plot_subtitle = NULL, 
         txt_size = 2.3, 
         box.padding = 0.5, 
         force = 2) 
  
# Brier scores for the unique time points -------
  
  insert_msg('Brier score for the unique time points')
  
  for(i in names(surv_summary$brier_scores)) {
    
    surv_plots$bs_plots[[i]] <- surv_summary$brier_scores[[i]] %>% 
      map(plot, cust_theme = globals$common_theme)
    
    surv_plots$bs_plots[[i]] <- 
      list(x = surv_plots$bs_plots[[i]], 
           y = globals$study_labels[names(surv_plots$bs_plots[[i]])]) %>% 
      pmap(function(x, y) x + 
             labs(title = y, 
                  subtitle = paste(globals$algo_labels[[i]], 'algorithm'), 
                  x = 'min/max scaled BCR-free survival') + 
             scale_color_manual(values = c(reference = 'gray60', 
                                           training = 'coral3'), 
                                labels = c(reference = 'random prediction', 
                                           training = 'model'), 
                                name = ''))
    
  }

# Survival in the score tertiles -------
  
  insert_msg('Kaplan-Meier plots for survival in the score tertiles')
  
  for(i in names(surv_summary$tertile_fits)) {
    
    surv_plots$km_plots[[i]] <- 
      list(fit = surv_summary$tertile_fits[[i]], 
           n_numbers = surv_summary$tertile_n[[i]], 
           p_value = surv_summary$tertile_test[[i]]$significance, 
           plot_title = globals$study_labels[names(surv_summary$tertile_fits[[i]])]) %>% 
      pmap(plot_tertile_km, 
           x_lab = 'min/max scaled BCR-free survival')
    
  }
  
# Variable importance -------
  
  insert_msg('Variable importance')
  
  surv_plots$importance_plots <- 
    list(data = surv_summary$importance, 
         imp_stat = c(rep('coef', 3), 
                      rep('delta_c_index', 2), 
                      rep('rel_influence', 3)), 
         labeller = c(rep(list(c(positive = 'top unfavorable', 
                                 negative = 'top favorable')), 
                          3), 
                      rep(list(c(positive = 'top improvement', 
                                 negative = 'top worsening')), 
                          5)), 
         plot_title = paste0(globals$algo_labels[names(surv_summary$importance)], 
                               ', variable importance'), 
         x_lab = c(list(expression('log HR'[Ridge]), 
                        expression('log HR'[ElasticNet]), 
                        expression('log HR'[LASSO])), 
                   rep(list(expression(Delta * ' C-index')), 
                       2), 
                   rep(list(expression(Delta * ' SSE')), 
                       3)), 
         n_top = c(rep(20, 5), rep(100, 3)), 
         palette = c(rep(list(c(positive = 'firebrick', 
                                negative = 'steelblue', 
                                ns = 'gray70')), 
                         5), 
                     rep(list(c(positive = 'steelblue', 
                                negative = 'steelblue', 
                                ns = 'gray70', 
                                clinical = 'aquamarine3', 
                                `collagen-related\ntranscript` = 'steelblue')), 
                         3))) %>% 
    pmap(plot_surv_importance, 
         form = 'bar', 
         flip = TRUE, 
         split_regulation = FALSE) %>% 
    map(~.x + 
          theme(plot.subtitle = element_blank()) + 
          geom_vline(xintercept = 0,
                     linetype = 'dashed') + 
          theme(legend.position = 'bottom'))
  
# Error structure: plots of square errors of prediction -------
  
  insert_msg('Plots of square errors of prediction')
  
  ## the square errors are available for Cox-like models
  
  surv_plots$square_plots <- list(ridge = ridge_surv, 
                                  elnet = elnet_surv, 
                                  lasso = lasso_surv, 
                                  gbm = gbm_surv, 
                                  gbm_clinic = surv_clin, 
                                  gbm_combi = surv_combi) %>% 
    map(~.x$calibration) %>% 
    future_map(map, 
        plot, 
        type = 'squares', 
        .options = furrr_options(seed = TRUE))
  
  ## plot titles and styling
  ## transposition in a more handy format
  
  for(i in names(surv_plots$square_plots)) {
    
    surv_plots$square_plots[[i]] <- 
      list(x = surv_plots$square_plots[[i]], 
           y = paste(globals$algo_labels[[i]], 
                     globals$study_labels[names(surv_plots$square_plots[[i]])], 
                     sep = ', ')) %>% 
      pmap(function(x, y) x %>% 
             map(~.x + 
                   labs(title = y) + 
                   globals$common_theme + 
                   theme(panel.grid.major.x = element_blank())))
    
  }
  
  surv_plots$square_plots <- surv_plots$square_plots %>% 
    map(transpose) %>%
    transpose
  
  for(i in names(surv_plots$square_plots$time)) {
    
    surv_plots$square_plots$time[[i]] <- 
      surv_plots$square_plots$time[[i]] %>% 
      map(~.x + 
            labs(x = 'min/max scaled BCR-free survival'))
    
  }
  
  surv_plots$square_plots$observation <- 
    surv_plots$square_plots$observation %>% 
    map(map, 
        ~.x + 
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank()) + 
          geom_hline(yintercept = 0.25, linetype = 'dashed'))
  
# END -------
  
  plan('sequential')
  
  rm(i)
  
  insert_tail()