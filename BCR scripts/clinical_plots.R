# Graphical comparison of performance of two GBM models: 
# 1) with the clinical predictors of BCR-free survival
# 2) with the clinical and transcriptional predictors of BCR-free survival 

  insert_head()
  
# container -----
  
  surv_clplots <- list()
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## models of interest
  
  surv_clplots$models <- c("gbm_clinic", "gbm", "gbm_combi")
  
  ## plotting labels and colors 
  
  surv_clplots$algo_labels <- 
    c(reference = 'random prediction', 
      gbm = paste0(globals$algo_labels["gbm"], ', expression'), 
      globals$algo_labels[c("gbm_clinic", "gbm_combi")])
  
  surv_clplots$algo_colors <- 
    c(reference = 'gray60', 
      globals$algo_colors[surv_clplots$models])
  
  surv_clplots$algo_labeller <- surv_clplots$algo_labels %>% 
    stri_replace(fixed = ', ', replacement = '\n') %>% 
    set_names(names(surv_clplots$algo_labels))
  
  ## global performance stats 
  
  surv_clplots$stats <- 
    surv_summary$stats[surv_clplots$models] %>% 
    compress(names_to = 'algorithm') %>% 
    mutate(cohort_lab = globals$study_labels[cohort], 
           cohort_lab = ifelse(algorithm == 'gbm_clinic', 
                               cohort_lab, NA), 
           algorithm = factor(algorithm, 
                              surv_clplots$models))
  
  ## Brier scores for the unique survival time points

  surv_clplots$brier_scores <- 
    surv_summary$brier_scores[surv_clplots$models] %>% 
    map2(., 
         list(c('time', 'reference', 'training'), 
              c('time', 'training'), 
              c('time', 'training')), 
         function(x, y) x %>% 
           map(select, all_of(y))) %>% 
    transpose %>% 
    map(reduce, full_join, by = c('time')) %>% 
    map(set_names, c('time', 'reference', surv_clplots$models)) %>% 
    map(pivot_longer, 
        cols = all_of(c('reference', surv_clplots$models)), 
        names_to = 'algorithm', 
        values_to = 'brier_score') %>% 
    map(mutate, 
        algorithm = factor(algorithm, 
                           c('reference', surv_clplots$models)))
  
  ## inference stats for the normalized predictor scores
  
  surv_clplots$z_inferfence <- surv_summary$z_inference %>% 
    filter(algorithm %in% surv_clplots$models) %>% 
    select(algorithm, cohort, hr, hr_se, hr_lower, hr_upper) %>% 
    mutate(plot_label = paste0(signif(hr, 2), ' [', 
                               signif(hr_lower, 2), ' to ', 
                               signif(hr_upper, 2), ']'), 
           algorithm = factor(algorithm, surv_clplots$models), 
           cohort = factor(cohort, rev(c('geo_pool', 'tcga', 'dkfz'))))
  
  ## ROC AUC and data for ROC plots
  
  surv_clplots$roc_obj <- surv_roc$roc_obj
  
  for(i in names(surv_clplots$roc_obj)) {
    
    for(j in names(surv_clplots$roc_obj[[i]])) {
      
      surv_clplots$roc_obj[[i]][[j]] <- 
        surv_clplots$roc_obj[[i]][[j]] %>% 
        map(filter, 
            algorithm %in% surv_clplots$models) %>% 
        map(mutate, 
            algorithm = factor(algorithm, surv_clplots$models))

    }
    
  }
  
  surv_clplots$auc_data <- surv_clplots$roc_obj %>% 
    map(map, ~.x$stats) %>% 
    transpose %>% 
    map(reduce, rbind) %>%
    compress(names_to = 'cohort') %>% 
    mutate(cohort = factor(cohort, c('geo_pool', 'tcga', 'dkfz')), 
           algorithm = factor(algorithm, surv_clplots$models))
  
# Plots of the global performance stats -------
  
  insert_msg('Plots of the globals performance stats')
  
  surv_clplots$stat_plot <- surv_clplots$stats %>% 
    ggplot(aes(x = c_index, 
               y = 1 - ibs_model, 
               fill = algorithm)) + 
    geom_vline(xintercept = 0.5, 
               linetype = 'dashed') + 
    geom_hline(yintercept = 0.75, 
               linetype = 'dashed') + 
    geom_line(aes(group = cohort)) + 
    geom_point(shape = 21, 
               size = 2, 
               color = 'gray30') + 
    geom_text_repel(aes(label = cohort_lab), 
                    size = 2.5) + 
    scale_fill_manual(values = surv_clplots$algo_colors, 
                      labels = surv_clplots$algo_labels, 
                      name = '') + 
    globals$common_theme + 
    labs(title = 'GBM model performance', 
         x = 'C-index', 
         y = '1 - IBS')
  
# Plots of Brier scores for the unique time points ------
  
  insert_msg('Plots of Brier scores for the unique time points')
  
  surv_clplots$brier_plots <- 
    list(x = surv_clplots$brier_scores, 
         y = globals$study_labels[names(surv_clplots$brier_scores)]) %>% 
    pmap(function(x, y) x %>% 
           filter(complete.cases(.)) %>% 
           ggplot(aes(x = time, 
                      y = brier_score, 
                      color = algorithm)) + 
           geom_line() + 
           scale_color_manual(values = surv_clplots$algo_colors, 
                              labels = surv_clplots$algo_labels,
                              name = '') + 
           globals$common_theme + 
           labs(title = y, 
                x = 'min/max scaled BCR-free survival', 
                y = 'Brier score'))
  
# Forest plots of the inference stats ------
  
  insert_msg('Forest plots for harzard ratios')
  
  surv_clplots$inference_plot <- surv_clplots$z_inferfence %>%
    ggplot(aes(x = hr, 
               y = cohort, 
               color = algorithm)) + 
    facet_grid(algorithm ~ ., 
               labeller = as_labeller(surv_clplots$algo_labeller))+ 
    geom_vline(xintercept = 1, 
               linetype = 'dashed') + 
    geom_errorbarh(aes(xmin = hr_lower, 
                       xmax = hr_upper), 
                   height = 0) + 
    geom_point(shape = 16, 
               size = 2) + 
    geom_text(aes(label = plot_label), 
              size = 2.5, 
              hjust = 0.3, 
              vjust = -1) + 
    scale_color_manual(values = surv_clplots$algo_colors, 
                       labels = surv_clplots$algo_labels, 
                       name = '') + 
    scale_y_discrete(labels = globals$study_labels) + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'GBM model inference', 
         x = 'normalized HR \u00B1 95%CI')
  
# Summary plots of AUC over time --------
  
  insert_msg('Summary plots of AUC values over time')
  
  surv_clplots$auc_summary <- surv_clplots$auc_data %>% 
    ggplot(aes(x = predict_time, 
               y = auc, 
               color = algorithm)) + 
    facet_grid(. ~ cohort,
               labeller = as_labeller(globals$study_labels)) + 
    geom_hline(yintercept = 0.5, linetype = 'dashed') + 
    geom_line() + 
    geom_point(shape = 16, 
               size = 2) + 
    scale_color_manual(values = surv_clplots$algo_colors, 
                       labels = surv_clplots$algo_labels, 
                       name = '') + 
    globals$common_theme + 
    labs(title = 'GBM model, BCR prediction', 
         x = 'time after diagnosis, months', 
         y = 'ROC AUC')
  
# ROC plots for selected survival time points -------
  
  insert_msg('ROC plots for selected time points')
  
  for(i in names(surv_clplots$roc_obj)) {
    
    surv_clplots$roc_plots[[i]] <- 
      list(roc_obj = surv_clplots$roc_obj[[i]], 
           plot_title_suffix = globals$study_labels[names(surv_clplots$roc_obj[[i]])]) %>% 
      pmap(plot_surv_roc_times, 
           palette = surv_clplots$algo_colors, 
           labels = surv_clplots$algo_labels, 
           txt_size = 2.25, 
           txt_x = 0.27)
    
  }

# END ------
  
  surv_clplots$stats <- NULL
  surv_clplots$brier_scores <- NULL
  surv_clplots$z_inferfence <- NULL
  surv_clplots$roc_obj <- NULL
  
  surv_clplots <- compact(surv_clplots)
  
  insert_tail()