# Figures for presentation of the revision results


# container ------

  pres <- list()

# Network plots: the pooled GEO cohort ---------
  
  insert_msg('Network analysis')
  
  ## network structure
  
  pres$net_plots <- 
    net$graph_plots$gene_group[c("geo_pool", "tcga", "dkfz")] %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    c(list(get_legend(net$graph_plots$gene_group[[1]]))) %>% 
    plot_grid(plotlist = .,
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'collagen_co_expression_network_plots', 
              ref_name = 'net_plots', 
              w = 210, 
              h = 210)
  
  ## vertex importance stats
  
  pres$net_stats <- net$stat_plots$geo_pool %>% 
    as_figure(label = 'collagen_coexpression_hub_genes', 
              ref_name = 'net_stats', 
              w = 180, 
              h = 150)

# Multi-parameter model of BCR, GBM ---------
  
  insert_msg('Multi-paramater modeling of survival GBM')
  
  ## performance stats, Brier scores and KM plots
  
  pres$bcr_gbm <- surv_plots$algorithm_stat_plots %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    c(list(get_legend(surv_plots$algorithm_stat_plots))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'bcr_free_survival_gbm_model', 
              ref_name = 'bcr_gbm', 
              w = 180, 
              h = 180)
  
  ## KM plots 
  
  pres$bcr_gbm_km <- surv_plots$km_plots$gbm %>% 
    map(~.x + theme(legend.title = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'bcr_free_survival_gbm_model_kaplan_meier', 
              ref_name = 'bcr_gbm_km', 
              w = 190, 
              h = 190)
  
  ## top most important variables
  
  pres$bcr_gbm_vimp <- surv_plots$importance_plots$gbm %>% 
    as_figure(label = 'bcr_free_survival_gbm_model_variable_importance', 
              ref_name = 'bcr_gbm_vimp', 
              w = 100, 
              h = 160) 

# Multi-parameter model of BCR, clinical and combined models ------
  
  insert_msg('GBM models: clinics and the combined predictor set')
  
  ## performance stats
  
  pres$bcr_clinic <- surv_clplots$stat_plot %>% 
    as_figure(label = 'bcr_free_survival_gbm_model_clinic', 
              ref_name = 'bcr_clinic', 
              w = 160, 
              h = 100)
  
  ## KM plots
  
  pres$bcr_clinic_km <- surv_plots$km_plots[c("gbm_clinic", "gbm_combi")] %>% 
    transpose %>% 
    unlist(recursive = FALSE) %>% 
    map(~.x + theme(legend.title = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'bcr_free_survival_gbm_model_clinic_kaplan_meier', 
              ref_name = 'bcr_clinic_km', 
              w = 190, 
              h = 270)
  
  ## top most important variables
  
  pres$bcr_gbm_vimp <- surv_plots$importance_plots$gbm_combi %>% 
    as_figure(label = 'bcr_free_survival_gbm_model_clinic_variable_importance', 
              ref_name = 'bcr_gbm_clinic_vimp', 
              w = 100, 
              h = 160) 
  
# OS model: RIDGE in the GSE16560 cohort -------
  
  insert_msg('Performance of the RIDGE model of OS')
  
  pres$os_ridge <- ridge_os$calibration %>% 
    map(plot, 
        palette = unname(globals$tertile_colors), 
        cust_theme = globals$common_theme, 
        show_cox = FALSE) %>% 
    map2(., paste('GSE16560,', c('training subset', 'test subset')), 
         ~.x + 
           theme(legend.position = 'bottom') + 
           labs(title = .y, 
                x = 'overall survival, months')) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr')
  
  pres$os_ridge <- pres$os_ridge %>% 
    as_figure(label = 'os_modeling_gse16560_best_ridge_model', 
              ref_name = 'os_ridge', 
              w = 180, 
              h = 100)
  
# Saving on the disc --------
  
  insert_msg('Saving on the disc')
  
  pres %>% 
    number_figures %>% 
    walk(pickle, 
         path = './paper/presentation', 
         device = cairo_pdf)
  
# END --------

  