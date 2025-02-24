# R code for building Supplementary Figures

  insert_head()
  
# container ------
  
  suppl_fig <- list()
  
# Co-expression networks of the collagen-related genes -------
  
  insert_msg('Co-expression of the collagen-related genes')
  
  suppl_fig$nets <- 
    map2(net$graph_plots$gene_group[c("geo_pool", "tcga", "dkfz")], 
         c("geo_pool", "tcga", "dkfz"), 
         ~.x + 
           labs(title = globals$study_labels[.y]) + 
           theme(legend.position = 'none')) %>% 
    c(list(get_legend(net$graph_plots$gene_group[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'co_expression_networks_collagen_related_genes', 
              ref_name = 'nets', 
              caption = paste('Co-expression networks of the collagen-related', 
                              'transcripts in PCA tissue.'), 
              w = 180, 
              h = 180)
  
# Co-expression networks, vertex importance ------
  
  insert_msg('Co-expression networks, vertex importance')
  
  suppl_fig$hubs <- net$stat_plots[c("geo_pool", "tcga", "dkfz")] %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    c(list(get_legend(net$stat_plots$geo_pool))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'co_expression_networks_gene_importance', 
              ref_name = 'hubs', 
              caption = paste('Vertex importance statistics for co-expression', 
                              'networks of the collagen-related transcripts in', 
                              'PCA.'), 
              w = 180, 
              h = 180) 
  
# KM plots for the shared standalone markers of BCR-free survival --------
  
  insert_msg('KM plots, shared standalone BCR predictors')
  
  ## top 3 unfavorable standalone markers
  ## top 3 favorable standalone markers
  
  suppl_fig[c('km_unfavorable', 'km_favorable')] <- 
    list(uni_plots$km_plots[c("COL1A1", "COL11A1", "COL11A2")], 
         uni_plots$km_plots[c("COL4A6", "LAMB3", "DST")]) %>% 
    map(unlist, recursive = FALSE) %>% 
    map(map, 
        ~.x + 
          theme(legend.position = 'bottom', 
                legend.text = element_text(size = 7))) %>%
    map(~plot_grid(plotlist = .x, 
                   ncol = 3, 
                   align = 'hv', 
                   axis = 'tbl'))
  
  ## figure object
  
  suppl_fig[c('km_unfavorable', 'km_favorable')] <- 
    suppl_fig[c('km_unfavorable', 'km_favorable')] %>% 
    list(x =., 
         label = c('top_3_unfavorable_standalone_markers_collagen_genes', 
                   'top_3_favorable_standalone_markers_collagen_genes'), 
         ref_name = names(.), 
         caption = paste('Top three', 
                         c('unfavorable', 'favorable'), 
                         'standalone collagen-related transcriptional', 
                         'markers of BCR-free survival in PCA.')) %>% 
    pmap(as_figure, 
         w = 180, 
         h = 220)
  
# GBM models: the cohort confounder --------
  
  insert_msg('GBM, the cohort confounder')
  
  suppl_fig$gbm_cohort <- 
    surv_cohort$forest_plot %>% 
    as_figure(label = 'bcr_survival_gbm_model_cohort_confounder_pooled_geo', 
              ref_name = 'gbm_cohort', 
              caption = paste('Investigation of the confounding study', 
                              'effect on prediction of BCR-free survival', 
                              'by the GBM model in the pooled GEO cohort.'), 
              w = 140, 
              h = 110)
  
# GBM models: clinical and combined ones --------
  
  insert_msg('GBM models, clinical and combined predictors')
  
  ## upper panel: the performance stats
  
  suppl_fig$clinic$upper <- 
    plot_grid(surv_clplots$stat_plot, 
              ncol = 2, 
              rel_widths = c(0.7, 0.3))
  
  ## middle panel: Brier scores
  
  suppl_fig$clinic$middle <- surv_clplots$brier_plots %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              nrow = 1, 
              align = 'hv', 
              axis = 'tblr') %>% 
    plot_grid(get_legend(surv_clplots$brier_plots[[1]] + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.85, 0.15))
  
  ## bottom panel: the most important predictors 
  ## for the GBM model with combined clinical and expression 
  ## predictors
  
  suppl_fig$clinic$bottom <- 
    surv_plots$importance_plots$gbm_combi + 
    guides(x = guide_axis(angle = 45)) + 
    theme(strip.background.y = element_blank(), 
          strip.text.y = element_blank(),
          axis.text.x = element_markdown(size = 6.5))  + 
    labs(fill = 'explanatory\nfactor')
  
  ## the entire figure
  
  suppl_fig$clinic <- 
    plot_grid(suppl_fig$clinic$upper, 
              suppl_fig$clinic$middle, 
              suppl_fig$clinic$bottom, 
              nrow = 3, 
              rel_heights = c(1, 1, 1.2), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'bcr_survival_modeling_gbm_clinic_expression', 
              ref_name = 'clinic',
              caption = paste('Modeling of BCR-free survival by GBM algorithm', 
                              'with clinical predictors and mRNA expression of', 
                              'the collagen-related genes.'), 
              w = 180, 
              h = 230)
  
# Kaplan-Meier plots for the combined GBM and clinical GBM model ------
  
  insert_msg('KM plots for the GBM models')
  
  suppl_fig$cl_plot <- 
    list(surv_plots$km_plots$gbm_clinic %>% 
           map(~.x + 
                 labs(title = paste0(.x$labels$title, 
                                    ', clinical'))),
         surv_plots$km_plots$gbm_combi %>% 
           map(~.x + 
                 labs(title = paste0(.x$labels$title, 
                                     ', clinical + expression'))))%>% 
    transpose %>% 
    unlist(recursive = FALSE) %>% 
    map(~.x + 
          theme(plot.subtitle = element_blank(), 
                legend.position = 'bottom',
                legend.title = element_blank(), 
                legend.text = element_text(size = 7))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv')
  
  suppl_fig$cl_plot <- suppl_fig$cl_plot %>% 
    as_figure(label = 'bcr_survival_gbm_clinic_expression_km_plots', 
              ref_name = 'cl_plot', 
              caption = paste('GBM modeling of BCR-free survival in PCA with', 
                              'clinical prediction and expression of', 
                              'collagen-related transcripts:', 
                              'survival in tertiles of predictor scores.'), 
              w = 180, 
              h = 220)
  
# Modeling of overall survival in the GSE16560 cohort -------
  
  insert_msg('Modeling of overall survival')
  
  ## upper panel: globals performance stats
  
  suppl_fig$os$upper <- os_plots$algorithm_stat_plots %>% 
    map(~.x + 
          scale_x_continuous(limits = c(0.5, 0.9)) + 
          scale_y_continuous(limits = c(0.75, 1)) + 
          theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr')

  ## middle panel: Brier scores for the unique time points
  
  suppl_fig$os$middle <- os_plots$bs_plots$ridge %>% 
    map(~.x + 
          scale_color_manual(values = c('gray60', 
                                        globals$algo_colors[["ridge"]]), 
                             labels = c('random prediction', 
                                        globals$algo_labels[["ridge"]])) + 
          theme(plot.subtitle = element_blank()))
  
  suppl_fig$os$middle <- suppl_fig$os$middle %>% 
    map(~.x + 
          labs(title = paste0(.x$labels$title, 
                              globals$algo_labels[["ridge"]], 
                              sep = ', ')) +
          theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv',
              axis = 'tblr') %>% 
    plot_grid(get_legend(suppl_fig$os$middle[[1]] + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.85, 0.1))
  
  ## bottom panel: KM plots
  
  suppl_fig$os$bottom <- os_plots$km_plots$ridge %>% 
    map(~.x + 
          labs(title = paste0(.x$labels$title, 
                              globals$algo_labels[["ridge"]], 
                              sep = ', ')) + 
          theme(plot.subtitle = element_blank(), 
                legend.position = 'bottom', 
                legend.title = element_blank(), 
                legend.text = element_text(size = 7))) %>% 
    plot_grid(plotlist = .,
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr')
  
  ## the entire figure
  
  suppl_fig$os <- 
    plot_grid(suppl_fig$os$upper, 
              suppl_fig$os$middle, 
              suppl_fig$os$bottom, 
              nrow = 3, 
              rel_heights = c(1, 1, 1.25), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'os_survival_modeling_gse16560_cohort_collagen_genes', 
              ref_name = 'os', 
              caption = paste('Modeling of overall survival in PCA with', 
                              'expression levels of the collagen-related', 
                              'transcripts as explanatory factors.'), 
              w = 180, 
              h = 230)
  
# Saving the figures on the disc ------
  
  insert_msg('Saving the figures on the disc')
  
  suppl_fig %>% 
    number_figures(prefix = 'supplementary_figure_') %>% 
    walk(pickle, 
         path = './paper/supplementary figures', 
         device = cairo_pdf) 
  
# END ------
  
  insert_tail()