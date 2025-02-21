# Builds figures for the main manuscript text

  insert_head()
  
# container ------
  
  fig <- list()
  
# univariable survival analysis, collagen-related genes -------
  
  insert_msg('Univariable survival analysis, transcriptome')
  
  ## Forest plots with HR +/- 95%CI and dot plots for C-indexes
  ## leaving some place for the protein data, if intended
  
  fig$uni_mrna <- 
    plot_grid(uni_plots$forest_plot + 
                labs(title = 'Standalone markers, BCR-free survival, mRNA', 
                     subtitle = 'model inference') + 
                theme(legend.position = 'none', 
                      axis.title.y = element_blank()), 
              uni_plots$c_index_plot + 
                labs(title = '', 
                     subtitle = 'performance statistics') + 
                theme(legend.position = 'none', 
                      axis.title.y = element_blank(), 
                      axis.text.y = element_blank()), 
              ncol = 2, 
              rel_widths = c(1.26, 1), 
              align = 'h', 
              axis = 'tblr') %>% 
    plot_grid(get_legend(uni_plots$c_index_plot + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.9, 0.1)) 
  
  fig$uni_mrna <- fig$uni_mrna %>% 
    as_figure(label = 'univariable_bcr_survival_collagen_mrna', 
              ref_name = 'uni_mrna', 
              caption = paste('Univariable analysis of BCR-free survival with', 
                              'collagen-related transcripts as candidate', 
                              'standalone prognostic factors.'), 
              w = 180, 
              h = 145)
  
# ML models of BCR-free survival, transcriptome --------
  
  insert_msg('BCR-free survivalML  models, transcriptome')
  
  ## showing the following: 
  ## the upper panel: global performance stats
  ## middle one: KM plots for survival in the score's tertiles
  ## bottom panel: variable importance
  
  fig$ml_mrna$upper <- surv_plots$algorithm_stat_plots %>% 
    map(~.x + 
          theme(plot.subtitle = element_blank(), 
                legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              axis = 'tblr') 
  
  fig$ml_mrna$middle <- surv_plots$km_plots$gbm %>% 
    map(~.x + 
          guides(color = guide_legend(nrow = 2)) + 
          theme(plot.subtitle = element_blank(), 
                legend.position = 'bottom', 
                legend.title = element_blank(), 
                legend.text = element_text(size = 7))) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              axis = 'tblr')
  
  fig$ml_mrna$bottom <- surv_plots$importance_plots$gbm + 
    guides(x = guide_axis(angle = 45)) + 
    theme(axis.text.x = element_markdown(size = 6.5), 
          strip.background.y = element_blank(), 
          strip.text.y = element_blank())
  
  ## the entire figure
  
  fig$ml_mrna <- 
    plot_grid(fig$ml_mrna$upper, 
              fig$ml_mrna$middle, 
              fig$ml_mrna$bottom, 
              nrow = 3, 
              rel_heights = c(1, 1.4, 1.2), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'ml_modeling_bcr_collagen_mrna', 
              ref_name = 'ml_mrna', 
              caption = paste('Machine learning modeling of BCR-free survival', 
                              'with expression levels of the collagen-related', 
                              'transcripts as explanatory factors.'), 
              w = 180, 
              h = 230)
  
# Saving the figures on the disc -------
  
  insert_msg('Saving the figures on the disc')
  
  fig %>% 
    number_figures %>% 
    walk(pickle, 
         path = './paper/figures', 
         device = cairo_pdf)
  
# END -----
  
  insert_tail()
  