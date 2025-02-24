# Testing for the confounding effects of the cohorts in GBM survival modeling 
#
# Three 'classical' Cox models of biochemical relapse-free survival are 
# constructed in the pooled GEO cohort:
# 
# 1) a cohort only model
#
# 2) a model with the GBM predictor score as the sole explanatory variable
#
# 3) an additive model with the cohort and GBM predictor scores as explanatory 
# variables
#
# We're checking inference for the GBM term in the univariable and in the 
# additive model. Contribution of the cohort and GBM score terms to the additive 
# model is checked by LRT (likelihood ratio test). Finally, we're comparing 
# fit stats such as integrated Brier scores and C-indexes between 
# these 3 models.

  insert_head()

# container ------
  
  surv_cohort <- list()

# analysis data -------
  
  insert_msg('Analysis data')
  
  ## cohort assignment, cohort confounder and the linear predictor scores
  
  surv_cohort$cohorts <- 
    geo_pool$clinic[, c("sample_id", "study")]

  surv_cohort$data <- 
    inner_join(gbm_surv$score_tbl$geo_pool, 
               surv_cohort$cohorts, 
               by = 'sample_id') %>% 
    mutate(gbm_score = zScores(gbm_score))
  
  ## numbers of cases and relapses
  
  surv_cohort$n_numbers <- 
    paste0('total: n = ', surv_cohort$data, 
           ', events: n = ', sum(surv_cohort$data$relapse))
  
# Differences in predictor score between the cohorts -------
  
  insert_msg('Differences in predictor score between the cohorts')
  
  ## descriptive stats
  
  surv_cohort$cohort_stats <- surv_cohort$data %>%
    fast_num_stats(split_fct = 'study', 
                   variables = 'gbm_score')

  ## normality and EOV in respect to the cohort splitting factor
  ## are violated, hence Kruskal-Wallis test
  
  surv_cohort$cohort_test <- 
    f_kruskal_test(surv_cohort$data$gbm_score, 
                   f = surv_cohort$data$study, 
                   as_data_frame = TRUE) %>% 
    re_adjust %>% 
    mutate(eff_size = paste('\u03B7\u00B2 =', signif(etasq, 2)), 
           plot_cap = paste(eff_size, significance)) %>% 
    as_tibble

  ## box plot of th predictor scores in the cohorts
  
  surv_cohort$cohort_plot <- surv_cohort$data %>% 
    mutate(study = globals$study_labels[as.character(study)], 
           study = unname(study), 
           study = factor(study, globals$study_labels)) %>% 
    plot_variable(variable = 'gbm_score', 
                  split_factor = 'study', 
                  type = 'box', 
                  point_hjitter = 0, 
                  plot_title = 'GBM predictor score', 
                  plot_subtitle = surv_cohort$cohort_test, x_lab = '', 
                  y_lab = 'GBM predictor, Z-score', 
                  cust_theme = globals$common_theme, 
                  x_n_labs = TRUE) + 
    scale_fill_brewer(palette = 'Set1')
  
# Construction of Cox models -------
  
  insert_msg('Construction of Cox models')
  
  surv_cohort$models <- 
    list(cohort = Surv(scaled_rfs_months, relapse) ~ study, 
         gbm_score = Surv(scaled_rfs_months, relapse) ~ gbm_score, 
         full = Surv(scaled_rfs_months, relapse) ~ study + gbm_score) %>% 
    map(~call2('coxph', 
               formula = .x, 
               data = surv_cohort$data, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map(as_coxex, surv_cohort$data)
  
# Model assumptions, fit stats and inference ------
  
  insert_msg('Model assumption, fit stats and inference')
  
  ## the proportional hazard assumption is clearly violated
  ## but as long as the effects off the GBM model can be reproduced 
  ## in the validation collectives this is not a big problem
  
  surv_cohort$assumptions <- surv_cohort$models %>% 
    map(summary, 'assumptions')
  
  surv_cohort$stats <- surv_cohort$models %>% 
    map(summary, 'fit')
  
  surv_cohort$inference <- surv_cohort$models %>% 
    map(summary, 'inference') %>% 
    map(mutate, 
        hr = exp(estimate), 
        hr_lower = exp(lower_ci), 
        hr_upper = exp(upper_ci))
  
  surv_cohort[c("stats", "inference")] <- 
    surv_cohort[c("stats", "inference")] %>% 
    map(compress, names_to = 'model_type') %>% 
    map(mutate, 
        model_type = factor(model_type, c('cohort', 'gbm_score', 'full')))
  
# Likelihood ration test (LRT) --------
  
  insert_msg('Likelihood ratio test')
  
  surv_cohort$lrt <- 
    list(cohort = surv_cohort$models[c("full", "gbm_score")], 
         gbm_score = surv_cohort$models[c("full", "cohort")]) %>% 
    map(map, as_coxph) %>% 
    map(~anova(.x[[1]], .x[[2]])) %>% 
    map(as_tibble) %>% 
    map(~.x[2, ]) %>% 
    compress(names_to = 'model_type')
  
  surv_cohort$lrt <- surv_cohort$lrt %>% 
    mutate(p_value = `Pr(>|Chi|)`) %>% 
    re_adjust(method = 'none') %>% 
    mutate(plot_cap = paste0('\u03C7\u00B2(', Df, ') = ', signif(Chisq, 2), 
                             '\n', significance), 
           plot_cap = paste(globals$cohort_model_labels[model_type], 
                            plot_cap, 
                            sep = '\n'))
  
# Plots of performance stats of the models -------
  
  insert_msg('Plot of performance stats')
  
  surv_cohort$stat_plot <- surv_cohort$stats %>% 
    plot_surv_stats(label_variable = 'model_type', 
                    color_variable = 'model_type', 
                    palette = unname(globals$cohort_model_colors), 
                    labels = globals$cohort_model_labels, 
                    plot_title = 'Study confounder, pooled GEO cohort', 
                    plot_subtitle = 'GBM survival modeling')
  
# Forest plot of the HRs -------
  
  insert_msg('Forest plot of the HR')

  ## LRT results to be shown in the Forest plot
  
  surv_cohort$forest_labs <- set_names(surv_cohort$lrt$plot_cap, 
                                       surv_cohort$lrt$model_type)
  
  surv_cohort$forest_labs['full'] <- 
    globals$cohort_model_labels["full"] %>% 
    stri_replace(fixed = ' + ', replacement = '\n+ ')
  
  ## LRT test results are shown in the plot facets
  
  surv_cohort$forest_plot <- surv_cohort$inference %>% 
    mutate(ax_lab = ifelse(variable == 'study', 
                           globals$study_labels[level], 
                           'GBM, Z-score'),
           ax_lab = ifelse(variable == 'study', 
                           paste(ax_lab, n, sep = '\nn = '), 
                           ax_lab), 
           hr_lab = paste0(signif(hr, 2), ' [96%CI: ', 
                           signif(hr_lower, 2), ' to ', 
                           signif(hr_upper, 2), ']')) %>% 
    ggplot(aes(x = hr, 
               y = ax_lab, 
               color = variable)) + 
    geom_vline(xintercept = 1, 
               linetype = 'dashed') + 
    geom_errorbarh(aes(xmin = hr_lower, 
                       xmax = hr_upper), 
                   height = 0) + 
    geom_point(shape = 16, 
               size = 2) + 
    geom_text(aes(label = hr_lab), 
              hjust = 0.2, 
              vjust = -1, 
              size = 2.5) +
    facet_grid(model_type ~ ., 
               scales = 'free', 
               space = 'free', 
               labeller = as_labeller(surv_cohort$forest_labs)) + 
    scale_color_manual(values = globals$cohort_model_colors) + 
    guides(color = 'none') + 
    globals$common_theme + 
    theme(axis.title.y = element_blank(), 
          strip.text.y = element_text(angle = 0, 
                                      hjust = 0)) + 
    labs(title = 'Study confounder, pooled GEO cohort', 
         subtitle = 'Inference and likelihood ratio test', 
         x = 'HR, 95% CI')
  
# END --------
  
  surv_cohort$data <- NULL
  surv_cohort$cohorts <- NULL
  
  surv_cohort <- compact(surv_cohort)
  
  insert_tail()