# Characteristic of the pooled GEO cohort used for survival modeling:
#
# 1) Comparison of demographic and clinical characteristic: factors that may
# impact on relapse and survival risk: PSA at diagnosis, tumor stage, node stage 
# and Gleason score. Age was available only for one cohort and hence excluded 
# from the analysis.
#
# 2) Comparison of overall relapse rates and follow-up time between the cohorts
# used to generate the pooled GEO collective
#
# 3) Comparison of relapse-free survival between the single cohorts.
#
# Clinical features are compared with Kruskal-Wallis test with eta-square effect 
# size statistic or Chi-square test with Cramer's V effect size statistic. 
# Survival is compared with Peto-Peto test. As a part of exploratory analysis, 
# no multiple tessting correction is applied to the p values.

  insert_head()
  
# container -------
  
  expl_pool <- list()
  
# analysis data -------
  
  insert_msg('Analysis data')
  
  ## cohort colors 
  
  expl_pool$study_colors <- brewer.pal(n = 4, name = 'Set1')

  ## variable lexicon
  
  expl_pool$lexicon <- globals$clinical_lexicon %>% 
    filter(variable %in% c('psa_diagnosis', 
                           'ct_stage', 
                           'pt_stage', 
                           'pn_stage', 
                           'gleason_sum', 
                           'gleason_simple', 
                           'relapse', 
                           'rfs_months')) %>% 
    mutate(variable = ifelse(variable == 'relapse', 
                             'relapse_factor', variable), 
           format = ifelse(variable == 'rfs_months', 'numeric', format), 
           test_type = ifelse(format == 'numeric', 
                              'kruskal_etasq', 'cramer_v'), 
           plot_type = ifelse(format == 'numeric', 
                              'box', 'stack'), 
           axis_lab = ifelse(format == 'factor', 
                             '% of cohort', 
                             ifelse(variable == 'rfs_months', 
                                    'months', 'ng/mL')))
  
  ## clinical data
  
  expl_pool$data <- 
    list("gse54460" = gse54460, 
         "gse70768" = gse70768, 
         "gse70769" = gse70769, 
         "gse220095" = gse220095) %>% 
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(select, sample_id, relapse, any_of(expl_pool$lexicon$variable)) %>% 
    map(mutate, 
        relapse_factor = car::recode(relapse, "0 = 'no'; 1 = 'yes'"), 
        relapse_factor = factor(relapse_factor, c('no', 'yes'))) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = factor(cohort, 
                           c('gse54460', 'gse70768', 'gse70769', 'gse220095')))
  
  ## cases with complete survival
  
  expl_pool$data <- expl_pool$data %>% 
    filter(!is.na(rfs_months), 
           !is.na(relapse))

# Descriptive stats --------
  
  insert_msg('Descriptive stats')
  
  expl_pool$stats <- expl_pool$data %>% 
    explore(variables = expl_pool$lexicon$variable, 
            split_factor = 'cohort', 
            what = 'table', 
            pub_styled = TRUE) %>% 
    format_desc
  
# Comparison of the clinical features --------
  
  insert_msg('Comparison of the clinical features')
  
  ## with a map loop, because some variables have only NAs in single cohorts
  
  expl_pool$test <-  
    map2_dfr(expl_pool$lexicon$variable, 
             expl_pool$lexicon$test_type, 
             function(x, y) expl_pool$data %>% 
               filter(!is.na(.data[[x]])) %>% 
               compare_variables(variables = x, 
                                 split_factor = 'cohort', 
                                 what = 'eff_size', 
                                 types = y, 
                                 exact = FALSE, 
                                 ci = FALSE, 
                                 pub_styled = TRUE)) %>% 
      mutate(plot_cap = paste(eff_size, significance, sep = ', '))

# Plots for the clinical features -------
  
  insert_msg('Plots for the clinical features')
  
  expl_pool$plots <- 
    list(variable = expl_pool$lexicon$variable, 
         plot_title = expl_pool$lexicon$label, 
         plot_subtitle = expl_pool$test$plot_cap, 
         y_lab = expl_pool$lexicon$axis_lab, 
         type = expl_pool$lexicon$plot_type) %>% 
    pmap(plot_variable, 
         expl_pool$data %>% 
           mutate(cohort = globals$study_labels[as.character(cohort)], 
                  cohort = unname(cohort), 
                  cohort = factor(cohort, levels = globals$study_labels), 
                  cohort = droplevels(cohort)), 
         split_factor = 'cohort', 
         scale = 'percent', 
         cust_theme = globals$common_theme, 
         x_n_labs = TRUE) %>% 
    set_names(expl_pool$lexicon$variable)
  
  ## styling, numeric plots
  
  expl_pool$plots[c("psa_diagnosis", "rfs_months")] <- 
    expl_pool$plots[c("psa_diagnosis", "rfs_months")] %>% 
    map2(c('log2', 'identity'), 
         ~.x + 
           scale_y_continuous(trans = .y) + 
           scale_fill_manual(values = unname(expl_pool$study_colors)))
  
  ## styling, stack plots
  
  expl_pool$plots[!names(expl_pool$plots) %in% c('psa_diagnosis', 'rfs_months')] <- 
    expl_pool$plots[!names(expl_pool$plots) %in% c('psa_diagnosis', 'rfs_months')] %>% 
    map(~.x + scale_fill_brewer(palette = 'Reds'))
  
# Surv-fit objects, median survival and comparison -------
  
  insert_msg('Survival stats and testing')
  
  expl_pool$surv_fit <- 
    survminer::surv_fit(Surv(rfs_months, relapse) ~ cohort, 
                        data = expl_pool$data)
  
  expl_pool$surv_stats <- expl_pool$surv_fit %>% 
    surv_median
  
  expl_pool$surv_test <- expl_pool$surv_fit %>% 
    surv_pvalue(method = 'S1') %>% 
    re_adjust(p_variable = 'pval', method = 'none')
  
# Kaplan-Meier plots -------
  
  insert_msg('Kaplan-Meier plot')
  
  ## plot labels: numbers of observations and events

  expl_pool$km_cap <- paste0('total: n = ', nrow(expl_pool$data), 
                             ', events: n = ', sum(expl_pool$data$relapse))
  
  expl_pool$legend_labs <- expl_pool$data %>% 
    blast(cohort) %>% 
    map2_chr(names(.), ., 
             ~paste0(globals$study_labels[.x], 
                     '\ntotal: n = ', nrow(.y), 
                     '\nevents: n = ', sum(.y$relapse)))
  
  ## plots
  
  expl_pool$km_plot <- 
    ggsurvplot(fit = expl_pool$surv_fit, 
               palette = expl_pool$study_colors, 
               pval = expl_pool$surv_test$significance, 
               pval.size = 2.75, 
               legend.labs = expl_pool$legend_labs, 
               legend.title = '')
  
  expl_pool$km_plot <- 
    expl_pool$km_plot$plot + 
    globals$common_theme + 
    labs(title = 'Pooled GEO collective, survival in the cohorts', 
         subtitle = expl_pool$km_cap, 
         x = 'Biochemical relapse-free survival, months')
  
# Common table with the results --------
  
  insert_msg('Common table with the results')
  
  expl_pool$result_tbl <- 
    left_join(expl_pool$stats, 
              expl_pool$test[c('variable', 'significance', 'eff_size')], 
              by = 'variable') %>% 
    format_summ_tbl(rm_n = FALSE) %>% 
    mutate(variable = exchange(variable, expl_pool$lexicon)) %>% 
    set_names(c('Variable', 
                globals$study_labels[levels(expl_pool$data$cohort)]), 
              'Significance', 'Effect size')
  
  ## clearing the cells without meaningfull data (NaN)
  
  expl_pool$result_tbl <- expl_pool$result_tbl %>% 
    map_dfc(~ifelse(stri_detect(.x, fixed = 'NaN'), '', .x))
  
# appending the result table with the N numbers ------
  
  insert_msg('Appending the result table with the N numbers')
  
  expl_pool$n_numbers <- expl_pool$data %>% 
    count(cohort) %>% 
    mutate(cohort = globals$study_labels[as.character(cohort)]) %>% 
    column_to_rownames('cohort') %>% 
    t %>% 
    as_tibble %>% 
    mutate(Variable = 'Cancer samples, N') %>% 
    relocate(Variable)
  
  expl_pool$result_tbl <- 
    full_rbind(expl_pool$n_numbers, 
               expl_pool$result_tbl)
  
# Appending the results with p values of Peto-Peto test for BCR-free survival --------
  
  insert_msg('Appending the results with BCR-free survival testing')

  expl_pool$result_tbl <- expl_pool$result_tbl %>% 
    mutate(Significance = ifelse(stri_detect(Variable, fixed = 'survival'), 
                                 expl_pool$surv_test$significance[[1]], 
                                 Significance), 
           `Effect size` = ifelse(stri_detect(Variable, fixed = 'survival'), 
                                 NA, `Effect size`))
    
# END -------
  
  expl_pool$data <- NULL
  
  expl_pool <- compact(expl_pool)
  
  insert_tail()