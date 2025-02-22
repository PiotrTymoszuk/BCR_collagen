# Comparison of the training and test portion of the GS16560 cohort: clinical 
# factors and overall survival. Only the tabular form: this table will be shown 
# in the revised manuscript. 
#
# Statistical significance: factors - chi-square test  with Cramer's V effect 
# size statistic, numerics: Mann-Whitney test with r effect size statistic, 
# survival: Peto-Peto test.

  insert_head()
  
# container -------
  
  expl_os <- list()
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## variable lexicon
  
  expl_os$lexicon <- globals$clinical_lexicon %>% 
    filter(variable %in% c('age', 
                           'gleason_sum', 
                           'gleason_simple', 
                           'death', 
                           'os_months')) %>% 
    mutate(variable = ifelse(variable == 'death', 
                             'death_factor', variable), 
           format = ifelse(variable == 'os_months', 'numeric', format), 
           test_type = ifelse(format == 'numeric', 
                              'kruskal_etasq', 'wilcoxon_r'), 
           plot_type = ifelse(format == 'numeric', 
                              'box', 'stack'), 
           axis_lab = ifelse(format == 'factor', 
                             '% of cohort', 
                             ifelse(variable == 'os_months', 
                                    'months', 'ng/mL')))
  
  ## analysis data
  
  expl_os$data <- gse16560$clinic %>% 
    filter(tissue_type == 'tumor') %>% 
    mutate(death_factor = car::recode(death, "0 = 'no'; 1 = 'yes'"), 
           death_factor = factor(death_factor, c('no', 'yes')), 
           observation = 1:nrow(.), 
           subset = ifelse(observation %in% gse16560$train_idx, 
                           'training', 'test'), 
           subset = factor(subset, c('training', 'test'))) %>% 
    select(sample_id, subset, death, all_of(expl_os$lexicon$variable))
  
# Descrptive stats -------
  
  insert_msg('Descriptive stats')
  
  expl_os$stats <- expl_os$data %>% 
    explore(split_factor = 'subset', 
            variables = expl_os$lexicon$variable, 
            what = 'table', 
            pub_styled = TRUE) %>% 
    format_desc
  
# Testing for differences -------
  
  insert_msg('Testing for differences between the modeling subsets')
  
  expl_os$test <- expl_os$data %>%
    compare_variables(variables = expl_os$lexicon$variable, 
                      split_factor = 'subset', 
                      what = 'eff_size', 
                      types = expl_os$lexicon$test_type, 
                      exact = FALSE, 
                      ci = FALSE, 
                      pub_styled = TRUE, 
                      adj_method = 'BH')
  
# comparison of the overall survival -------
  
  insert_msg('Comparison of the overall survival')
  
  expl_os$surv_fit <- 
    survminer::surv_fit(Surv(os_months, death) ~ subset, 
                        data = expl_os$data)
  
  expl_os$surv_stats <- expl_os$surv_fit %>% 
    surv_median
  
  expl_os$surv_test <- expl_os$surv_fit %>% 
    surv_pvalue(method = 'S1') %>% 
    re_adjust(p_variable = 'pval', method = 'none')
  
# Result table -------
  
  insert_msg('Result table')

  expl_os$result_tbl <- 
    left_join(expl_os$stats, 
              expl_os$test[c('variable', 'significance', 'eff_size')], 
              by = 'variable') %>% 
    format_summ_tbl(rm_n = FALSE) %>% 
    mutate(variable = exchange(variable, expl_os$lexicon)) %>% 
    set_names(c('Variable', 
                'Training', 'Test', 
                'Significance', 'Effect size'))
    
# N numbers in the training and test subsets, appending the results table ------
  
  insert_msg('N numbers')
  
  expl_os$n_numbers <- expl_os$data %>% 
    count(subset) %>% 
    column_to_rownames('subset') %>% 
    t %>% 
    as_tibble %>% 
    mutate(Variable = 'Cancer samples, N') %>% 
    relocate(Variable) %>% 
    set_names(c('Variable', 'Training', 'Test'))
  
  expl_os$result_tbl <- 
    full_rbind(expl_os$n_numbers, 
               expl_os$result_tbl)
  
# Appending the result table with the survival p values ------
  
  insert_msg('Survival p value to the result table')
  
  expl_os$result_tbl <- expl_os$result_tbl %>% 
    mutate(Significance = ifelse(stri_detect(Variable, fixed = 'survival'), 
                                 expl_os$surv_test$significance, 
                                 Significance), 
           `Effect size` = ifelse(stri_detect(Variable, fixed = 'survival'), 
                                  NA, `Effect size`))
  
# END ------
  
  expl_os$data <- NULL
  
  expl_os <- compact(expl_os)
  
  insert_tail()