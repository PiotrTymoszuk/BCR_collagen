# characteristic of the investigates cohorts

  insert_head()
  
# container ------
  
  cohorts <- list()
  
# globals ------
  
  insert_msg('Globals')
  
  ## variables
  
  cohorts$lexicon <- globals$clinical_lexicon

  ## data set-specific variable lists
  
  cohorts$data <- globals$study_exprs %>%
    eval %>% 
    map(~.x$clinic) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(select, any_of(cohorts$lexicon$variable))

  ## analysis tables
  ## wrangling for a common variable format
  ## stages are reduced to the main ones
  
  cohorts$data <- cohorts$data %>% 
    map(safely_mutate, 
        death = ifelse(death == 1, 'yes', 'no'),
        death = factor(death)) %>% 
    map(safely_mutate, 
        relapse = ifelse(relapse == 1, 'yes', 'no'), 
        relapse = factor(relapse))
  
  cohorts$var_list <- cohorts$data %>% 
    map(names) %>% 
    map(~cohorts$lexicon$variable[cohorts$lexicon$variable %in% .x])

# Descriptive statistic ------
  
  insert_msg('Descriptive stats')
  
  cohorts$stats <- 
    map2(cohorts$data, 
         cohorts$var_list, 
         ~explore(data = .x, 
                  variables = .y, 
                  what = 'table', 
                  pub_styled = TRUE)) %>% 
    format_desc() %>% 
    set_names(c('Variable', 
                globals$study_labels[names(cohorts$data)])) %>% 
    format_summ_tbl(rm_n = FALSE) %>% 
    mutate(Variable = factor(Variable, cohorts$lexicon$variable)) %>% 
    arrange(Variable) %>% 
    mutate(Variable = exchange(as.character(Variable), 
                               cohorts$lexicon))
  
# appending the stat table with the total N numbers ogf cancer samples ------
  
  insert_msg('N numbers')
  
  cohorts$n_nunbers <- cohorts$data %>% 
    map_dbl(nrow) %>% 
    compress(names_to = 'cohort', 
             values_to = 'n') %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    column_to_rownames('cohort') %>% 
    t %>% 
    as_tibble %>% 
    mutate(Variable = 'Cancer samples, N') %>% 
    relocate(Variable)
  
  cohorts$stats <- 
    full_rbind(cohorts$n_nunbers, 
               cohorts$stats)
  
# END -----
  
  cohorts <- cohorts[c("lexicon", "stats")]
  
  insert_tail()