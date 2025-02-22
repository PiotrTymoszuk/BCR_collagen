# This script provides project specific tools ----

# tools ------

  library(tidyverse)
  library(cowplot)
  library(furrr)
  library(caret)
  library(doParallel)
  library(survminer)
  library(exda)
  library(furrr)
  library(stringi)
  library(rlang)
  library(ggrepel)
  library(caret)
  library(figur)

# data import ------

  annotate_raw_symbol <- function(data) {
    
    data %>%
      filter(!is.na(entrez_id),
             entrez_id != '') %>%
      mutate(gene_symbol = mapIds(org.Hs.eg.db,
                                  column = 'SYMBOL',
                                  keys = entrez_id,
                                  keytype = 'ENTREZID')) %>%
      mutate(gene_symbol = unlist(gene_symbol)) %>%
      filter(gene_symbol != '') %>%
      filter(complete.cases(.))
    
  }

# Result formatting --------
  
  format_desc <- function(desc_data, merge_fun = full_join) {
    
    levs <- names(desc_data)
    
    desc_data %>% 
      reduce(merge_fun, by = 'variable') %>% 
      set_names(c('variable', levs))
    
  }
  
  format_summ_tbl <- function(data, 
                              rm_n = TRUE, 
                              rm_mean = TRUE, 
                              rm_complete = TRUE) {
    
    ## formats a summary table with descriptive stats
    
    data <- data %>% 
      map_dfc(stri_replace, regex = 'no:.*\\nyes:\\s{1}', replacement = '') %>% 
      map_dfc(stri_replace, regex = '\\nno:.*$', replacement = '\n') %>% 
      map_dfc(stri_replace_all, fixed = '% (', replacement = '% (n = ') %>% 
      map_dfc(stri_replace, fixed = 'Median =', replacement = 'median:') %>% 
      map_dfc(stri_replace, fixed = 'Mean =', replacement = 'mean:') %>% 
      map_dfc(stri_replace, fixed = 'Range', replacement = 'range') %>% 
      map_dfc(stri_replace, fixed = 'Complete', replacement = 'complete') %>% 
      map_dfc(stri_replace, fixed = ' [', replacement = '\n[')
    
    if(rm_n) {
      
      data <- data %>% 
        map_dfc(stri_replace, regex = '\\ncomplete.*$', replacement = '')
      
    }
    
    if(rm_mean) {
      
      data <- data %>% 
        map_dfc(stri_replace, 
                regex = 'mean.*\\nmedian:\\s{1}', 
                replacement = '') %>% 
        map_dfc(stri_replace, fixed = '\n[', replacement = ' [')
      
    }
    
    if(rm_complete) {
      
      data <- data %>% 
        map_dfc(stri_replace, fixed = 'complete: ', replacement = '')
      
    }
    
    data
    
  }
  
# ML modeling of survival --------
  
  cut_tertiles <- function(x) {
    
    terts <- quantile(x, c(1/3, 2/3), na.rm = TRUE)
    
    cuts <- c(-Inf, terts, Inf)
    
    cut(x, cuts, c('low', 'int', 'high'))
    
  }

  plot_surv_stats <- function(stats, 
                              plot_title = 'Modeling of BCR-free survival', 
                              plot_subtitle = 'Ridge Cox algorithm', 
                              x_lab = 'C-index', 
                              y_lab = '1 - IBS', 
                              palette = globals$model_colors, 
                              labels = globals$study_labels, 
                              label_variable = 'cohort', 
                              txt_size = 2.75, 
                              color_variable = 'dataset', ...) {
    
    ## a scatter plot of 1 - IBS and C-index
    ## cohort type is color-coded

    stats[[label_variable]] <- 
      labels[stats[[label_variable]]]

    stats %>% 
      ggplot(aes(x = c_index, 
                 y = 1 - ibs_model, 
                 fill = .data[[color_variable]], 
                 color = .data[[color_variable]])) + 
      geom_vline(xintercept = 0.5, 
                 linetype = 'dashed') +
      geom_hline(yintercept = 0.75, 
                 linetype = 'dashed') + 
      geom_point(size = 2, 
                 shape = 21, 
                 color = 'black') + 
      geom_text_repel(aes(label = .data[[label_variable]]), 
                      size = txt_size, 
                      show.legend = FALSE, ...) + 
      scale_fill_manual(values = palette, 
                        name = '') + 
      scale_color_manual(values = palette, 
                         name = '') + 
      globals$common_theme + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           x = x_lab, 
           y = y_lab)
    
  }
  
  plot_tertile_km <- function(fit, 
                              n_numbers, 
                              p_value, 
                              palette = globals$tertile_colors, 
                              plot_title = NULL, 
                              x_lab = 'min/max scaled BCR-free survival',
                              legend_title = 'Score tertile', 
                              txt_size = 2.75) {
    
    ## Kaplan-Meier plots for tertiles of a transcriptional score
    
    plot_subtitle <- paste0('total: n = ', sum(n_numbers$n_total), 
                            ', events: n = ', sum(n_numbers$n_events))
    
    scale_labs <- n_numbers %>% 
      mutate(scale_lab = paste0(score_cuts, '\ntotal: n = ', 
                                n_total, '\nevents: n = ', 
                                n_events, '\n'))
    
    scale_labs <- set_names(scale_labs$scale_lab, 
                            as.character(scale_labs$score_cuts))
    
    
    km_plot <- ggsurvplot(fit = fit, 
                          pval = p_value, 
                          pval.size = txt_size, 
                          pval.coord = c(0.01, 0.05))
    
    km_plot$plot + 
      scale_color_manual(values = unname(palette),
                         labels = unname(scale_labs), 
                         name = legend_title) +
      globals$common_theme + 
      theme(legend.position = 'bottom') + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           x = x_lab)
    
  }
  
  plot_surv_importance <- function(data, 
                                   imp_stat = 'coef', 
                                   n_top = NULL, 
                                   form = c('bar', 'bubble'), 
                                   labeller = c(positive = 'unfavorable', 
                                                negative = 'favorable', 
                                                ns = 'ns'), 
                                   plot_title = 'Transcriptional Collagen Score, variable importance', 
                                   plot_subtitle = 'Ridge Cox regression', 
                                   size_title = 'abs(log HR)', 
                                   max_size = 5,
                                   x_lab = expression('log HR'[Ridge]), 
                                   line_color = 'black', 
                                   palette = c(positive = 'firebrick', 
                                               negative = 'steelblue', 
                                               ns = 'gray70'), 
                                   gene_variables = globals$genes, 
                                   flip = FALSE, 
                                   split_regulation = TRUE) {
    
    ## a bar or bubble plot of variable importance stats
    
    form <- match.arg(form, c('bar', 'bubble'))
    
    gene_regex <- gene_variables %>% 
      sort(decreasing = TRUE) %>% 
      paste(collapse = '|')
    
    data <- data %>% 
      mutate(regulation = ifelse(.data[[imp_stat]] > 0, 'positive', 
                                 ifelse(.data[[imp_stat]]< 0, 'negative', NA)),
             regulation = factor(regulation, c('positive', 'negative')), 
             variable = ifelse(stri_detect(variable, regex = gene_regex), 
                               html_italic(variable), 
                               ifelse(variable == 'KLK3', 
                                      html_italic(variable), 
                                      exchange(variable, 
                                               globals$clinical_lexicon))), 
             variable = stri_replace(variable, 
                                     fixed = '_sq', 
                                     replacement = '<sup>2</sup>'))
    
    ## handling the special case: clinical and expression values
    
    if(any(stri_detect(data$variable, fixed = 'stage'))) {
      
      data <- data %>% 
        mutate(regulation = ifelse(stri_detect(variable, regex = gene_regex), 
                                   'collagen-related\ntranscript', 'clinical'), 
               regulation = factor(regulation, 
                                   c('collagen-related\ntranscript', 'clinical'))) 
      
    }
    
    if(!is.null(n_top)) {
      
      data <- data %>% 
        blast(regulation) %>% 
        map_dfr(top_n, n = n_top, abs(.data[[imp_stat]]))
      
    }
    
    data <- data %>% 
      filter(!is.na(regulation), 
             .data[[imp_stat]] != 0)
    
    if(form == 'bar') {
      
      if(!flip) {
        
        imp_plot <- data %>% 
          ggplot(aes(x = .data[[imp_stat]], 
                     y = reorder(variable, .data[[imp_stat]]), 
                     fill = regulation))
        
      } else {
        
        imp_plot <- data %>% 
          ggplot(aes(y = .data[[imp_stat]], 
                     x = reorder(variable, -.data[[imp_stat]]), 
                     fill = regulation))
        
      }
      
      imp_plot <- imp_plot + 
        geom_bar(stat = 'identity', 
                 color = line_color) 
      
    } else {
      
      if(!flip) {
        
        imp_plot <- data %>% 
          ggplot(aes(x = .data[[imp_stat]], 
                     y = reorder(variable, .data[[imp_stat]]), 
                     fill = regulation, 
                     size = abs(.data[[imp_stat]])))
        
      } else {
        
        imp_plot <- data %>% 
          ggplot(aes(y = .data[[imp_stat]], 
                     x = reorder(variable, -.data[[imp_stat]]), 
                     fill = regulation, 
                     size = abs(.data[[imp_stat]])))
        
      }
      
      imp_plot <- imp_plot + 
        geom_point(shape = 21) + 
        scale_size_area(max_size = max_size, 
                        name = size_title)
      
    }
    
    if(split_regulation) {
      
      imp_plot <- imp_plot + 
        facet_grid(regulation ~ ., 
                   scales = 'free', 
                   space = 'free', 
                   labeller = as_labeller(labeller))
      
    }
    
    imp_plot <- imp_plot + 
      scale_fill_manual(values = palette) + 
      globals$common_theme + 
      labs(title = plot_title, 
           subtitle = plot_subtitle)
    
    if(!flip) {
      
      imp_plot <- imp_plot + 
        theme(axis.title.y = element_blank(), 
              axis.text.y = element_markdown()) + 
        labs(x = x_lab)
      
    } else {
      
      imp_plot <- imp_plot + 
        theme(axis.title.x = element_blank(), 
              axis.text.x = element_markdown()) + 
        labs(y = x_lab)
      
    }
    
  }
  
# Univariable survival modeling -------
  
  hr_from_cut <- function(survcut_obj) {
    
    ## fits a Cox PH model (high versus low expression strata) to expression 
    # data stratified by the optimal cutoff
    ## used to identify favorable and unfavorable markers.

    ## construction of the Cox PH model ------
    
    cox_data <- survcut_obj %>% 
      model.frame %>% 
      select(rfs_months, relapse, ends_with('_strata'))
    
    cox_model <- call2('coxph', 
                       Surv(rfs_months, relapse) ~ ., 
                       data = cox_data, 
                       x = TRUE, 
                       y = TRUE) %>% 
      eval %>% 
      as_coxex(cox_data)
    
    ## fit stats --------
    
    stats <- cox_model %>% 
      summary('fit') %>%
      transmute(c_index = c_index, 
                c_lower = lower_ci, 
                c_upper = upper_ci, 
                ibs_refence = ibs_reference, 
                ibs_model = ibs_model)
    
    ## inference ------
    
    inference <- cox_model %>% 
      summary('inference')
    
    variable <- inference$variable
    
    inference <- inference %>% 
      transmute(gene_symbol = stri_replace(variable, 
                                           regex = '_.*', 
                                           replacement = ''), 
                z_stat = stat, 
                hr = exp(estimate), 
                se = exp(se), 
                hr_lower = exp(lower_ci), 
                hr_upper = exp(upper_ci), 
                p_value = p_value)
    
    ## cutoff value, numbers of observations and relapses -------
    
    best_cutoff <- cutoff(survcut_obj)[[1]]
 
    counts <- list()
    
    counts$total <- cox_data %>% 
      count(.data[[variable]]) %>% 
      mutate(strata = c('n_total_low', 'n_total_high'))
    
    counts$events <- cox_data %>% 
      filter(relapse == 1) %>% 
      count(.data[[variable]]) %>% 
      mutate(strata = c('n_events_low', 'n_events_high'))
  
    counts <- counts %>% 
      map(select, strata, n) %>% 
      map_dfr(column_to_rownames, 'strata') %>% 
      t %>% 
      as.data.frame

    ## the output -------
    
    cbind(counts, inference, stats) %>% 
      mutate(cutoff = best_cutoff) %>% 
      as_tibble %>% 
      relocate(gene_symbol, cutoff)

  } 
  
# Markdown functions ------

  my_word <- function(...) {
    
    form <- word_document2(number_sections = FALSE, 
                           reference_docx = 'ms_template.docx')
    
    form$pandoc$lua_filters <- c(form$pandoc$lua_filters, 
                                 'scholarly-metadata.lua', 
                                 'author-info-blocks.lua')
    
    form
    
  }
  
  format_strata_n <- function(count_data) {
    
    map2_chr(count_data[[1]], count_data[[2]], 
             paste, sep = ': n = ') %>% 
      paste(collapse = ', ')
    
  }
  
# END -----