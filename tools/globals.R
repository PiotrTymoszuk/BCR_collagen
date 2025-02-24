# Project globals 

  insert_head()

# container --------

  globals <- list()
  
# genes of interest --------

  ## classification of the genes by their function
  
  globals$genes_lexicon <-
    c('ADAMTS2' = 'ECM processing', 
      'ALDH18A1' = 'proline metabolism', 
      'BMP1' = 'ECM processing', 
      'CD151' = 'adhesion', 
      'COL11A1' = 'ECM component', 
      'COL11A2' = 'ECM component', 
      'COL14A1' = 'ECM component', 
      'COL15A1' = 'ECM component', 
      'COL16A1' = 'ECM component', 
      'COL17A1' = 'ECM component', 
      'COL18A1' = 'ECM component', 
      'COL19A1' = 'ECM component', 
      'COL1A1' = 'ECM component', 
      'COL1A2' = 'ECM component', 
      'COL2A1' = 'ECM component', 
      'COL3A1' = 'ECM component', 
      'COL4A1' = 'ECM component', 
      'COL4A2' = 'ECM component', 
      'COL4A3' = 'ECM component', 
      'COL4A5' = 'ECM component', 
      'COL4A6' = 'ECM component', 
      'COL5A1' = 'ECM component', 
      'COL5A2' = 'ECM component', 
      'COL6A1' = 'ECM component', 
      'COL6A2' = 'ECM component', 
      'COL6A3' = 'ECM component', 
      'COL7A1' = 'ECM component', 
      'COL9A1' = 'ECM component', 
      'COL9A2' = 'ECM component', 
      'COL9A3' = 'ECM component', 
      'CTSS' = 'ECM processing', 
      'DST' = 'adhesion', 
      'ITGA6' = 'adhesion', 
      'ITGB4' = 'adhesion', 
      'LAMA3' = 'ECM component', 
      'LAMB3' = 'ECM component', 
      'LAMC2' = 'ECM component', 
      'LOX' = 'collagen modification', 
      'LOXL1' = 'collagen modification',
      'LOXL2' = 'collagen modification', 
      'MMP13' = 'ECM processing', 
      'MMP7' = 'ECM processing', 
      'MMP9' = 'ECM processing', 
      'P4HA1' = 'collagen modification', 
      'P4HA2' = 'collagen modification', 
      'P4HB' = 'collagen modification', 
      'PCOLCE' = 'ECM processing', 
      'PCOLCE2' = 'ECM processing', 
      'PEPD' = 'proline metabolism', 
      'PLOD1' = 'collagen modification', 
      'PLOD2' = 'collagen modification', 
      'PLOD3' = 'collagen modification', 
      'PPIB' = 'collagen modification', 
      'PYCR1' = 'proline metabolism', 
      'SERPINH1' = 'ECM processing') %>% 
    compress(names_to = 'gene_symbol', 
             values_to = 'gene_group') %>% 
    mutate(gene_group = factor(gene_group, 
                               c('proline metabolism', 
                                 'collagen modification', 
                                 'ECM component', 
                                 'ECM processing', 
                                 'adhesion')))
  
  globals$genes <- globals$genes_lexicon$gene_symbol
  
# clinical variables -------
  
  globals$clinical_lexicon <- 
    c('age' = 'Age at diagnosis, years', 
      'psa_diagnosis' = 'PSA at diagnosis', 
      'ct_stage' = 'Clinical tumor stage', 
      'pt_stage' = 'pT stage', 
      'pn_stage' = 'pN stage', 
      'pm_stage' = 'pM stage', 
      'gleason_sum' = 'Gleason score', 
      'gleason_simple' = 'ISUP grade', 
      'surgical_margins' = 'Surgical margins', 
      'ece' = 'Extracapsular extension', 
      'death' = 'Death', 
      'os_months' = 'Overall survival, months', 
      'relapse' = 'Biochemical relapse', 
      'rfs_months' = 'Biochemical relapse-free survival, months') %>% 
    compress(names_to = 'variable', 
             values_to = 'label') %>% 
    mutate(format = ifelse(variable %in% c('age', 
                                           'psa_diagnosis', 
                                           'gleason', 
                                           'vitality_fup', 
                                           'relapse_fup'), 
                           'numeric', 'factor'))
  
# graphics --------
  
  ## theme for most plots
  
  globals$common_text <- element_text(size = 8, 
                                      face = 'plain', 
                                      color = 'black')
  
  globals$common_margin <- ggplot2::margin(t = 3, 
                                           l = 3, 
                                           r = 3, 
                                           unit = 'mm')
  
  globals$common_theme <- theme_classic() + 
    theme(axis.text = globals$common_text, 
          axis.title = globals$common_text, 
          plot.title = element_text(size = 8, 
                                    face = 'bold'), 
          plot.subtitle = globals$common_text, 
          plot.tag = element_text(size = 8, 
                                  face = 'plain', 
                                  color = 'black', 
                                  hjust = 0, 
                                  vjust = 1), 
          plot.tag.position = 'bottom', 
          legend.text = globals$common_text, 
          legend.title = globals$common_text, 
          strip.text = globals$common_text,
          strip.background = element_rect(fill = 'white'), 
          plot.margin = globals$common_margin, 
          panel.grid.major = element_line(color = 'gray90'))
  
  ## theme for network plots
  
  globals$net_theme <- theme_void() + 
    theme(plot.title = element_text(size = 8, 
                                    face = 'bold'), 
          plot.subtitle = globals$common_text, 
          legend.text = globals$common_text, 
          legend.title = globals$common_text,
          plot.margin = globals$common_margin)
  
# Study labels and colors --------
  
  ## GSE116918 is excluded from the analysis

  globals$study_labels <- c('gse16560' = 'GSE16560', 
                            'gse54460' = 'GSE54460', 
                            'gse70768' = 'GSE70768', 
                            'gse70769' = 'GSE70769', 
                             'gse220095' = 'GSE220095', 
                            'tcga' = 'TCGA', 
                            'dkfz' = 'DKFZ', 
                            'geo_pool' = 'pooled GEO')

  globals$study_arrays <- c('gse16560' = TRUE, 
                            'gse54460' = FALSE, 
                            'gse70768' = TRUE, 
                            'gse70769' = TRUE, 
                            'gse220095' = FALSE, 
                            'tcga' = FALSE, 
                            'dkfz' = FALSE)
  
  ## a ready to use expression for calling a cohort list
  
  globals$analysis_studies <- names(globals$study_labels)
  
  globals$study_exprs <- globals$analysis_studies %>% 
    map_chr(~paste(.x, .x, sep = ' = ')) %>% 
    paste(collapse = ', ') %>% 
    paste0('list(', ., ')') %>% 
    parse_expr
  
  ## colors for the collectives used in BCR and OS modeling
  
  globals$bcr_study_colors <- 
    c('geo_pool' = 'indianred3', 
      'tcga' = 'steelblue', 
      'dkfz' = 'aquamarine3')
  
# labels and colors for algorithms -------
  
  insert_msg('Colors and labels')
  
  globals$tertile_colors <- c(low = 'darkolivegreen', 
                                   int = 'steelblue', 
                                   high = 'firebrick')
  
  globals$model_colors <-
    c(test = 'steelblue', 
      training = 'indianred3')
  
  globals$algo_labels <- 
    c(ridge = 'RIDGE Cox', 
      elnet = 'Elastic Net Cox', 
      lasso = 'LASSO Cox', 
      svm = 'SVM', 
      rf = 'Random Forest', 
      gbm = 'GBM', 
      gbm_clinic = 'GBM, clinic', 
      gbm_combi = 'GBM, clinic + expression')
  
  globals$algo_colors <- 
    c(ridge = 'plum4', 
      elnet = 'steelblue2', 
      lasso = 'indianred3', 
      svm = 'gray60', 
      rf = 'darkolivegreen', 
      gbm = 'orangered2', 
      gbm_clinic = 'aquamarine3', 
      gbm_combi = 'firebrick4')
  
  ## comparison of the collagen models with survival models
  ## with clinical risk factors
  
  globals$type_colors <- 
    c(clinic = 'gray60', 
      collagen = 'orangered2', 
      full = 'firebrick4')
  
  globals$type_labels <- 
    c(clinic = 'clinical factors', 
      collagen = 'collagen genes', 
      full = 'collagen genes + clinical factors')
  
  ## testing for the confounding effect of the cohort
  ## in the TCGA BLCA cohort
  
  globals$cohort_model_colors <- 
    c(cohort = 'darkseagreen4', 
      gbm_score = 'orangered2', 
      full = 'firebrick4')
  
  globals$cohort_model_labels <- 
    c(cohort = 'study', 
      gbm_score = 'GBM predictor', 
      full = 'GBM predictor + study')

# END -----
  
  insert_tail()