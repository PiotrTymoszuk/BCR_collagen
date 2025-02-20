# Import of the study data for GSE220095

  insert_head()

# container ------

  gse220095 <- list()

# reading from GEO ------

  insert_msg('Reading from GEO')

  gse220095$raw <- getGEO(GEO = 'GSE220095', destdir = './data/GSE220095')
  
# Clinical information --------
  
  insert_msg('Clinical information')

  gse220095$clinic <- pData(gse220095$raw[[1]]) %>% 
    as_tibble
  
  gse220095$clinic <- gse220095$clinic %>% 
    transmute(sample_id = geo_accession, 
              patient_id = stri_extract(title, regex = 'RIB.*$'), 
              tissue_type = factor('tumor', c('normal', 'tumor')), 
              relapse = car::recode(`biochemical relapse_(bcr):ch1`, 
                                    "'yes' = 1; 'no' = 0"), 
              relapse = as.numeric(relapse), 
              rfs_months = as.numeric(`time to_bcr_(month),_censored_at_last_follow-up,_if_no_bcr:ch1`), 
              gleason_sum = stri_extract_all(`gleason score:ch1`, 
                                             regex = '\\d{1}'), 
              gleason_major = map_chr(gleason_sum, function(x) if(all(is.na(x))) NA else x[[2]]), 
              gleason_minor = map_chr(gleason_sum, function(x) if(all(is.na(x))) NA else x[[3]]), 
              gleason_sum = map_chr(gleason_sum, function(x) if(all(is.na(x))) NA else x[[1]]), 
              gleason_sum = as.numeric(gleason_sum),
              gleason_simple = cut(gleason_sum, 
                                   c(-Inf, 6, 7, Inf), 
                                   c('ISUP1', 'ISUP2', 'ISUP3+')), 
              pn_stage = stri_extract(`pathological n_stage:ch1`, 
                                      regex = 'N\\d{1}'), 
              pn_stage = factor(pn_stage, paste0('N', 0:5)), 
              pt_stage = stri_extract(`pathological t_stage:ch1`, 
                                      regex = 'T\\d{1}'),
              pt_stage = factor(pt_stage, paste0('T', 1:4)), 
              psa_diagnosis = as.numeric(`total psa_(ng/ml,_pre-biopsy):ch1`), 
              purity = as.numeric(`tumor cell_content_(%):ch1`))
  
  gse220095$clinic[c('gleason_sum', 'gleason_major', 'gleason_minor')] <- 
    gse220095$clinic[c('gleason_sum', 'gleason_major', 'gleason_minor')] %>% 
    map_dfc(factor)
  
  gse220095$clinic <- gse220095$clinic %>% 
    map_dfc(function(x) if(is.factor(x)) droplevels(x) else x)
    
# Expression and annotation --------
  
  insert_msg('Annotation and expression')
  
  gse220095$expression <- 
    read_csv('./data/GSE220095/GSE220095_VST_counts.csv')
  
  ## annotation
  
  gse220095$annotation <- gse220095$expression %>% 
    transmute(probe_id = as.character(1:nrow(.)), 
              ensembl_id = Ensembl_ID) %>% 
    mutate(entrez_id = mapIds(org.Hs.eg.db, 
                              keys = ensembl_id, 
                              keytype = 'ENSEMBL', 
                              column = 'ENTREZID')) %>% 
    annotate_raw_symbol
  
  ## expression: provided already as log2 counts
  
  gse220095$expression <- gse220095$expression %>% 
    select(-Ensembl_ID) %>% 
    as.matrix
  
  rownames(gse220095$expression) <- 
    as.character(1:nrow(gse220095$expression))
  
  gse220095$expression <- gse220095$expression %>% 
    integrate_expression(gse220095$annotation) %>% 
    mutate(patient_id = sample_id) %>% 
    select(- sample_id) %>% 
    left_join(gse220095$clinic[c('sample_id', 'patient_id', 'tissue_type')], ., 
              by = 'patient_id')
  
# Caching the results -------
  
  insert_msg('Caching the results')
  
  gse220095 <- gse220095[c("clinic", "annotation", "expression")]
  
  save(gse220095, file = './data/gse220095.RData')
  
# END -----
  
  insert_tail()