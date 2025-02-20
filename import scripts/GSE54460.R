# Import of the study data for GSE54460

  insert_head()
  
# container ------
  
  gse54460 <- list()
  
# reading from GEO ------
  
  insert_msg('Reading from GEO')
  
  gse54460$raw <- getGEO(GEO = 'GSE54460', destdir = './data/GSE54460')
  
# Clinical information --------
  
  insert_msg('Clinical information')
  
  gse54460$clinic <- pData(gse54460$raw[[1]]) %>% 
    as_tibble
  
  gse54460$clinic <- gse54460$clinic %>% 
    transmute(sample_id = geo_accession,
              tissue_type = factor('tumor', c('normal', 'tumor')), 
              gleason_major = as.numeric(stri_extract(`gleason score:ch1`, regex = '^\\d{1}')), 
              gleason_sum = as.numeric(stri_extract(`gleason score:ch1`, regex = '\\d{1}$')), 
              gleason_minor = gleason_sum - gleason_major, 
              gleason_simple = cut(gleason_sum, 
                                   c(-Inf, 6, 7, Inf), 
                                   c('ISUP1', 'ISUP2', 'ISUP3+')), 
              gleason_sum = factor(gleason_sum), 
              gleason_minor = factor(gleason_minor), 
              gleason_major = factor(gleason_major), 
              relapse = as.numeric(`bcr:ch1`), 
              rfs_months = as.numeric(`months to bcr:ch1`), 
              rfs_months = ifelse(is.na(rfs_months), 
                                  as.numeric(`months total f/u:ch1`), 
                                  rfs_months), 
              race = ifelse(`race:ch1` == 'NA', NA, `race:ch1`), 
              race = factor(race), 
              surgical_margins = factor(tolower(`sms:ch1`), 
                                        c('negative', 'positive')), 
              psa_diagnosis = as.numeric(`prepsa:ch1`), 
              pt_stage = stri_extract(`pstage:ch1`, regex = 'T\\d{1}'), 
              pt_stage = factor(pt_stage, c('T1', 'T2', 'T3', 'T4')))
  
  ## appending with a sample mapping scheme corrected by hand: 
  ## there was a quite mess up...
  
  gse54460$sample_map <- read_xlsx('./data/GSE54460/sample_map.xlsx')
  
  gse54460$clinic <- 
    left_join(gse54460$clinic, gse54460$sample_map, by = 'sample_id')
  
# Expression and annotation --------
  
  insert_msg('Expression and annotation')
  
  gse54460$expression <-
    read_tsv('./data/GSE54460/GSE54460_FPKM-genes-TopHat2-106samples-12-4-13.txt')

  gse54460$expression <- gse54460$expression[-14:-1, ] %>% 
    map_dfc(parse_guess)
  
  ## annotation of the genes
  
  gse54460$annotation <- 
    tibble(probe_id = as.character(1:nrow(gse54460$expression)), 
           entrez_id = as.character(gse54460$expression[[1]])) %>% 
    annotate_raw_symbol
  
  ## expression: they seem to be already log2 with a single outlier
  
  gse54460$expression <- gse54460$expression[, -2:-1] %>% 
    map_dfc(as.numeric) %>% 
    as.matrix
  
  rownames(gse54460$expression) <- 
    as.character(1:nrow(gse54460$expression))
  
  gse54460$expression <- gse54460$expression %>% 
    integrate_expression(gse54460$annotation) %>% 
    mutate(assay_id = sample_id) %>% 
    select(- sample_id) %>% 
    left_join(gse54460$clinic[c('sample_id', 'patient_id', 'assay_id', 'tissue_type')], ., 
              by = 'assay_id') %>% 
    select(-assay_id)
  
# Caching the results -------
  
  insert_msg('Caching the results')
  
  gse54460 <- gse54460[c("clinic", "annotation", "expression")]
  
  save(gse54460, file = './data/gse54460.RData')
  
# END -----
  
  insert_tail()
  
  