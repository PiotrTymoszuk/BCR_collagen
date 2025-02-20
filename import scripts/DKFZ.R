# Import of the DKFZ data set 

  insert_head()
  
# container ------
  
  dkfz <- list()
  
# clinical information --------
  
  insert_msg('Clinical information')
  
  dkfz$clinic[c('patient', 'sample')] <- 
    c('./data/DKFZ/data_clinical_patient.txt', 
      './data/DKFZ/data_clinical_sample.txt') %>% 
    map(read_tsv, skip = 4)
  
  ## patient's data
  
  dkfz$clinic$patient <- dkfz$clinic$patient %>% 
    transmute(patient_id = PATIENT_ID, 
              age = as.numeric(AGE), 
              pt_stage = stri_extract(STAGE, regex = 'T\\d{1}'), 
              pt_stage = factor(pt_stage, c('T1', 'T2', 'T3', 'T4')),
              psa_diagnosis = as.numeric(PREOPERATIVE_PSA), 
              treatment = factor(INITIAL_TREATMENT), 
              relapse = as.numeric(BCR_STATUS), 
              rfs_months = as.numeric(TIME_FROM_SURGERY_TO_BCR_LASTFU))
  
  ## sample
  
  dkfz$clinic$sample <- dkfz$clinic$sample %>% 
    transmute(patient_id = PATIENT_ID, 
              sample_id = SAMPLE_ID, 
              tissue_type = factor('tumor', c('normal', 'tumor')), 
              purity = as.numeric(MEDIAN_PURITY), 
              gleason_major = as.numeric(stri_extract(GLEASON_SCORE, regex = '^\\d{1}')), 
              gleason_minor = as.numeric(stri_extract(GLEASON_SCORE, regex = '\\d{1}$')), 
              gleason_sum = gleason_major + gleason_minor, 
              gleason_simple = cut(gleason_sum, 
                                   c(-Inf, 6, 7, Inf), 
                                   c('ISUP1', 'ISUP2', 'ISUP3+')),
              tmb = as.numeric(TMB_NONSYNONYMOUS))
  
  dkfz$clinic$sample[c('gleason_sum', 'gleason_minor', 'gleason_major')] <- 
    dkfz$clinic$sample[c('gleason_sum', 'gleason_minor', 'gleason_major')] %>% 
    map_dfc(factor)
  
  dkfz$clinic <- reduce(dkfz$clinic, inner_join, by = 'patient_id')
  
# Expression and annotation --------
  
  insert_msg('Expression and annotation')
  
  dkfz$expression <- read_tsv('./data/DKFZ/data_mrna_seq_rpkm.txt')
  
  ## annotation via symbol and alias
  
  dkfz$annotation <- dkfz$expression %>% 
    transmute(probe_id = as.character(1:nrow(.)), 
              gene_symbol = Hugo_Symbol) %>% 
    filter(!is.na(gene_symbol),
           gene_symbol != '') %>% 
    mutate(entrez_id = mapIds(org.Hs.eg.db, 
                              keys = gene_symbol, 
                              keytype = 'SYMBOL', 
                              column = 'ENTREZID')) %>% 
    mutate(entrez_id = ifelse(is.na(entrez_id), 
                              mapIds(org.Hs.eg.db, 
                                     keys = gene_symbol, 
                                     keytype = 'ALIAS', 
                                     column = 'ENTREZID'), 
                              entrez_id)) %>% 
    annotate_raw_symbol
  
  ## expression: log2 and agggregation of duplicated genes
  ## by arithmetic means
  
  dkfz$expression <- dkfz$expression %>% 
    select(-Hugo_Symbol, -Entrez_Gene_Id) %>% 
    as.matrix
  
  rownames(dkfz$expression) <- as.character(1:nrow(dkfz$expression))
  
  dkfz$expression <- log2(dkfz$expression + 1)

  dkfz$expression <- dkfz$expression %>% 
    integrate_expression(dkfz$annotation) %>% 
    inner_join(dkfz$clinic[c('sample_id', 'patient_id', 'tissue_type')], ., 
               by = 'sample_id')

# Mutations -------
  
  insert_msg('Mutations')
  
  plan('multisession')
  
  ## non-silent ones, assignment to a gene required
  
  dkfz$mutations <- read_tsv('./data/DKFZ/data_mutations.txt')
  
  dkfz$mutations <- dkfz$mutations %>% 
    filter(!Consequence %in% c('synonymous_variant', 
                               'intron_variant', 
                               'intergenic_variant', 
                               'intron_variant,non_coding_transcript_variant', 
                               'downstream_gene_variant', 
                               'upstream_gene_variant'), 
           !is.na(Entrez_Gene_Id), 
           !is.na(Hugo_Symbol)) %>%
    extract_mutations
  
  plan('sequential')
  
# Common samples with clinical and transcriptome -------
  
  insert_msg('Common samples')
  
  dkfz$common_samples <- dkfz[c("clinic", "expression")] %>% 
    map(~.x$sample_id) %>% 
    reduce(intersect)
  
  dkfz[c("clinic", "expression", "mutations")] <- 
    dkfz[c("clinic", "expression", "mutations")] %>% 
    map(filter, sample_id %in% dkfz$common_samples)
  
# caching the results -------
  
  insert_msg('Caching the results')
  
  dkfz <- dkfz[c("clinic", "annotation", "expression", "mutations")]
  
  save(dkfz, file = './data/dkfz.RData')
  
# END -------
  
  insert_tail()