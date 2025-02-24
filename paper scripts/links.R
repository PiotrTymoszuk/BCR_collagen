# Handling of common links and biobiography in the Rmarkdown project

  insert_head()

# container ------

  proj_links <- list()

# Links to development packages --------

  insert_msg('Links to development packages')

  proj_links <-
    list('ExDA' = 'https://github.com/PiotrTymoszuk/ExDA',
         'trafo' = 'https://github.com/PiotrTymoszuk/trafo',
         'figur' = 'https://github.com/PiotrTymoszuk/figur',
         'clustTools' = 'https://github.com/PiotrTymoszuk/clustTools',
         'microViz' = 'https://github.com/PiotrTymoszuk/microViz', 
         'coxExtensions' = 'https://github.com/PiotrTymoszuk/coxExtensions', 
         'kmOptimizer' = 'https://github.com/PiotrTymoszuk/kmOptimizer', 
         'fastTest' = 'https://github.com/PiotrTymoszuk/fastTest', 
         'graphExtra' = 'https://github.com/PiotrTymoszuk/graphExtra', 
         'htGLMNET' = 'https://github.com/PiotrTymoszuk/htGLMNET') %>%
    compress(names_to = 'obj_name',
             values_to = 'x') %>%
    mutate(ref_name = paste0('_', obj_name, '_'))

  proj_links <- proj_links[c('ref_name', 'x')] %>%
    pmap(mdlink) %>%
    set_names(proj_links$obj_name)
  
# reading the bibliography -------
  
  insert_msg('Rending the bibliography')
  
  coll_bib <- read_bib('./paper/markdown/coll_biblio.bib')

# END ------

  insert_tail()
