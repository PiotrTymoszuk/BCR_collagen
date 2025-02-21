# Renders the manuscript and the supplementary material

  insert_head()
  
# reading the bibliography -------
  
  insert_msg('Raeding the bibliography')
  
  coll_bib <- read_bib('./paper/markdown/coll_biblio.bib')

# analysis report -----
  
  insert_msg('Rendering manuscript parts')
  
  render('./paper/markdown/manuscript_figures_tables.Rmd', 
         output_format = word_document2(number_sections = FALSE, 
                                        reference_docx = 'ms_template.docx'), 
         output_dir = './paper')
  
  render('./paper/markdown/manuscript_supplement.Rmd', 
         output_format = word_document2(number_sections = FALSE, 
                                        reference_docx = 'ms_template.docx'), 
         output_dir = './paper')
  
# END -----
  
  insert_tail()