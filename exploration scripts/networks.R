# co-expression networks of the collagen-related genes. 
# The networks are built for correlations of ComBat-adjusted log2 expression 
# levels with spearman's rho > 0.5. 
#
# Communities are identified by Leiden's algorithm.

  insert_head()
  
# container -------
  
  net <- list()
  
# analysis globals: variables, their attributes and ComBat expression -------
  
  insert_msg('Analysis globals')
  
  ## variables and their attribute (functional group)
  
  net$variables <- globals$genes
  
  net$attr_tbl <- globals$genes_lexicon %>% 
    set_names(c('variable', 'gene_group'))
  
  ## ComBat-adjusted log2 expression levels
  
  net$data <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$combat) %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, all_of(net$variables))
  
  ## attribute lexicon as a list of tuples used later for plotting
  
  net$attr_lexicon <- 
    list(community_id = c('community_id', 'gene communities'), 
         gene_group = c('gene_group', 'gene function'), 
         degree = c('degree', 'degree'), 
         hub_score = c('hub_score', 'hub score'), 
         betweennes = c('betweenness', 'betweenness')) 
  
# Graph objects -------
  
  insert_msg('Graph objects')
  
  ## isolated nodes are removed
  ## setting the attributes
  
  net$graph_obj <- net$data %>% 
    map(as_iGraph, 
        feature_type = 'columns', 
        fun = cor, 
        cutoff = 0.5, 
        weighted = TRUE, 
        diag = FALSE) %>% 
    map(prune_degree, cutoff = 0) %>% 
    map(set_vertex_attributes, 
        net$attr_tbl)
  
# Communities ---------
  
  insert_msg('Communities')
  
  set.seed(12345)
  
  ## communities with single genes (if present are lumped)
  
  net$communities <- net$graph_obj %>% 
    map(cluster_leiden, 
        objective_function = 'modularity', 
        resolution_parameter = 1,  
        n_iterations = 100) %>% 
    map(comm_lump, 
        cutoff = 1, 
        other_name = 'other')
  
  ## appending the graph nodes with the community 
  ## assignment information 
  
  net$graph_obj <- 
    map2(net$graph_obj, 
         net$communities, 
         add_communities)
  
# Vertex importance stats ------
  
  insert_msg('Vertex importance stats')
  
  ## and adding them to the node attributes
  
  net$stats <- net$graph_obj %>% 
    map(summary) %>% 
    map(select, -index)
  
  net$graph_obj <- 
    map2(net$graph_obj, 
         net$stats, 
         set_vertex_attributes)
  
  net$stats <- net$graph_obj %>% 
    map(get_vertex_attributes)
  
  ## calculating of the min/max-normalized betweenness
  
  net$stats <- net$stats %>% 
    map(mutate, 
        norm_betweenness = minMax(betweenness))
  
# Plots of the the vertex importance stats --------
  
  insert_msg('Plots of the globals vertex importance stats')
  
  net$stat_plots <- 
    list(x = net$stats, 
         y = globals$study_labels[names(net$stats)]) %>% 
    pmap(function(x, y) x %>% 
           ggplot(aes(x = degree, 
                      y = hub_score, 
                      size = norm_betweenness, 
                      fill = gene_group)) + 
           geom_point(shape = 21) + 
           geom_text_repel(aes(label = name),
                           size = 2.3, 
                           fontface = 'italic') + 
           scale_size_area(max_size = 4.5, 
                           name = 'min/max scaled\nbetweenness') + 
           scale_fill_brewer(palette = 'Set1', 
                             drop = FALSE, 
                             name = 'gene function') + 
           globals$common_theme + 
           labs(title = y, 
                 x = 'degree', 
                 y = 'hub score'))
  
# Network plots ---------
  
  insert_msg('Network plots')
  
  ## color codes for community_id, functional gene group, 
  ## degree, hub score, and betweenness
  
  for(i in net$attr_lexicon) {
    
    net$graph_plots[[i[[1]]]] <- 
      list(x = net$graph_obj, 
           plot_title = paste(globals$study_labels[names(net$graph_obj)], 
                              i[[2]], sep = ', ')) %>% 
      pmap(plot, 
           layout = layout.fruchterman.reingold, 
           vertex_fill_variable = i[[1]], 
           vertex_color_variable = NULL, 
           vertex_color = 'black', 
           edge_alpha = 0.4, 
           label_vertices = TRUE, 
           vertex_txt_color_variable = i[[1]], 
           vertex_txt_size = 2.5, 
           vertex_txt_face = 'italic', 
           cust_theme = globals$net_theme, 
           seed = 1234)
    
  }
  
  ## fill scales
  
  net$graph_plots$community_id <- net$graph_plots$community_id %>% 
    map(~.x + 
          scale_fill_brewer(palette = 'Set2', 
                            labels = function(x) paste0('#', x), 
                            name = 'gene community') + 
          scale_color_brewer(palette = 'Set2', 
                             labels = function(x) paste0('#', x), 
                             name = 'gene community'))
  
  net$graph_plots$gene_group <- net$graph_plots$gene_group %>% 
    map(~.x + 
          scale_fill_brewer(palette = 'Set1', 
                            drop = FALSE, 
                            name = 'gene function') + 
          scale_color_brewer(palette = 'Set1', 
                             drop = FALSE, 
                             name = 'gene function'))
  
  net$graph_plots[c("degree", "hub_score", "betweenness")] <- 
    net$graph_plots[c("degree", "hub_score", "betweenness")] %>% 
    map(map, 
        ~.x + 
          scale_fill_gradient(low = 'steelblue', 
                              high = 'firebrick') + 
          scale_color_gradient(low = 'steelblue', 
                               high = 'firebrick'))

# END ------
  
  insert_tail()