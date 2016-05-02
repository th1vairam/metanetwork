modPreservStats <- function(ref.mod.name, ref.g, test.g, ref.mod, test.mod){
  
  # Filter genes that intersect reference and test modules
  ref.mod = filter(ref.mod, EnsembleID %in% test.mod$EnsembleID)
  test.mod = filter(test.mod, EnsembleID %in% ref.mod$EnsembleID)

  # Collect members of reference modules
  ref.mod.mem = ref.mod$EnsembleID[ref.mod$moduleLabel == ref.mod.name]

  # Get reference module size
  ref.mod.sz = length(ref.mod.mem)

  # Extract sub networks from reference and test
  ref.sg = igraph::induced_subgraph(ref.g, which(V(ref.g)$name %in% ref.mod.mem))
  test.sg = igraph::induced_subgraph(test.g, which(V(test.g)$name %in% ref.mod.mem))

  ## Calculate modularity of test network
  tmp = test.mod$moduleNumber+1
  names(tmp) = test.mod$EnsembleID
  mod.test.sg = igraph::modularity(test.sg, tmp[intersect(names(tmp), ref.mod.mem)])

  ### Density preservation stats
  # Correlation between adjacency matrices
  ref.adj = igraph::get.adjacency(ref.sg)
  test.adj = igraph::get.adjacency(test.sg) %>% as.matrix
  test.adj = test.adj[rownames(ref.adj), colnames(ref.adj)]
  cor.adj = WGCNA::bicor(as.vector(ref.adj), as.vector(test.adj),  use = 'pairwise.complete.obs') %>%
    as.numeric

  # Module density in test network
  mean.adj = igraph::graph.density(test.sg)

  # Fold change in adjacency
  change.adj = mean.adj/igraph::graph.density(ref.sg)

  ### Connectivity preservation stats
  # Degree of nodes (k.all)
  test.allDegree = igraph::degree(test.g)

  # Degree of all nodes in sub-netowrk (k.in)
  test.degree = igraph::degree(test.sg)
  ref.degree = igraph::degree(ref.sg)
  ref.degree = ref.degree[names(test.degree)]
  meankIM = mean(test.degree, na.rm=T)

  # Correlation between 
  cor.kIM = WGCNA::bicor(as.vector(ref.degree), as.vector(test.degree), use = 'pairwise.complete.obs') %>%
    as.numeric

  # Calculate summary stats of intra modular connectivity
  kIM2kALL = test.degree/test.allDegree[names(test.degree)]
  meankIM2kALL = mean(kIM2kALL, na.rm=T)
  mediankIM2kALL = median(kIM2kALL, na.rm=T)
  sdkIM2kALL = sd(kIM2kALL, na.rm = T)

  return(data.frame(
    module.name = ref.mod.name,
    size = ref.mod.sz,
    modularity.test = mod.test.sg,
    cor.adj = cor.adj,
    mean.adj = mean.adj,
    change.adj = change.adj,
    meankIM = meankIM,
    cor.kIM = cor.kIM,
    meankIM2kALL = meankIM2kALL,
    mediankIM2kALL = mediankIM2kALL,
    sdkIM2kALL = sdkIM2kALL))
}