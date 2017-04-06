score.nodes <- function(g, G, h=3, mode = 'all'){
  
  # Get background genes from graph
  background.genes = igraph::V(g)$name
  G = G[intersect(background.genes, names(G))]
  G.not = rep(0, length(setdiff(background.genes, names(G))))
  names(G.not) = setdiff(background.genes, names(G))
  G = c(G, G.not)
  
  # Calculate node degree for identifying global drivers
  node.degree = igraph::degree(g)
  mean.node.degree = mean(node.degree, na.rm = T)
  stddev.node.degree = sd(node.degree, na.rm = T)
  
  # Find H-Layer neighborhood
  neighbor.nodes = lapply(1:h, function(hi,sg){
    igraph::ego(sg, order = hi, nodes = igraph::V(sg), mode = mode, mindist = hi)
  }, g)
  
  # Find summed weight of current node based on h-layer neighbors
  node.scores = foreach::foreach (x = 1:length(neighbor.nodes), .combine = cbind) %dopar% {
    score = sapply(neighbor.nodes[[x]], function(y, G){
      sum(G[names(y)] * (1/node.degree[names(y)]))
    }, G)
  } 
  node.scores = rowSums(node.scores * t(matrix(rep(1/c(1:h), length(G)), h, length(G))), na.rm = T)
  node.scores = node.scores + G
  names(node.scores) = background.genes
  
  return(node.scores)
}