source('nested.kmeans.R')
nested.kmeans.all <- function (g, single.size = 20, sig.method = "pval.perm", k.max = Inf,
                               n.skip = 20, n.perm = 100, module.pvalue = 0.05, 
                               d.func = function(x) {x}, alpha.range = seq(0.01, 10, 0.01)){
  # Remove singleton nodes
  el = igraph::as_data_frame(g, what = 'edges')
  g1 = graph.data.frame(el, directed = F)
  gc()
  
  singletons = setdiff(V(g)$name, V(g1)$name)
  
  # Get connected componenets
  con.comp = get.connectedModule(V(g1)$name, g1)
  
  # Remove singletongs and small motifs
  len.con.comp = sapply(con.comp, length)
  singletons = c(singletons, 
                 con.comp[len.con.comp <= single.size] %>%
                   unlist)
  con.comp = con.comp[len.con.comp > single.size]
  
  output.mod <- vector("list", length(con.comp))
  for (i in 1:length(con.comp)) {
    gg <- igraph::induced.subgraph(g1, con.comp[[i]])
    output.mod[[i]] <- nested.kmeans(sub.g = gg, sig.method = sig.method, 
                                     n.singleton = single.size, k.max = k.max, n.skip = n.skip, 
                                     module.pvalue = module.pvalue, d.func = d.func, alpha.range = alpha.range)
  }
  singletons <- c(singletons)
  modules <- list()
  module.relation <- matrix(0, nrow = 0, ncol = 2)
  module.alpha <- c()
  module.pvalue <- c()
  module.compRel <- c()
  for (i in 1:length(output.mod)) {
    new.mod <- output.mod[[i]]$modules
    names(new.mod) <- paste("comp", i, "_", 1:length(new.mod), 
                            sep = "")
    modules <- c(modules, new.mod)
    mod.rel <- output.mod[[i]]$module.relation
    mod.rel <- cbind(names(new.mod)[mod.rel[, 1]], names(new.mod)[mod.rel[, 2]])
    module.relation <- rbind(module.relation, mod.rel)
    singletons <- c(singletons, unlist(output.mod[[i]]$singletons))
    module.alpha <- c(module.alpha, output.mod[[i]]$module.alpha)
    module.pvalue <- c(module.pvalue, output.mod[[i]]$module.pvalue)
    module.compRel <- c(module.compRel, output.mod[[i]]$module.compRel)
  }
  j <- which(sapply(modules, length) > single.size)
  modules <- modules[j]
  module.alpha <- module.alpha[j]
  module.pvalue <- module.pvalue[j]
  module.compRel <- module.compRel[j]
  names(singletons) <- NULL
  integer.relation <- cbind(match(module.relation[, 1], names(modules)), 
                            match(module.relation[, 2], names(modules)))
  output <- list(modules = modules, singletons = singletons, 
                 module.relation = integer.relation, module.alpha = module.alpha, 
                 module.pvalue = module.pvalue, module.compRel = module.compRel)
  return(output)
}