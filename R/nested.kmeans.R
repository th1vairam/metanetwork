nested.kmeans <- function (sub.g, k.max = Inf, n.skip = 20, n.perm = 100, module.pvalue = 0.05, 
                           sig.method = c("pval.perm", "pval.norm"), d.func = function(x) {x}, 
                           n.singleton = 20, alpha.range = seq(0.01, 10, 0.01))
{
  if (!is.connected(sub.g))
    stop("input network is not connected.")
  cat("Calculating distance metric and similarity...\n")
  
  D <- get.network.metric(sub.g, d.func = d.func)
  S <- get.LPI(sub.g)
  
  modules.keep <- list(V(sub.g)$name)
  modules.singleton <- list()
  
  do.test <- 1
  do.update <- TRUE
  alpha.defined <- c(Inf)
  pvalue.calc <- c(NA)
  comp.score <- c(NA)
  subset.relation <- matrix(0, ncol = 2, nrow = 0)
  n.iter <- 1
  while (do.update) {
    test.i <- which(do.test == 1)
    cat(paste("iteration:", n.iter, "\n", sep = ""))
    cat(paste("- #. tested:", length(test.i), "\n", sep = ""))
    for (i in 1:length(test.i)) {
      do.test[test.i[i]] <- 0
      alph <- alpha.range[which(alpha.range <= alpha.defined[test.i[i]])]
      D.rand <- get.random.D(nv = length(modules.keep[[test.i[i]]]), 
                             w = E(sub.g)[modules.keep[[test.i[i]]] %--% modules.keep[[test.i[i]]]]$weight, 
                             d.func = d.func, name.vec = V(sub.g)$name)
      comp.p <- compact.indiv(modules.keep[[test.i[i]]], 
                              D, alpha = alpha.range)
      sub.net <- igraph::induced.subgraph(sub.g, modules.keep[[test.i[i]]])
      cat("- ")
      pam.out <- pam.cleanSplit(sub.net, D, k.max, n.skip)
      pam.out <- do.call(c, lapply(pam.out, get.connectedModule, 
                                   sub.net))
      i.single <- which(sapply(pam.out, length) <= n.singleton)
      modules.singleton <- c(modules.singleton, pam.out[i.single])
      pam.out <- pam.out[setdiff(1:length(pam.out), i.single)]
      cat(paste("- #. of split:", length(pam.out), "\n", 
                sep = ""))
      if (length(pam.out) > 1) {
        cat("- assess improvements over compactness\n")
        comp.c <- lapply(pam.out, function(v, d, alph) compact.indiv(v, 
                                                                     d, alph), d = D, alph = alpha.range)
        comp.rel <- lapply(comp.c, function(x, y) x/y, 
                           y = comp.p)
        def.alpha <- rep(NA, length(comp.rel))
        comp.sig <- rep(NA, length(comp.rel))
        comp.sc <- rep(NA, length(comp.rel))
        for (j in 1:length(comp.rel)) {
          if (any(comp.rel[[j]] < 1)) {
            alpha.i <- max(which(comp.rel[[j]] < 1))
            def.alpha[j] <- min(c(alpha.range[alpha.i], 
                                  alpha.defined[test.i[i]]))
            stat.out <- evaluate.significance(score.o = comp.rel[[j]][alpha.i], 
                                              nv = length(pam.out[[j]]), D.rand = D.rand, 
                                              alpha = alph[alpha.i], n.perm = n.perm)
            if (sig.method == "pval.perm") {
              pval = stat.out$pval
            }
            else {
              pval = stat.out$pval.norm
            }
            comp.sig[j] <- pval
            comp.sc[j] <- comp.rel[[j]][alpha.i]
          }
        }
        subset.relation <- rbind(subset.relation, cbind(rep(test.i[i], 
                                                            length(pam.out)), (length(modules.keep) + 1):(length(modules.keep) + 
                                                                                                            length(pam.out))))
        modules.keep <- c(modules.keep, pam.out)
        if (all(is.na(def.alpha) & comp.sig >= module.pvalue)) {
          do.test <- c(do.test, rep(0, length(pam.out)))
        }
        else {
          update.test <- rep(0, length(pam.out))
          update.test[which(!is.na(def.alpha) & comp.sig < 
                              module.pvalue)] <- 1
          do.test <- c(do.test, update.test)
          def.alpha[is.na(def.alpha)] <- alpha.defined[test.i[i]]
        }
        alpha.defined <- c(alpha.defined, def.alpha)
        pvalue.calc <- c(pvalue.calc, comp.sig)
        comp.score <- c(comp.score, comp.sc)
      }
    }
    if (all(do.test == 0)) 
      do.update <- FALSE
    n.iter <- n.iter + 1
  }
  alpha.defined[is.na(alpha.defined)] <- 0
  names(modules.keep) <- 1:length(modules.keep)
  output <- list(modules = modules.keep, singletons = modules.singleton, 
                 module.relation = subset.relation, module.alpha = alpha.defined, 
                 module.pvalue = pvalue.calc, module.compRel = comp.score)
  return(output)
}