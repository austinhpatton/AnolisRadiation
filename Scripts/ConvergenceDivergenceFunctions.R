# This file contains all functions necessaty to conduct the 
# convergence/divergence analyses for the island/mainland anolis
# study (Patton et al.)

# This function will calculate phylogenetically standardized trait 
# dissimilary. It expects SimTraitList to be a list of matrices,
# each representing a single set of simulated traits. Column names (traits)
# should be the same as in empricial data, and each row should be named
# according to the species for which those simulated trait
# values correspond. "EmpiricalTraits" should be an analogous data.frame,
# but of empirical trait values for each species in the tree. 
# THIS ASSUMES ALL TRAITS ARE CONTINUOUS - euclidean distances. 


GetTraitDists <- 
  function(sims = 'SimTraitList' , obs = 'EmpiricalTraits'){
    sim.trait.arr <- array(unlist(sims), dim = c(nrow(sims[[1]]), ncol(sims[[1]]), length(sims)))
    sim.trait.dists <- list()
    
    for(i in 1:length(sims)){
      sim.trait.dists[[i]] <- as.matrix(dist(x = sim.trait.arr[,,i]))
      if(i %% length(sims) / 100 == 0){print(paste0(i, 'Simulations Done'))}
    }
    
    # Now, calculate observed pairwise distances of trait values. 
    obs.dists <- as.matrix(dist(x = obs))
    
    # Now put together into list
    traitDists <- list('simulated' = sim.trait.dists, 
                       'observed' = obs.dists)
    return(traitDists)
  }

# This function calculates the Phylogenetically Standardized Trait Dissimilarity (PSTD)
# metric first introduced in Mazel et al. 2017. This uses as input the trait diststances 
# obtained from the GetTraitDists function above

GetPSTD <- 
  function(traitDists = 'traitDists'){
    
    # Calculate the mean
    mean.dists <- apply(simplify2array(traitDists$simulated), 1:2, mean)
    
    # Calculate the median
    median.dists <- apply(simplify2array(traitDists$simulated), 1:2, median)
    
    # Lastly the standard deviation
    sd.dists <- apply(simplify2array(traitDists$simulated), 1:2, sd)
    
    # Now, calculate observed pairwise distances of trait values. 
    obs.dists <- as.matrix(dist(x = traitDists$observed))
    
    # Now calculate the phylogenetically standardized trait distance
    # of Mazel et al. 2017 (PSTD)
    # (Dist Observed - mean(Dist Neutral))/StDev(Dist Neutral)
    pstd.mean <- (obs.dists - mean.dists)/sd.dists
    pstd.mean[upper.tri(pstd.mean, diag = T)] <- NA
    pstd.median <- (obs.dists - median.dists)/sd.dists
    pstd.median[upper.tri(pstd.median, diag = T)] <- NA
    
    pstd <- list('mean' = pstd.mean, 'median' = pstd.median)
    return(pstd)
  }


# Now prepare to identify actual convergent/divergent
# species pairs. We will correspond PSTD values and species
# pairs to regions

# This function requires as input the result of GetPSTD, as well as
# a dataframe that corresponds species to regions (or to be more general,
# any categoriacal variable of interest). The latter should be formatted 
# such that species are in the first row, region in the second. 

PreparePSTD <- 
  function(pstd = 'PSTD', spp.cats = 'SpeciesCategories'){
    dist.list <- list()
    
    for(m in ls(pstd)){
      # Convert to lists of distances
      # We can then reduce down to top and bottop extremes
      # and replace species names with group names
      # Thus can get at which types of convergences are most common
      dist.list[[m]] <- matrixConvert(as.dist(pstd[[m]]), colname = c("Spp1", "Spp2", "PSTD"))
      # dist.list[[m]] <- dist2list(as.dist(pstd[[m]])) %>% na.omit(.)
      # dist.list[[m]] <- dist.list[[m]][which(dist.list[[m]]$value != 0),]
      #dist.list[as.numeric(dist.list$col) > as.numeric(dist.list$row),]
      
      # Now correspond species ID with region in the pairwise distance dataframe
      dist.list[[m]]$Spp1.Category <- spp.cats[,2][match(dist.list[[m]]$Spp1, 
                                                         spp.cats[,1])]
      dist.list[[m]]$Spp2.Category <- spp.cats[,2][match(dist.list[[m]]$Spp2, 
                                                         spp.cats[,1])]
      
      # Remove rows with NAs due to secondary colonization of islands/mainland. 
      dist.list[[m]] <- na.omit(dist.list[[m]])
      colnames(dist.list[[m]]) <- c('Spp1', 'Spp2', 'PSTD', 
                                    'Spp1.Category', 
                                    'Spp2.Category')
    }
    
    return(dist.list)
  }


# This function, given the formatted PSTD values, a threshold (value < 1 
# to alpha) for identification of putatively convergent taxa, and the 
# complete phylogeny, will identify putative convergence. 
# Because this method cannot distinguish between stasis and true convergence,
# we will set a minimum nodal distance, or a minimum number of nodes 
# intervening putatively convergent taxa. 
GetConvergent <-
  function(pstd.dists = 'PreparedPSTD', threshold = 'Quantile', 
           tree = 'Phylogeny', minDist = 5,
           dist.type = c('PatristicDist', 'NodeDist')){
    bott <- list()
    
    # Get the distances here so you don't have to calculate it each time within the function
    patr.distances <- distTips(tree, tips = "all", method = c("patristic"), 
                               useC = TRUE)
    node.distances <- distTips(tree, tips = "all", method = c("nNodes"), 
                               useC = TRUE)
    
    patr.distances <- matrixConvert(as.dist(patr.distances), 
                                    colname = c("Spp1", "Spp2", "PatristicDist"))
    node.distances <- matrixConvert(as.dist(node.distances), 
                                    colname = c("Spp1", "Spp2", "NodeDist"))
    dists <- cbind(patr.distances, node.distances$NodeDist)
    colnames(dists) <- c('Spp1', 'Spp2', 'PatristicDist', 'NodeDist')
    
    for(m in ls(pstd.dists)){
      
      upper <- 1-threshold
      lower <- 0+threshold
      
      quantiles <- quantile(pstd.dists[[m]]$PSTD, probs = c(lower, upper))
      
      bott[[m]] <- pstd.dists[[m]][pstd.dists[[m]]$PSTD <= quantiles[[1]],]
      
      print(paste0('There are ', nrow(bott[[m]]),
                   ' putative convergences'))
      
      bott[[m]] <- merge(dists, bott[[m]], by = c('Spp1', 'Spp2'))
      colnames(bott[[m]]) <- c('Spp1', 'Spp2', 'PatristicDist', 'NodeDist', 
                               'PSTD', 'Spp1.Category', 'Spp2.Category')
      bott[[m]] <- bott[[m]][which(bott[[m]][[dist.type]] >= minDist),]
      
    }
    return(bott)
  }

GetDivergent <-
  function(pstd.dists = 'PreparedPSTD', threshold = 'Quantile', 
           tree = 'Phylogeny'){
    top <- list()
    
    for(m in ls(pstd.dists)){
      
      upper <- 1-threshold
      lower <- 0+threshold
      
      quantiles <- quantile(pstd.dists[[m]]$PSTD, probs = c(lower, upper))
      
      top[[m]] <- pstd.dists[[m]][pstd.dists[[m]]$PSTD >= quantiles[[2]],]
      
      print(paste0('There are ', nrow(top[[m]]),
                   ' putative divergences'))
      
    }
    
    return(top)
  }


# Now we need to take some more intermediate steps, namely 
# 1) identifying the expected prevalences for each region/category
# 2) determine the proportional species richness for each region/category
# 3) determine the number of convergences expected given the cutoff we used, and
# 4) count the number of convergences/divergences in each group
# 5) determine the number of convergences/divergences observed, 
#    relative to the number expected.

SummarizeConvDiv <- 
  function(conv = 'Convergences', div = 'Divergences',
           cats = 'Categories', dist.list = 'PreparePSTD.Out', 
           threshold = 0.01, phy = phy){
    
    res <- list()
    
    for(m in ls(conv)){
      
      spp.cats <- data.frame('Species' = c(c(as.character(conv[[m]]$Spp1), 
                                             as.character(conv[[m]]$Spp2)),
                                           c(as.character(div[[m]]$Spp1), 
                                             as.character(div[[m]]$Spp2))),
                             'Category' = c(c(as.character(conv[[m]]$Spp1.Category), 
                                              as.character(conv[[m]]$Spp2.Category)),
                                            c(as.character(div[[m]]$Spp1.Category), 
                                              as.character(div[[m]]$Spp2.Category))))
      spp.cats <- distinct(spp.cats)
      
      res[[m]] <- 
        data.frame('ExpPrev' = rep(NA, length(cats)*2),
                   'NumExpConv' = rep(NA, length(cats)*2), 
                   'NumObsConv' = rep(NA, length(cats)*2),
                   'ObsPrev' = rep(NA, length(cats)*2),
                   'ObsOverExpPrev' = rep(NA, length(cats)*2),
                   'ObsMinusExpPrev' = rep(NA, length(cats)*2),
                   'BinomP' = rep(NA, length(cats)*2),
                   'FDR.q' = rep(NA, length(cats)*2),
                   'Category' = rep(NA, length(cats)*2),
                   'Type' = rep(NA, length(cats)*2))
      
      
      # First convergences 
      for(i in 1:(nrow(res[[m]])/2)){
        cat <- cats[i]
        res[[m]][i,1] <- (sum(c(as.character(dist.list[[m]]$Spp1.Category),
                                as.character(dist.list[[m]]$Spp2.Category)) == cat)) /
          length(c(dist.list[[m]]$Spp1.Category, dist.list[[m]]$Spp2.Category))
        res[[m]][i,2] <- res[[m]][i,1] * (nrow(conv[[m]])*2)
        res[[m]][i,3] <- (length(which(conv[[m]]$Spp1.Category == cat)) +
                            (length(which(conv[[m]]$Spp2.Category == cat)))) 
        res[[m]][i,4] <- res[[m]][i,3] / (nrow(conv[[m]])*2)        
        res[[m]][i,5] <- res[[m]][i,4] / res[[m]][i,1]
        res[[m]][i,6] <- ((res[[m]][i,3] / (nrow(conv[[m]]))*2) - 
                            res[[m]][i,1])
        res[[m]][i,7] <- binom.test(x = res[[m]][i,3], n = (nrow(conv[[m]])*2), p = res[[m]][i,1])$p.value
      }
      
      # Then divergences
      for(i in ((nrow(res[[m]])/2)+1):nrow(res[[m]])){
        cat <- cats[i-nrow(res[[m]])/2]
        res[[m]][i,1] <- (sum(c(as.character(dist.list[[m]]$Spp1.Category),
                                as.character(dist.list[[m]]$Spp2.Category)) == cat)) /
          length(c(dist.list[[m]]$Spp1.Category, dist.list[[m]]$Spp2.Category))
        res[[m]][i,2] <- res[[m]][i,1] * (nrow(conv[[m]])*2) 
        res[[m]][i,3] <- (length(which(div[[m]]$Spp1.Category == cat)) +
                            (length(which(div[[m]]$Spp2.Category == cat)))) 
        res[[m]][i,4] <- res[[m]][i,3] / (nrow(div[[m]])*2)        
        res[[m]][i,5] <- res[[m]][i,4] / res[[m]][i,1]
        res[[m]][i,6] <- ((res[[m]][i,3] / (nrow(div[[m]]))*2) - 
                            res[[m]][i,1])
        res[[m]][i,7] <- binom.test(x = res[[m]][i,3], n = (nrow(conv[[m]])*2), p = res[[m]][i,1])$p.value
      }
      
      res[[m]][1:length(cats),8] <- p.adjust(res[[m]][1:length(cats),7], method = 'fdr')
      res[[m]][(length(cats)+1):(length(cats)*2),8] <- p.adjust(res[[m]][(length(cats)+1):(length(cats)*2),7], method = 'fdr')
      res[[m]][,9] <- rep(cats, 2)
      
      res[[m]][,10] <- c(rep('Convergence', length(cats)),
                         rep('Divergence', length(cats)))
    }
    
    return(res)
  }

# This function prepares the convergence/divergence results for plotting/analyzing 
# pairwise convegence/divergences. 

PrepareHeatmap <-
  function(Convergent = 'ConvergenceResults', 
           Divergent = 'DivergenceResults', 
           pstdResults = 'pstdResults',
           pair.types = 'PossiblePairs',
           pair.levels = 'PairLevels',
           distType = 'PatristicDist',
           threshold = 0.01,
           maxDist = 20){
    
    # Get the distances here so you don't have to calculate it each time within the function
    patr.distances <- distTips(tree, tips = "all", method = c("patristic"), 
                               useC = TRUE)
    node.distances <- distTips(tree, tips = "all", method = c("nNodes"), 
                               useC = TRUE)
    
    patr.distances <- matrixConvert(as.dist(patr.distances), 
                                    colname = c("Spp1", "Spp2", "PatristicDist"))
    node.distances <- matrixConvert(as.dist(node.distances), 
                                    colname = c("Spp1", "Spp2", "NodeDist"))
    dists <- cbind(patr.distances, node.distances$NodeDist)
    colnames(dists) <- c('Spp1', 'Spp2', 'PatristicDist', 'NodeDist')
    
    for(m in ls(pstdResults)){
      Convergent[[m]]$pair <- 
        as.factor(paste(Convergent[[m]]$Spp1.Category, 
                        Convergent[[m]]$Spp2.Category, 
                        sep = "--"))
      Divergent[[m]]$pair <- 
        as.factor(paste(Divergent[[m]]$Spp1.Category, 
                        Divergent[[m]]$Spp2.Category, 
                        sep = "--"))
      pstdResults[[m]]$pair <- 
        as.factor(paste(pstdResults[[m]]$Spp1.Category, 
                        pstdResults[[m]]$Spp2.Category, 
                        sep = "--"))
      
      for(i in 1:length(pair.types)){
        Convergent[[m]]$pair <- as.factor(sub(pair.types[[i]][1], 
                                              pair.types[[i]][2], 
                                              x = Convergent[[m]]$pair))
        Divergent[[m]]$pair <- as.factor(sub(pair.types[[i]][1], 
                                             pair.types[[i]][2], 
                                             x = Divergent[[m]]$pair))
        pstdResults[[m]]$pair <- as.factor(sub(pair.types[[i]][1], 
                                               pair.types[[i]][2], 
                                               x = pstdResults[[m]]$pair))
      }
      
      Convergent[[m]]$pair <- factor(Convergent[[m]]$pair,
                                     levels = c(pair.levels))
      Divergent[[m]]$pair <- factor(Divergent[[m]]$pair,
                                    levels = c(pair.levels))
      pstdResults[[m]]$pair <- factor(pstdResults[[m]]$pair,
                                      levels = c(pair.levels))
      
      combined <- merge(pstdResults[[m]], dists, by = c('Spp1', 'Spp2'))
      pair.exp.freqs <- c()
      for(i in 1:length(summary(pstdResults[[m]]$pair))){ 
        pair <- names(summary(pstdResults[[m]]$pair))[i]
        pair <- factor(pair, levels = c(pair.levels))
        pair.exp.freqs[i] <- summary(pstdResults[[m]]$pair)[i] / 
          nrow(pstdResults[[m]])
        names(pair.exp.freqs)[i] <- as.character(pair)
      }
    }
    
    Obs <- 
      data.frame('Pair' = pair.levels,
                 'ConvMed' = summary(Convergent$median$pair),
                 'ConvMedFreq' = (summary(Convergent$median$pair) / 
                                    sum(summary(Convergent$median$pair))),
                 'ConvMean' = summary(Convergent$mean$pair),
                 'ConvMeanFreq' = (summary(Convergent$mean$pair) / 
                                     sum(summary(Convergent$mean$pair))),
                 'DivMed' = summary(Divergent$median$pair),
                 'DivMedFreq' = (summary(Divergent$median$pair) / 
                                   sum(summary(Divergent$median$pair))),
                 'DivMean' = summary(Divergent$mean$pair),
                 'DivMeanFreq' = (summary(Divergent$mean$pair) / 
                                   sum(summary(Divergent$mean$pair)))
      )
    
    Exp <- 
      data.frame(
        # It doesn't matter whether mean or median is chosen - these are constants
        # given your dataset 
        'Pair' = names(summary(pstdResults$mean$pair)),
        'ExpFreq' = pair.exp.freqs,
        'ExpNumConv' = pair.exp.freqs * nrow(Convergent$mean),
        'ExpNumDiv' = pair.exp.freqs * nrow(Divergent$mean)
      )

    # Chi <-
    #   data.frame('Pair' = Obs$Pair,
    #              'ConvMed.Chi' = ((Obs$ConvMedFreq - pair.exp.freqs)^2) / pair.exp.freqs,
    #              'ConvMean.Chi' = ((Obs$ConvMeanFreq - pair.exp.freqs)^2) / pair.exp.freqs,
    #              'DivMed.Chi' = ((Obs$DivMedFreq - pair.exp.freqs)^2) / pair.exp.freqs,
    #              'DivMean.Chi' = ((Obs$DivMeanFreq - pair.exp.freqs)^2) / pair.exp.freqs
    #   )   
    # 
    Binom <- 
      data.frame('Pair' = Obs$Pair,
                 'ConvMed.P' = rep(NA, length(Obs$Pair)),
                 'ConvMean.P' = rep(NA, length(Obs$Pair)),
                 'DivMed.P' = rep(NA, length(Obs$Pair)),
                 'DivMean.P' = rep(NA, length(Obs$Pair)))  
    
    for(i in 1:nrow(Obs)){
      Binom[i,2] <- binom.test(x = Obs[i,2], n = nrow(Convergent$median), p = Exp[i,2], 
                               alternative = 'two.sided')$p.value
      Binom[i,3] <- binom.test(x = Obs[i,4], n = nrow(Convergent$mean), p = Exp[i,2], 
                               alternative = 'two.sided')$p.value
      Binom[i,4] <- binom.test(x = Obs[i,6], n = nrow(Divergent$median), p = Exp[i,2], 
                               alternative = 'two.sided')$p.value
      Binom[i,5] <- binom.test(x = Obs[i,8], n = nrow(Divergent$mean), p = Exp[i,2], 
                               alternative = 'two.sided')$p.value
    }
    
    Fdr <- 
      data.frame('Pair' = Binom$Pair,
                 'ConvMed.q' = rep(NA, length(Binom$Pair)),
                 'ConvMean.q' = rep(NA, length(Binom$Pair)),
                 'DivMed.q' = rep(NA, length(Binom$Pair)),
                 'DivMean.q' = rep(NA, length(Binom$Pair)))  
    Fdr$ConvMed.q <- p.adjust(Binom$ConvMed.P, method = 'fdr')
    Fdr$ConvMean.q <- p.adjust(Binom$ConvMean.P, method = 'fdr')
    Fdr$DivMed.q <- p.adjust(Binom$DivMed.P, method = 'fdr')
    Fdr$DivMean.q <- p.adjust(Binom$DivMean.P, method = 'fdr')
    
    ObsExp <-
      data.frame('Pair' = Obs$Pair,
                 'ExpFreq' = pair.exp.freqs,
                 'ConvMed.ObsExpFreqs' = Obs$ConvMedFreq / pair.exp.freqs,
                 'ConvMean.ObsExpFreqs' = Obs$ConvMeanFreq / pair.exp.freqs,
                 'DivMed.ObsExpFreqs' = Obs$DivMedFreq / pair.exp.freqs,
                 'DivMean.ObsExpFreqs' = Obs$DivMeanFreq / pair.exp.freqs)
    
    ObsExp <- merge(ObsExp, Binom)
    ObsExp <- merge(ObsExp, Fdr)
    
    HeatmapData <- list('Observed' = Obs,
                        'Expected' = Exp,
                        'ObsExp' = ObsExp)
    
    tmp <- unlist(strsplit(as.character(HeatmapData$ObsExp$Pair), split = '--'))
    
    HeatmapData$ObsExp$Region1 <- tmp[seq(from = 1, to = length(tmp)-1, by = 2)]
    HeatmapData$ObsExp$Region2 <- tmp[seq(from = 2, to = length(tmp), by = 2)]
    
    
    return(HeatmapData)
  }
