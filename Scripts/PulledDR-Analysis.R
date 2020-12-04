library(castor)
library(ape)
library(dispRity)
library(reshape2)
library(matrixStats)

setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/')
all.spp <- 
  read.csv('./FinalData/Male-ALL-DATA-Trait-Ecol-Regions-FINAL.csv')

drop.spp <- 
  c('agassizi', 'bicaorum', 'lineatus', 'roatanensis', 
    'concolor', 'pinchoti', 'utilensis', 'townsendi')
###########################
# Variable rates function #
###########################

fit_pdr_variable_grid <- function(tree, rho=1, starting_grid_size, age0=0,
                                  ntry_fit, ntry_search,
                                  nboot, ntry_boot, nthreads, 
                                  max_time, to_return="sim"){
  max_age <- get_tree_span(tree)$max_distance
  n_grid <- starting_grid_size
  LTT <- count_lineages_through_time(tree, Ntimes=1000)
  age_grid_i <- castor:::get_inhomogeneous_grid_1D(Xstart=0.000001, Xend=max_age, 
                                                   Ngrid=n_grid, densityX=rev(max_age-LTT$times),
                                                   densityY=sqrt(rev(LTT$lineages)))
  age_grid_i[1] <- 0
  f_i <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid_i, age0=age0, Ntrials = ntry_search,
                             min_PDR = -5, max_PDR=5, Nthreads = nthreads, condition = "crown",
                             max_model_runtime = max_time)
  age_grid_j <- castor:::get_inhomogeneous_grid_1D(Xstart=0.000001, Xend=max_age, 
                                                   Ngrid=n_grid+1, densityX=rev(max_age-LTT$times),
                                                   densityY=sqrt(rev(LTT$lineages)))
  age_grid_j[1] <- 0
  f_j <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid_j, age0=age0, Ntrials = ntry_search,
                             min_PDR = -5, max_PDR=5, Nthreads = nthreads, condition = "crown",
                             max_model_runtime = max_time)
  aic_i <- f_i$AIC
  aic_j <- f_j$AIC
  ## set threshold of AIC 
  while (aic_j <= aic_i){
    aic_i <- aic_j
    n_grid <- n_grid + 1
    age_grid_j <- castor:::get_inhomogeneous_grid_1D(Xstart=0.000001, Xend=max_age, 
                                                     Ngrid=n_grid+1, densityX=rev(max_age-LTT$times),
                                                     densityY=sqrt(rev(LTT$lineages)))
    age_grid_j[1] <- 0
    f_j <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid_j, age0=age0, Ntrials = ntry_search,
                               min_PDR = -5, max_PDR=5, Nthreads = nthreads, condition = "crown",
                               max_model_runtime = max_time)
    aic_j <- f_j$AIC
  }
  ## use the estimated age grid 
  age_grid <- castor:::get_inhomogeneous_grid_1D(Xstart=0.000001, Xend=max_age, 
                                                 Ngrid=n_grid, densityX=rev(max_age-LTT$times),
                                                 densityY=sqrt(rev(LTT$lineages)))
  age_grid[1] <- 0
  if (is.null(nboot)) { ## if nboot is null just return the model fitted
    f <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid, age0=age0, Ntrials = ntry_fit,
                             min_PDR = -5, max_PDR=5, Nthreads = nthreads, condition = "crown",
                             Nbootstraps = 0, Ntrials_per_bootstrap = ntry_boot,
                             max_model_runtime = max_time)
    if (to_return == "fit"){
      return(list(fit=f)) ## return just the fitted object
    } else {
      ntip <- length(tree$tip.label)
      pdr_mle <- simulate_deterministic_hbd(LTT0=ntip, oldest_age=max_age, 
                                            age_grid = f$age_grid, lambda0=f$fitted_rholambda0/rho,
                                            mu=0, PDR=f$fitted_PDR, splines_degree = 1)
      return(list(fit = f, pdr_mle = pdr_mle)) ## return both the fit and imputed pdr
    }
  } else { ## bootstrap the data
    f_boot <- fit_hbd_pdr_on_grid(tree, age_grid = age_grid, age0=age0, Ntrials = ntry_fit,
                                  min_PDR = -5, max_PDR=5, Nthreads = nthreads, condition = "crown",
                                  Nbootstraps = nboot, Ntrials_per_bootstrap = ntry_boot,
                                  max_model_runtime = max_time)
    if (to_return == "fit"){
      return(list(fit=f_boot)) ## return just the fitted object
    } else {
      ntip <- length(tree$tip.label)
      pdr_mle <- simulate_deterministic_hbd(LTT0=ntip, oldest_age=max_age, rho=rho,
                                            age_grid = f_boot$age_grid, lambda0=f_boot$fitted_rholambda0/rho,
                                            mu=0, PDR=f_boot$fitted_PDR, splines_degree = 1)
      pdr_lower <- simulate_deterministic_hbd(LTT0=ntip, oldest_age=max_age,rho=rho,
                                              age_grid = f_boot$age_grid, lambda0=f_boot$CI95lower$rholambda0/rho,
                                              mu=0, PDR=f_boot$CI95lower$PDR, splines_degree = 1)
      pdr_upper <- simulate_deterministic_hbd(LTT0=ntip, oldest_age=max_age, rho=rho,
                                              age_grid = f_boot$age_grid, lambda0=f_boot$CI95upper$rholambda0/rho,
                                              mu=0, PDR=f_boot$CI95upper$PDR, splines_degree = 1)
      ## return the fitted object, the imputed pdr and the imputed CIs
      return(list(fit = f_boot, pdr_mle = pdr_mle, pdr_lower = pdr_lower, pdr_upper = pdr_upper)) 
    }
  } 
}
## Misc functions for running analyses
## starting grid size function
sgs_f <- function(tree)
  1 + round(log(length(tree$tip.label)))
## function for setting the maximum time
maxt_f <- function(tree)
  max(1,length(tree$tip.label)/ 5e4)

############################################################################################################

#I2
i2 <- read.tree('./FinalData/Island_Clade2_TimeTree.tree')
i2 <- drop.tip(i2, drop.spp)

grid <- seq.int(from = 0, to = max(node.depth.edgelength(i2)), length = 25)
i2.rho <- length(i2$tip.label) / nrow(all.spp[which(all.spp$Group == 'I2'),])

##############

ntry_fit <- 10
nthreads <- 10
nboot <- 0
ntry_boot <- 0 
ntry_search <- 10

## Starting parameters
sgs <- sgs_f(i2)
maxt <- maxt_f(i2)

##  Fit variable rate models (with bootstraps) 
vrm.i2  <- fit_pdr_variable_grid(i2, rho=i2.rho, starting_grid_size=sgs, age0=0,
                              ntry_fit=ntry_fit, ntry_search=ntry_search,
                              nboot=nboot, ntry_boot=ntry_boot, nthreads=nthreads,
                              max_time = maxt)
#M1
m1 <- read.tree('./FinalData/Mainland_Clade1_TimeTree.tree')
m1 <- drop.tip(m1, drop.spp)

grid <- seq.int(from = 5, to = max(node.depth.edgelength(m1)), length = 5)
m1.rho <- length(m1$tip.label) / nrow(all.spp[which(all.spp$Group == 'M1'),])

##############

ntry_fit <- 10
nthreads <- 5
nboot <- 0
ntry_boot <- 0 
ntry_search <- 10

## Starting parameters
sgs <- sgs_f(m1)
maxt <- maxt_f(m1)

##  Fit variable rate models (with bootstraps) 
vrm.m1  <- fit_pdr_variable_grid(m1, rho=m1.rho, starting_grid_size=sgs, age0=0,
                              ntry_fit=ntry_fit, ntry_search=ntry_search,
                              nboot=nboot, ntry_boot=ntry_boot, nthreads=nthreads,
                              max_time = maxt)
#M2
m2 <- read.tree('./FinalData/Mainland_Clade2_TimeTree.tree')
m2 <- drop.tip(m2, drop.spp)

grid <- round(seq.int(from = 5, to = max(node.depth.edgelength(m2)), length = 25), 4)
m2.rho <- length(m2$tip.label) / nrow(all.spp[which(all.spp$Group == 'M2'),])

##############

ntry_fit <- 10
nthreads <- 10
nboot <- 0
ntry_boot <- 0 
ntry_search <- 10

## Starting parameters
sgs <- sgs_f(m2)
maxt <- maxt_f(m2)

##  Fit variable rate models (with bootstraps) 
vrm.m2  <- fit_pdr_variable_grid(m2, rho=m2.rho, starting_grid_size=sgs, age0=0,
                                 ntry_fit=ntry_fit, ntry_search=ntry_search,
                                 nboot=nboot, ntry_boot=ntry_boot, nthreads=nthreads,
                                 max_time = maxt)


#M2a
m2a <- read.tree('./FinalData/Mainland_Clade2a_TimeTree.tree')
m2a <- drop.tip(m2a, drop.spp)

grid <- round(seq.int(from = 5, to = max(node.depth.edgelength(m2a)), length = 25), 4)
m2a.rho <- length(m2a$tip.label) / nrow(all.spp[which(all.spp$FineGroup == 'M2a'),])

##############

ntry_fit <- 10
nthreads <- 10
nboot <- 0
ntry_boot <- 0 
ntry_search <- 10

## Starting parameters
sgs <- sgs_f(m2a)
maxt <- maxt_f(m2a)

##  Fit variable rate models (with bootstraps) 
vrm.m2a  <- fit_pdr_variable_grid(m2a, rho=m2a.rho, starting_grid_size=sgs, age0=0,
                                  ntry_fit=ntry_fit, ntry_search=ntry_search,
                                  nboot=nboot, ntry_boot=ntry_boot, nthreads=nthreads,
                                  max_time = maxt)


#M2b
m2b <- read.tree('./FinalData/Mainland_Clade2b_TimeTree.tree')
m2b <- drop.tip(m2b, drop.spp)

grid <- round(seq.int(from = 5, to = max(node.depth.edgelength(m2b)), length = 25), 4)
m2b.rho <- length(m2b$tip.label) / nrow(all.spp[which(all.spp$FineGroup == 'M2b'),])

##############

ntry_fit <- 10
nthreads <- 10
nboot <- 0
ntry_boot <- 0 
ntry_search <- 10

## Starting parameters
sgs <- sgs_f(m2b)
maxt <- maxt_f(m2b)

##  Fit variable rate models (with bootstraps) 
vrm.m2b  <- fit_pdr_variable_grid(m2b, rho=m2b.rho, starting_grid_size=sgs, age0=0,
                                  ntry_fit=ntry_fit, ntry_search=ntry_search,
                                  nboot=nboot, ntry_boot=ntry_boot, nthreads=nthreads,
                                  max_time = maxt)


#All Species
anolis <- read.tree('./FinalData/Anole_TimeCal_NoOut.tree')
anolis <- drop.tip(anolis, drop.spp)

grid <- round(seq.int(from = 5, to = max(node.depth.edgelength(anolis)), length = 25), 4)
anolis.rho <- length(anolis$tip.label) / nrow(all.spp)

##############

ntry_fit <- 10
nthreads <- 10
nboot <- 0
ntry_boot <- 0 
ntry_search <- 10

## Starting parameters
sgs <- sgs_f(anolis)
maxt <- maxt_f(anolis)

##  Fit variable rate models (with bootstraps) 
vrm.anolis  <- fit_pdr_variable_grid(anolis, rho=anolis.rho, starting_grid_size=sgs, age0=0,
                                     ntry_fit=ntry_fit, ntry_search=ntry_search,
                                     nboot=nboot, ntry_boot=ntry_boot, nthreads=nthreads,
                                     max_time = maxt)




# Now simulate sets of trees to plot as the bootstrap distribution
i2.sim.pdr.trees <-
  generate_tree_hbd_reverse(Ntips = length(i2$tip),
                            age_grid = vrm.i2$fit$age_grid,
                            PDR = vrm.i2$fit$fitted_PDR,
                            crown_age = max(tree.age(i2)[,1]),
                            rholambda0 = vrm.i2$fit$fitted_rholambda0,
                            Ntrees = 1000,
                            rho = i2.rho)

m1.sim.pdr.trees <-
  generate_tree_hbd_reverse(Ntips = length(m1$tip),
                            age_grid = vrm.m1$fit$age_grid,
                            PDR = vrm.m1$fit$fitted_PDR,
                            crown_age = max(tree.age(m1)[,1]),
                            rholambda0 = vrm.m1$fit$fitted_rholambda0,
                            Ntrees = 1000,
                            rho = m1.rho)

m2.sim.pdr.trees <-
  generate_tree_hbd_reverse(Ntips = length(m2$tip),
                            age_grid = round(vrm.m2$fit$age_grid, 3),
                            PDR = vrm.m2$fit$fitted_PDR,
                            crown_age = max(tree.age(m2)[,1]),
                            rholambda0 = vrm.m2$fit$fitted_rholambda0,
                            Ntrees = 1000,
                            rho = m2.rho)

m2a.sim.pdr.trees <-
  generate_tree_hbd_reverse(Ntips = length(m2a$tip),
                            age_grid = round(vrm.m2a$fit$age_grid, 3),
                            PDR = vrm.m2a$fit$fitted_PDR,
                            crown_age = max(tree.age(m2a)[,1]),
                            rholambda0 = vrm.m2a$fit$fitted_rholambda0,
                            Ntrees = 1000,
                            rho = m2a.rho)
m2b.sim.pdr.trees <-
  generate_tree_hbd_reverse(Ntips = length(m2b$tip),
                            age_grid = round(vrm.m2b$fit$age_grid, 3),
                            PDR = vrm.m2b$fit$fitted_PDR,
                            crown_age = max(tree.age(m2b)[,1]),
                            rholambda0 = vrm.m2b$fit$fitted_rholambda0,
                            Ntrees = 1000,
                            rho = m2b.rho)
anolis.sim.pdr.trees <-
  generate_tree_hbd_reverse(Ntips = length(anolis$tip),
                            age_grid = round(vrm.anolis$fit$age_grid, 3),
                            PDR = vrm.anolis$fit$fitted_PDR,
                            crown_age = max(tree.age(anolis)[,1]),
                            rholambda0 = vrm.anolis$fit$fitted_rholambda0,
                            Ntrees = 1000,
                            rho = anolis.rho)


# Now estimate the time variable PDR for every simulated M1 tree
m1.sim.pdr.fit <- data.frame('Time' = vrm.m1$fit$age_grid)
for(phy in 1:length(m1.sim.pdr.trees$trees)){
  print(paste0('Fitting PDR for Simulated Tree #', phy))
  maxt <- maxt_f(m1.sim.pdr.trees$trees[[phy]])

  fit <-
    fit_hbd_pdr_on_grid(m1.sim.pdr.trees$trees[[phy]],
                        fixed_rholambda0 = vrm.m1$fit$fitted_rholambda0,
                        age_grid = vrm.m1$fit$age_grid, age0 = 0,
                        max_model_runtime = maxt,
                        Nthreads = nthreads,
                        Ntrials = ntry_fit,
                        Nbootstraps = 0)

  if(fit$success == 'TRUE'){

    m1.sim.pdr.fit[,phy+1] <- fit$fitted_PDR
  }else(m2.sim.pdr.fit[,phy+1] <- rep(NA, length(vrm.m2$fit$age_grid)))

  rm(fit)

  colnames(m1.sim.pdr.fit)[phy+1] <- paste0('Sim', phy)
}

#write.csv(m1.sim.pdr.fit, 'M1-PDR-Sim-Fits.csv', row.names = F)
m1.sim.pdr.fit <- read.csv('M1-PDR-Sim-Fits.csv')
m1.sim.pdr.fit <- m1.sim.pdr.fit[,seq(1, ncol(m1.sim.pdr.fit), 3)]

m1.sim.pdr.fit <- m1.sim.pdr.fit[,which(m1.sim.pdr.fit[5,] > -10)]
m1.good.sims <- c()
for(i in 2:ncol(m1.sim.pdr.fit)){
  m1.good.sims[i] <- any(dist(m1.sim.pdr.fit[,i]) > 10)
}
if(sum(na.omit(m1.good.sims)) > 0){
  m1.sim.pdr.fit <- m1.sim.pdr.fit[,-which(m1.good.sims == T)]
}

m1.sim.pdr.fit$Time <- as.factor(m1.sim.pdr.fit$Time)
m1.sim.pdr.fit <- melt(m1.sim.pdr.fit)
m1.sim.pdr.fit$Time <- as.numeric(as.character(m1.sim.pdr.fit$Time))

# And every I2 tree
i2.sim.pdr.fit <- data.frame('Time' = vrm.i2$fit$age_grid)
for(phy in 1:length(i2.sim.pdr.trees$trees)){
  print(paste0('Fitting PDR for Simulated Tree #', phy))
  maxt <- maxt_f(i2.sim.pdr.trees$trees[[phy]])
  
  fit <- 
    fit_hbd_pdr_on_grid(i2.sim.pdr.trees$trees[[phy]], 
                        fixed_rholambda0 = vrm.i2$fit$fitted_rholambda0,
                        age_grid = vrm.i2$fit$age_grid, age0 = 0,
                        max_model_runtime = maxt,
                        Nthreads = nthreads,
                        Ntrials = ntry_fit,
                        Nbootstraps = 0)
  
  if(fit$success == 'TRUE'){
    i2.sim.pdr.fit[,phy+1] <- fit$fitted_PDR
  }else(i2.sim.pdr.fit[,phy+1] <- rep(NA, length(vrm.i2$fit$age_grid)))
  
  rm(fit)
  
  colnames(i2.sim.pdr.fit)[phy+1] <- paste0('Sim', phy)
}

#write.csv(i2.sim.pdr.fit, 'I2-PDR-Sim-Fits.csv', row.names = F)
i2.sim.pdr.fit <- read.csv('I2-PDR-Sim-Fits.csv')
i2.sim.pdr.fit <- i2.sim.pdr.fit[,seq(1, ncol(i2.sim.pdr.fit), 3)]

i2.sim.pdr.fit <- i2.sim.pdr.fit[,which(i2.sim.pdr.fit[5,] > -10)]
i2.good.sims <- c()
for(i in 2:ncol(i2.sim.pdr.fit)){
  i2.good.sims[i] <- any(dist(i2.sim.pdr.fit[,i]) > 10)
}
if(sum(na.omit(i2.good.sims)) > 0){
  i2.sim.pdr.fit <- i2.sim.pdr.fit[,-which(i2.good.sims == T)]
}

i2.sim.pdr.fit$Time <- as.factor(i2.sim.pdr.fit$Time)
i2.sim.pdr.fit <- melt(i2.sim.pdr.fit)
i2.sim.pdr.fit$Time <- as.numeric(as.character(i2.sim.pdr.fit$Time))


# And every M2 tree
m2.sim.pdr.fit <- data.frame('Time' = vrm.m2$fit$age_grid)
for(phy in 1:length(m2.sim.pdr.trees$trees)){
  print(paste0('Fitting PDR for Simulated Tree #', phy))
  maxt <- maxt_f(m2.sim.pdr.trees$trees[[phy]])
  
  fit <- 
    fit_hbd_pdr_on_grid(m2.sim.pdr.trees$trees[[phy]], 
                        fixed_rholambda0 = vrm.m2$fit$fitted_rholambda0,
                        age_grid = round(vrm.m2$fit$age_grid, 3), age0 = 0,
                        max_model_runtime = maxt,
                        Nthreads = nthreads,
                        Ntrials = ntry_fit,
                        Nbootstraps = 0)
  
  if(fit$success == 'TRUE'){
    m2.sim.pdr.fit[,phy+1] <- fit$fitted_PDR
  }else(m2.sim.pdr.fit[,phy+1] <- rep(NA, length(vrm.m2$fit$age_grid)))
  
  rm(fit)
  
  colnames(m2.sim.pdr.fit)[phy+1] <- paste0('Sim', phy)
}

#write.csv(m2.sim.pdr.fit, 'M2-PDR-Sim-Fits.csv', row.names = F)
m2.sim.pdr.fit <- read.csv('M2-PDR-Sim-Fits.csv')
m2.sim.pdr.fit <- m2.sim.pdr.fit[,seq(1, ncol(m2.sim.pdr.fit), 3)]

m2.sim.pdr.fit <- m2.sim.pdr.fit[,which(m2.sim.pdr.fit[5,] > -10)]
m2.good.sims <- c()
for(i in 2:ncol(m2.sim.pdr.fit)){
  m2.good.sims[i] <- any(dist(m2.sim.pdr.fit[,i]) > 10)
}
if(sum(na.omit(m2.good.sims)) > 0){
  m2.sim.pdr.fit <- m2.sim.pdr.fit[,-which(m2.good.sims == T)]
}

m2.sim.pdr.fit$Time <- as.factor(m2.sim.pdr.fit$Time)
m2.sim.pdr.fit <- melt(m2.sim.pdr.fit, measure.vars = colnames(m2.sim.pdr.fit)[-1])
m2.sim.pdr.fit$Time <- as.numeric(as.character(m2.sim.pdr.fit$Time))

# And every M2a tree
m2a.sim.pdr.fit <- data.frame('Time' = vrm.m2a$fit$age_grid)
for(phy in 1:length(m2a.sim.pdr.trees$trees)){
  print(paste0('Fitting PDR for Simulated Tree #', phy))
  maxt <- maxt_f(m2a.sim.pdr.trees$trees[[phy]])
  
  fit <- 
    fit_hbd_pdr_on_grid(m2a.sim.pdr.trees$trees[[phy]], 
                        fixed_rholambda0 = vrm.m2a$fit$fitted_rholambda0,
                        age_grid = round(vrm.m2a$fit$age_grid, 3), age0 = 0,
                        max_model_runtime = maxt,
                        Nthreads = nthreads,
                        Ntrials = ntry_fit,
                        Nbootstraps = 0)
  
  if(fit$success == 'TRUE'){
    m2a.sim.pdr.fit[,phy+1] <- fit$fitted_PDR
  }else(m2a.sim.pdr.fit[,phy+1] <- rep(NA, length(vrm.m2a$fit$age_grid)))
  
  rm(fit)
  
  colnames(m2a.sim.pdr.fit)[phy+1] <- paste0('Sim', phy)
}

#write.csv(m2a.sim.pdr.fit, 'M2a-PDR-Sim-Fits.csv', row.names = F)
m2a.sim.pdr.fit <- read.csv('M2a-PDR-Sim-Fits.csv')
m2a.sim.pdr.fit <- m2a.sim.pdr.fit[,seq(1, ncol(m2a.sim.pdr.fit), 3)]

m2a.sim.pdr.fit <- m2a.sim.pdr.fit[,which(m2a.sim.pdr.fit[5,] > -10)]
m2a.good.sims <- c()
for(i in 2:ncol(m2a.sim.pdr.fit)){
  m2a.good.sims[i] <- any(dist(m2a.sim.pdr.fit[,i]) > 10)
}
if(sum(na.omit(m2a.good.sims)) > 0){
  m2a.sim.pdr.fit <- m2a.sim.pdr.fit[,-which(m2a.good.sims == T)]
}

m2a.quants <- rowQuantiles(as.matrix(m2a.sim.pdr.fit), probs = c(0.025, 0.975), na.rm = T)

m2a.sim.pdr.fit$Time <- as.factor(m2a.sim.pdr.fit$Time)
m2a.sim.pdr.fit <- melt(m2a.sim.pdr.fit, measure.vars = colnames(m2a.sim.pdr.fit)[-1])
m2a.sim.pdr.fit$Time <- as.numeric(as.character(m2a.sim.pdr.fit$Time))


# And every M2b tree
m2b.sim.pdr.fit <- data.frame('Time' = vrm.m2b$fit$age_grid)
for(phy in 1:length(m2b.sim.pdr.trees$trees)){
  print(paste0('Fitting PDR for Simulated Tree #', phy))
  maxt <- maxt_f(m2b.sim.pdr.trees$trees[[phy]])
  
  fit <- 
    fit_hbd_pdr_on_grid(m2b.sim.pdr.trees$trees[[phy]], 
                        fixed_rholambda0 = vrm.m2b$fit$fitted_rholambda0,
                        age_grid = round(vrm.m2b$fit$age_grid, 3), age0 = 0,
                        max_model_runtime = maxt,
                        Nthreads = nthreads,
                        Ntrials = ntry_fit,
                        Nbootstraps = 0)
  
  if(fit$success == 'TRUE'){
    m2b.sim.pdr.fit[,phy+1] <- fit$fitted_PDR
  }else(m2b.sim.pdr.fit[,phy+1] <- rep(NA, length(vrm.m2b$fit$age_grid)))
  
  rm(fit)
  
  colnames(m2b.sim.pdr.fit)[phy+1] <- paste0('Sim', phy)
}

#write.csv(m2b.sim.pdr.fit, 'M2b-PDR-Sim-Fits.csv', row.names = F)
m2b.sim.pdr.fit <- read.csv('M2b-PDR-Sim-Fits.csv')
m2b.sim.pdr.fit <- m2b.sim.pdr.fit[,seq(1, ncol(m2b.sim.pdr.fit), 3)]

m2b.sim.pdr.fit <- m2b.sim.pdr.fit[,which(m2b.sim.pdr.fit[5,] > -10)]
m2b.good.sims <- c()
for(i in 2:ncol(m2b.sim.pdr.fit)){
  m2b.good.sims[i] <- any(dist(m2b.sim.pdr.fit[,i]) > 10)
}
if(sum(na.omit(m2b.good.sims)) > 0){
  m2b.sim.pdr.fit <- m2b.sim.pdr.fit[,-which(m2b.good.sims == T)]
}

m2b.quants <- rowQuantiles(as.matrix(m2b.sim.pdr.fit), probs = c(0.025, 0.975), na.rm = T)

m2b.sim.pdr.fit$Time <- as.factor(m2b.sim.pdr.fit$Time)
m2b.sim.pdr.fit <- melt(m2b.sim.pdr.fit, measure.vars = colnames(m2b.sim.pdr.fit)[-1])
m2b.sim.pdr.fit$Time <- as.numeric(as.character(m2b.sim.pdr.fit$Time))


# And finally all anolis
anolis.sim.pdr.fit <- data.frame('Time' = vrm.anolis$fit$age_grid)
for(phy in 1:length(anolis.sim.pdr.trees$trees)){
  print(paste0('Fitting PDR for Simulated Tree #', phy))
  maxt <- maxt_f(anolis.sim.pdr.trees$trees[[phy]])
  
  fit <- 
    fit_hbd_pdr_on_grid(anolis.sim.pdr.trees$trees[[phy]], 
                        fixed_rholambda0 = vrm.anolis$fit$fitted_rholambda0,
                        age_grid = round(vrm.anolis$fit$age_grid, 3), age0 = 0,
                        max_model_runtime = maxt,
                        Nthreads = nthreads,
                        Ntrials = ntry_fit,
                        Nbootstraps = 0)
  
  if(fit$success == 'TRUE'){
    anolis.sim.pdr.fit[,phy+1] <- fit$fitted_PDR
  }else(anolis.sim.pdr.fit[,phy+1] <- rep(NA, length(vrm.anolis$fit$age_grid)))
  
  rm(fit)
  
  colnames(anolis.sim.pdr.fit)[phy+1] <- paste0('Sim', phy)
}

#write.csv(anolis.sim.pdr.fit, 'Anolis-PDR-Sim-Fits.csv', row.names = F)
anolis.sim.pdr.fit <- read.csv('Anolis-PDR-Sim-Fits.csv')
anolis.sim.pdr.fit <- anolis.sim.pdr.fit[,seq(1, ncol(anolis.sim.pdr.fit), 3)]

anolis.sim.pdr.fit <- anolis.sim.pdr.fit[,which(anolis.sim.pdr.fit[5,] > -10)]
anolis.good.sims <- c()
for(i in 2:ncol(anolis.sim.pdr.fit)){
  anolis.good.sims[i] <- any(dist(anolis.sim.pdr.fit[,i]) > 10)
}
if(sum(na.omit(anolis.good.sims)) > 0){
  anolis.sim.pdr.fit <- anolis.sim.pdr.fit[,-which(anolis.good.sims == T)]
}

anolis.sim.pdr.fit$Time <- as.factor(anolis.sim.pdr.fit$Time)
anolis.sim.pdr.fit <- melt(anolis.sim.pdr.fit, measure.vars = colnames(anolis.sim.pdr.fit)[-1])
anolis.sim.pdr.fit$Time <- as.numeric(as.character(anolis.sim.pdr.fit$Time))

# 
# cols <- c("black", "#d72027", "#add8a4", "#2b83bb")
# 
# combined.sim <- rbind(i1.sim.pdr.fit, m1.sim.pdr.fit, i2.sim.pdr.fit)
# combined.sim$Group <- c(rep('I1', nrow(i1.sim.pdr.fit)),
#                         rep('M1', nrow(m1.sim.pdr.fit)),
#                         rep('I2', nrow(i2.sim.pdr.fit)))

cols <- c('All Anolis' = 'black',  'M1' = '#add8a4', 'I2' = '#d72027', 'M2' = '#2b83bb')


c("#af8dc3", "#b2182b", "#ef8a62", 
  "#2166ac","#67a9cf")
main.cols <- c('M2' = '#2b83bb', 'M2a' = '#8da0cb', 'M2b' = '#fc8d62')

# base_breaks <- function(n = 10){
#   function(x) {
#     axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
#   }
# }

anolis.sim.pdr.fit <- anolis.sim.pdr.fit[which(anolis.sim.pdr.fit$Time != 0),]
m1.sim.pdr.fit <- m1.sim.pdr.fit[which(m1.sim.pdr.fit$Time != 0),]
m2.sim.pdr.fit <- m2.sim.pdr.fit[which(m2.sim.pdr.fit$Time != 0),]
i2.sim.pdr.fit <- i2.sim.pdr.fit[which(i2.sim.pdr.fit$Time != 0),]
m2a.sim.pdr.fit <- m2a.sim.pdr.fit[which(m2a.sim.pdr.fit$Time != 0),]
m2b.sim.pdr.fit <- m2b.sim.pdr.fit[which(m2b.sim.pdr.fit$Time != 0),]

smooth <- 
  ggplot() +
  # geom_line(alpha = 0.025, data = i1.sim.pdr.fit, aes(x = Time, y = value, group = variable),
  #           color = '#faae61') +
  # geom_line(alpha = 0.025, data = anolis.sim.pdr.fit, aes(x = Time, y = value, group = variable, color = 'All Anolis'),
  #           span = 0.5, se = F, size = 0.75, stat = 'smooth') +
  geom_point(alpha = 0.025, data = anolis.sim.pdr.fit, aes(x = Time, y = value, group = variable, color = 'All Anolis'),
             size = 1.2) +
  # geom_line(alpha = 0.05, data = m1.sim.pdr.fit, aes(x = Time, y = value, group = variable, color = 'M1'),
  #           span = 0.5, se = F, size = 0.75, stat = 'smooth') +
  geom_point(alpha = 0.05, data = m1.sim.pdr.fit, aes(x = Time, y = value, group = variable, color = 'M1'),
             size = 1.2) +
  # geom_line(alpha = 0.025, data = i2.sim.pdr.fit, aes(x = Time, y = value, group = variable, color = 'I2'),
  #           span = 0.5, se = F, size = 0.75, stat = 'smooth') +
  geom_point(alpha = 0.025, data = i2.sim.pdr.fit, aes(x = Time, y = value, group = variable, color = 'I2'),
             size = 1.2) +
  # geom_line(alpha = 0.025, data = m2.sim.pdr.fit, aes(x = Time, y = value, group = variable, color = 'M2'),
  #           span = 0.5, se = F, size = 0.75, stat = 'smooth') +
  geom_point(alpha = 0.025, data = m2.sim.pdr.fit, aes(x = Time, y = value, group = variable, color = 'M2'),
            size = 1.2) +
  geom_line(aes(x = vrm.m1$fit$age_grid[-1], y = vrm.m1$fit$fitted_PDR[-1]),
            color = 'black', size = 2.5, span = 0.5, se = F, stat = 'smooth') +
  geom_line(aes(x = vrm.anolis$fit$age_grid[-1], y = vrm.anolis$fit$fitted_PDR[-1], color = 'All Anolis'), size = 2.5, 
            span = 0.5, se = F, stat = 'smooth') +
  geom_line(aes(x = vrm.m1$fit$age_grid[-1], y = vrm.m1$fit$fitted_PDR[-1], color = 'M1'), size = 2, 
            span = 0.5, se = F, stat = 'smooth') +
  geom_line(aes(x = vrm.i2$fit$age_grid[-1], y = vrm.i2$fit$fitted_PDR[-1]),
            color = 'black', size = 2.5, span = 0.5, se = F, stat = 'smooth') +
  geom_line(aes(x = vrm.i2$fit$age_grid[-1], y = vrm.i2$fit$fitted_PDR[-1], color = 'I2'), size = 2, 
            span = 0.5, se = F, stat = 'smooth') +
  geom_line(aes(x = vrm.m2$fit$age_grid[-1], y = vrm.m2$fit$fitted_PDR[-1]),
            color = 'black', size = 2.5, span = 0.5, se = F, stat = 'smooth') +
  geom_line(aes(x = vrm.m2$fit$age_grid[-1], y = vrm.m2$fit$fitted_PDR[-1], color = 'M2'), size = 2, 
            span = 0.5, se = F, stat = 'smooth') +
  coord_cartesian(ylim = c(-0.5, 1.5),
                  xlim = c(50, 0)) +
  labs(x = 'Mybp', y = 'Pulled DR', color = 'Group') +
  scale_color_manual(values = cols) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.8)) + 
  scale_x_reverse()

main.smooth <- 
  ggplot() +
  # geom_line(alpha = 0.025, data = m2.sim.pdr.fit, aes(x = Time, y = value, group = variable, color = 'M2'),
  #           span = 0.5, se = F, size = 0.75, stat = 'smooth') +
  # geom_line(alpha = 0.05, data = m2a.sim.pdr.fit, aes(x = Time, y = value, group = variable, color = 'M2a'),
  #           span = 0.5, se = F, size = 0.75, stat = 'smooth') +
  # geom_line(alpha = 0.05, data = m2b.sim.pdr.fit, aes(x = Time, y = value, group = variable, color = 'M2b'),
  #           span = 0.5, se = F, size = 0.75, stat = 'smooth') +
  geom_point(alpha = 0.05, data = m2a.sim.pdr.fit, aes(x = Time, y = value, group = variable, color = 'M2a'),
             size = 1.2) +
  geom_point(alpha = 0.05, data = m2b.sim.pdr.fit, aes(x = Time, y = value, group = variable, color = 'M2b'),
             size = 1.2) +
  # geom_line(aes(x = vrm.m2$fit$age_grid, y = vrm.m2$fit$fitted_PDR),
  #           color = 'black', size = 2.5, span = 0.5, se = F, stat = 'smooth') +
  # geom_line(aes(x = vrm.m2$fit$age_grid, y = vrm.m2$fit$fitted_PDR, color = 'M2'), size = 2, 
  #           span = 0.5, se = F, stat = 'smooth') +
  geom_line(aes(x = vrm.m2a$fit$age_grid[-1], y = vrm.m2a$fit$fitted_PDR[-1]),
            color = 'black', size = 2.5, span = 0.5, se = F, stat = 'smooth') +
  geom_line(aes(x = vrm.m2a$fit$age_grid[-1], y = vrm.m2a$fit$fitted_PDR[-1], color = 'M2a'), size = 2, 
            span = 0.5, se = F, stat = 'smooth') +
  # geom_line(aes(x = vrm.m2a$fit$age_grid, y = m2a.quants[,1]), size = 1.5, color = 'black',
  #           span = 0.5, se = F, stat = 'smooth') +
  # geom_line(aes(x = vrm.m2a$fit$age_grid, y = m2a.quants[,2]), size = 1.5, color = 'black',
  #           span = 0.5, se = F, stat = 'smooth') +
  # geom_line(aes(x = vrm.m2a$fit$age_grid, y = m2a.quants[,1], color = 'M2a'), size = 1, 
  #           span = 0.5, se = F, stat = 'smooth') +
  # geom_line(aes(x = vrm.m2a$fit$age_grid, y = m2a.quants[,2], color = 'M2a'), size = 1, 
  #           span = 0.5, se = F, stat = 'smooth') +
  geom_line(aes(x = vrm.m2b$fit$age_grid[-1], y = vrm.m2b$fit$fitted_PDR[-1]),
            color = 'black', size = 2.5, span = 0.5, se = F, stat = 'smooth') +
  geom_line(aes(x = vrm.m2b$fit$age_grid[-1], y = vrm.m2b$fit$fitted_PDR[-1], color = 'M2b'), size = 2, 
            span = 0.5, se = F, stat = 'smooth') +
  # geom_line(aes(x = vrm.m2b$fit$age_grid, y = m2b.quants[,1]), size = 1.5, color = 'black',
  #           span = 0.5, se = F, stat = 'smooth') +
  # geom_line(aes(x = vrm.m2b$fit$age_grid, y = m2b.quants[,2]), size = 1.5, color = 'black',
  #           span = 0.5, se = F, stat = 'smooth') +
  # geom_line(aes(x = vrm.m2b$fit$age_grid, y = m2b.quants[,1], color = 'M2b'), size = 1, 
  #           span = 0.5, se = F, stat = 'smooth') +
  # geom_line(aes(x = vrm.m2b$fit$age_grid, y = m2b.quants[,2], color = 'M2b'), size = 1, 
  #           span = 0.5, se = F, stat = 'smooth') +
  coord_cartesian(ylim = c(-0.5, 1.5),
                  xlim = c(50, 0)) +
  labs(x = 'Mybp', y = 'Pulled DR', color = 'Group') +
  scale_color_manual(values = main.cols) +
  theme_bw() +
  theme(legend.position = c(0.1, 0.85))  + 
  scale_x_reverse()


pdf('Smoothed-Truncated-PDR-Anolis.pdf', width = 12, height = 6)
smooth | main.smooth
dev.off()

             