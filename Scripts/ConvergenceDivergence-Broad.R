##################################################
# Testing for convergence a la Mazel et al. 2017 #
##################################################
######### Fit models of trait evolution ##########
library(treeplyr)
library(purrr)
library(geiger)
library(mvMORPH)
library(plyr)
library(RRphylo)
library(adephylo)
library(spaa)
library(otuSummary)

setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/Convergence/')
source('../ConvergenceDivergenceFunctions.R')

m.data <- read.csv('../FinalData/Male_Transf-Scaled-AllTraitsEcol-FINAL.csv', 
                   header = T)
# f.data <- read.csv('../FinalData/Female_Transf-Scaled-AllTraitsEcol-FINAL.csv', 
#                    header = T)

# Remove body height as this leads to a substantial reduction in 
# the number of species included. For these convergence analyses,
# we want to ensure we've got as complete a dataset as possible to
# make sure that we've got as many of the 'classic' island ecomorphs
# as possible. We will also remove agasizzi, since this species is a 
# bit obscure for a number of reasons (trees)
m.data <- m.data[,c(1,3:10,12:13,23:24)]
m.data <- m.data[-which(m.data$Species %in% c('agassizi', 'bicaorum', 'lineatus', 'roatanensis', 'medemi')),]

# Now remove missing data
m.data <- na.omit(m.data)
# f.data <- na.omit(f.data)

# Read in the time-calibrated phylogeny
tree <- read.tree('../FinalData/Anole_TimeCal_NoOut.tree')

# And correspond species' trait data to their position in the phylogeny.
m.full.td <- make.treedata(tree, m.data)
m.full.data <- as.data.frame(m.full.td$dat)
rownames(m.full.data) <- m.full.td$phy$tip.label

bm.fit.m <- mvBM(m.full.td$phy, m.full.data[,-c(11:12)], model = 'BM1', method = 'rpf', diagnostic = T)
saveRDS(object = bm.fit.m, file = 'mvBM1-Males_AllTraits.rds')

ou.fit.m <- mvOU(m.full, m.full.data[,-c(11:12
                                         )], model = 'OU1', method = 'rpf', diagnostic = T)
saveRDS(object = ou.fit.m, file = 'mvOU1-Males_AllTraits.rds')


# Read in basic model fitting results. We will simulate under these models
setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/Convergence/')
mv.bm.m <- readRDS('mvBM1-Males_AllTraits.rds')
mv.ou.m <- readRDS('mvOU1-Males_AllTraits.rds')

model.comp <- 
  data.frame('Sex' = rep('Males', 2), 
             'Model' = c('Multivariate BM', 'Multivariate OU'),
             'AICc' = c(mv.bm.m$AICc, mv.ou.m$AICc),
             'Delta AICc' = c(0, abs(mv.bm.m$AICc - mv.ou.m$AICc)))

write.table(model.comp, file = 'TraitEvolution-ModelComp.txt', sep = '\t')

######### Simulate Trait Data and Calculte PSTD ##########
m.full <- m.full.td$phy

bm.sim.m <- mvSIM(m.full, nsim = 20000, model = 'BM1', param = mv.bm.m)

saveRDS(object = bm.sim.m, file = 'mvBM-Sim-Males.rds')

# Now calculate pairwise distances for each simulation.
library(StatMatch)

bm.sim.m <- readRDS(file = 'mvBM-Sim-Males.rds')

m.traitDists <- GetTraitDists(sims = bm.sim.m, obs = m.full.data)

m.pstd <- GetPSTD(m.traitDists) 

write.csv(m.pstd$mean, 'Male_mvBM-Final_PSTD_Mean_Values.csv')
write.csv(m.pstd$median, 'Male_mvBM-Final_PSTD_Median_Values.csv')

# saveRDS(m.pstd, 'Male_mvBM_PSTD.rds')
# 
######### Identify Convergent and Divergent Taxa - Island/Mainland #########
setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/Convergence/')

# Get species-region correspondance
m1 <- read.tree('../FinalData/Mainland_Clade1_TimeTree.tree')
m2 <- read.tree('../FinalData/Mainland_Clade2_TimeTree.tree')
i1 <- read.tree('../FinalData/Island_Clade1_TimeTree.tree')
i2 <- read.tree('../FinalData/Island_Clade2_TimeTree.tree')

spp <- c(i1$tip.label, i2$tip.label, m1$tip.label, m2$tip.label)
regions <- c(rep('I1', length(i1$tip.label)), rep('I2', length(i2$tip.label)),
             rep('M1', length(m1$tip.label)), rep('M2', length(m2$tip.label)))
spp.regions <- data.frame('Species' = spp, 'Region' = regions)

spp.regions.m <- spp.regions[which(spp.regions$Species %in% m.full$tip.label),]

m.pstd <- readRDS('Male_mvBM-Final_PSTD.rds')

# Now prepare to identify actual convergent/divergent
# species pairs. We will correspond PSTD values and species
# pairs to regions

m.pstd.results <- PreparePSTD(pstd = m.pstd, spp.cats = spp.regions.m)

m.percentile <- ecdf(x = m.pstd.results$mean$PSTD)

m.quantiles <- c()

for(i in 1:length(m.pstd.results$mean$PSTD)){
  x <- m.pstd.results$mean$PSTD[i]
  m.quantiles[i] <- m.percentile(x)
}

m.pstd.results$mean$Quantile <- m.quantiles

m.bott5.pstd <- m.pstd.results$mean[which(m.pstd.results$mean$Quantile <= 0.05),]
m.top5.pstd <- m.pstd.results$mean[which(m.pstd.results$mean$Quantile >= 0.95),]

write.table(m.bott5.pstd, 'Male_Bottom5percent_PSTD.txt', sep = '\t', row.names = F, quote = F)
write.table(m.top5.pstd, 'Male_Top5percent_PSTD.txt', sep = '\t', row.names = F, quote = F)

# And now get your convergent and divergent taxa
m.bott <- GetConvergent(pstd.dists = m.pstd.results, 
                        threshold = 0.01, 
                        tree = tree, minPatrDist = 20,
                        minNodeDist = 2)
m.top <- GetDivergent(pstd.dists = m.pstd.results, 
                      threshold = 0.01, 
                      tree = tree)

#Test for abundance of convergences among allopatric vs sympatric species.
setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/Convergence/')

m.pstd <- readRDS('Male_mvBM-Final_PSTD.rds')
m.pstd.results <- PreparePSTD(pstd = m.pstd, spp.cats = spp.regions.m)

m.pstd.all <- GetConvergent(pstd.dists = m.pstd.results, 
                            threshold = 1, 
                            tree = tree, minPatrDist = 0,
                            minNodeDist = 0)


m.m1.pstd <- 
  m.pstd.all$mean[which(m.pstd.all$mean$Spp1.Category == 'M1' &
                            m.pstd.all$mean$Spp2.Category == 'M1' ),]
m.m2.pstd <- 
  m.pstd.all$mean[which(m.pstd.all$mean$Spp2.Category == 'M2' &
                            m.pstd.all$mean$Spp2.Category == 'M2' ),]

m.i2.pstd <- 
  m.pstd.all$mean[which(m.pstd.all$mean$Spp2.Category == 'I2' &
                            m.pstd.all$mean$Spp2.Category == 'I2' ),]

m.m1.pstd$Within <- rep('M1', nrow(m.m1.pstd))
m.m2.pstd$Within <- rep('M2', nrow(m.m2.pstd))
m.i2.pstd$Within <- rep('I2', nrow(m.i2.pstd))

all.pstd <- rbind(m.i2.pstd, m.m1.pstd, m.m2.pstd)
mainland.pstd <- rbind(m.m1.pstd, m.m2.pstd)
mainland.conv.pstd <- mainland.pstd[which(mainland.pstd$PSTD <= 0),]

mainland.conv.pstd$Spp1.NumM1.Sympat <- NA
mainland.conv.pstd$Spp1.NumM2a.Sympat <- NA
mainland.conv.pstd$Spp1.NumM2b.Sympat <- NA
mainland.conv.pstd$Spp2.NumM1.Sympat <- NA
mainland.conv.pstd$Spp2.NumM2a.Sympat <- NA
mainland.conv.pstd$Spp2.NumM2b.Sympat <- NA

data <- read.csv('../FinalData/Mainland-Species-Sympatry-Summary.csv', header = T, row.names = 1)
rownames(data) <- data$Species
data <- data[,-1]

colnames(data) <- c('Lineage',  "Sympatric W/ M1?", "Total Number Symp Spp",
                    "NumM1.Sympat", "NumM2a.Sympat", "NumM2b.Sympat",
                    "Prop Sympatric M1", "Prop Sympatric M2a", 
                    "Prop Sympatric M2b")


data$Species <- rownames(data)
mainland.conv.pstd$Spp1 <- as.character(mainland.conv.pstd$Spp1)
mainland.conv.pstd$Spp2 <- as.character(mainland.conv.pstd$Spp2)

for(i in 1:nrow(mainland.conv.pstd)){
  mainland.conv.pstd$Spp1.NumM1.Sympat[i] <- 
    data[which(data$Species == mainland.conv.pstd$Spp1[i]),]$NumM1.Sympat
  
  mainland.conv.pstd$Spp1.NumM2a.Sympat[i] <- 
    data[which(data$Species == mainland.conv.pstd$Spp1[i]),]$NumM2a.Sympat
  
  mainland.conv.pstd$Spp1.NumM2b.Sympat[i] <- 
    data[which(data$Species == mainland.conv.pstd$Spp1[i]),]$NumM2b.Sympat
  
  ##
  
  mainland.conv.pstd$Spp2.NumM1.Sympat[i] <- 
    data[which(data$Species == mainland.conv.pstd$Spp2[i]),]$NumM1.Sympat
  
  mainland.conv.pstd$Spp2.NumM2a.Sympat[i] <- 
    data[which(data$Species == mainland.conv.pstd$Spp2[i]),]$NumM2a.Sympat
  
  mainland.conv.pstd$Spp2.NumM2b.Sympat[i] <- 
    data[which(data$Species == mainland.conv.pstd$Spp2[i]),]$NumM2b.Sympat
}

for(i in 1:nrow(mainland.conv.pstd)){
  mainland.conv.pstd$MeanNumM1Sympat[i] <- mean(mainland.conv.pstd$Spp1.NumM1.Sympat[i], 
                                                mainland.conv.pstd$Spp2.NumM1.Sympat[i])
  mainland.conv.pstd$MeanNumM2aSympat[i] <- mean(mainland.conv.pstd$Spp1.NumM2a.Sympat[i], 
                                                 mainland.conv.pstd$Spp2.NumM2a.Sympat[i])
  mainland.conv.pstd$MeanNumM2bSympat[i] <- mean(mainland.conv.pstd$Spp1.NumM2b.Sympat[i], 
                                                 mainland.conv.pstd$Spp2.NumM2b.Sympat[i])
  mainland.conv.pstd$MeanNumM2Sympat[i] <- (mainland.conv.pstd$MeanNumM2aSympat[i] + 
                                              mainland.conv.pstd$MeanNumM2bSympat[i])
}

# Determine if spp1 is sympatric with spp2
reg.occup <- read.csv('../FinalData/PoeDat_AP.csv', header = T) %>% .[,-c(2:6,14:21)]
mland.spp <- spp.regions.m[which(spp.regions.m$Region %in% c('M1', 'M2')),]
mland.patry <- merge(mland.spp, reg.occup, by = 'Species')
rownames(mland.patry) <- mland.patry$Species
mland.patry <- mland.patry[,-1]

nearctic <- rownames(mland.patry[which(mland.patry$Nearctic == 1),])
uca <- rownames(mland.patry[which(mland.patry$UpperCentralAmerica == 1),])
lca <- rownames(mland.patry[which(mland.patry$LowerCentralAmerica == 1),])
choco <- rownames(mland.patry[which(mland.patry$Choco == 1),])
carib <- rownames(mland.patry[which(mland.patry$CaribbeanCoast == 1),])
andes <- rownames(mland.patry[which(mland.patry$Andes == 1),])
amaz <- rownames(mland.patry[which(mland.patry$Amazonia == 1),])

mainland.conv.pstd$Patry <- NA
for(i in 1:nrow(mainland.conv.pstd)){
  # go down the list of PSTD and test if the two species share a region. 
  spp1 <- mainland.conv.pstd$Spp1[i]
  spp2 <- mainland.conv.pstd$Spp2[i]
  
  symp <- 
    sum(spp1 %in% nearctic & spp2 %in% nearctic, 
        spp1 %in% uca & spp2 %in% uca,
        spp1 %in% lca & spp2 %in% lca,
        spp1 %in% choco & spp2 %in% choco,
        spp1 %in% carib & spp2 %in% carib,
        spp1 %in% andes & spp2 %in% andes,
        spp1 %in% amaz & spp2 %in% amaz)
  
  mainland.conv.pstd$Patry[i] <- 
    if_else(symp == 0, 'Allopatric', 'Sympatric')
}

# Look for all - island/mainland

spp.dists <- read.csv('../FinalData/PoeDat_AP.csv')
spp.dists <- spp.dists[which(spp.dists$Species %in% c(i2$tip.label, m2$tip.label)),]

# Reduce down to only those species involved in convergences
# spp.dists <- spp.dists[which(spp.dists$Species %in% 
#                                c(as.character(m.conv$Spp1), as.character(m.conv$Spp2))),]


amazo <- as.character(spp.dists[which(spp.dists$Amazonia == 1),1])
andes <- as.character(spp.dists[which(spp.dists$Andes == 1),1])
baham <- as.character(spp.dists[which(spp.dists$Bahamas == 1),1])
caribCo <- as.character(spp.dists[which(spp.dists$CaribbeanCoast == 1),1])
cuba <- as.character(spp.dists[which(spp.dists$CubaCayman == 1),1])
hisp <- as.character(spp.dists[which(spp.dists$Hispaniola == 1),1])
jama <- as.character(spp.dists[which(spp.dists$Jamaica == 1),1])
lca <- as.character(spp.dists[which(spp.dists$LowerCentralAmerica == 1),1])
near <- as.character(spp.dists[which(spp.dists$Nearctic == 1),1])
pr <- as.character(spp.dists[which(spp.dists$PuertoRico == 1),1])
smIsl <- as.character(spp.dists[which(spp.dists$SmallIslands == 1),1])
uca <- as.character(spp.dists[which(spp.dists$UpperCentralAmerica == 1),1])
GrAntil <- c(baham, cuba, hisp, jama, pr)


m.conv <- read.csv('Male_mvBM_Mean_Convergences.csv')
m.conv$Patry <- NA
for(i in 1:nrow(m.conv)){
  # go down the list of PSTD and test if the two species share a region. 
  spp1 <- m.conv$Spp1[i]
  spp2 <- m.conv$Spp2[i]
  
  symp <- 
    sum(spp1 %in% nearctic & spp2 %in% nearctic, 
        spp1 %in% uca & spp2 %in% uca,
        spp1 %in% lca & spp2 %in% lca,
        spp1 %in% choco & spp2 %in% choco,
        spp1 %in% carib & spp2 %in% carib,
        spp1 %in% andes & spp2 %in% andes,
        spp1 %in% amaz & spp2 %in% amaz)
  
  
  m.conv$Patry[i] <- 
    ifelse(symp == 0, 'Allopatric', 'Sympatric')
}

m.pstd.results$mean$Patry <- NA
for(i in 1:nrow(m.pstd.results$mean)){
  # go down the list of PSTD and test if the two species share a region. 
  spp1 <- m.pstd.results$mean$Spp1[i]
  spp2 <- m.pstd.results$mean$Spp2[i]
  
  symp <- 
    sum(spp1 %in% nearctic & spp2 %in% nearctic, 
        spp1 %in% uca & spp2 %in% uca,
        spp1 %in% lca & spp2 %in% lca,
        spp1 %in% choco & spp2 %in% choco,
        spp1 %in% caribCo & spp2 %in% caribCo,
        spp1 %in% andes & spp2 %in% andes,
        spp1 %in% amaz & spp2 %in% amaz,
        spp1 %in% baham & spp2 %in% baham,
        spp1 %in% cuba & spp2 %in% cuba,
        spp1 %in% hisp & spp2 %in% hisp,
        spp1 %in% jama & spp2 %in% jama,
        spp1 %in% pr & spp2 %in% pr,
        spp1 %in% smIsl & spp2 %in% smIsl)
  
  m.pstd.results$mean$Patry[i] <- 
    ifelse(symp == 0, 'Allopatric', 'Sympatric')
}

i2.pstd.results <- m.pstd.results$mean[which(m.pstd.results$mean$Spp1.Category == 'I2' & 
                                           m.pstd.results$mean$Spp2.Category == 'I2'),]
m1.pstd.results <- m.pstd.results$mean[which(m.pstd.results$mean$Spp1.Category == 'M1' & 
                                           m.pstd.results$mean$Spp2.Category == 'M1'),]
m2.pstd.results <- m.pstd.results$mean[which(m.pstd.results$mean$Spp1.Category == 'M2' & 
                                           m.pstd.results$mean$Spp2.Category == 'M2'),]

expected.allo <- c(sum(i2.pstd.results$Patry == 'Allopatric'), 
                   sum(m2.pstd.results$Patry == 'Allopatric'))
expected.prob.allo <- c(sum(i2.pstd.results$Patry == 'Allopatric') / nrow(i2.pstd.results), 
                        sum(m2.pstd.results$Patry == 'Allopatric') / nrow(m2.pstd.results))

expected.symp <- c(sum(i2.pstd.results$Patry == 'Sympatric'), 
                   sum(m2.pstd.results$Patry == 'Sympatric'))
expected.prob.symp <- c(sum(i2.pstd.results$Patry == 'Sympatric') / nrow(i2.pstd.results), 
                        sum(m2.pstd.results$Patry == 'Sympatric') / nrow(m2.pstd.results))

expected.prop.patry <- expected.allo / expected.symp
expected.probs.patry <- expected.prob.allo / expected.prob.symp

  
#m.conv <- read.csv('Male_mvBM_Mean_Convergences.csv')

obs.allo <- c(sum(m.conv[which(m.conv$Spp1.Category == 'I2' &
                                 m.conv$Spp2.Category == 'I2'),]$Patry == 'Allopatric'), 
              sum(m.conv[which(m.conv$Spp1.Category == 'M2' &
                                 m.conv$Spp2.Category == 'M2'),]$Patry == 'Allopatric'))

obs.symp <- c(sum(m.conv[which(m.conv$Spp1.Category == 'I2' &
                                 m.conv$Spp2.Category == 'I2'),]$Patry == 'Sympatric'), 
              sum(m.conv[which(m.conv$Spp1.Category == 'M2' &
                                 m.conv$Spp2.Category == 'M2'),]$Patry == 'Sympatric'))

exp.mat <- matrix(c(expected.symp*0.01, expected.allo*0.01), ncol = 2, 
                  dimnames = list(c('I2', 'M2'), c('Sympatric', 'Allopatric')))
exp.prob <- matrix(c(expected.prob.symp, expected.prob.allo), ncol = 2, 
                   dimnames = list(c('I2', 'M2'), c('Sympatric', 'Allopatric')))

obs.mat <- matrix(c(obs.symp, obs.allo), ncol = 2, 
                  dimnames = list(c('I2', 'M2'), c('Sympatric', 'Allopatric')))



capture.output(chisq.test(obs.mat[1,], p = exp.prob[1,], simulate.p.value = T), 
               file = 'NumConvergences_AllopatrySympatry_I2-ChiSq.csv')
capture.output(chisq.test(obs.mat[2,], p = exp.prob[2,], simulate.p.value = T), 
               file = 'NumConvergences_AllopatrySympatry_M2-ChiSq.csv')

capture.output(chisq.test(obs.mat, p = exp.prob, simulate.p.value = T), 
               file = 'NumConvergences_AllopatrySympatry_I2-M2-ChiSq.csv')

capture.output(fisher.test(obs.mat), 
               file = 'NumConvergences_WithinAmongRegions-FishersExact.csv')

obs.symp <- c(sum(i2.pstd.results$Patry == 'Sympatric'), 
              sum(m1.pstd.results$Patry == 'Sympatric'), 
              sum(m2.pstd.results$Patry == 'Sympatric'))
mainland.top.conv <- mainland.conv.pstd[which(mainland.conv.pstd$PSTD < -1.58),]


cols <- c("#d72027", "#add8a4", "#2b83bb")

pstd.patry <- 
  ggplot(data = mainland.conv.pstd, 
         aes(x = Patry, y = PSTD)) +
  geom_boxplot(aes(fill = Within), outlier.colour = NA) +
  geom_point(aes(fill = Within, shape = Patry), 
             position = position_jitterdodge(jitter.width = 0.1),
             pch = 21, size = 3, alpha = 0.5) +
  scale_color_manual(values = cols[-1]) +
  scale_fill_manual(values = cols[-1]) + 
  xlab('Distribution of Species Pair') +
  theme_bw()

dist.patry <- 
  ggplot(data = mainland.conv.pstd, 
         aes(x = Patry, y = PatristicDist)) +
  geom_boxplot(aes(fill = Within), outlier.colour = NA) +
  geom_point(aes(fill = Within, shape = Patry), 
             position = position_jitterdodge(jitter.width = 0.1),
             pch = 21, size = 3, alpha = 0.5) +
  scale_color_manual(values = cols[-1]) +
  scale_fill_manual(values = cols[-1]) + 
  xlab('Distribution of Species Pair') +
  theme_bw()

pstd.patrist.ratio.patry <- 
  ggplot(data = mainland.conv.pstd, 
         aes(x = Patry, y = PSTD/PatristicDist)) +
  geom_boxplot(aes(fill = Within), outlier.colour = NA) +
  geom_point(aes(fill = Within, shape = Patry), 
             position = position_jitterdodge(jitter.width = 0.1),
             pch = 21, size = 3, alpha = 0.5) +
  scale_color_manual(values = cols[-1]) +
  scale_fill_manual(values = cols[-1]) + 
  xlab('Distribution of Species Pair') +
  theme_bw()


sig.conv.patry.pstd <- read.csv('Male_mvBM-Final_Mean_Convergences.csv')

sig.conv.patry.pstd$Patry <- NA
for(i in 1:nrow(sig.conv.patry.pstd)){
  # go down the list of PSTD and test if the two species share a region. 
  spp1 <- sig.conv.patry.pstd$Spp1[i]
  spp2 <- sig.conv.patry.pstd$Spp2[i]
  
  symp <- 
    sum(spp1 %in% nearctic & spp2 %in% nearctic, 
        spp1 %in% uca & spp2 %in% uca,
        spp1 %in% lca & spp2 %in% lca,
        spp1 %in% choco & spp2 %in% choco,
        spp1 %in% caribCo & spp2 %in% caribCo,
        spp1 %in% andes & spp2 %in% andes,
        spp1 %in% amaz & spp2 %in% amaz,
        spp1 %in% baham & spp2 %in% baham,
        spp1 %in% cuba & spp2 %in% cuba,
        spp1 %in% hisp & spp2 %in% hisp,
        spp1 %in% jama & spp2 %in% jama,
        spp1 %in% pr & spp2 %in% pr,
        spp1 %in% smIsl & spp2 %in% smIsl)
  
  sig.conv.patry.pstd$Patry[i] <- 
    ifelse(symp == 0, 'Allopatric', 'Sympatric')
}

sig.conv.patry.pstd$Within <- sig.conv.patry.pstd$Spp1.Category

cols <- c("#d72027", "#add8a4", "#2b83bb")

pstd.patry <- 
  ggplot(data = sig.conv.patry.pstd, 
         aes(x = Patry, y = PSTD)) +
  geom_boxplot(aes(fill = Within), outlier.colour = NA) +
  geom_point(aes(fill = Within, shape = Patry), 
             position = position_jitterdodge(jitter.width = 0.1),
             pch = 21, size = 3, alpha = 0.5) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) + 
  xlab('Distribution of Species Pair') +
  theme_bw()

dist.patry <- 
  ggplot(data = sig.conv.patry.pstd, 
         aes(x = Patry, y = PatristicDist)) +
  geom_boxplot(aes(fill = Within), outlier.colour = NA) +
  geom_point(aes(fill = Within, shape = Patry), 
             position = position_jitterdodge(jitter.width = 0.1),
             pch = 21, size = 3, alpha = 0.5) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) + 
  xlab('Distribution of Species Pair') +
  theme_bw()

pstd.patrist.ratio.patry <- 
  ggplot(data = sig.conv.patry.pstd, 
         aes(x = Patry, y = PSTD/PatristicDist)) +
  geom_boxplot(aes(fill = Within), outlier.colour = NA) +
  geom_point(aes(fill = Within, shape = Patry), 
             position = position_jitterdodge(jitter.width = 0.1),
             pch = 21, size = 3, alpha = 0.5) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) + 
  xlab('Distribution of Species Pair') +
  theme_bw()



pdf('PSTD-PatristicDist-Mainland-Sympat-Allopat.pdf', width = 16, height = 8)
pstd.patry + dist.patry
dev.off()

ggplot(data = mainland.conv.pstd, 
       aes(x = PatristicDist, y = PSTD, 
           group = Within)) +
  geom_point(aes(fill = Within, shape = Patry), size = 3) +
  stat_smooth(aes(group = Within),
              method = 'lm', formula = y ~ poly(x, 2),
              color = "black", size = 3, se = F) +
  stat_smooth(aes(color = Within),
              method = 'lm', formula = y ~ poly(x, 2), 
              size = 2, se = F) +
  stat_smooth(aes(group = Patry, linetype = Patry),
              method = 'lm', formula = y ~ poly(x, 2),
              color = "black", size = 3, se = F) +
  stat_smooth(aes(color = Patry, linetype = Patry),
              method = 'lm', formula = y ~ poly(x, 2), 
              size = 2, se = F) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) + 
  scale_shape_manual(values = c(21, 22)) +
  theme_bw()


pdf('PSTD-Patristic-M1-M2-I2.pdf', width = 8, height = 7)
ggplot(data = all.pstd[which(all.pstd$PSTD <= 0),], 
       aes(x = PatristicDist, y = PSTD, group = Within)) +
  geom_point(aes(fill = Within), pch = 21, size = 3) +
  stat_smooth(aes(group = Within),
              method = 'lm', formula = y ~ poly(x, 2),
              color = "black", size = 3, se = F) +
  stat_smooth(aes(color = Within),
              method = 'lm', formula = y ~ poly(x, 2), 
              size = 2, se = F) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) + 
  theme_bw()
dev.off()  
  

cols <- c("#d72027", "#add8a4", "#2b83bb")
ggplot() +
  # geom_point(data = m.i2.pstd, aes(x = PatristicDist, y = PSTD), 
  #            fill = "#d72027", pch = 21, size = 3) +
  # stat_smooth(data = m.i2.pstd, aes(x = PatristicDist, y = PSTD), 
  #             method = 'lm', formula = y ~ poly(x, 2), 
  #             color = "black", size = 2) +
  # stat_smooth(data = m.i2.pstd, aes(x = PatristicDist, y = PSTD), 
  #             method = 'lm', formula = y ~ poly(x, 2), 
  #             color = "#d72027", size = 1) +
  geom_point(data = m.m2.pstd, aes(x = PatristicDist, y = PSTD), 
             fill = "#2b83bb", pch = 21, size = 3) +
  stat_smooth(data = m.m2.pstd, aes(x = PatristicDist, y = PSTD), 
              method = 'lm', formula = y ~ poly(x, 2), 
              color = "black", size = 2) +
  stat_smooth(data = m.m2.pstd, aes(x = PatristicDist, y = PSTD), 
              method = 'lm', formula = y ~ poly(x, 2), 
              color = "#2b83bb", size = 1) +
  geom_point(data = m.m1.pstd, aes(x = PatristicDist, y = PSTD), 
             fill = "#add8a4", pch = 21, size = 3) +
  stat_smooth(data = m.m1.pstd, aes(x = PatristicDist, y = PSTD), 
              method = 'lm', formula = y ~ poly(x, 2), 
              color = "black", size = 2) +
  stat_smooth(data = m.m1.pstd, aes(x = PatristicDist, y = PSTD), 
              method = 'lm', formula = y ~ poly(x, 2), 
              color = "#add8a4", size = 1)



ggplot() +
  geom_point(data = m.i2.pstd, aes(x = PatristicDist, y = PSTD), 
             fill = "#d72027", pch = 21, size = 3) +
  stat_smooth(data = m.i2.pstd, aes(x = PatristicDist, y = PSTD), 
              method = 'lm', formula = y ~ poly(x, 2), 
              color = "black", size = 2) +
  stat_smooth(data = m.i2.pstd, aes(x = PatristicDist, y = PSTD), 
              method = 'lm', formula = y ~ poly(x, 2), 
              color = "#d72027", size = 1) +
  geom_point(data = m.m2.pstd, aes(x = PatristicDist, y = PSTD), 
             fill = "#2b83bb", pch = 21, size = 3) +
  stat_smooth(data = m.m2.pstd, aes(x = PatristicDist, y = PSTD), 
              method = 'lm', formula = y ~ poly(x, 2), 
              color = "black", size = 2) +
  stat_smooth(data = m.m2.pstd, aes(x = PatristicDist, y = PSTD), 
              method = 'lm', formula = y ~ poly(x, 2), 
              color = "#2b83bb", size = 1) +
  geom_point(data = m.m1.pstd, aes(x = PatristicDist, y = PSTD), 
             fill = "#add8a4", pch = 21, size = 3) +
  stat_smooth(data = m.m1.pstd, aes(x = PatristicDist, y = PSTD), 
              method = 'lm', formula = y ~ poly(x, 2), 
              color = "black", size = 2) +
  stat_smooth(data = m.m1.pstd, aes(x = PatristicDist, y = PSTD), 
              method = 'lm', formula = y ~ poly(x, 2), 
              color = "#add8a4", size = 1)



   
# Now plot the distances among convergent taxa to get a sense of how 
# robust results are to stasis. 

pdf('ConvergenceDistances_MeanPSTD_BestMods.pdf', width = 10, height = 6)
par(mfrow = c(1,2))
hist(m.bott$mean$PatristicDist, breaks = 10, xlab = 'Patristic Distance Among Male Convergence', main = '')
hist(m.bott$mean$NodeDist, breaks = 20, xlim=c(0,35), xlab = '# Intervening Nodes Among Male Convergence', main = '')
dev.off()

write.csv(m.top$mean, 'Male_mvBM_Mean_Divergences.csv', row.names = F)
write.csv(m.bott$mean, 'Male_mvBM_Mean_Convergences.csv', row.names = F)

# Calculate proportion of prevalences in extremes
library(dplyr)

m.summs <- SummarizeConvDiv(conv = m.bott, div = m.top, 
                            cats = c('I1', 'M1', 'I2', 'M2'),
                            dist.list = m.pstd.results,
                            threshold = 0.01, phy = m.full)

capture.output(m.summs, file = 'ConvDiv_mvBM_ByGroup_Scaled_Males.txt')

library(ggplot2)
library(plotrix)

pdf('ConvDiv_mvBM-Mean_ByGroup_Scaled_Males.pdf', width=7, height=6)
ggplot(data = m.summs$mean, aes(x = Category, y = ObsOverExpPrev)) + 
  geom_bar(stat = "identity", aes(fill = Type), 
           position=position_dodge()) + 
  scale_fill_manual(values = c('#67a9cf', '#ef8a62')) +
  geom_hline(yintercept = 1, lty = 2) +
  theme_bw()
dev.off()

################ Now plot as a heatmap #####################

library(igraph)
library(scales)
library(plyr)
library(tidyr)
library(stringr)

# First we need to define the possible sets of pairs. Because the order can flip, 
# we need to account for this. Likewise, we need to define factor 
# levels for the pairs

pair.types <- list(c('I1-I2', 'I2-I1'), c('I1-M1', 'M1-I1'), 
                   c('I1-M2', 'M2-I1'), c('I2-M2', 'M2-I2'), 
                   c('M1-I2', 'I2-M1'), c('M1-M2', 'M2-M1'))
pair.levels <- c('I1-I1', 'I2-I1', 'M1-I1', 'M2-I1', 'I2-I2',
                 'M2-I2', 'I2-M1', 'M1-M1', 'M2-M1', 'M2-M2')


m.heat <-
  PrepareHeatmap(Convergent = m.bott, Divergent = m.top,
                 pstdResults = m.pstd.results, 
                 pair.types = pair.types,
                 pair.levels = pair.levels,
                 threshold = 0.01)

# Prop Test for M2 and I2
#m.heat <- readRDS(file = 'PairwiseConvDiv-Males-Broad.rds')
m2.succ <- m.heat$Observed$ConvMed[10]
i2.succ <- m.heat$Observed$ConvMed[5]
m2.fail <- sum(m.heat$Observed$ConvMed)-m2.succ
i2.fail <- sum(m.heat$Observed$ConvMed)-i2.succ

prop.test(matrix(c(i2.succ, m2.succ, i2.fail, m2.fail), ncol = 2), alternative = 'greater')

# Rearrange slightly for plotting
m.heat$ObsExp[4,c(15,16)] <- c('M1', 'I2')

saveRDS(m.heat, file = 'PairwiseConvDiv-Males-Broad.rds')
capture.output(m.heat, file = 'PairwiseConvDiv-Males-Broad.txt')

pdf('Male_Convergences-Mean_Heatmap.pdf', width = 7, height = 7)
ggplot(data = m.heat$ObsExp, aes(x = Region1, y = Region2, 
                          fill = log(ConvMean.ObsExpFreqs))) +
  geom_tile() +
  scale_fill_distiller(type = 'div', palette = 7, limits = c(-2.5, 1.5)) +
  guides(fill=guide_legend(title="log(Obs/Exp Convergences)")) +
  theme_bw()  +
  theme(legend.position = c(0.025, 0.975), 
        legend.justification = c(0, 1)) 
dev.off()

########
# Now, relaxing to 2.5%
m.bott <- GetConvergent(pstd.dists = m.pstd.results, 
                        threshold = 0.025, 
                        tree = tree, minPatrDist = 20,
                        minNodeDist = 2)
m.heat <-
  PrepareHeatmap(Convergent = m.bott, Divergent = m.top,
                 pstdResults = m.pstd.results, 
                 pair.types = pair.types,
                 pair.levels = pair.levels,
                 threshold = 0.025)

# Rearrange slightly for plotting
m.heat$ObsExp[4,c(15,16)] <- c('M1', 'I2')

saveRDS(m.heat, file = 'PairwiseConvDiv-0.025-Males-Broad.rds')
capture.output(m.heat, file = 'PairwiseConvDiv-0.025-Males-Broad.txt')

pdf('Male_Convergences-Mean_0.025-Heatmap.pdf', width = 7, height = 7)
ggplot(data = m.heat$ObsExp, aes(x = Region1, y = Region2, 
                                 fill = log(ConvMean.ObsExpFreqs))) +
  geom_tile() +
  scale_fill_distiller(type = 'div', palette = 7, limits = c(-2, 2)) +
  guides(fill=guide_legend(title="log(Obs/Exp Convergences)")) +
  theme_bw()  +
  theme(legend.position = c(0.025, 0.975), 
        legend.justification = c(0, 1)) 
dev.off()

########
# Now, relaxing to 2.5%
m.bott <- GetConvergent(pstd.dists = m.pstd.results, 
                        threshold = 0.05, 
                        tree = tree, minPatrDist = 20,
                        minNodeDist = 2)
m.heat <-
  PrepareHeatmap(Convergent = m.bott, Divergent = m.top,
                 pstdResults = m.pstd.results, 
                 pair.types = pair.types,
                 pair.levels = pair.levels,
                 threshold = 0.05)

# Rearrange slightly for plotting
m.heat$ObsExp[4,c(15,16)] <- c('M1', 'I2')

saveRDS(m.heat, file = 'PairwiseConvDiv-0.05-Males-Broad.rds')
capture.output(m.heat, file = 'PairwiseConvDiv-0.05-Males-Broad.txt')

pdf('Male_Convergences-Mean_0.05-Heatmap.pdf', width = 7, height = 7)
ggplot(data = m.heat$ObsExp, aes(x = Region1, y = Region2, 
                                 fill = log(ConvMean.ObsExpFreqs))) +
  geom_tile() +
  scale_fill_distiller(type = 'div', palette = 7, limits = c(-2, 2)) +
  guides(fill=guide_legend(title="log(Obs/Exp Convergences)")) +
  theme_bw()  +
  theme(legend.position = c(0.025, 0.975), 
        legend.justification = c(0, 1)) 
dev.off()

######### Co-Phylo showing convergences ######### 
m.bott <- GetConvergent(pstd.dists = m.pstd.results, 
                        threshold = 0.01, 
                        tree = tree, minPatrDist = 20,
                        minNodeDist = 2)

m.bott$mean[,1] <- as.character(m.bott$mean[,1])
m.bott$mean[,2] <- as.character(m.bott$mean[,2])
for(i in 1:nrow(m.bott$mean)){
  if(i %% 2 == 0){
    a <- as.character(m.bott$mean[i,1])
    b <- as.character(m.bott$mean[i,2])
    m.bott$mean[i,1] <- b
    m.bott$mean[i,2] <- a 
  } else {
    m.bott$mean[i,1] <- as.character(m.bott$mean[i,1])
    m.bott$mean[i,2] <- as.character(m.bott$mean[i,2])
  }
}

conv.trees <- cophylo(m.full.td$phy, m.full.td$phy, assoc = m.bott$mean[,1:2], 
                      rotate = F)
pdf('Convergence-mvBM-0.1_CoPhylo-Males.pdf', width=13, height=7)
plot.cophylo(conv.trees, fsize = 0.4, link.lwd = 1, link.col='darkgrey', node.size = 0)
dev.off()



# Now define edge widths
m.top$mean$pair <- as.factor(paste(m.top$mean$Spp1.Category, m.top$mean$Spp2.Category, sep = "-"))
m.bott$mean$pair <- as.factor(paste(m.bott$mean$Spp1.Category, m.bott$mean$Spp2.Category, sep = "-"))

m.bott$mean$spp.pair <- as.factor(paste(m.bott$mean$Spp1, m.bott$mean$Spp2, sep = "-"))

m.top.pair.counts <- count(m.top$mean$pair)
m.bott.pair.counts <- count(m.bott$mean$pair)

m.bott.spp.pair.counts <- count(m.bott$mean$spp.pair)

m.top.pair.counts <- separate(m.top.pair.counts, x, c('from', 'to'), sep = '-')
m.bott.pair.counts <- separate(m.bott.pair.counts, x, c('from', 'to'), sep = '-')

m.bott.spp.pair.counts <- separate(m.bott.spp.pair.counts, x, c('from', 'to'), sep = '-')
f.bott.spp.pair.counts <- separate(f.bott.spp.pair.counts, x, c('from', 'to'), sep = '-')

m.group.richness <- c(length(which(spp.regions.m$Region == 'I1')), 
                      length(which(spp.regions.m$Region == 'M1')),
                      length(which(spp.regions.m$Region == 'I2')),
                      length(which(spp.regions.m$Region == 'M2')))

m.top.pair.counts$Group1.Rich <- as.integer(mapvalues(m.top.pair.counts$from,
                                                      from = c('I1', 'M1', 'I2', 'M2'),
                                                      to = m.group.richness))
m.top.pair.counts$Group2.Rich <- as.integer(mapvalues(m.top.pair.counts$to,
                                                      from = c('I1', 'M1', 'I2', 'M2'),
                                                      to = m.group.richness))
m.bott.pair.counts$Group1.Rich <- as.integer(mapvalues(m.bott.pair.counts$from,
                                                       from = c('I1', 'M1', 'I2', 'M2'),
                                                       to = m.group.richness))
m.bott.pair.counts$Group2.Rich <- as.integer(mapvalues(m.bott.pair.counts$to,
                                                       from = c('I1', 'M1', 'I2', 'M2'),
                                                       to = m.group.richness))


# Calculate edge weights. 
m.top.pair.counts$weight <- ifelse(m.top.pair.counts$Group1.Rich != m.top.pair.counts$Group2.Rich,
                                   (m.top.pair.counts$freq / (m.top.pair.counts$Group1.Rich * m.top.pair.counts$Group2.Rich)),
                                   (m.top.pair.counts$freq / m.top.pair.counts$Group1.Rich))
m.bott.pair.counts$weight <- ifelse(m.bott.pair.counts$Group1.Rich != m.bott.pair.counts$Group2.Rich,
                                    (m.bott.pair.counts$freq / (m.bott.pair.counts$Group1.Rich * m.bott.pair.counts$Group2.Rich)),
                                    (m.bott.pair.counts$freq / m.bott.pair.counts$Group1.Rich))

m.div.net <- graph_from_data_frame(d = m.top$mean[,c(4:5,3)], directed = F)
m.cov.net <- graph_from_data_frame(d = m.bott$mean[,c(6:7,3)], directed = F)

m.cov.spp.net <- graph_from_data_frame(d = m.bott$mean[,c(1:2,3)], directed = F)

m.div.edges <- igraph::as_data_frame(m.div.net, what="edges")
m.div.links <- m.top.pair.counts[,c(1,2,3)]
m.div.links$weight <- rescale(m.div.links$freq, to = c(2,8))
m.div.links <- m.div.links[order(m.div.links$from, m.div.links$to),]
m.div.links$type <- rep('Divergence', nrow(m.div.links))
m.div.links$color <- rep('#ef8a62', nrow(m.div.links))
rownames(m.div.links) <- NULL
m.div.nodes <- igraph::as_data_frame(m.div.net, what="vertices")
m.div.net <- graph_from_data_frame(d=m.div.links, vertices=m.div.nodes, directed=F)
E(m.div.net)$width <- E(m.div.net)$weight
V(m.div.net)$color <- 'grey87'

m.cov.edges <- igraph::as_data_frame(m.cov.net, what="edges")
m.cov.links <- m.bott.pair.counts[,c(1,2,3)]
m.cov.links$weight <- rescale(m.cov.links$freq, to = c(2,8))
m.cov.links <- m.cov.links[order(m.cov.links$from, m.cov.links$to),]
m.cov.links$type <- rep('Convergence', nrow(m.cov.links))
m.cov.links$color <- rep('#67a9cf', nrow(m.cov.links))
rownames(m.cov.links) <- NULL
m.cov.nodes <- igraph::as_data_frame(m.cov.net, what="vertices")
m.cov.net <- graph_from_data_frame(d=m.cov.links, vertices=m.cov.nodes, directed=F)
E(m.cov.net)$width <- E(m.cov.net)$weight
V(m.cov.net)$color <- 'grey87'

m.cov.spp.net <- graph_from_data_frame(d = m.bott[,c(1:2)], directed = F)
m.cov.spp.edges <- igraph::as_data_frame(m.cov.spp.net, what="edges")
m.cov.spp.links <- m.bott.spp.pair.counts[,c(1,2)]
m.cov.spp.links$type <- rep('Convergence', nrow(m.cov.spp.links))
m.cov.spp.links$color <- rep('#67a9cf', nrow(m.cov.spp.links))
rownames(m.cov.spp.links) <- NULL
m.cov.spp.nodes <- igraph::as_data_frame(m.cov.spp.net, what="vertices")
m.cov.spp.net <- graph_from_data_frame(d=m.cov.spp.links, vertices=m.cov.spp.nodes, directed=F)
V(m.cov.spp.net)$color <- 'grey87'

# We can even set the network layout:
#graph_attr(net, "layout") <- layout_with_lgl

pdf('Male_NoHeight_mvBM_Convergence_Network.pdf', width=9, height=8)
plot(m.cov.net, edge.curved=0.2, 
     edge.color = E(m.cov.net)$color,
     vertex.size=20, vertex.label.color="black") 
dev.off()

pdf('Male_mvBM_Species_Convergence_Network.pdf', width=9, height=8)
plot(m.cov.spp.net, edge.curved=0.3, 
     edge.color = E(m.cov.spp.net)$color, cex = 0.5,
     vertex.size=1, vertex.label.color="black") 
dev.off()

pdf('Male_NoHeight_mvBM_Divergence_Network.pdf', width=9, height=8)
plot(m.div.net, edge.curved=0.2, 
     edge.color = E(m.div.net)$color,
     vertex.size=20, vertex.label.color="black") 
dev.off()

