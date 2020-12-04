library(tidyverse)
library(reshape)
library(wesanderson)
library(ape)
library(otuSummary)
library(adephylo)

setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/Convergence/')
source('../ConvergenceDivergenceFunctions.R')

# Get species-region correspondance
tree <- read.tree('../FinalData/Anole_TimeCal_NoOut.tree')
m1 <- read.tree('../FinalData/Mainland_Clade1_TimeTree.tree')
m2 <- read.tree('../FinalData/Mainland_Clade2_TimeTree.tree')
i1 <- read.tree('../FinalData/Island_Clade1_TimeTree.tree')
i2 <- read.tree('../FinalData/Island_Clade2_TimeTree.tree')

spp <- c(i1$tip.label, i2$tip.label, m1$tip.label, m2$tip.label)
regions <- c(rep('I1', length(i1$tip.label)), rep('I2', length(i2$tip.label)),
             rep('M1', length(m1$tip.label)), rep('M2', length(m2$tip.label)))
spp.regions <- data.frame('Species' = spp, 'Region' = regions)

spp.regions.m <- spp.regions[which(spp.regions$Species %in% tree$tip.label),]

all.pstd <- readRDS('Male_mvBM_PSTD.rds')

all.pstd <- PreparePSTD(pstd = all.pstd, spp.cats = spp.regions.m)

all.pstd <- GetConvergent(pstd.dists = all.pstd, 
                          threshold = 1, 
                          tree = tree, minDist = 0,
                          dist.type = 'PatristicDist')

spp.dists <- read.csv('../FinalData/PoeDat_AP.csv')
spp.dists <- spp.dists[which(spp.dists$Species %in% c(i2$tip.label, m2$tip.label)),]

# Reduce down to only those species involved in convergences
# spp.dists <- spp.dists[which(spp.dists$Species %in% 
#                                c(as.character(m.conv$Spp1), as.character(m.conv$Spp2))),]

amaz <- as.character(spp.dists[which(spp.dists$Amazonia == 1),1])
andes <- as.character(spp.dists[which(spp.dists$Andes == 1),1])
baham <- as.character(spp.dists[which(spp.dists$Bahamas == 1),1])
caribCo <- as.character(spp.dists[which(spp.dists$CaribbeanCoast == 1),1])
cuba <- as.character(spp.dists[which(spp.dists$CubaCayman == 1),1])
hisp <- as.character(spp.dists[which(spp.dists$Hispaniola == 1),1])
jama <- as.character(spp.dists[which(spp.dists$Jamaica == 1),1])
lca <- as.character(spp.dists[which(spp.dists$LowerCentralAmerica == 1),1])
choco <- as.character(spp.dists[which(spp.dists$Choco == 1),1])
near <- as.character(spp.dists[which(spp.dists$Nearctic == 1),1])
pr <- as.character(spp.dists[which(spp.dists$PuertoRico == 1),1])
smIsl <- as.character(spp.dists[which(spp.dists$SmallIslands == 1),1])
uca <- as.character(spp.dists[which(spp.dists$UpperCentralAmerica == 1),1])
GrAntil <- c(baham, cuba, hisp, jama, pr)


all.pstd$median$Patry <- NA
for(i in 1:nrow(all.pstd$median)){
  # go down the list of PSTD and test if the two species share a region. 
  spp1 <- all.pstd$median$Spp1[i]
  spp2 <- all.pstd$median$Spp2[i]
  
  symp <- 
    sum(spp1 %in% near & spp2 %in% near, 
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
  
  all.pstd$median$Patry[i] <- 
    ifelse(symp == 0, 'Allopatric', 'Sympatric')
}


conv <- all.pstd$median[which(all.pstd$median$PSTD <= quantile(all.pstd$median$PSTD, 0.01)[[1]]),]

# Test of overabundance of allopatric convergence
i2.pstd.results <- all.pstd$median[which(all.pstd$median$Spp1.Category == 'I2' &
                                     all.pstd$median$Spp2.Category == 'I2'),]
m1.pstd.results <- all.pstd$median[which(all.pstd$median$Spp1.Category == 'M1' &
                                           all.pstd$median$Spp2.Category == 'M1'),]
m2.pstd.results <- all.pstd$median[which(all.pstd$median$Spp1.Category == 'M2' &
                                           all.pstd$median$Spp2.Category == 'M2'),]


expected.allo <- c(sum(i2.pstd.results$Patry == 'Allopatric'), 
                   sum(m2.pstd.results$Patry == 'Allopatric'))
expected.prob.allo <- c(sum(i2.pstd.results$Patry == 'Allopatric') / nrow(i2.pstd.results), 
                        sum(m2.pstd.results$Patry == 'Allopatric') / nrow(m2.pstd.results))

expected.symp <- c(sum(i2.pstd.results$Patry == 'Sympatric'), 
                   sum(m2.pstd.results$Patry == 'Sympatric'))
expected.prob.symp <- c(sum(i2.pstd.results$Patry == 'Sympatric') / nrow(i2.pstd.results), 
                        sum(m2.pstd.results$Patry == 'Sympatric') / nrow(m2.pstd.results))

# expected.prop.patry <- c(expected.allo,expected.symp)
# expected.probs.patry <- c(expected.prob.allo,expected.prob.symp)
# 

#m.conv <- read.csv('Male_mvBM_Median_Convergences.csv')

obs.allo <- c(sum(conv[which(conv$Spp1.Category == 'I2' &
                                              conv$Spp2.Category == 'I2'),]$Patry == 'Allopatric'), 
              sum(conv[which(conv$Spp1.Category == 'M2' &
                                              conv$Spp2.Category == 'M2'),]$Patry == 'Allopatric'))

obs.symp <- c(sum(conv[which(conv$Spp1.Category == 'I2' &
                               conv$Spp2.Category == 'I2'),]$Patry == 'Sympatric'), 
              sum(conv[which(conv$Spp1.Category == 'M2' &
                               conv$Spp2.Category == 'M2'),]$Patry == 'Sympatric'))

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

conv$Within <- conv$Spp1.Category

pstd.patry <- 
  ggplot(data = conv, 
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
  ggplot(data = conv, 
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
  ggplot(data = conv, 
         aes(x = Patry, y = PSTD/PatristicDist)) +
  geom_boxplot(aes(fill = Within), outlier.colour = NA) +
  geom_point(aes(fill = Within, shape = Patry), 
             position = position_jitterdodge(jitter.width = 0.1),
             pch = 21, size = 3, alpha = 0.5) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) + 
  xlab('Distribution of Species Pair') +
  theme_bw()



############################################
# Look at which regions are more involved in 
# allopatric vs sympatric convergences in M2
m2.conv <- conv[which(conv$Spp1 %in% m2$tip.label),]

m2.spp.reg <- spp.dists[,-c(2:6,14:21)]

# Now we want to merge the species region matrix
# First Species 1
colnames(m2.spp.reg)[1] <- 'Spp1'
m2.conv <- merge(m2.conv, m2.spp.reg, by = 'Spp1')

# Then Species 2 - regions will get renamed as '"'region.x' and 'region.y'
colnames(m2.spp.reg)[1] <- 'Spp2'
m2.conv <- merge(m2.conv, m2.spp.reg, by = 'Spp2')

#amaz.conv <- m2.conv[,]

region.pairs <- outer(colnames(m2.spp.reg)[-1], colnames(m2.spp.reg)[-1], paste, sep = '-')

region.count.mat <- matrix(ncol = 7, nrow = 7, 0)
region.count.mat[upper.tri(region.count.mat)] <- NA
colnames(region.count.mat) <- c(colnames(m2.spp.reg[-1]))
rownames(region.count.mat) <- c(colnames(m2.spp.reg[-1]))

upper <- region.pairs[upper.tri(region.pairs)]
lower <- region.pairs[lower.tri(region.pairs)]
within <- diag(region.pairs)




region.pairs <- list(c(within[1], within[1]), c(within[2], within[2]),
                     c(within[3], within[3]), c(within[4], within[4]),
                     c(within[5], within[5]), c(within[6], within[6]),
                     c(within[7], within[7]),
                     c(lower[1], upper[1]), c(lower[2], upper[2]),
                     c(lower[3], upper[4]), c(lower[4], upper[7]),
                     c(lower[5], upper[11]), c(lower[6], upper[16]),
                     c(lower[7], upper[3]), c(lower[8], upper[5]),
                     c(lower[9], upper[8]), c(lower[10], upper[12]),
                     c(lower[11], upper[17]), c(lower[12], upper[6]),
                     c(lower[13], upper[9]), c(lower[14], upper[13]),
                     c(lower[15], upper[18]), c(lower[16], upper[10]),
                     c(lower[17], upper[14]), c(lower[18], upper[19]),
                     c(lower[19], upper[15]), c(lower[20], upper[20]),
                     c(lower[21], upper[21]))
region.pair.levels <- c(within, upper)



for(i in 1:nrow(m2.conv)){
  # For that row get species ids
  spp1 <- as.character(m2.conv$Spp1[i])
  spp2 <- as.character(m2.conv$Spp2[i])
  
  # # Get the row ID
  # r1 <- which(m2.conv$Species == spp1)
  # r2 <- which(m2.conv$Species == spp2)
  
  # Now summarize how many times each region is involved in a sympatric or allopatric convergence
  
  # Here, if the two species are both present in the same region (Sum = 2), add one to respective diagonal 
  if(m2.conv[i,]$Nearctic.x + m2.conv[i,]$Nearctic.y == 2) region.count.mat[1,1] <- region.count.mat[1,1] + 1
  if(m2.conv[i,]$UpperCentralAmerica.x + m2.conv[i,]$UpperCentralAmerica.y == 2) region.count.mat[2,2] <- region.count.mat[2,2] + 1
  if(m2.conv[i,]$LowerCentralAmerica.x + m2.conv[i,]$LowerCentralAmerica.y == 2) region.count.mat[3,3] <- region.count.mat[3,3] + 1
  if(m2.conv[i,]$Choco.x + m2.conv[i,]$Choco.y == 2) region.count.mat[4,4] <- region.count.mat[4,4] + 1
  if(m2.conv[i,]$CaribbeanCoast.x + m2.conv[i,]$CaribbeanCoast.y == 2) region.count.mat[5,5] <- region.count.mat[5,5] + 1
  if(m2.conv[i,]$Andes.x + m2.conv[i,]$Andes.y == 2) region.count.mat[6,6] <- region.count.mat[6,6] + 1
  if(m2.conv[i,]$Amazonia.x + m2.conv[i,]$Amazonia.y == 2) region.count.mat[7,7] <- region.count.mat[7,7] + 1
  
  # Now, test if the two species are present in any of the pairwise region comparisons. If sum >= 2, add 1 to that cell
  if(sum(m2.conv[i,]$Nearctic.x + m2.conv[i,]$UpperCentralAmerica.y) + 
     sum(m2.conv[i,]$Nearctic.y + m2.conv[i,]$UpperCentralAmerica.x) >= 2) region.count.mat[2,1] <- region.count.mat[2,1] + 1
  if(sum(m2.conv[i,]$Nearctic.x + m2.conv[i,]$LowerCentralAmerica.y) + 
     sum(m2.conv[i,]$Nearctic.y + m2.conv[i,]$LowerCentralAmerica.x) >= 2) region.count.mat[3,1] <- region.count.mat[3,1] + 1
  if(sum(m2.conv[i,]$Nearctic.x + m2.conv[i,]$Choco.y) + 
     sum(m2.conv[i,]$Nearctic.y + m2.conv[i,]$Choco.x) >= 2) region.count.mat[4,1] <- region.count.mat[4,1] + 1
  if(sum(m2.conv[i,]$Nearctic.x + m2.conv[i,]$CaribbeanCoast.y) + 
     sum(m2.conv[i,]$Nearctic.y + m2.conv[i,]$CaribbeanCoast.x) >=2) region.count.mat[5,1] <- region.count.mat[5,1] + 1
  if(sum(m2.conv[i,]$Nearctic.x + m2.conv[i,]$Andes.y) + 
     sum(m2.conv[i,]$Nearctic.y + m2.conv[i,]$Andes.x) >=2) region.count.mat[6,1] <- region.count.mat[6,1] + 1
  if(sum(m2.conv[i,]$Nearctic.x + m2.conv[i,]$Amazonia.y) + 
     sum(m2.conv[i,]$Nearctic.y + m2.conv[i,]$Amazonia.x) >=2) region.count.mat[7,1] <- region.count.mat[7,1] + 1
     
  if(sum(m2.conv[i,]$UpperCentralAmerica.x + m2.conv[i,]$LowerCentralAmerica.y) + 
     sum(m2.conv[i,]$UpperCentralAmerica.y + m2.conv[i,]$LowerCentralAmerica.x) >=2) region.count.mat[3,2] <- region.count.mat[3,2] + 1
  if(sum(m2.conv[i,]$UpperCentralAmerica.x + m2.conv[i,]$Choco.y) + 
     sum(m2.conv[i,]$UpperCentralAmerica.y + m2.conv[i,]$Choco.x) >=2) region.count.mat[4,2] <- region.count.mat[4,2] + 1
  if(sum(m2.conv[i,]$UpperCentralAmerica.x + m2.conv[i,]$CaribbeanCoast.y) + 
     sum(m2.conv[i,]$UpperCentralAmerica.y + m2.conv[i,]$CaribbeanCoast.x) >=2) region.count.mat[5,2] <- region.count.mat[5,2] + 1
  if(sum(m2.conv[i,]$UpperCentralAmerica.x + m2.conv[i,]$Andes.y) + 
     sum(m2.conv[i,]$UpperCentralAmerica.y + m2.conv[i,]$Andes.x) >=2) region.count.mat[6,2] <- region.count.mat[6,2] + 1
  if(sum(m2.conv[i,]$UpperCentralAmerica.x + m2.conv[i,]$Amazonia.y) + 
     sum(m2.conv[i,]$UpperCentralAmerica.y + m2.conv[i,]$Amazonia.x) >=2) region.count.mat[7,2] <- region.count.mat[7,2] + 1
  
  if(sum(m2.conv[i,]$LowerCentralAmerica.x + m2.conv[i,]$Choco.y) + 
     sum(m2.conv[i,]$LowerCentralAmerica.y + m2.conv[i,]$Choco.x) >=2) region.count.mat[4,3] <- region.count.mat[4,3] + 1
  if(sum(m2.conv[i,]$LowerCentralAmerica.x + m2.conv[i,]$CaribbeanCoast.y) + 
     sum(m2.conv[i,]$LowerCentralAmerica.y + m2.conv[i,]$CaribbeanCoast.x) >=2) region.count.mat[5,3] <- region.count.mat[5,3] + 1
  if(sum(m2.conv[i,]$LowerCentralAmerica.x + m2.conv[i,]$Andes.y) + 
     sum(m2.conv[i,]$LowerCentralAmerica.y + m2.conv[i,]$Andes.x) >=2) region.count.mat[6,3] <- region.count.mat[6,3] + 1
  if(sum(m2.conv[i,]$LowerCentralAmerica.x + m2.conv[i,]$Amazonia.y) + 
     sum(m2.conv[i,]$LowerCentralAmerica.y + m2.conv[i,]$Amazonia.x) >=2) region.count.mat[7,3] <- region.count.mat[7,3] + 1 
  
  if(sum(m2.conv[i,]$Choco.x + m2.conv[i,]$CaribbeanCoast.y) + 
     sum(m2.conv[i,]$Choco.y + m2.conv[i,]$CaribbeanCoast.x) >=2) region.count.mat[5,4] <- region.count.mat[5,4] + 1
  if(sum(m2.conv[i,]$Choco.x + m2.conv[i,]$Andes.y) + 
     sum(m2.conv[i,]$Choco.y + m2.conv[i,]$Andes.x) >=2) region.count.mat[6,4] <- region.count.mat[6,4] + 1
  if(sum(m2.conv[i,]$Choco.x + m2.conv[i,]$Amazonia.y) + 
     sum(m2.conv[i,]$Choco.y + m2.conv[i,]$Amazonia.x) >=2) region.count.mat[7,4] <- region.count.mat[7,4] + 1 
  
  if(sum(m2.conv[i,]$CaribbeanCoast.x + m2.conv[i,]$Andes.y) + 
     sum(m2.conv[i,]$CaribbeanCoast.y + m2.conv[i,]$Andes.x) >=2) region.count.mat[6,5] <- region.count.mat[6,5] + 1
  if(sum(m2.conv[i,]$CaribbeanCoast.x + m2.conv[i,]$Amazonia.y) + 
     sum(m2.conv[i,]$CaribbeanCoast.y + m2.conv[i,]$Amazonia.x) >=2) region.count.mat[7,5] <- region.count.mat[7,5] + 1 
  
  if(sum(m2.conv[i,]$Andes.x + m2.conv[i,]$Amazonia.y) + 
     sum(m2.conv[i,]$Andes.y + m2.conv[i,]$Amazonia.x) >=2) region.count.mat[7,6] <- region.count.mat[7,6] + 1 
}


region.counts <- melt(region.count.mat) %>% na.omit()
colnames(region.counts) <- c('Region1', 'Region2', 'Count')

region.counts$Region1 <- factor(region.counts$Region1, 
                                levels = c('Amazonia', 'Andes', 'Choco', 'CaribbeanCoast',
                                           'LowerCentralAmerica', 'UpperCentralAmerica', 'Nearctic'))
region.counts$Region2 <- factor(region.counts$Region2, 
                                levels = c('Amazonia', 'Andes', 'Choco', 'CaribbeanCoast',
                                           'LowerCentralAmerica', 'UpperCentralAmerica', 'Nearctic'))

#region.counts[2,c(1:2)] <- region.counts[2,c(2,1)]
region.counts[20,c(1:2)] <- region.counts[20,c(2,1)]

pdf('ConvergenceCounts-PerRegion-M2.pdf', width = 8, height = 7)
ggplot(data = region.counts, aes(x = Region2, y = Region1, fill = Count)) +
  scale_fill_gradient(low = '#ffffb2', high = '#b10026') +
  geom_tile(color = 'black') + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
# 
# 
#   # Get which regions each species is in
#   region1 <- colnames(mland.spp.reg)[which(mland.spp.reg[r1,] == 1)]
#   region2 <- colnames(mland.spp.reg)[which(mland.spp.reg[r2,] == 1)]
#   
#   
#   # Now tally up counts of convergences for region-pairs
#   
#   for(r in 1:nrow(region.conv.count)){
#     if(m2.conv$Patry[i] == 'Allopatric'){
#       if(region.conv.count[r,1] %in% region1){
#         region.conv.count[r,2] <-  region.conv.count[r,2] + 1
#       }
#       if(region.conv.count[r,1] %in% region2){
#         region.conv.count[r,2] <-  region.conv.count[r,2] + 1
#       }
#     }else{
#       if(region.conv.count[r,1] %in% region1){
#         region.conv.count[r,3] <-  region.conv.count[r,3] + 1
#       }
#       if(region.conv.count[r,1] %in% region2){
#         region.conv.count[r,3] <-  region.conv.count[r,3] + 1
#       }
#     }
#   }
# }
# 
# region.conv <- data.frame('Region' = rep(region.conv.count[,1], 2),
#                           'Type' = c(rep('Allopatric', 7), 
#                                      rep('Sympatric', 7)),
#                           'Count' = c(region.conv.count[,2], 
#                                       region.conv.count[,3]))
# 
# region.conv$Region <- factor(region.conv$Region, 
#                              levels = c('Amazonia', 'Andes', 'Choco', 
#                                         'CaribbeanCoast', 'LowerCentralAmerica', 
#                                         'UpperCentralAmerica', 'Nearctic'))
# 
# pdf('M2-Region-ConvergenceCounts.pdf', width = 8, height = 8)
# ggplot(data = region.conv, aes(x = Region, y = Count)) +
#   geom_point(aes(fill = Type), pch = 21, size = 3,
#              position = position_jitterdodge(jitter.width = 0.25)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# dev.off()
# 


