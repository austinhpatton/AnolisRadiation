########### Morphospace Comparison - I1, I2, M1, M2 #####################
library(treeplyr)
library(purrr)
library(phytools)
library(geiger)
library(gridExtra)
library(ggpubr)

setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/')

m.data <- read.csv('./FinalData/Male-ALL-DATA-Trait-Ecol-Regions-FINAL.csv', 
                   header = T)

tree <- read.tree('./FinalData/Anole_TimeCal_NoOut.tree')
m1 <- read.tree('./FinalData/M1.tree')
m2 <- read.tree('./FinalData/M2.tree')
i1 <- read.tree('./FinalData/I1.tree')
i2 <- read.tree('./FinalData/I2.tree')

m.data <- m.data[,c(1,3:10,12:13,23,41)]
m.data <- m.data[-which(m.data$Species %in% c('agassizi', 'bicaorum', 'lineatus', 'roatanensis', 
                                           'concolor', 'pinchoti', 'utilensis', 'townsendi')),]

# Now remove missing data
m.data <- na.omit(m.data[,-c(12:13)])
rownames(m.data) <- m.data$Species
m.data <- m.data[,-1]

m.full.td <- make.treedata(tree, m.data[,-c(12:13)])
m.full.data <- as.data.frame(m.full.td$dat)
rownames(m.full.data) <- m.full.td$phy$tip.label

m.full <- m.full.td$phy

# Now run the pPCA
m.ppca <- phyl.pca(tree, m.full.data, mode = 'corr')

# Extract the scores
m.scores <- m.ppca$S 

i1.scores <- m.scores[which(rownames(m.scores) %in% i1$tip.label),]
m1.scores <- m.scores[which(rownames(m.scores) %in% m1$tip.label),]
i2.scores <- m.scores[which(rownames(m.scores) %in% i2$tip.label),]
m2.scores <- m.scores[which(rownames(m.scores) %in% m2$tip.label),]

# Make clade specific trees
m.i1 <- make.treedata(i1, m.data)$phy
m.m1 <- make.treedata(m1, m.data)$phy
m.i2 <- make.treedata(i2, m.data)$phy
m.m2 <- make.treedata(m2, m.data)$phy

# And create dataframes corresponding scores to regions
m.scores <- rbind(i1.scores, m1.scores, i2.scores, m2.scores)
m.scores <- as.data.frame(m.scores, row.names = row.names(m.scores))
m.scores$Region <- c(rep('I1', nrow(i1.scores)),
                     rep('M1', nrow(m1.scores)),
                     rep('I2', nrow(i2.scores)), 
                     rep('M2', nrow(m2.scores)))


# Now we need to attribute species to regions. Each external node will be a pie chart, with colors equal to biogeographic regions
spp.dists <- read.csv('./FinalData/PoeDat_AP.csv')
spp.dists <- spp.dists[which(spp.dists$Species %in% c(m1$tip.label, m2$tip.label, 
                                                      i1$tip.label, i2$tip.label)),c(1,4,7:19)]

for(i in 1:nrow(spp.dists)){
  spp.dists[i,-1] <- spp.dists[i,-1]/sum(spp.dists[i,-1])
}
m.scores$Species <- rownames(m.scores)

m.scores <- merge(m.scores, spp.dists, by = 'Species')

m1.regs <- as.matrix(m.scores[which(m.scores$Species %in% m1$tip.label),14:20])
rownames(m1.regs) <- m.scores[which(m.scores$Species %in% m1$tip.label),]$Species
m2.regs <- as.matrix(m.scores[which(m.scores$Species %in% m2$tip.label),14:20])
rownames(m2.regs) <- m.scores[which(m.scores$Species %in% m2$tip.label),]$Species

i1.regs <- as.matrix(m.scores[which(m.scores$Species %in% i1$tip.label),c(13,21:26)])
rownames(i1.regs) <- m.scores[which(m.scores$Species %in% i1$tip.label),]$Species
i2.regs <- as.matrix(m.scores[which(m.scores$Species %in% i2$tip.label),c(13,21:26)])
rownames(i2.regs) <- m.scores[which(m.scores$Species %in% i2$tip.label),]$Species



colors <- c('#8dd3c7', '#ffffb3', '#bebada', '#fb8072', 
            '#80b1d3', '#fdb462', '#b3de69')

setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/')
pdf('Mainland-Phylomorphospace-pPCs1-3.pdf', width = 6, height = 10)
par(mfrow = c(3,2))

phylomorphospace(tree = m.m1, X = data.frame(m1.scores[,1:2]), 
                 xlim = c(-27,24), ylim = c(-27,24), label = F)
tiplabels(pie=m1.regs, piecol=colors, cex=1.5) 
legend('topleft',
       legend = colnames(m2.regs), 
       title = "Region",
       fill = colors,
       cex = 1,
       bty = "n")

phylomorphospace(tree = m.m2, X = data.frame(m2.scores[,1:2]), 
                 xlim = c(-27,24), ylim = c(-27,24), label = F)
tiplabels(pie=m2.regs, piecol=colors, cex=1.5)


phylomorphospace(tree = m.m1, X = data.frame(m1.scores[,2:3]), 
                 xlim = c(-27,24), ylim = c(-27,24), label = F)
tiplabels(pie=m1.regs, piecol=colors, cex=1.5)


phylomorphospace(tree = m.m2, X = data.frame(m2.scores[,2:3]), 
                 xlim = c(-27,24), ylim = c(-27,24), label = F)
tiplabels(pie=m2.regs, piecol=colors, cex=1.5)

phylomorphospace(tree = m.m1, X = data.frame(m1.scores[,c(1,3)]), 
                 xlim = c(-27,24), ylim = c(-27,24), label = F)
tiplabels(pie=m1.regs, piecol=colors, cex=1.5)

phylomorphospace(tree = m.m2, X = data.frame(m2.scores[,c(1,3)]), 
                 xlim = c(-27,24), ylim = c(-27,24), label = F)
tiplabels(pie=m2.regs, piecol=colors, cex=1.5)

dev.off()


setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/')
pdf('Mainland-Phylomorphospace-pPCs1-2.pdf', width = 10, height = 10)
par(mfrow = c(2,2))

phylomorphospace(tree = m.m1, X = data.frame(m1.scores[,1:2]), 
                 xlim = c(-27,24), ylim = c(-27,24), label = F)
tiplabels(pie=m1.regs, piecol=colors, cex=1) 
legend('topleft',
       legend = colnames(m2.regs), 
       title = "Region",
       fill = colors,
       cex = 1,
       bty = "n")

phylomorphospace(tree = m.m2, X = data.frame(m2.scores[,1:2]), 
                 xlim = c(-27,24), ylim = c(-27,24), label = F)
tiplabels(pie=m2.regs, piecol=colors, cex=1)


phylomorphospace(tree = m.m1, X = data.frame(m1.scores[,2:3]), 
                 xlim = c(-27,24), ylim = c(-27,24), label = F)
tiplabels(pie=m1.regs, piecol=colors, cex=1)


phylomorphospace(tree = m.m2, X = data.frame(m2.scores[,2:3]), 
                 xlim = c(-27,24), ylim = c(-27,24), label = F)
tiplabels(pie=m2.regs, piecol=colors, cex=1)

dev.off()


pdf('Mainland-Phylomorphospace-pPCs3-5.pdf', width = 10, height = 9)
par(mfrow = c(2,3))

phylomorphospace(tree = m.m1, X = data.frame(m.m1.scores[,3:4]), 
                 xlim = c(-27,24), ylim = c(-27,24))
tiplabels(pie=m1.regs, piecol=colors, cex=1)
legend('topleft',
       legend = colnames(m2.regs), 
       title = "Region",
       fill = colors,
       cex = 1,
       bty = "n")

phylomorphospace(tree = m.m1, X = data.frame(m.m1.scores[,4:5]), 
                 xlim = c(-27,24), ylim = c(-27,24))
tiplabels(pie=m1.regs, piecol=colors, cex=1)

phylomorphospace(tree = m.m1, X = data.frame(m.m1.scores[,c(3,5)]), 
                 xlim = c(-27,24), ylim = c(-27,24))
tiplabels(pie=m1.regs, piecol=colors, cex=1)


phylomorphospace(tree = m.m2, X = data.frame(m.m2.scores[,3:4]), 
                 xlim = c(-27,24), ylim = c(-27,24))
tiplabels(pie=m2.regs, piecol=colors, cex=1)

phylomorphospace(tree = m.m2, X = data.frame(m.m2.scores[,4:5]), 
                 xlim = c(-27,24), ylim = c(-27,24))
tiplabels(pie=m2.regs, piecol=colors, cex=1)

phylomorphospace(tree = m.m2, X = data.frame(m.m2.scores[,c(3,5)]), 
                 xlim = c(-27,24), ylim = c(-27,24))
tiplabels(pie=m2.regs, piecol=colors, cex=1)

dev.off()

# not including carolinensis here
i2.pmorph <- drop.tip(m.i2, 'carolinensis')
##### Now the islands
pdf('Island-Phylomorphospace-pPCs1-3.pdf', width = 10, height = 9)
par(mfrow = c(2,3))

phylomorphospace(tree = m.i1, X = data.frame(m.i1.scores[,1:2]), 
                 xlim = c(-1, 1), ylim = c(-1, 1))
tiplabels(pie=i1.regs, piecol=colors, cex=1) 
legend('topleft',
       legend = colnames(i1.regs), 
       title = "Region",
       fill = colors,
       cex = 1,
       bty = "n")

phylomorphospace(tree = m.i1, X = data.frame(m.i1.scores[,2:3]), 
                 xlim = c(-1, 1), ylim = c(-1, 1))
tiplabels(pie=i1.regs, piecol=colors, cex=1)

phylomorphospace(tree = m.i1, X = data.frame(m.i1.scores[,c(1,3)]), 
                 xlim = c(-1, 1), ylim = c(-1, 1))
tiplabels(pie=i1.regs, piecol=colors, cex=1)


# Again, making sure to remove carolinensis
phylomorphospace(tree = i2.pmorph, X = data.frame(m.i2.scores[-29,1:2]), 
                 xlim = c(-1, 1), ylim = c(-1, 1))
tiplabels(pie=i2.regs[-12,], piecol=colors, cex=1)

phylomorphospace(tree = i2.pmorph, X = data.frame(m.i2.scores[-29,2:3]), 
                 xlim = c(-1, 1), ylim = c(-1, 1))
tiplabels(pie=i2.regs[-12,], piecol=colors, cex=1)

phylomorphospace(tree = i2.pmorph, X = data.frame(m.i2.scores[-29,c(1,3)]), 
                 xlim = c(-1, 1), ylim = c(-1, 1))
tiplabels(pie=i2.regs[-12,], piecol=colors, cex=1)

dev.off()


pdf('Island-Phylomorphospace-pPCs3-5.pdf', width = 10, height = 9)
par(mfrow = c(2,3))

phylomorphospace(tree = m.i1, X = data.frame(m.i1.scores[,3:4]), 
                 xlim = c(-1, 1), ylim = c(-1, 1))
tiplabels(pie=i1.regs, piecol=colors, cex=1)
legend('topleft',
       legend = colnames(i2.regs), 
       title = "Region",
       fill = colors,
       cex = 1,
       bty = "n")

phylomorphospace(tree = m.i1, X = data.frame(m.i1.scores[,4:5]), 
                 xlim = c(-1, 1), ylim = c(-1, 1))
tiplabels(pie=i1.regs, piecol=colors, cex=1)

phylomorphospace(tree = m.i1, X = data.frame(m.i1.scores[,c(3,5)]), 
                 xlim = c(-1, 1), ylim = c(-1, 1))
tiplabels(pie=i1.regs, piecol=colors, cex=1)


phylomorphospace(tree = i2.pmorph, X = data.frame(m.i2.scores[-29,3:4]), 
                 xlim = c(-1, 1), ylim = c(-1, 1))
tiplabels(pie=i2.regs[-12,], piecol=colors, cex=1)

phylomorphospace(tree = i2.pmorph, X = data.frame(m.i2.scores[-29,4:5]), 
                 xlim = c(-1, 1), ylim = c(-1, 1))
tiplabels(pie=i2.regs[-12,], piecol=colors, cex=1)

phylomorphospace(tree = i2.pmorph, X = data.frame(m.i2.scores[-29,c(3,5)]), 
                 xlim = c(-1, 1), ylim = c(-1, 1))
tiplabels(pie=i2.regs[-12,], piecol=colors, cex=1)
dev.off()


pheno.phy <- drop.tip(m.full.td$phy, 
                      m.full.td$phy$tip.label[-which(m.full.td$phy$tip.label %in% m.scores$Species)])

regs <- m.scores$Region
names(regs) <- m.scores$Species
tmp <-m.scores$PC2
names(tmp) <- m.scores$Species

pheno.phy.regs <- make.simmap(pheno.phy, regs, model = 'SYM')
saveRDS(pheno.phy.regs, file = 'Male-Traits-Clade-Simmap.rds')

colors <- c('#faae61', '#d72027', '#add8a4', '#2b83bb')

x <- getStates(pheno.phy.regs,"tips")
cols <- setNames(colors[1:length(unique(x))],sort(unique(x)))

pdf('PC1-Phenogram.pdf', height = 10, width = 5)
tmp <- m.scores[,2]
names(tmp) <- m.scores$Species
phenogram(pheno.phy.regs, tmp, ylim = range(m.scores[,2]),
          fsize = 0.5, colors = cols, ftype = 'i',
          ylab = 'PC1')
add.simmap.legend(colors=cols,x=0.5*par()$usr[1],
                  y=0.9*par()$usr[4],prompt=FALSE)
dev.off()

pdf('PC2-Phenogram.pdf', height = 10, width = 5)
tmp <- m.scores[,3]
names(tmp) <- m.scores$Species
phenogram(pheno.phy.regs, tmp, ylim = range(m.scores[,3]),
          fsize = 0.5, colors = cols, ftype = 'i',
          ylab = 'PC2')
add.simmap.legend(colors=cols,x=0.5*par()$usr[1],
                  y=0.9*par()$usr[4],prompt=FALSE)
dev.off()

pdf('PC3-Phenogram.pdf', height = 10, width = 5)
tmp <- m.scores[,4]
names(tmp) <- m.scores$Species
phenogram(pheno.phy.regs, tmp, ylim = range(m.scores[,4]),
          fsize = 0.5, colors = cols, ftype = 'i',
          ylab = 'PC3')
add.simmap.legend(colors=cols,x=0.5*par()$usr[1],
                  y=0.9*par()$usr[4],prompt=FALSE)
dev.off()

pdf('PC4-Phenogram.pdf', height = 10, width = 5)
tmp <- m.scores[,5]
names(tmp) <- m.scores$Species
phenogram(pheno.phy.regs, tmp, ylim = range(m.scores[,5]),
          fsize = 0.5, colors = cols, ftype = 'i',
          ylab = 'PC4')
add.simmap.legend(colors=cols,x=0.5*par()$usr[1],
                  y=0.9*par()$usr[4],prompt=FALSE)
dev.off()

pdf('PC5-Phenogram.pdf', height = 10, width = 5)
tmp <- m.scores[,6]
names(tmp) <- m.scores$Species
phenogram(pheno.phy.regs, tmp, ylim = range(m.scores[,6]),
          fsize = 0.5, colors = cols, ftype = 'i',
          ylab = 'PC5')
add.simmap.legend(colors=cols,x=0.5*par()$usr[1],
                  y=0.9*par()$usr[4],prompt=FALSE)
dev.off()

for(i in 3:6){
  tmp <- m.scores[,i]
  names(tmp) <- m.scores$Species
  
  phenogram(pheno.phy.regs, tmp, ylim = range(m.scores[,i]),
            fsize = 0.5, colors = cols, ftype = 'i', add  = T)
}
dev.off()





# Now calculate the overlap among M1, M2, I1, and I2
library(hypervolume)
library(alphahull)

row.names(m.scores) <- m.scores$Species
m.scores <- m.scores[,-1]

m.regions <- read.csv('./FinalData/Male-ALL-DATA-Trait-Ecol-Regions-FINAL.csv', row.names = 1)
m2.regions <- m.regions[which(m.regions$Group != 'M1'),]
m2.symp <- rownames(m2.regions[which(m2.regions$MainlandPatry == 'M2-Sympatric'),])
m2.allo <- rownames(m2.regions[which(m2.regions$MainlandPatry == 'M2-Allopatric'),])
m2.symp.scores <- m.scores[which(rownames(m.scores) %in% m2.symp),]
m2.allo.scores <- m.scores[which(rownames(m.scores) %in% m2.allo),]


m.m1.hv.svm = hypervolume_svm(m1.scores[,1:3], name = 'M1', 
                              svm.gamma = 0.6, svm.nu = 0.005)
m.m2.hv.svm = hypervolume_svm(m2.scores[,1:3], name = 'M2', 
                              svm.gamma = 0.6, svm.nu = 0.005)
m.i2.hv.svm = hypervolume_svm(i2.scores[,1:3], name = 'I2', 
                              svm.gamma = 0.6, svm.nu = 0.005)

get_volume(m.m1.hv.svm)
get_volume(m.m2.hv.svm)
get_volume(m.i2.hv.svm)

m1m2.hv.svm.set <- hypervolume_set(m.m1.hv.svm, m.m2.hv.svm, check.memory=F)
m1i2.hv.svm.set <- hypervolume_set(m.m1.hv.svm, m.i2.hv.svm, check.memory=F)
m2i2.hv.svm.set <- hypervolume_set(m.m2.hv.svm, m.i2.hv.svm, check.memory=F)

hypervolume_overlap_statistics(m1m2.hv.svm.set)
hypervolume_overlap_statistics(m1i2.hv.svm.set)
hypervolume_overlap_statistics(m2i2.hv.svm.set)

hypervolume_distance(hv1 = m.m1.hv.svm, hv2 = m.m2.hv.svm )
hypervolume_distance(hv1 = m.m1.hv.svm, hv2 = m.i2.hv.svm )
hypervolume_distance(hv1 = m.m2.hv.svm, hv2 = m.i2.hv.svm )

colors <- c('#faae61', '#d72027', '#add8a4', '#2b83bb')

hv.cols <- colors[c(2,4,3)]
hv.list <- hypervolume_join(m.i2.hv.svm, m.m2.hv.svm, m.m1.hv.svm)

pdf('Morphospace-Hypervolumes-CorrMat-PCs1-3.pdf',width = 6, height = 6)
plot(hv.list, colors = hv.cols, cex.centroid = 2.5, contour.lwd = 2, 
     contour.kde.level = 0.01, point.alpha.min = 0.01,
     limits = c(-30, 40),
     cex.data = 1, cex.random = 0.4, 
     num.points.max.random = 1500,
     cex.legend = 2)
dev.off()


m.m2.symp.hv.svm = hypervolume_svm(m2.symp.scores[,1:3], name = 'M2-Sympatric', 
                              svm.gamma = 0.6, svm.nu = 0.005)
m.m2.allo.hv.svm = hypervolume_svm(m2.allo.scores[,1:3], name = 'M2-Allopatric', 
                              svm.gamma = 0.6, svm.nu = 0.005)

m1m2.symp.hv.svm.set <- hypervolume_set(m.m1.hv.svm, m.m2.symp.hv.svm, check.memory=F)
m1m2.allo.hv.svm.set <- hypervolume_set(m.m1.hv.svm, m.m2.allo.hv.svm, check.memory=F)
m2.symp.m2.allo.hv.svm.set <- hypervolume_set(m.m2.symp.hv.svm, m.m2.allo.hv.svm, check.memory=F)

mland.list <- hypervolume_join(m.m1.hv.svm, m.m2.symp.hv.svm, m.m2.allo.hv.svm)

hypervolume_overlap_statistics(m1m2.symp.hv.svm.set)
hypervolume_overlap_statistics(m1m2.allo.hv.svm.set)
hypervolume_overlap_statistics(m2.symp.m2.allo.hv.svm.set)

hypervolume_distance(hv1 = m.m1.hv.svm, hv2 = m.m2.symp.hv.svm )
hypervolume_distance(hv1 = m.m1.hv.svm, hv2 = m.m2.allo.hv.svm )
hypervolume_distance(hv1 = m.m2.symp.hv.svm, hv2 = m.m2.allo.hv.svm )

mland.cols <- c('#add8a4', '#ffc866', '#b2182b')

pdf('Mainland-Morphospace-Hypervolumes-CorrMat-PCs1-3.pdf',width = 6, height = 6)
plot(mland.list, colors = mland.cols, cex.centroid = 2.5, contour.lwd = 2, 
     contour.kde.level = 0.01, point.alpha.min = 0.01,
     limits = c(-26, 20),
     cex.data = 1, cex.random = 0.4, 
     num.points.max.random = 1500,
     cex.legend = 2)
dev.off()


# Now we need to attribute species to regions. Each external node will be a pie chart, with colors equal to biogeographic regions
spp.dists <- read.csv('./FinalData/PoeDat_AP.csv')
spp.dists <- spp.dists[which(spp.dists$Species %in% c(m1$tip.label, m2$tip.label, 
                                                      i1$tip.label, i2$tip.label)),c(1,4,7:19)]

for(i in 1:nrow(spp.dists)){
  spp.dists[i,-1] <- spp.dists[i,-1]/sum(spp.dists[i,-1])
}
m.scores$Species <- rownames(m.scores)

m.scores <- merge(m.scores, spp.dists, by = 'Species')

m2.sa.scores <- m.scores[which(m.scores$Region == 'M2' &
                                 rowSums(m.scores[17:20]) >= 1),1:11]
m2.rest.scores <- m.scores[which(m.scores$Region == 'M2' &
                                 rowSums(m.scores[17:20]) != 1),1:11]

m2.sa.rest <- rbind(m2.sa.scores, m2.rest.scores)
m2.sa.rest$Dist <- c(rep('SA', nrow(m2.sa.scores)),
                     rep('Rest', nrow(m2.rest.scores)))
m2.sa.rest.manova <- manova(as.matrix(m2.sa.rest[,2:11]) ~ Dist, m2.sa.rest)
summary(m2.sa.rest.manova)


m.m2.sa.hv.svm = hypervolume_svm(m2.sa.scores[,1:3], name = 'M2-SA', 
                                   svm.gamma = 0.6, svm.nu = 0.005)
m.m2.rest.hv.svm = hypervolume_svm(m2.rest.scores[,1:3], name = 'M2-Rest', 
                                   svm.gamma = 0.6, svm.nu = 0.005)

m1m2.sa.hv.svm.set <- hypervolume_set(m.m1.hv.svm, m.m2.sa.hv.svm, check.memory=F)
m1m2.rest.hv.svm.set <- hypervolume_set(m.m1.hv.svm, m.m2.rest.hv.svm, check.memory=F)
m2.sa.m2.rest.hv.svm.set <- hypervolume_set(m.m2.sa.hv.svm, m.m2.rest.hv.svm, check.memory=F)

mland.list.2 <- hypervolume_join(m.m1.hv.svm, m.m2.sa.hv.svm, m.m2.rest.hv.svm)

hypervolume_overlap_statistics(m1m2.sa.hv.svm.set)
hypervolume_overlap_statistics(m1m2.rest.hv.svm.set)
hypervolume_overlap_statistics(m2.sa.m2.rest.hv.svm.set)

hypervolume_distance(hv1 = m.m1.hv.svm, hv2 = m.m2.sa.hv.svm )
hypervolume_distance(hv1 = m.m1.hv.svm, hv2 = m.m2.rest.hv.svm )
hypervolume_distance(hv1 = m.m2.sa.hv.svm, hv2 = m.m2.rest.hv.svm )

mland.cols <- c('#add8a4', '#ffc866', '#b2182b')

pdf('Mainland-M2-SA-Rest-Hypervolumes-CorrMat-PCs1-3.pdf',width = 6, height = 6)
plot(mland.list.2, colors = mland.cols, cex.centroid = 2.5, contour.lwd = 2, 
     contour.kde.level = 0.01, point.alpha.min = 0.01,
     limits = c(-26, 20),
     cex.data = 1, cex.random = 0.4, 
     num.points.max.random = 1500,
     cex.legend = 2)
dev.off()





# Now determine what is loading on each PC and how much variance
# is explained by each PC
m.ppca.perVar <- round(diag(m.ppca$Eval)/sum(diag(m.ppca$Eval)),4) * 100

barplot(m.ppca.perVar, xlab = 'Variable', 
        ylab = 'Percent Variance Explained')

pdf('Male_CorrMat-pPCA_Loadings.pdf')
par(mfrow = c(2, 3))
biplot(m.ppca, col = c('grey','red'), choices = c(1,2), cex = 0.65)
biplot(m.ppca, col = c('grey','red'), choices = c(2,3), cex = 0.65)
biplot(m.ppca, col = c('grey','red'), choices = c(3,4), cex = 0.65)
biplot(m.ppca, col = c('grey','red'), choices = c(2,4), cex = 0.65)
biplot(m.ppca, col = c('grey','red'), choices = c(4,5), cex = 0.65)
barplot(m.ppca.perVar, xlab = 'Variable', 
        ylab = 'Percent Variance Explained')
dev.off()

pdf('./Morphospace/Female_pPCA_Loadings.pdf')
par(mfrow = c(2, 3))
biplot(f.ppca, col = c('white','red'), choices = c(1,2), cex = 0.65)
biplot(f.ppca, col = c('white','red'), choices = c(2,3), cex = 0.65)
biplot(f.ppca, col = c('white','red'), choices = c(3,4), cex = 0.65)
biplot(f.ppca, col = c('white','red'), choices = c(2,4), cex = 0.65)
biplot(f.ppca, col = c('white','red'), choices = c(4,5), cex = 0.65)
barplot(f.ppca.perVar, xlab = 'Variable', 
        ylab = 'Percent Variance Explained')
dev.off()

