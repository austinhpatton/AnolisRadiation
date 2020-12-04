# Individual trait boxplots for M1/M2/I2, and M1/M2-Allo/M2-Symp
setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/')

library(ggpubr)


m.dat <- read.csv('./FinalData/Male-ALL-DATA-Trait-Ecol-Regions-FINAL.csv')
m.dat <- m.dat[-which(is.na(m.dat$Group)),]
m.dat <- m.dat[-which(m.dat$Species %in% c('agassizi', 'bicaorum', 'lineatus', 'roatanensis', 
                                           'concolor', 'pinchoti', 'utilensis', 'townsendi')),]

m.dat <- m.dat[,c(1,3:10,12:13,42,44,23,41)]
m.dat <- m.dat[which(m.dat$Group != 'I1'),]
colnames(m.dat)[12:13] <- c('PH', 'PD')

# We need to rename the greater antilles radiation from its old name (I2) to its new name (I1)
m.dat$Group <- as.factor(gsub(m.dat$Group, pattern = 'I2', replacement = 'I1'))

comparisons <- list(c('I1', 'M2'),
                    c('M1', 'M2'),
                    c('I1', 'M1'))

colors <- c("#d72027", "#add8a4", "#2b83bb")

# Okay, now we're just going to plot trait by trait
ph.m <- 
  ggboxplot(m.dat, 
            x = "Group", y = "PH", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 4,
                     label.y = c(max(na.omit(m.dat$PH)*0.99), 
                                 max(na.omit(m.dat$PH))*0.95, 
                                 max(na.omit(m.dat$PH))*0.90)) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab('Ln(Perch Height)') 

t.test(
  m.eco[which(m.eco$Group == 'M1'),2],
  m.eco[which(m.eco$Group == 'M2'),2]
)
cohensD( m.eco[which(m.eco$Group == 'M1'),2],
         m.eco[which(m.eco$Group == 'M2'),2], 
         method = 'unequal')

pd.m <- 
  ggboxplot(m.dat, 
            x = "Group", y = "PD", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 4,
                     label.y = c(max(na.omit(m.dat$PD)*1.05), 
                                 max(na.omit(m.dat$PD))*0.99, 
                                 max(na.omit(m.dat$PD))*0.93)) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab('Ln(Perch Diameter)')

morph.ecol.legend <- 
  ggplot(data = m.dat, aes(y = PD, x = Group, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

morph.ecol.legend <- get_legend(morph.ecol.legend)



# Now plot just the ecological variables
library(cowplot)
pdf('Habitat-Use-Cladewise.pdf', width = 9, height = 9)
eco.broad.p <- 
  plot_grid(ph.m, pd.m, ncol = 2, align = 'v', 
            labels = c('A.', 'B.', 'C.', 'D.', morph.ecol.legend), label_size = 16)
plot_grid(eco.broad.p, morph.ecol.legend, ncol = 1, rel_heights = c(1, .05))
dev.off()

# Okay, and do the same for the comparison between M1 and sympatric/allopatric M2
comparisons <- list(c('M2-Allopatric', 'M2-Sympatric'),
                 c('M1', 'M2-Sympatric'),
                 c('M1', 'M2-Allopatric'))

main.cols <- c("#add8a4", "#b2182b", "#ffc866")
main.dat <- m.dat[which(m.dat$MainlandPatry != 'I2'),]
main.dat$MainlandPatry <- droplevels(main.dat$MainlandPatry)

ph.m <- 
  ggboxplot(main.dat, 
            x = "MainlandPatry", y = "PH", fill = "MainlandPatry", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = MainlandPatry), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 4,
                     label.y = c(max(na.omit(m.dat$PH)*0.99), 
                                 max(na.omit(m.dat$PH))*0.95, 
                                 max(na.omit(m.dat$PH))*0.90)) +
  scale_fill_manual(values=main.cols) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab('Ln(Perch Height)') 


pd.m <- 
  ggboxplot(main.dat, 
            x = "MainlandPatry", y = "PD", fill = "MainlandPatry", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = MainlandPatry), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 4,
                     label.y = c(max(na.omit(m.dat$PD)*1.13), 
                                 max(na.omit(m.dat$PD))*1.05, 
                                 max(na.omit(m.dat$PD))*0.98)) +
  scale_fill_manual(values=main.cols) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab('Ln(Perch Diameter)')


morph.ecol.legend <- 
  ggplot(data = main.dat, aes(y = PD, x = MainlandPatry, fill = MainlandPatry)) +
  geom_boxplot() +
  scale_fill_manual(values=main.cols) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

morph.ecol.legend <- get_legend(morph.ecol.legend)

pdf('Habitat-Use-Mainland-Males.pdf', width = 9, height = 9)
eco.main.p <- 
  plot_grid(ph.m, pd.m, ncol = 2, align = 'v', 
            labels = c('A.', 'B.', 'C.', 'D.', morph.ecol.legend), label_size = 16)
plot_grid(eco.main.p, morph.ecol.legend, ncol = 1, rel_heights = c(1, .05))
dev.off()

################################################
# Now traits

# Figure five comparison of significantly different traits between M1/M2 when including I2
setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/')

################ Simple comparison among groups ###################
library(MASS)
library(reshape)


m.data <- read.csv('./FinalData/Male-ALL-DATA-Trait-Ecol-Regions-FINAL.csv')
m.data <- m.data[-which(is.na(m.data$Group)),]
m.data <- m.data[-which(m.data$Species %in% c('agassizi', 'bicaorum', 'lineatus', 'roatanensis', 
                                              'concolor', 'pinchoti', 'utilensis', 'townsendi')),]

m.traits <- m.data[,-c(2,11,14:22,24:45)]
rownames(m.traits) <- m.traits$Species
m.traits <- m.traits[,-1]
m.traits <- m.traits[which(m.traits$Group != 'I1'),]

# m.group <- as.factor(m.full.data$Lineage)
# names(m.group) <- m.full.data$Species

# m.traits$Lineage <- m.group
m.traits$Group <- factor(m.traits$Group,  
                         levels = c('I2','M1', 'M2'))

m.traits.svl <- m.traits[,c(1,11)]
m.traits.NoSVL <- m.traits[,-1]

box.leg <- ggplot(melt(m.traits)) +
  geom_boxplot(aes(x = variable, y = value, fill = Group))
box.leg <- get_legend(box.leg)
box.leg <- as_ggplot(box.leg)


comparisons <- list(c('M1', 'M2'),
                    c('I2', 'M2'),
                    c('I2', 'M1'))

colors <- c('#d72027', '#add8a4', '#2b83bb')

svl.p <- 
  ggboxplot(m.traits[,c(1,11)], 
            x = "Group", y = "SVL", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative SVL')

flimb.len.p <- 
  ggboxplot(m.traits[,c(8,11)], 
            x = "Group", y = "Forelimb.length", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Fore Limb Length')

h.width.p <- 
  ggboxplot(m.traits[,c(4,11)], 
            x = "Group", y = "Head.width", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Head Width')

h.depth.p <- 
  ggboxplot(m.traits[,c(5,11)], 
            x = "Group", y = "Head.depth", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Head Depth')

hlimb.len.p <- 
  ggboxplot(m.traits[,c(6,11)], 
            x = "Group", y = "Hind.limb.length", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Hind Limb Length')


m.eco <- m.data[,c(23,42:45)]
m.eco <- m.eco[which(m.eco$Group != 'I1'),]
ph.p <- 
  ggboxplot(m.eco[,1:2], 
            x = "Group", y = "PH.Ln.Unscaled", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6,
                     label.y = c(max(na.omit(m.dat$PH)*0.99), 
                                 max(na.omit(m.dat$PH))*0.95, 
                                 max(na.omit(m.dat$PH))*0.90)) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Ln(Perch Height)')
t.test(
  m.eco[,c(1:2)][which(m.eco[,c(1:2)]$Group == 'M1'),2],
  m.eco[,c(1:2)][which(m.eco[,c(1:2)]$Group == 'M2'),2]
)
cohensD( m.eco[,c(1:2)][which(m.eco[,c(1:2)]$Group == 'M1'),2],
         m.eco[,c(1:2)][which(m.eco[,c(1:2)]$Group == 'M2'),2], 
         method = 'unequal')

legend <- 
  ggboxplot(m.eco[,1:2], 
            x = "Group", y = "PH.Ln.Unscaled", fill = "Group", 
            outlier.shape = NA) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(legend.position = "bottom") +
  ylab('Ln(Perch Height)')

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  legend + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# Now plot them all together!
library(cowplot)

plots.all <- 
  plot_grid(svl.p, flimb.len.p, hlimb.len.p, h.width.p, h.depth.p, ph.p, 
            ncol = 3, align = 'v', 
            labels = c('A.', 'B.', 'C.', 'D.', 'E.', 'F.', legend), label_size = 28)

pdf('SignifDiff-Traits-PH-AmongM1M2.pdf', width = 14.5, height = 10)
plot_grid(plots.all, legend, ncol = 1, rel_heights = c(1, .05))
dev.off()

# And just make the rest of them for completeness
tail.p <- 
  ggboxplot(m.traits[,c(2,11)], 
            x = "Group", y = "Tail", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Tail Length')

h.length.p <- 
  ggboxplot(m.traits[,c(3,11)], 
            x = "Group", y = "Head.length", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Head Length')

longtoe.p <- 
  ggboxplot(m.traits[,c(7,11)], 
            x = "Group", y = "Longest.toe", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Length of Longest Toe')

flam.p <- 
  ggboxplot(m.traits[,c(9,11)], 
            x = "Group", y = "ForeLam", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Forefoot Lamella Count')

hlam.p <- 
  ggboxplot(m.traits[,c(10,11)], 
            x = "Group", y = "HindLam", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Hindfoot Lamella Count')



plots.rest <- 
  plot_grid(tail.p, h.length.p, longtoe.p, flam.p, hlam.p, legend, ncol = 3, align = 'v', 
            labels = c('A.', 'B.', 'C.', 'D.', 'E.'), label_size = 28)

pdf('RemainingTraitBoxplots-CladeWise.pdf', width = 14.5, height = 12)
plot_grid(plots.rest, rel_widths = c(3, .3))
dev.off()

###################################################################################################
# Now between the mainland groups                                                                 #
###################################################################################################
m.traits <- m.data[,-c(2,11,14:40, 42:45)]
colnames(m.traits)[12] <- 'Group'
rownames(m.traits) <- m.traits$Species
m.traits <- m.traits[,-1]
m.traits <- m.traits[which(m.traits$Group != 'I1'),]


box.leg <- ggplot(melt(m.traits)) +
  geom_boxplot(aes(x = variable, y = value, fill = Group))
box.leg <- get_legend(box.leg)
box.leg <- as_ggplot(box.leg)

comparisons <- list(c('M2-Allopatric', 'M2-Sympatric'),
                    c('M1', 'M2-Sympatric'),
                    c('M1', 'M2-Allopatric'))
m.traits <- m.traits[-which(m.traits$Group %in% c('I1', 'I2')),]
colors <- c("#add8a4", "#b2182b", "#ffc866")

svl.p <- 
  ggboxplot(m.traits[,c(1,11)], 
            x = "Group", y = "SVL", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative SVL')

tail.p <- 
  ggboxplot(m.traits[,c(2,11)], 
            x = "Group", y = "Tail", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Tail Length')

h.width.p <- 
  ggboxplot(m.traits[,c(4,11)], 
            x = "Group", y = "Head.width", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Head Width')

h.depth.p <- 
  ggboxplot(m.traits[,c(5,11)], 
            x = "Group", y = "Head.depth", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Head Depth')

hlimb.len.p <- 
  ggboxplot(m.traits[,c(6,11)], 
            x = "Group", y = "Hind.limb.length", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Hind Limb Length')

flimb.len.p <- 
  ggboxplot(m.traits[,c(8,11)], 
            x = "Group", y = "Forelimb.length", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Fore Limb Length')

h.length.p <- 
  ggboxplot(m.traits[,c(3,11)], 
            x = "Group", y = "Head.length", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Head Length')

longtoe.p <- 
  ggboxplot(m.traits[,c(7,11)], 
            x = "Group", y = "Longest.toe", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Length of Longest Toe')

flam.p <- 
  ggboxplot(m.traits[,c(9,11)], 
            x = "Group", y = "ForeLam", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Forefoot Lamella Count')

hlam.p <- 
  ggboxplot(m.traits[,c(10,11)], 
            x = "Group", y = "HindLam", fill = "Group", 
            outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1), 
             pch = 21, aes(fill = Group), size = 3, stroke = 1) +
  stat_compare_means(comparisons = comparisons,
                     method = 't.test',
                     size = 6) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Relative Hindfoot Lamella Count')

legend <- 
  ggboxplot(m.traits[,c(1,11)], 
            x = "Group", y = "SVL", fill = "Group", 
            outlier.shape = NA) +
  scale_fill_manual(values=colors) +
  theme_bw(base_size = 24) +
  theme(legend.position = "bottom") +
  ylab('Ln(Perch Height)')

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  legend + theme(legend.box.margin = margin(0, 0, 0, 12))
)


main.plots <- 
  plot_grid(svl.p, flimb.len.p, hlimb.len.p, h.width.p, h.depth.p, 
            tail.p, h.length.p, longtoe.p, flam.p, hlam.p, ncol = 5, align = 'V', 
            labels = c('A.', 'B.', 'C.', 'D.', 'E.', 
                       'F.', 'G.', 'H.', 'I.', 'J.'), label_size = 28)

pdf('Morph-Traits-MainlandGroups.pdf', width = 20, height = 12)
plot_grid(main.plots, legend, ncol = 1, rel_heights = c(1, .05))
dev.off()


