# Boxplots of traits and ecological variables for M1 and M2, broken down by broad geographic area. 
########### See if M1 in expansion range is different than M1 in ancestral #####################
library(treeplyr)
library(purrr)
library(phytools)
library(geiger)
library(gridExtra)
library(ggpubr)
library(reshape)
library(QsRutils)

setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/')

m.data <- read.csv('FinalData/Male-ALL-DATA-Trait-Ecol-Regions-FINAL.csv', 
                   header = T)
m.data <- m.data[-which(m.data$Species %in% c('agassizi', 'bicaorum', 'lineatus', 'roatanensis', 
                                              'concolor', 'pinchoti', 'utilensis', 'townsendi')),]

# FirstReduce to just M1
m1.traits <- m.data[which(m.data$MainlandPatry == 'M1'),c(1,3:10, 12:13,41)]

# Get ecology
m1.eco <- m.data[which(m.data$MainlandPatry == 'M1'),c(1,42:45)]

# get regions 
m1.regs <- m.data[which(m.data$MainlandPatry == 'M1'),c(1,26:41)]

# And sort into expansion (lower CA) and ancestral (SA)
m1.exp.spp <- as.character(m1.regs$Species[which(m1.regs$Lower.Central.America == 1)])
m1.anc.spp <- as.character(m1.regs$Species[which(m1.regs$SouthAmerica == 1)])

# Now reduce data
exp.data <- m1.traits[which(m1.traits$Species %in% m1.exp.spp),-c(12:13)]
exp.eco <- m1.eco[which(m1.eco$Species %in% m1.exp.spp),-c(6:7)]
anc.data <- m1.traits[which(m1.traits$Species %in% m1.anc.spp),-c(12:13)]
anc.eco <- m1.eco[which(m1.eco$Species %in% m1.anc.spp),-c(6:7)]

exp.dat <- merge(exp.data, exp.eco, by = 'Species', all = T)
exp.dat$Group <- 'M1-LowerCA'
anc.dat <- merge(anc.data, anc.eco, by = 'Species', all = T)
anc.dat$Group <- 'M1-SouthAmerica'

m1.dat <- rbind(exp.dat, anc.dat)
m1.dat$Group <- factor(m1.dat$Group, levels = c('M1-SouthAmerica', 'M1-LowerCA'))

comparisons <- list(c("M1-SouthAmerica", "M1-LowerCA"))

m1.dat.svl <- melt(m1.dat[,c(1:2,16)], id.vars = c('Species', 'Group'))
m1.dat.no.svl <- melt(m1.dat[,c(1,3:11,16)], id.vars = c('Species', 'Group'))
m1.dat.eco <- melt(m1.dat[,c(1,12:16)], id.vars = c('Species', 'Group'))


# Now look at M2. We'll combine the above with these later. 
m2.traits <- m.data[which(m.data$Group == 'M2'),c(1,3:10, 12:13,41)]

# Get ecology
m2.eco <- m.data[which(m.data$Group == 'M2'),c(1,42:45)]

# get regions 
m2.regs <- m.data[which(m.data$Group == 'M2'),c(1,26:41)]

m2.uca.spp <- as.character(m2.regs$Species[which(rowSums(m2.regs[,3:4]) != 0)])
m2.lca.spp <- as.character(m2.regs$Species[which(m2.regs$Lower.Central.America == 1)])
m2.sa.spp <- as.character(m2.regs$Species[which(m2.regs$SouthAmerica == 1)])


# Now reduce data
m2.uca.data <- m2.traits[which(m2.traits$Species %in% m2.uca.spp),-c(12:13)]
m2.uca.eco <- m2.eco[which(m2.eco$Species %in% m2.uca.spp),-7]
m2.lca.data <- m2.traits[which(m2.traits$Species %in% m2.lca.spp),-c(12:13)]
m2.lca.eco <- m2.eco[which(m2.eco$Species %in% m2.lca.spp),-7]
m2.sa.data <- m2.traits[which(m2.traits$Species %in% m2.sa.spp),-c(12:13)]
m2.sa.eco <- m2.eco[which(m2.eco$Species %in% m2.sa.spp),-7]

m2.uca.dat <- merge(m2.uca.data, m2.uca.eco, by = 'Species', all = T)
m2.uca.dat$Group <- 'M2-UpperCA-Nearctic'
m2.lca.dat <- merge(m2.lca.data, m2.lca.eco, by = 'Species', all = T)
m2.lca.dat$Group <- 'M2-LowerCA'
m2.sa.dat <- merge(m2.sa.data, m2.sa.eco, by = 'Species', all = T)
m2.sa.dat$Group <- 'M2-SouthAmerica'

m2.dat <- rbind(m2.uca.dat, m2.lca.dat, m2.sa.dat)
m2.dat$Group <- factor(m2.dat$Group, levels = c('M2-SouthAmerica',
                                                'M2-LowerCA',
                                                'M2-UpperCA-Nearctic'))

comparisons <- list(c("M2-UpperCA-Nearctic", "M2-LowerCA"),
                    c("M2-UpperCA-Nearctic", "M2-SouthAmerica"),
                    c("M2-LowerCA", "M2-SouthAmerica"))

m2.dat.svl <- melt(m2.dat[,c(1:2,16)], id.vars = c('Species', 'Group'))
m2.dat.no.svl <- melt(m2.dat[,c(1,3:11,16)], id.vars = c('Species', 'Group'))
m2.dat.eco <- melt(m2.dat[,c(1,12:16)], id.vars = c('Species', 'Group'))

m1m2.dat <- rbind(m1.dat, m2.dat)
m1m2.dat.svl <- rbind(m1.dat.svl, m2.dat.svl)
m1m2.dat.no.svl <- rbind(m1.dat.no.svl, m2.dat.no.svl)
m1m2.dat.eco <- rbind(m1.dat.eco, m2.dat.eco)

comparisons <- list(c("M1-SouthAmerica", "M2-UpperCA-Nearctic"),
                    c("M2-UpperCA-Nearctic","M2-LowerCA"),
                    c("M1-SouthAmerica", "M2-LowerCA"),
                    c("M2-UpperCA-Nearctic","M2-SouthAmerica"),
                    c("M1-SouthAmerica", "M2-SouthAmerica"),
                    c("M2-UpperCA-Nearctic","M1-LowerCA"),
                    c("M1-SouthAmerica", "M1-LowerCA"),
                    c("M2-SouthAmerica", "M2-LowerCA"),
                    c("M1-LowerCA", "M2-SouthAmerica"),
                    c("M1-LowerCA", "M2-LowerCA"))

pvals <- pairwise.t.test(x = m1m2.dat$SVL, m1m2.dat$Group, p.adjust.method = 'fdr', pool.sd = F)
pval.table <- as.data.frame(pvals$p.value)

ltrs <- make_letter_assignments(pvals)
svl.p <- 
  ggplot(m1m2.dat, aes(x = Group, y = SVL, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(), pch = 21, size = 3, stroke = 1, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Ln(SVL)')

# simple test of ecological filtering - are the lower CA M1 species larger on 
# average than would be expected from a random sample?
obs <- mean(na.omit(m1m2.dat$SVL[which(m1m2.dat$Group == 'M1-LowerCA')]))
m1.sa.dat <- c(na.omit(m1m2.dat[which(m1m2.dat$Group == 'M1-SouthAmerica'), 'SVL']))
exp <- c()
for(i in 1:1000){
  exp[i] <- mean(sample(x = m1.sa.dat, size = 4, replace = F))
}
p <- 1-(sum(obs > exp)/1000)
exp <- data.frame('Expected' = exp)
probs <- c(0, 0.025, 0.5, 0.975, 1)
quantiles <- quantile(exp$Expected, prob=probs)
dens <- density(exp$Expected)
dens <- data.frame('Expected' = dens$x, 
                   'Density' = dens$y)
dens$quant <- factor(findInterval(dens$Expected, quantiles))

ecol.filt.p <- 
  ggplot(data = dens, aes(x = Expected, y = Density, fill = quant)) + 
  geom_ribbon(aes(ymin=0, ymax=Density, fill=quant), 
              outline.type = 'lower', color = 'black') + 
  scale_fill_brewer(guide="none") +
  geom_vline(xintercept = obs, color = 'red', lty = 2) +
  annotate(geom = 'text', label = paste0('p = ', round(p, 3)), 
           y = max(dens$Density), x = min(dens$Expected), hjust = 0,
           size = 6) + 
  theme_bw(base_size = 16) + 
  xlab('Expected SVL')
  
ggsave(ecol.filt.p, filename = 'ObsVsExp-LowerCA-M1-SVL.pdf')

x <- rownames(ltrs$LetterMatrix)
y <- graphics::boxplot(data = m1m2.dat, SVL ~ Group, plot = FALSE)$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame('Group' = x, y, 'Letters' = ltrs$Letters)

lmts <- get_plot_limits(svl.p)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
svl.p <- 
  svl.p + geom_text(data = ltr_df, aes(x=x, y=y, label=Letters), nudge_y = y.nudge, size = 9)

#################
# Forelimb length
pvals <- pairwise.t.test(x = m1m2.dat$Forelimb.length, m1m2.dat$Group, p.adjust.method = 'fdr', pool.sd = F)
pval.table <- rbind(pval.table, pvals$p.value)

ltrs <- make_letter_assignments(pvals)
flimb.p <- 
  ggplot(m1m2.dat, aes(x = Group, y = Forelimb.length, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(), pch = 21, size = 3, stroke = 1, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Rel. Forelimb Length')


x <- rownames(ltrs$LetterMatrix)
y <- graphics::boxplot(data = m1m2.dat, Forelimb.length ~ Group, plot = FALSE)$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame('Group' = x, y, 'Letters' = ltrs$Letters)

lmts <- get_plot_limits(flimb.p)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
flimb.p <- 
  flimb.p + geom_text(data = ltr_df, aes(x=x, y=y, label=Letters), nudge_y = y.nudge, size = 9)

#################
# Hindlimb length
pvals <- pairwise.t.test(x = m1m2.dat$Hind.limb.length, m1m2.dat$Group, p.adjust.method = 'fdr', pool.sd = F)
pval.table <- rbind(pval.table, pvals$p.value)

ltrs <- make_letter_assignments(pvals)
hlimb.p <- 
  ggplot(m1m2.dat, aes(x = Group, y = Hind.limb.length, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(), pch = 21, size = 3, stroke = 1, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Rel. Hindlimb Length')


x <- rownames(ltrs$LetterMatrix)
y <- graphics::boxplot(data = m1m2.dat, Hind.limb.length ~ Group, plot = FALSE)$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame('Group' = x, y, 'Letters' = ltrs$Letters)

lmts <- get_plot_limits(hlimb.p)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
hlimb.p <- 
  hlimb.p + geom_text(data = ltr_df, aes(x=x, y=y, label=Letters), nudge_y = y.nudge, size = 9)

#################
# Head Width
pvals <- pairwise.t.test(x = m1m2.dat$Head.width, m1m2.dat$Group, p.adjust.method = 'fdr', pool.sd = F)
pval.table <- rbind(pval.table, pvals$p.value)

ltrs <- make_letter_assignments(pvals)
hwidth.p <- 
  ggplot(m1m2.dat, aes(x = Group, y = Head.width, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(), pch = 21, size = 3, stroke = 1, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Rel. Head Width')


x <- rownames(ltrs$LetterMatrix)
y <- graphics::boxplot(data = m1m2.dat, Head.width ~ Group, plot = FALSE)$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame('Group' = x, y, 'Letters' = ltrs$Letters)

lmts <- get_plot_limits(hwidth.p)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
hwidth.p <- 
  hwidth.p + geom_text(data = ltr_df, aes(x=x, y=y, label=Letters), nudge_y = y.nudge, size = 9)

#################
# Head Depth
pvals <- pairwise.t.test(x = m1m2.dat$Head.depth, m1m2.dat$Group, p.adjust.method = 'fdr', pool.sd = F)
pval.table <- rbind(pval.table, pvals$p.value)

ltrs <- make_letter_assignments(pvals)
hdepth.p <- 
  ggplot(m1m2.dat, aes(x = Group, y = Head.depth, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(), pch = 21, size = 3, stroke = 1, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Rel. Head Depth')


x <- rownames(ltrs$LetterMatrix)
y <- graphics::boxplot(data = m1m2.dat, Head.depth ~ Group, plot = FALSE)$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame('Group' = x, y, 'Letters' = ltrs$Letters)

lmts <- get_plot_limits(hdepth.p)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
hdepth.p <- 
  hdepth.p + geom_text(data = ltr_df, aes(x=x, y=y, label=Letters), nudge_y = y.nudge, size = 9)

#################
# Tail Length
pvals <- pairwise.t.test(x = m1m2.dat$Tail, m1m2.dat$Group, p.adjust.method = 'fdr', pool.sd = F)
pval.table <- rbind(pval.table, pvals$p.value)

ltrs <- make_letter_assignments(pvals)
tail.p <- 
  ggplot(m1m2.dat, aes(x = Group, y = Tail, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(), pch = 21, size = 3, stroke = 1, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Rel. Tail Length')


x <- rownames(ltrs$LetterMatrix)
y <- graphics::boxplot(data = m1m2.dat, Tail ~ Group, plot = FALSE)$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame('Group' = x, y, 'Letters' = ltrs$Letters)

lmts <- get_plot_limits(tail.p)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
tail.p <- 
  tail.p + geom_text(data = ltr_df, aes(x=x, y=y, label=Letters), nudge_y = y.nudge, size = 9)

#################
# Head Length
pvals <- pairwise.t.test(x = m1m2.dat$Head.length, m1m2.dat$Group, p.adjust.method = 'fdr', pool.sd = F)
pval.table <- rbind(pval.table, pvals$p.value)

ltrs <- make_letter_assignments(pvals)
hlength.p <- 
  ggplot(m1m2.dat, aes(x = Group, y = Head.length, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(), pch = 21, size = 3, stroke = 1, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Rel. Head Length')


x <- rownames(ltrs$LetterMatrix)
y <- graphics::boxplot(data = m1m2.dat, Head.length ~ Group, plot = FALSE)$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame('Group' = x, y, 'Letters' = ltrs$Letters)

lmts <- get_plot_limits(hlength.p)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
hlength.p <- 
  hlength.p + geom_text(data = ltr_df, aes(x=x, y=y, label=Letters), nudge_y = y.nudge, size = 9)

#################
# Longest Toe
pvals <- pairwise.t.test(x = m1m2.dat$Longest.toe, m1m2.dat$Group, p.adjust.method = 'fdr', pool.sd = F)
pval.table <- rbind(pval.table, pvals$p.value)

ltrs <- make_letter_assignments(pvals)
toe.p <- 
  ggplot(m1m2.dat, aes(x = Group, y = Longest.toe, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(), pch = 21, size = 3, stroke = 1, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Rel. Length Longest Toe')


x <- rownames(ltrs$LetterMatrix)
y <- graphics::boxplot(data = m1m2.dat, Longest.toe ~ Group, plot = FALSE)$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame('Group' = x, y, 'Letters' = ltrs$Letters)

lmts <- get_plot_limits(toe.p)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
toe.p <- 
  toe.p + geom_text(data = ltr_df, aes(x=x, y=y, label=Letters), nudge_y = y.nudge, size = 9)

#################
# Forefoot Lamella count
pvals <- pairwise.t.test(x = m1m2.dat$ForeLam, m1m2.dat$Group, p.adjust.method = 'fdr', pool.sd = F)
pval.table <- rbind(pval.table, pvals$p.value)

ltrs <- make_letter_assignments(pvals)
flam.p <- 
  ggplot(m1m2.dat, aes(x = Group, y = ForeLam, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(), pch = 21, size = 3, stroke = 1, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Rel. Forefoot Lamella Count')


x <- rownames(ltrs$LetterMatrix)
y <- graphics::boxplot(data = m1m2.dat, ForeLam ~ Group, plot = FALSE)$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame('Group' = x, y, 'Letters' = ltrs$Letters)

lmts <- get_plot_limits(flam.p)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
flam.p <- 
  flam.p + geom_text(data = ltr_df, aes(x=x, y=y, label=Letters), nudge_y = y.nudge, size = 9)

#################
# Hindfoot Lamella count
pvals <- pairwise.t.test(x = m1m2.dat$HindLam, m1m2.dat$Group, p.adjust.method = 'fdr', pool.sd = F)
pval.table <- rbind(pval.table, pvals$p.value)

ltrs <- make_letter_assignments(pvals)
hlam.p <- 
  ggplot(m1m2.dat, aes(x = Group, y = HindLam, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(), pch = 21, size = 3, stroke = 1, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Rel. Hindfoot Lamella Count')


x <- rownames(ltrs$LetterMatrix)
y <- graphics::boxplot(data = m1m2.dat, HindLam ~ Group, plot = FALSE)$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame('Group' = x, y, 'Letters' = ltrs$Letters)

lmts <- get_plot_limits(hlam.p)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
hlam.p <- 
  hlam.p + geom_text(data = ltr_df, aes(x=x, y=y, label=Letters), nudge_y = y.nudge, size = 9)

pval.table <- 
  cbind(pval.table,
        'Trait' = c(rep('SVL', 4), rep('Forelimb', 4), rep('Hindlimb', 4), 
                   rep('HeadWidth', 4), rep('HeadDepth', 4), rep('TailLength', 4), 
                   rep('HeadLength', 4), rep('LongestToe', 4), rep('ForeLam', 4), 
                   rep('HindLam', 4)))
write.table(pval.table, 'Mainland-Traits-ByRegion-M1-M2.tsv', sep = '\t', row.names = T, col.names = T, quote = F)

# Get a legend
legend <- 
  ggplot(m1m2.dat, aes(x = Group, y = HindLam, fill = Group)) +
  geom_boxplot() +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
legend <- get_legend(legend)

library(cowplot)

trait.main.p <- 
  plot_grid(svl.p, flimb.p, hlimb.p, hwidth.p, hdepth.p,
            tail.p, hlength.p, toe.p, flam.p, hlam.p, ncol = 5, align = 'v', 
            labels = c('A.', 'B.', 'C.', 'D.', 'E.', 
                       'F.', 'G.', 'H.', 'I.', 'J.'), label_size = 16)
pdf('Traits-Mainland-ByRegion-Males.pdf', width = 18, height = 10)
plot_grid(trait.main.p, legend, ncol = 1, rel_heights = c(1, .05))
dev.off()


############################################################################
# Now the ecological variables

#################
# Perch height
pvals <- pairwise.t.test(x = m1m2.dat$PH.Ln.Unscaled, m1m2.dat$Group, p.adjust.method = 'fdr', pool.sd = F)
pval.table <- as.data.frame(pvals$p.value)

ltrs <- make_letter_assignments(pvals)
ph.p <- 
  ggplot(m1m2.dat, aes(x = Group, y = PH.Ln.Unscaled, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(), pch = 21, size = 3, stroke = 1, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Ln(Perch Height)')


x <- rownames(ltrs$LetterMatrix)
y <- graphics::boxplot(data = m1m2.dat, PH.Ln.Unscaled ~ Group, plot = FALSE)$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame('Group' = x, y, 'Letters' = ltrs$Letters)

lmts <- get_plot_limits(ph.p)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
ph.p <- 
  ph.p + geom_text(data = ltr_df, aes(x=x, y=y, label=Letters), nudge_y = y.nudge, size = 9)

#################
# Perch Diameter
pvals <- pairwise.t.test(x = m1m2.dat$PD.Ln.Unscaled, m1m2.dat$Group, p.adjust.method = 'fdr', pool.sd = F)
pval.table <- rbind(pval.table, pvals$p.value)

ltrs <- make_letter_assignments(pvals)
pd.p <- 
  ggplot(m1m2.dat, aes(x = Group, y = PD.Ln.Unscaled, fill = Group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(), pch = 21, size = 3, stroke = 1, alpha = 0.5) +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  ylab('Ln(Perch Diameter)')


x <- rownames(ltrs$LetterMatrix)
y <- graphics::boxplot(data = m1m2.dat, PD.Ln.Unscaled ~ Group, plot = FALSE)$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame('Group' = x, y, 'Letters' = ltrs$Letters)

lmts <- get_plot_limits(pd.p)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
pd.p <- 
  pd.p + geom_text(data = ltr_df, aes(x=x, y=y, label=Letters), nudge_y = y.nudge, size = 9)

pval.table <- 
  cbind(pval.table,
        'Trait' = c(rep('PH', 4), rep('PD', 4)))
write.table(pval.table, 'Mainland-Ecol-ByRegion-M1-M2.tsv', sep = '\t', row.names = T, col.names = T, quote = F)

# Get a legend
legend <- 
  ggplot(m1m2.dat, aes(x = Group, y = PD.Ln.Unscaled, fill = Group)) +
  geom_boxplot() +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
legend <- get_legend(legend)


eco.main.p <- 
  plot_grid(ph.p, pd.p, ncol = 2, align = 'v', 
            labels = c('A.', 'B.', label_size = 16))

pdf('Habitat-Use-Mainland-ByRegion-Males.pdf', width = 10, height = 7)
plot_grid(eco.main.p, legend, ncol = 1, rel_heights = c(1, .05))
dev.off()

