rm(list = ls(all.names = TRUE))
setwd('~/Library/CloudStorage/Box-Box/Huang_lab/manuscripts/pQTL/analysis')
library('plyr')
library('dplyr')
library('readxl')
library('ggrepel')
library('ggplot2')
library('formattable')
library('ComplexHeatmap')

color.palette =  setNames(c("#8DD3C7","#FB8072","#BEBADA"),
                          c('missense', 'truncating', 'synonymous'))
color.palette1 =  setNames(c("#FAD2D9", "#9EDDF9"),
                         c('mutated', 'nonmutated'))

pan.coh = c('BRCA', 'CRC', 'CCRCC', 'LUAD', 'OV', 'UCEC')
mutations = c('missense', 'truncating', 'synonymous') 

################################################################################
# stacked bar plots (Fig 2A)
################################################################################
rna = read.csv('../data/DNA_RNA_regression/Table.DNA.RNA.regression.linearLIMMA.RNAVsMut.csv')
pro = read.csv('../data/DNA_Pro_regression/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.csv')

dat = rna
# dat = pro
p <- ggplot(dat, aes(x=cancer, fill=mutation)) 
p <- p + ggtitle('mRNA') + theme_bw(base_size=12)
# p <- p + ggtitle('Protein') + theme_bw(base_size=12)
p <- p + geom_bar(position='stack', alpha = .7, stat = 'count', show.legend = FALSE)
p <- p + scale_fill_manual(values = color.palette)
p <- p + theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.text.x = element_text(colour="black",angle = 45, vjust=0.5))
print(p)
ggsave(file=paste0('../doc/RNA_barplot.png'), w=2.5, h=2.5, unit='in')
# ggsave(file=paste0('../doc/Pro_barplot.png'), w=2.5, h=2.5, unit='in')

dat = rna
dat = dat[dat$FDR < 0.05, ] # eQTLs
eQTL_counts = dat %>% group_by(cancer, mutation) %>% dplyr ::summarise(count_muts = n())
rna_counts = rna %>% group_by(cancer, mutation) %>% dplyr ::summarise(count_muts = n())
counts = merge(eQTL_counts, rna_counts, by=c('cancer', 'mutation'))
eQTL_counts$percentage = counts$count_muts.x / counts$count_muts.y
p <- ggplot(eQTL_counts, aes(x=cancer, y=percentage, fill=mutation)) 
p <- p + ggtitle('eQTL') + theme_bw(base_size=12)
p <- p + theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_bar(position='fill', alpha = .7, stat = 'identity', show.legend = FALSE)
p <- p + scale_fill_manual(values = color.palette)
p <- p + theme(axis.text.x = element_text(colour="black",angle = 45, vjust=0.5))
print(p)
ggsave(file=paste0('../doc/eQTLs_barplot.png'), w=2.5, h=2.5, unit='in')

dat = pro
dat = dat[dat$FDR < 0.05, ] # pQTLs
pQTL_counts = dat %>% group_by(cancer, mutation) %>% dplyr ::summarise(count_muts = n())
pro_counts = pro %>% group_by(cancer, mutation) %>% dplyr ::summarise(count_muts = n())
counts = merge(pQTL_counts, pro_counts, by=c('cancer', 'mutation'))
pQTL_counts$percentage = counts$count_muts.x / counts$count_muts.y
p <- ggplot(pQTL_counts, aes(x=cancer, y=percentage, fill=mutation)) 
p <- p + ggtitle('pQTL') + theme_bw(base_size=12)
p <- p + theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_bar(position='fill', alpha = .7, stat = 'identity', show.legend = FALSE)
p <- p + scale_fill_manual(values = color.palette)
p <- p + theme(axis.text.x = element_text(colour="black",angle = 45, vjust=0.5))
print(p)
ggsave(file=paste0('../doc/pQTLs_barplot.png'), w=2.5, h=2.5, unit='in')

################################################################################
# volcano plots (Fig 2B)
################################################################################
# res <- read.csv('../data/DNA_RNA_regression/Table.DNA.RNA.regression.linearLIMMA.RNAVsMut.csv')
res <- read.csv('../data/DNA_Pro_regression/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.csv')

top = 5
thr.lfc = 1
thr.FDR = 0.05

cancer = 'CCRCC'
for (cancer in pan.coh) {
  res_can = res[res$cancer == cancer,]
  gg.text = res_can[(abs(res_can$logFC) > thr.lfc) & (res_can$FDR < thr.FDR), ]
  gg.text = head(unique(gg.text[order(-log(gg.text$FDR), decreasing = T),]), top)
  p <- ggplot(res_can) 
  p <- p + geom_vline(xintercept = -thr.lfc, alpha = .3) + 
    geom_vline(xintercept = thr.lfc, alpha = .3)
  p <- p + geom_hline(yintercept = -log10(thr.FDR), alpha = .3)
  p <- p + ggtitle(cancer) + theme_bw()
  p <- p + theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p <- p + geom_point(res_can, mapping = aes(logFC, -log10(FDR), col=mutation), alpha = .7, show.legend = FALSE)
  p <- p + geom_label_repel(data=gg.text, mapping = aes(logFC, -log10(FDR), label=Gene, col=mutation), 
                            min.segment.length = 0, alpha = .9, max.overlaps=50, size=4, show.legend = FALSE)
  p <- p + scale_color_manual(values = color.palette)
  p <- p + xlab('FC') + ylab('-Log10(FDR)')
  p <- p + theme(axis.text = element_text(colour="black", vjust=0.5, size=11))
  # ggsave(file=paste0('../doc/volcano\ plot/', cancer, 'eQTL.png'), w=2.5, h=2.5)
  ggsave(file=paste0('../doc/volcano\ plot/', cancer, 'pQTL.png'), w=2.5, h=2.5)
  print(p)
  dev.off()
}

################################################################################
# scatter plots (Fig 3A and Fig 4A)
################################################################################
dat <- read_excel('../data/results/concordant.xlsx')
p <- ggplot(dat) 
p <- p + ggtitle('concordant QTLs') + theme_bw(base_size=12)
p <- p + geom_point(dat, mapping = aes(RNAlogFC, PrologFC, col=mutation), alpha = .7)
p <- p + geom_label_repel(data=dat, mapping = aes(RNAlogFC, PrologFC, label=Gene, col=mutation),
                          min.segment.length = 0, alpha = .9, max.overlaps=50, size=4, show.legend = FALSE)
p <- p + scale_color_manual(values = color.palette)
p <- p + xlab('RNAlogFC') + ylab('PrologFC')
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_abline(intercept = 0, slope = 1, color='grey', size=0.8)
print(p)
ggsave(file='../doc/concordantQTLs.png', w=6, h=4.5)

dat <- read_excel('../data/results/discordant.xlsx')
dat <- dat[dat$overlap == TRUE,]
p <- ggplot(dat) 
p <- p + ggtitle('discordant QTLs') + theme_bw(base_size=12)
p <- p + geom_point(dat, mapping = aes(RNAlogFC, PrologFC, col=mutation), alpha = .7)
p <- p + scale_color_manual(values = color.palette)
p <- p + xlab('RNAlogFC') + ylab('PrologFC')
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_abline(intercept = 0, slope = 1, color='grey', size=0.8)
ggsave(file='../doc/discordantQTLs.png', w=4, h=3)

################################################################################
# violin plots (Fig 3B and Fig 4B)
################################################################################
data.level = 'Pro'
data.level = 'RNA'

mutation = 'truncating'
for (mutation in mutations) {
  # exp <- read.csv(paste0('../data/wildVSmut/concordant/', data.level, '_exp_wildVsMut.csv'))
  # dir <- read_excel('../data/results/concordant.xlsx')
  exp <- read.csv(paste0('../data/wildVSmut/discordant/', data.level, '_exp_wildVsMut.csv'))
  dir <- read_excel('../data/results/discordant.xlsx')
  
  exp <- exp[exp$mutation==mutation,]
  exp[exp$isMut == 'mutant', ]$isMut = 'mutated'
  exp[exp$isMut == 'wild', ]$isMut = 'nonmutated'
  genes = unique(exp$gene)
  gene = 'ATM'
  for (gene in genes) {
    if (sum(exp[exp$gene==gene,]$isMut == 'mutated') < 3) {
      next
    } else {
      p <- ggplot(exp[exp$gene==gene,], aes(x=isMut, y=expression, fill=isMut, color=isMut))
      p <- p + facet_grid(.~cancer, scale = "free", space = "free", drop=T)
      p <- p + geom_dotplot(dotsize=3, binwidth=.01, binaxis= "y", color='black', stackdir ="centerwhole",alpha=0.5, show.legend = FALSE)
      p <- p + geom_violin(alpha=0.6, show.legend = FALSE)
      p <- p + theme_bw(base_size=21)
      p <- p + theme(axis.text.x = element_text(colour="black",angle = 45, vjust=0.5, size=22))
      p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      p <- p + scale_color_manual(values = color.palette1) + scale_fill_manual(values = color.palette1)
      
      gg.text = exp[(exp$gene==gene) & (exp$isMut == 'mutated'),]
      if (data.level == 'Pro') {
        if (dir[(dir$Gene==gene) & (dir$mutation == mutation),]$PrologFC > 0) {
          gg.text = head(gg.text[order(gg.text$expression, decreasing = T), ], 3)
        } else {
          gg.text = head(gg.text[order(gg.text$expression, decreasing = F), ], 3)
        }
      } else {
        if (dir[(dir$Gene==gene) & (dir$mutation == mutation),]$RNAlogFC > 0) {
          gg.text = head(gg.text[order(gg.text$expression, decreasing = T), ], 3)
        } else {
          gg.text = head(gg.text[order(gg.text$expression, decreasing = F), ], 3)
        }
      }
      
      p <- p + geom_label_repel(gg.text, mapping = aes(x=isMut, y=expression, label=hgvsp, col='white'),
                                min.segment.length = 0, alpha = .9, max.overlaps=50, size=6, show.legend = FALSE)
      p <- p + xlab('mutation status')
      p <- p + ggtitle(paste(cancer, gene, mutation))
      print(p)
      # ggsave(file=paste0('../doc/MutVsWild/concordant/', mutation, '/', data.level, '/', gene, 'ExpMutVsWild.png'))
      ggsave(file=paste0('../doc/MutVsWild/discordant/', mutation, '/', data.level, '/', gene, 'ExpMutVsWild.png'))
    }
  }
}

################################################################################
# pie chart of psQTLs (Fig 4C)
################################################################################
dat <- read_excel('../data/results/discordant.xlsx')
dat <- dat[dat$overlap == TRUE,]
p <- ggplot(dat, aes(x="", fill=mutation)) +
  geom_bar(stat="count", width=1, alpha = .7) + 
  scale_fill_manual(values = color.palette) + 
  coord_polar("y", start=0) + 
  geom_text(stat='count', aes(label=percent(..count../55)), size=4, position = position_stack(vjust = 0.4)) +
  theme_void(base_size=12)
print(p)
ggsave(file=paste0('../doc/psQTLs_mut_count.png'), w=4, h=3)

################################################################################
# heat map (Fig 4D)
################################################################################
sig_pro <- read_excel('../data/results/discordant.xlsx')
sig_pro_rna = read.csv('../data/sig_pro_rna.csv')
sig_rna <- sig_pro_rna[(sig_pro_rna$ProFDR >= 0.05) & (sig_pro_rna$RNAFDR < 0.05),]

mutation = 'truncating'
# get heatmap data
for (mutation in mutations) {
  spro <- sig_pro[sig_pro$mutation == mutation, ]
  srna <- sig_rna[sig_rna$mutation == mutation, ]
  dat <- read.csv(paste0('../data/heatmapData/heatmapData', mutation, '.csv'))
  rownames(dat) <- dat$X
  dat <- dat[-c(1)]
  p <- Heatmap(as.matrix(dat),  
               column_title = mutation,
               row_split = 1:nrow(dat), row_title = NULL,
               column_split = 1:ncol(dat), 
               border = NA, cluster_rows = F, cluster_columns = F,
               layer_fun = function(j, i, x, y, width, height, fill) {
                 gene = rownames(dat)[i]
                 cancer = colnames(dat)[j]
                 if (dim(spro[(spro$Gene==gene) & (spro$cancer==cancer), ])[1] != 0) {
                   grid.rect(x, y, gp = gpar(lwd = 2.5, fill = "transparent", col = 'brown'))
                 } #else if ((dim(srna[(srna$Gene==gene) & (srna$cancer==cancer), ])[1] != 0)) {
                 #grid.rect(x, y, gp = gpar(lwd = 2.5, fill = "transparent", col = 'black'))
                 #}
               }
  )
  png(paste0('../doc/heatmap/', mutation, 'heatmap.png'), w=4, h=4, unit='in', res=100)
  print(p)
  dev.off()
}

