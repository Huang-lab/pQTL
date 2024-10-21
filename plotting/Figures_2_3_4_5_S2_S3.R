# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#     Figures 2-5     [Figs.R]
#     for Figure S2, please see the end of this script (line# 148+)
#     for Figure S3, please see the end of this script (line# 400+)
#
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|

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
color.palette2 =  setNames(c("#ED2891", "#9EDDF9", "#6E7BA2", "#DAF1FC","#00A99D", "#FBE3C7"),
                           c('BRCA', 'CRC', 'CCRCC', 'LUAD', 'OV', 'UCEC'))

pan.coh = c('BRCA', 'CRC', 'CCRCC', 'LUAD', 'OV', 'UCEC')
mutations = c('missense', 'truncating', 'synonymous') 

################################################################################
# stacked bar plots (Fig 2A)
################################################################################
rna = read.csv('../data/DNA_RNA_regression/Table.DNA.RNA.regression.linearLIMMA.RNAVsMut.csv')
pro = read.csv('../data/DNA_Pro_regression/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.csv')

# dat = rna
dat = pro
p <- ggplot(dat, aes(x=cancer, fill=mutation)) 
# p <- p + ggtitle('mRNA') + theme_bw(base_size=12)
p <- p + ggtitle('Protein') + theme_bw(base_size=12)
p <- p + geom_bar(position='stack', alpha = .7, stat = 'count', show.legend = FALSE)
p <- p + scale_fill_manual(values = color.palette)
p <- p + theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.text.x = element_text(colour="black",angle = 45, vjust=0.5))
print(p)
# ggsave(file=paste0('../doc/RNA_barplot.png'), w=2.5, h=2.5, unit='in')
ggsave(file=paste0('../doc/Pro_barplot.png'), w=2.5, h=2.5, unit='in')

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
res <- read.csv('../data/DNA_RNA_regression/Table.DNA.RNA.regression.linearLIMMA.RNAVsMut.csv')
# res <- read.csv('../data/DNA_Pro_regression/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.csv')

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
  ggsave(file=paste0('../doc/volcano\ plot/', cancer, 'eQTL.png'), w=2.5, h=2.5)
  # ggsave(file=paste0('../doc/volcano\ plot/', cancer, 'pQTL.png'), w=2.5, h=2.5)
  print(p)
  dev.off()
}

################################################################################
# scatter plots (Fig 3A and Fig 4A)
################################################################################
dat <- read_excel('../data/results/concordant.xlsx')
p <- ggplot(dat) 
top = 3
gg.text = head(unique(dat[order(dat$PrologFC, decreasing = T),]), top)
gg.text = rbind(gg.text, head(unique(dat[order(dat$PrologFC, decreasing = F),]), top))

p <- p + ggtitle('concordant QTLs') + theme_bw(base_size=14) 
p <- p + geom_point(dat, mapping = aes(RNAlogFC, PrologFC, col=mutation), alpha = .7)
p <- p + geom_label_repel(data=gg.text, mapping = aes(RNAlogFC, PrologFC, label=paste(cancer, Gene), col=mutation),
                          min.segment.length = 0, alpha = .9, max.overlaps=50, size=4, show.legend = FALSE)
p <- p + scale_color_manual(values = color.palette)
p <- p + xlab('RNAlogFC') + ylab('PrologFC')
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_abline(intercept = 0, slope = 1, color='grey', size=0.8)
print(p)
ggsave(file='../doc/concordantQTLs.png', w=5, h=3.75)

dat <- read_excel('../data/results/discordant.xlsx')
dat <- dat[dat$overlap == TRUE,]
top = 3
gg.text = head(unique(dat[order(dat$PrologFC, decreasing = T),]), top)
gg.text = rbind(gg.text, head(unique(dat[order(dat$PrologFC, decreasing = F),]), top))
p <- ggplot(dat) 
p <- p + ggtitle('discordant QTLs') + theme_bw(base_size=14)
p <- p + geom_point(dat, mapping = aes(RNAlogFC, PrologFC, col=mutation), alpha = .7)
p <- p + geom_label_repel(data=gg.text, mapping = aes(RNAlogFC, PrologFC, label=paste(cancer, Gene), col=mutation),
                          min.segment.length = 0, alpha = .9, max.overlaps=50, size=4, show.legend = FALSE)
p <- p + scale_color_manual(values = color.palette)
p <- p + xlab('RNAlogFC') + ylab('PrologFC')
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + geom_abline(intercept = 0, slope = 1, color='grey', size=0.8)
ggsave(file='../doc/discordantQTLs.png', w=5, h=3.75)

################################################################################
# violin plots (Fig 3B, Fig 4C, Fig S2)
################################################################################
data.level = 'Pro'
data.level = 'RNA'

mutation = 'missense'
for (mutation in mutations) {
  # exp <- read.csv(paste0('../data/wildVSmut/concordant/', data.level, '_exp_wildVsMut.csv')) #S2
  # dir <- read_excel('../data/results/concordant.xlsx')
  exp <- read.csv(paste0('../data/wildVSmut/discordant/', data.level, '_exp_wildVsMut.csv'))
  dir <- read_excel('../data/results/discordant.xlsx')
  
  exp <- exp[exp$mutation==mutation,]
  exp[exp$isMut == 'mutant', ]$isMut = 'mutated'
  exp[exp$isMut == 'wild', ]$isMut = 'nonmutated'
  genes = unique(exp$gene)
  gene = 'TP53'
  for (gene in genes) {
    if (sum(exp[exp$gene==gene,]$isMut == 'mutated') < 3) {
      next
    } else {
      p <- ggplot(exp[exp$gene==gene,], aes(x=isMut, y=expression, fill=isMut, color=isMut))
      p <- p + facet_grid(.~cancer, scale = "free", space = "free", drop=T)
      p <- p + geom_dotplot(dotsize=3, binwidth=.01, binaxis= "y", color='black', stackdir ="centerwhole", alpha=0.5, show.legend = FALSE)
      p <- p + geom_violin(alpha=0.6, show.legend = FALSE)
      p <- p + theme_bw(base_size=18)
      p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=16))
      p <- p + scale_color_manual(values = color.palette1) + scale_fill_manual(values = color.palette1)
      
      all.text = exp[(exp$gene==gene) & (exp$isMut == 'mutated'),]
      cancers = unique(exp[exp$gene==gene,]$cancer)
      gg.text = c()
      
      for (cancer in cancers){
        text = all.text[all.text$cancer == cancer,]
        if (data.level == 'Pro') {
          if (dir[(dir$Gene==gene) & (dir$mutation == mutation) & (dir$cancer == cancer),]$PrologFC > 0) {
            gg.text = rbind(gg.text, head(text[order(text$expression, decreasing = T), ], 3))
          } else {
            gg.text = rbind(gg.text, head(text[order(text$expression, decreasing = F), ], 3))
          }
        } else {
          if (dir[(dir$Gene==gene) & (dir$mutation == mutation) & (dir$cancer == cancer),]$RNAlogFC > 0) {
            gg.text = rbind(gg.text, head(text[order(text$expression, decreasing = T), ], 3))
          } else {
            gg.text = rbind(gg.text, head(text[order(text$expression, decreasing = F), ], 3))
          }
        }
        
      }
      
      p <- p + geom_label_repel(gg.text, mapping = aes(x=isMut, y=expression, label=hgvsp, col='white'),
                                min.segment.length = 0, alpha = .9, max.overlaps=50, size=6, show.legend = FALSE)
      p <- p + ggtitle(paste(gene, mutation))
      print(p)
      # ggsave(file=paste0('../doc/MutVsWild/concordant/', mutation, '/', data.level, '/', gene, 'ExpMutVsWild.png'), w=4, h=6)
      ggsave(file=paste0('../doc/MutVsWild/discordant/', mutation, '/', data.level, '/', gene, 'ExpMutVsWild.png'),  w=6, h=9)
    }
  }
}

################################################################################
# pie chart of psQTLs for Fig 4B
################################################################################
dat <- read_excel('../data/results/discordant.xlsx')
dat <- dat[dat$overlap == TRUE,]
len <- dim(dat[dat$overlap == TRUE,])[1]
p <- ggplot(dat, aes(x="", fill=mutation)) +
  geom_bar(stat="count", width=1, alpha = .7) + 
  scale_fill_manual(values = color.palette) + 
  coord_polar("y", start=0) + 
  geom_text(stat='count', aes(label=percent(..count../len)), size=4, position = position_stack(vjust = 0.4)) +
  theme_void(base_size=12)
print(p)
ggsave(file=paste0('../doc/psQTLs_mut_count.png'), w=4, h=3)

################################################################################
# heat map (Fig 4B)
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

################################################################################
# volcano plots: Fig 2B and Fig S3 (overlap with DEP)
################################################################################
# res <- read.csv('../data/DNA_RNA_regression/Table.DNA.RNA.regression.linearLIMMA.RNAVsMut.csv')
res <- read.csv('../data/DNA_Pro_regression/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.csv')
concordant_overlap <- read_excel('../data/DEP/concordant_overlap.xlsx')
# discordant_overlap <- read_excel('../data/DEP/discordant_overlap.xlsx')

top = 5
thr.lfc = 1
thr.FDR = 0.05

text = merge(res, concordant_overlap, by.x=c('Gene', 'cancer', 'mutation'), by.y=c('common_gene', 'cancer', 'mutation'))

cancer = 'CCRCC'
for (cancer in pan.coh) {
  res_can = res[res$cancer == cancer,]
  gg.text = text[text$cancer == cancer, ]
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
  # ggsave(file=paste0('../doc/overlap_volcano_plot/', cancer, 'discordant_eQTL.png'), w=2.5, h=2.5)
  ggsave(file=paste0('../doc/overlap_volcano_plot/', cancer, 'concordant_pQTL.png'), w=2.5, h=2.5)
  print(p)
  dev.off()
}

################################################################################
# volcano plots (DEP volcano)
################################################################################
concordant_overlap <- read_excel('../data/DEP/concordant_overlap.xlsx')
discordant_overlap <- read_excel('../data/DEP/discordant_overlap.xlsx')

top = 5
thr.lfc = 1
thr.FDR = 0.05

cancer = 'CCRCC'
for (cancer in pan.coh) {
  res_can <- read_excel(paste0('../data/DEP/', cancer, '.DEP.paired.notcorrected.xlsx'))
  res_can <- read_excel(paste0('../data/DEP_RNA/', cancer, '.DEG.xlsx'))
  if (dim(res_can)[1] == 1)
    next
  overlap_can <- merge(concordant_overlap[concordant_overlap$cancer == cancer, ], res_can, by.x='common_gene', by.y='...1')
  p <- ggplot(res_can) 
  p <- p + geom_vline(xintercept = -thr.lfc, alpha = .3) + 
    geom_vline(xintercept = thr.lfc, alpha = .3)
  p <- p + geom_hline(yintercept = -log10(thr.FDR), alpha = .3)
  p <- p + ggtitle(cancer) + theme_bw()
  p <- p + theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p <- p + geom_point(res_can, mapping = aes(logFC, -log10(adj.P.Val)), alpha = .7, show.legend = FALSE)
  if (dim(overlap_can)[1] != 0) {
    gg.text = overlap_can[(abs(overlap_can$logFC) > thr.lfc) & (overlap_can$adj.P.Val < thr.FDR), ]
    p <- p + geom_label_repel(data=gg.text, mapping = aes(logFC, -log10(adj.P.Val), label=common_gene), 
                              min.segment.length = 0, alpha = .9, max.overlaps=50, size=4, show.legend = FALSE)
  }
  p <- p + scale_color_manual(values = color.palette)
  p <- p + xlab('FC') + ylab('-Log10(FDR)')
  p <- p + theme(axis.text = element_text(colour="black", vjust=0.5, size=11))
  ggsave(file=paste0('../doc/DEP_volcano_plot/', cancer, 'discordant_eQTL.png'), w=2.5, h=2.5)
  #ggsave(file=paste0('../doc/DEP_volcano_plot/', cancer, 'concordant_pQTL.png'), w=2.5, h=2.5)
  print(p)
  dev.off()
}

################################################################################
# TP53 expression (Fig 5A)
################################################################################
pro <- read.csv('../data/tp53_expression/pro_exp_wildVsMut_hgvsp.csv')
rna <- read.csv('../data/tp53_expression/rna_exp_wildVsMut_hgvsp.csv')

norm_pro = c()
for (cancer in pan.coh) {
  tmp = pro[pro$cancer == cancer,]
  exp = tmp$expression
  norm_exp = (exp - min(exp)) / (max(exp) - min(exp))
  tmp$norm_exp = norm_exp
  norm_pro = rbind(norm_pro, tmp)
}

norm_rna = c()
for (cancer in pan.coh) {
  tmp = rna[rna$cancer == cancer,]
  exp = tmp$expression
  norm_exp = (exp - min(exp)) / (max(exp) - min(exp))
  tmp$norm_exp = norm_exp
  norm_rna = rbind(norm_rna, tmp)
}

pro = norm_pro
rna = norm_rna

kotler1 <- read.csv('../data/tp53_kotler1.csv')
kotler2 <- read.csv('../data/tp53_kotler2.csv')

pro_mut <- pro[pro$isMut == 'mutant',]
rna_mut <- rna[rna$isMut == 'mutant',]

tp53_mut = merge(pro_mut, rna_mut, by=c('sample_id', 'mutation', 'cancer', 'hgvsp'))
tp53_mut_k1 = merge(tp53_mut, kotler1, by='hgvsp', all.x=TRUE)
tp53_mut_k2 = merge(tp53_mut, kotler2, by='hgvsp', all.x=TRUE)

p <- ggplot(tp53_mut) 
p <- p + ggtitle('TP53 Mutation') + theme_bw(base_size=14) 
p <- p + geom_point(tp53_mut, mapping = aes(norm_exp.y, norm_exp.x, col=mutation), alpha = .7, show.legend = FALSE)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + scale_color_manual(values = color.palette)
p <- p + xlab('mRNA(FPKM)') + ylab('Protein(log ratio)')
print(p)
ggsave(file=paste0('../doc/tp53_mutation_expression_norm.png'), w=2.5, h=2.5)

p <- ggplot(tp53_mut) 
p <- p + ggtitle('TP53 Mutation') + theme_bw(base_size=14) 
p <- p + geom_point(tp53_mut, mapping = aes(norm_exp.y, norm_exp.x, col=cancer), alpha = .7)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + scale_color_manual(values = color.palette2)
p <- p + xlab('mRNA(FPKM)') + ylab('Protein(log ratio)')
print(p)
ggsave(file=paste0('../doc/tp53_mutation_expression1_norm.png'), w=3.5, h=2.5)

# 
p <- ggplot(tp53_mut) 
p <- p + ggtitle('TP53 Mutation') + theme_bw(base_size=14) 
p <- p + geom_point(tp53_mut_k1, mapping = aes(norm_exp.y, norm_exp.x, col=Enrichment_invivo>-1))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab('mRNA(FPKM)') + ylab('Protein(log ratio)')
ggsave(file=paste0('../doc/tp53_mutation_expression2_norm.png'), w=5, h=2.5)







# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#     Figure 5B (overlap_with_kotler.R)   
#
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|

rm(list = ls(all.names = TRUE))
setwd('~/Library/CloudStorage/Box-Box/Huang_lab/manuscripts/pQTL/analysis')
library('readxl')
library('stringr')
library('ggplot2')
library('cowplot')
library('ggrepel')
library('gsubfn')

pan.coh = c('BRCA', 'CRC', 'CCRCC', 'LUAD', 'OV', 'UCEC')


pro <- read.csv('../data/tp53_expression/pro_exp_wildVsMut_hgvsp.csv') #  all tp53
rna <- read.csv('../data/tp53_expression/rna_exp_wildVsMut_hgvsp.csv') #  all tp53
# pro <- read.csv(paste0('../data/wildVSmut/discordant/', 'Pro', '_exp_wildVsMut.csv'))
# rna <- read.csv(paste0('../data/wildVSmut/discordant/', 'RNA', '_exp_wildVsMut.csv'))
# pro_c <- read.csv(paste0('../data/wildVSmut/concordant/', 'Pro', '_exp_wildVsMut.csv'))
# rna_c <- read.csv(paste0('../data/wildVSmut/concordant/', 'RNA', '_exp_wildVsMut.csv'))

norm_pro = c()
for (cancer in pan.coh) {
  tmp = pro[pro$cancer == cancer,]
  exp = tmp$expression
  norm_exp = (exp - min(exp)) / (max(exp) - min(exp))
  tmp$norm_exp = norm_exp
  norm_pro = rbind(norm_pro, tmp)
}

norm_rna = c()
for (cancer in pan.coh) {
  tmp = rna[rna$cancer == cancer,]
  exp = tmp$expression
  norm_exp = (exp - min(exp)) / (max(exp) - min(exp))
  tmp$norm_exp = norm_exp
  norm_rna = rbind(norm_rna, tmp)
}

pro = norm_pro
rna = norm_rna

pro_wt_median = median(pro[pro$isMut == 'wild', ]$norm_exp)
rna_wt_median = median(rna[rna$isMut == 'wild', ]$norm_exp)

pro_tp53 <- pro[(pro$gene == 'TP53') & (pro$isMut == 'mutant'), ]
rna_tp53 <- rna[(rna$gene == 'TP53') & (rna$isMut == 'mutant'), ]

tp53_df <- merge(pro_tp53, rna_tp53, by='sample_id')
tp53_df <- tp53_df[,c('sample_id', 'mutation.x', 'cancer.x', 'hgvsp.x', 'norm_exp.x', 'norm_exp.y')]
colnames(tp53_df) <- c('sample_id', 'mutation', 'cancer', 'hgvsp', 'expression.pro', 'expression.rna')
tp53_df$Codon_num = strapply(tp53_df$hgvsp, "(\\d+).*", as.numeric, simplify = c)

kotler1 <- read_excel('../../../Huang_lab_data/TP53_Kotler_MolCell2018/1-s2.0-S1097276518304544-mmc7.xlsx') # Enrichment_invivo
kotler2 <- read_excel('../../../Huang_lab_data/TP53_Kotler_MolCell2018/1-s2.0-S1097276518304544-mmc4.xlsx') # RFS and IARC
hgvsp <- paste0('p.', str_sub(kotler1$AA_change, 1, 1), kotler1$Codon_num, str_sub(kotler1$AA_change, 3))
kotler1$hgvsp = hgvsp
kotler1$mutation = 'missense'
kotler1[grepl(".S", kotler1$AA_change) | grepl("\\*", kotler1$AA_change), ]$mutation = 'truncating'
hgvsp <- paste0('p.', str_sub(kotler2$AA_change, 1, 1), kotler2$Codon_num, str_sub(kotler2$AA_change, 3))
kotler2$hgvsp = hgvsp
kotler2$mutation = 'missense'
kotler2[grepl(".S", kotler2$AA_change) | grepl("\\*", kotler2$AA_change), ]$mutation = 'truncating'

tp53_kotler1 <- merge(tp53_df, kotler1, by='hgvsp')
write.csv(tp53_kotler1, '../data/tp53_kotler1.csv')
tp53_kotler2 <- merge(tp53_df, kotler2, by='hgvsp')
write.csv(tp53_kotler2, '../data/tp53_kotler2.csv')




# top %20 (RFS)
top = 3

tp53 = tp53_df[(tp53_df$mutation == 'missense'),]

tp53_top = head(tp53[order(-tp53$expression.pro),], round(dim(tp53)[1] / 5))
kotler2 = kotler2[kotler2$mutation == 'missense',]
tp53_rest = tp53[!tp53$hgvsp %in% tp53_top$hgvsp,]

top_overlap = merge(tp53_top, kotler2, by='hgvsp')
top_overlap$category = 'category 1'
top_overlap = top_overlap[, c('RFS_H1299', 'category', 'hgvsp')]
gg.text = head(unique(top_overlap[order(top_overlap$RFS_H1299, decreasing = T),]), top)

rest_overlap = merge(tp53_rest, kotler2, by='hgvsp')
rest_overlap$category = 'category 2'
rest_overlap = rest_overlap[, c('RFS_H1299', 'category', 'hgvsp')]

kotler2_rest = kotler2[!kotler2$hgvsp %in% tp53$hgvsp,]
kotler2_rest$category = "category 3"
kotler2_rest = kotler2_rest[, c('RFS_H1299', 'category', 'hgvsp')]

dat = rbind(kotler2_rest, rest_overlap, top_overlap)

p <- ggplot(dat, aes(x=category, y=RFS_H1299, fill="#DAF1FC")) + ylab('RFS') + theme_bw(base_size=14)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p <- p + geom_dotplot(dotsize=3, binwidth=.01, binaxis= "y", color="#DAF1FC", stackdir ="centerwhole",alpha=0.5, show.legend = FALSE)
p <- p + geom_violin(alpha=0.6, show.legend = FALSE) 
p <- p + scale_color_manual(values = "#DAF1FC") + scale_fill_manual(values = "#DAF1FC")
p <- p + geom_label_repel(data=gg.text, mapping = aes(category, RFS_H1299, label=hgvsp),
                          min.segment.length = 0, alpha = .9, max.overlaps=50, size=4, show.legend = FALSE)
p <- p + geom_hline(yintercept = -1, alpha = .3)
p <- p + theme(axis.text.x = element_text(colour="black",angle = 45, vjust=0.5))

print(p)

wilcox.test(top_overlap$RFS_H1299, rest_overlap$RFS_H1299, alternative="greater")
wilcox.test(top_overlap$RFS_H1299, kotler2_rest$RFS_H1299, alternative="greater")
wilcox.test(rest_overlap$RFS_H1299, kotler2_rest$RFS_H1299)

ggsave(paste0(file='../doc/validation_with_kotler/RFS_violin_pro.png'), w=4, h=4)


# top %20 (enrichment in vivo)
top = 3

tp53 = tp53_df[(tp53_df$mutation == 'missense'),]

tp53_top = head(tp53[order(-tp53$expression.pro),], round(dim(tp53)[1] / 5))
kotler1 = kotler1[kotler1$mutation == 'missense',]
tp53_rest = tp53[!tp53$hgvsp %in% tp53_top$hgvsp,]

enrichment_median = median(kotler1$Enrichment_invivo)

top_overlap = merge(tp53_top, kotler1, by=c('hgvsp', 'mutation'))
top_overlap$category = 'category 1'
top_overlap = top_overlap[, c('Enrichment_invivo', 'category', 'hgvsp')]
gg.text = head(unique(top_overlap[order(top_overlap$Enrichment_invivo, decreasing = T),]), top)

rest_overlap = merge(tp53_rest, kotler1, by=c('hgvsp', 'mutation'))
rest_overlap$category = 'category 2'
rest_overlap = rest_overlap[, c('Enrichment_invivo', 'category', 'hgvsp')]

kotler1_rest = kotler1[!kotler1$hgvsp %in% tp53$hgvsp,]
kotler1_rest$category = "category 3"
kotler1_rest = kotler1_rest[, c('Enrichment_invivo', 'category', 'hgvsp')]

dat = rbind(kotler1_rest, rest_overlap, top_overlap)

p <- ggplot(dat, aes(x=category, y=Enrichment_invivo, fill="#DAF1FC")) + ylab('Enrichment score in vivo') + theme_bw(base_size=14)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p <- p + geom_dotplot(dotsize=3, binwidth=.01, binaxis= "y", color="#DAF1FC", stackdir ="centerwhole", alpha=0.5, show.legend = FALSE)
p <- p + geom_violin(alpha=0.6, show.legend = FALSE) 
p <- p + scale_color_manual(values = "#DAF1FC") + scale_fill_manual(values = "#DAF1FC")
p <- p + geom_label_repel(data=gg.text, mapping = aes(category, Enrichment_invivo, label=hgvsp),
                          min.segment.length = 0, alpha = .9, max.overlaps=50, size=4, show.legend = FALSE)
p <- p + geom_hline(yintercept = enrichment_median, alpha = .3)
p <- p + theme(axis.text.x = element_text(colour="black",angle = 45, vjust=0.5))

print(p)

wilcox.test(top_overlap$Enrichment_invivo, rest_overlap$Enrichment_invivo, alternative="greater")
wilcox.test(top_overlap$Enrichment_invivo, kotler1_rest$Enrichment_invivo, alternative="greater")
wilcox.test(rest_overlap$Enrichment_invivo, kotler1_rest$Enrichment_invivo)

ggsave(paste0(file='../doc/validation_with_kotler/Enrichment_in_vivo_violin_pro.png'), w=4, h=4)

# top %20 (IARC)
top = 3

tp53 = tp53_df[(tp53_df$mutation == 'missense'),]

tp53_top = head(tp53[order(-tp53$expression.rna),], round(dim(tp53)[1] / 5))
kotler2 = kotler2[kotler2$mutation == 'missense',]
tp53_rest = tp53[!tp53$hgvsp %in% tp53_top$hgvsp,]

top_overlap = merge(tp53_top, kotler2, by='hgvsp')
top_overlap$category = 'category 1'
top_overlap = top_overlap[, c('IARC_occurrences', 'category', 'hgvsp')]
gg.text = head(unique(top_overlap[order(top_overlap$IARC_occurrences, decreasing = T),]), top)

rest_overlap = merge(tp53_rest, kotler2, by='hgvsp')
rest_overlap$category = 'category 2'
rest_overlap = rest_overlap[, c('IARC_occurrences', 'category', 'hgvsp')]

kotler2_rest = kotler2[!kotler2$hgvsp %in% tp53$hgvsp,]
kotler2_rest$category = "category 3"
kotler2_rest = kotler2_rest[, c('IARC_occurrences', 'category', 'hgvsp')]

dat = rbind(kotler2_rest, rest_overlap, top_overlap)

p <- ggplot(dat, aes(x=category, y=IARC_occurrences, fill="#DAF1FC")) + ylab('IARC_occurrences') + theme_bw(base_size=14)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p <- p + geom_dotplot(dotsize=3, binwidth=.01, binaxis= "y", color="#DAF1FC", stackdir ="centerwhole",alpha=0.5, show.legend = FALSE)
p <- p + geom_violin(alpha=0.6, show.legend = FALSE) 
p <- p + scale_color_manual(values = "#DAF1FC") + scale_fill_manual(values = "#DAF1FC")
p <- p + geom_label_repel(data=gg.text, mapping = aes(category, IARC_occurrences, label=hgvsp),
                          min.segment.length = 0, alpha = .9, max.overlaps=50, size=4, show.legend = FALSE)
p <- p + geom_hline(yintercept = -1, alpha = .3)
p <- p + theme(axis.text.x = element_text(colour="black",angle = 45, vjust=0.5))

print(p)

wilcox.test(top_overlap$IARC_occurrences, rest_overlap$IARC_occurrences, alternative="greater")
wilcox.test(top_overlap$IARC_occurrences, kotler2_rest$IARC_occurrences, alternative="greater")
wilcox.test(rest_overlap$IARC_occurrences, kotler2_rest$IARC_occurrences)

ggsave(paste0(file='../doc/validation_with_kotler/IARC_violin_.png'), w=4, h=4)

