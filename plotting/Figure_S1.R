# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#     Differential expression of gene-cancer seQTL/spQTL pairs in retrospective CPTAC cohorts: BRCA, CRC, and OV (Figure S1B) [retro_volcano.R]  
#
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|

rm(list = ls(all.names = TRUE))
setwd('~/Library/CloudStorage/Box-Box/Huang_lab/manuscripts/pQTL/analysis')
library('ggplot2')
library('ggrepel')
library('readxl')

color.palette =  setNames(c("#8DD3C7","#FB8072","#BEBADA"),
                          c('missense', 'truncating', 'synonymous'))

top = 5
thr.lfc = 1
thr.FDR = 0.05

mut.types = c('missense', 'truncating', 'synonymous')
pan.coh = c('CRC', 'BRCA', 'OV')

cancer = 'OV'
mut = 'missense'

for (cancer in pan.coh) {
  dat = c()
  for (mut in mut.types) {
    df = read_excel(paste0('../data/RNA_regression_retro_cptac/Table.DNA.RNA.regression.linearLIMMA.RNAVsMut.', cancer, '.', mut, '.xlsx'))
    # df = read_excel(paste0('../data/Pro_regression_retro_cptac/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.', cancer, '.', mut, '.xlsx'))
    df$mutation = mut
    dat = rbind(dat, df)
  }
    gg.text = dat[(abs(dat$logFC) > thr.lfc) & (dat$FDR < thr.FDR), ]
    gg.text = head(unique(gg.text[order(-log(gg.text$FDR), decreasing = T),]), top)
    p <- ggplot(dat) 
    p <- p + geom_vline(xintercept = -thr.lfc, alpha = .3) + 
      geom_vline(xintercept = thr.lfc, alpha = .3)
    p <- p + geom_hline(yintercept = -log10(thr.FDR), alpha = .3)
    p <- p + ggtitle(cancer) + theme_bw()
    p <- p + theme(legend.position="bottom", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + geom_point(dat, mapping = aes(logFC, -log10(FDR), col=mutation), alpha = .7, show.legend = FALSE)
    p <- p + geom_label_repel(data=gg.text, mapping = aes(logFC, -log10(FDR), label=Gene, col=mutation), 
                              min.segment.length = 0, alpha = .9, max.overlaps=50, show.legend = FALSE)
    p <- p + scale_color_manual(values = color.palette)
    p <- p + xlab('FC') + ylab('-Log10(FDR)')
    p <- p + theme(axis.text = element_text(colour="black", vjust=0.5, size=11))
  
    ggsave(file=paste0('../doc/retro_volcano/', cancer, 'eQTL.png'), w=2.5,h=2.5)
    # ggsave(file=paste0('../doc/retro_volcano/', cancer, 'pQTL.png'), w=2.5,h=2.5)
    print(p)
    dev.off()
}
