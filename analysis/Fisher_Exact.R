rm(list = ls(all.names = TRUE))
setwd('~/Library/CloudStorage/Box-Box/Huang_lab/manuscripts/pQTL/analysis')
library('dplyr')

cancers = c('BRCA', 'CRC', 'CCRCC', 'LUAD', 'OV', 'UCEC')

pro = read.csv('../data/DNA_Pro_regression/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.csv')

dat = pro
dat = dat[dat$FDR < 0.05, ] # pQTLs

tru = pro[pro$mutation == 'truncating',]
not_tru = pro[!pro$mutation == 'truncating',]
ori_tru_counts = tru %>% group_by(cancer) %>% dplyr ::summarise(count_muts = n())
ori_not_tru_counts = not_tru %>% group_by(cancer) %>% dplyr ::summarise(count_muts = n())

tru = dat[dat$mutation == 'truncating',]
not_tru = dat[!dat$mutation == 'truncating',]
qtl_tru_counts = tru %>% group_by(cancer) %>% dplyr ::summarise(count_muts = n())
qtl_not_tru_counts = not_tru %>% group_by(cancer) %>% dplyr ::summarise(count_muts = n())
# qtl_not_tru_counts = rbind(qtl_not_tru_counts, c('CCRCC', 0))

can = 'CCRCC'
can = 'BRCA'
for (can in cancers) {
  print(can)
  table = matrix(c(ori_tru_counts[ori_tru_counts$cancer == can,]$count_muts, 
                   qtl_tru_counts[qtl_tru_counts$cancer == can,]$count_muts, 
                   ori_not_tru_counts[ori_not_tru_counts$cancer == can,]$count_muts,
                   qtl_not_tru_counts[qtl_not_tru_counts$cancer == can,]$count_muts),
                   byrow=TRUE, nrow=2, 
                   dimnames = list(c('tru', 'not_tru'),
                                 c('original_cohort', 'QTLs')))

  print(fisher.test(t(table)))
}


