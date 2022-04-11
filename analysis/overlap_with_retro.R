rm(list = ls(all.names = TRUE))
setwd('~/Library/CloudStorage/Box-Box/Huang_lab/manuscripts/pQTL/analysis')
library('readxl')
library('ggplot2')

mut.types = c('missense', 'truncating', 'synonymous')
cancers = c('CRC', 'OV', 'BRCA')
mut = 'missense'
cancer = 'OV'

# add the places of OV to data frame
places = c('JHU', 'PNNL')
for (mut in mut.types) {
  df = data.frame()
  for (place in places) {
  dat = read_excel(paste0('../data/Pro_regression_retro_cptac/Table.DNA.PRO.',place,'.regression.linearLIMMA.ProVsMut.OV.',mut,'.xlsx'))
  dat$institution = place
  df = rbind(df, dat)
  xlsx::write.xlsx(df, paste0('../data/Pro_regression_retro_cptac/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.OV.',mut,'.xlsx'))
  }
}

for (mut in mut.types) {
  # cur <- read_excel(paste0('../data/DNA_PRO_regression/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.',cancer,'.',mut,'.xlsx'))
  # retro <- read_excel(paste0('../data/Pro_regression_retro_cptac/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.',cancer,'.',mut,'.xlsx'))
  
  cur <- read_excel(paste0('../data/DNA_RNA_regression/Table.DNA.RNA.regression.linearLIMMA.RNAVsMut.',cancer,'.',mut,'.xlsx'))
  retro <- read_excel(paste0('../data/RNA_regression_retro_cptac/Table.DNA.RNA.regression.linearLIMMA.RNAVsMut.',cancer,'.',mut,'.xlsx'))

  cur_sig <- cur[cur$P.value < 0.05,]
  cur_sig <- cur_sig[, c('Gene', 'logFC', 'FDR')]
  retro_sig <- retro[retro$P.value < 0.05,]
  retro_sig <- retro_sig[, c('Gene', 'logFC', 'FDR')]
  
  overlap_sig <- merge(cur_sig, retro_sig, by='Gene')
  
  colnames(overlap_sig) <- sub(".y",".retro", colnames(overlap_sig))
  
  overlap_sig$mut.typ <- mut
  overlap_sig$cancer <- cancer
  
  # xlsx::write.xlsx(overlap_sig, paste0('../data/overlap_with_retro/pQTLs/', cancer, mut, 'p_value.05_overlap_pQTLs.xlsx'))
  xlsx::write.xlsx(overlap_sig, paste0('../data/overlap_with_retro/eQTLs/', cancer, mut, 'p_value.05_overlap_eQTLs.xlsx'))
}

# put all results in one file
file_names <- list.files('../data/overlap_with_retro/eQTLs', pattern = 'eQTLs.xlsx', full.names = TRUE)
df = data.frame()
for (file_name in file_names) {
  dat = read_excel(file_name)
  df = rbind(df, dat)
}
xlsx::write.xlsx(df, paste0('../data/overlap_with_retro/eQTLs/overlap_eQTLs.xlsx'))