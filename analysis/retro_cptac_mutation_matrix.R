rm(list = ls(all.names = TRUE))
setwd('~/Library/CloudStorage/Box-Box/Huang_lab/manuscripts/pQTL/analysis')
library(limma)
library(readxl)

cancers = c('OV', 'CRC', 'BRCA')
cancer = 'OV'

for (cancer in cancers) {
  tcga = read.table('../../../Huang_lab_data/TCGA_PanCanAtlas_2018/MC3_Ellrott_CellSys2018/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample.gz', header = T)
  cptac_pro = read.table(paste0('../../../Huang_lab_data/CPTAC2retrospective/', cancer,'/', cancer, '_PRO_formatted.txt'), quote = "", sep = '\t', header = T)
  cptac_rna = read.table(paste0('../../../Huang_lab_data/CPTAC2retrospective/', cancer,'/', cancer, '_mRNA_formatted.txt'), quote = "", sep = '\t', header = T)
  colnames(cptac_rna) <- sub("X","", colnames(cptac_rna))
  colnames(cptac_pro) <- sub("X","", colnames(cptac_pro))
  cptac_sam = c(colnames(cptac_pro), colnames(cptac_rna))
  
  format_mut = function(col){
    col = substr(col, 1, 16)
    col = sub("TCGA-","", col)
    col = gsub("-",".", col)
    return(col)
  }
  
  tcga$cptacID = format_mut(tcga$Tumor_Sample_Barcode)
  tcga_cptac <- tcga[tcga$cptacID %in% cptac_sam,]
  rownames(tcga_cptac) <- 1:nrow(tcga_cptac)
  
  mutation <- data.frame(matrix(0, length(unique(tcga_cptac$Hugo_Symbol)), length(unique(tcga_cptac$cptacID))))
  rownames(mutation) = unique(tcga_cptac$Hugo_Symbol)
  colnames(mutation) = unique(tcga_cptac$cptacID)
  
  # mut.type = 'Missense_Mutation'; mut = 'missense'
  mut.type = 'Silent'; mut = 'synonymous'
  # mut.type = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Nonsense_Mutation","Splice_Site","Start_Codon_Del","Start_Codon_Del"); mut = 'truncating'
  df = mutation
  for (i in 1:nrow(tcga_cptac)) {
    gene = tcga_cptac[i,]$Hugo_Symbol
    cptacID = tcga_cptac[i,]$cptacID
    if (tcga_cptac[i,]$Variant_Classification %in% mut.type)
      df[gene, cptacID] = 1
  }
  
  write.csv(df, paste0('../data/retro_cptac/', cancer, '_', mut, '.csv'))
}

