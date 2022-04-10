rm(list = ls(all.names = TRUE))
setwd('~/Library/CloudStorage/Box-Box/Huang_lab/manuscripts/pQTL/analysis')
library(limma)
library(readxl)

cancer = 'OV'
mut.types = c('missense', 'truncating', 'synonymous')
mut.type = 'missense'

for (mut.type in mut.types) {
  rna = read.table(paste0('../../../Huang_lab_data/CPTAC2retrospective/', cancer, '/', cancer, '_mRNA_formatted_normalized.txt'))
  cli = read.table(paste0('../../../Huang_lab_data/CPTAC2retrospective/', cancer, '/', cancer, '_clinical.txt'), quote = "", sep = '\t', header = T)
  mut = read.csv(paste0('../data/retro_cptac_mutation/', cancer, '_', mut.type, '.csv'), header = T)
  
  rownames(cli) <- cli$X
  rownames(mut) <- mut$X
  colnames(mut) <- sub("X","", colnames(mut))
  colnames(cli) <- sub("X","", colnames(cli))
  colnames(rna) <- sub("X","", colnames(rna))
  
  sam = intersect(colnames(rna), colnames(cli))
  sam = intersect(sam, colnames(mut))
  
  gen = intersect(rownames(mut), rownames(rna))
  
  mut = mut[gen, sam]
  rna = rna[gen, sam]
  cli = cli[, sam]

  # confounding covariates
  des.con = c()
  if (sum(!is.na(cli['DAYS_TO_BIRTH',])) > 1) {
      des.con = cbind(des.con, age = - as.numeric(cli['DAYS_TO_BIRTH',])/365)
  }
  # if (sum(!is.na(cli['DAYS_TO_BIRTH',])) > 1) {
  #   des.con = cbind(des.con, age = - as.numeric(cli['birth_days_to',]) / 365)
  # }
  if (sum(!is.na(cli['GENDER',])) > 1) {
    des.con = cbind(des.con, sex = factor(as.character(cli['GENDER',])))
  }
  if (sum(!is.na(cli['ETHNICITY',])) > 1) {
      des.con = cbind(des.con, eth = factor(as.character(cli['ETHNICITY',])))
  }
  # if (sum(!is.na(cli['ethnicity',])) > 1) {
  #   des.con = cbind(des.con, eth = factor(as.character(cli['ethnicity',])))
  # }
  
  daf = rna
  # regression
  res = c(); genes = c(); transcriptome = c()
  gene = 'A2M'
  for (gene in gen) {
    mutation = as.numeric(mut[gene,])
    isna = is.na(mutation) | rowSums(is.na(des.con))>0
    mutation = mutation[!isna]
    des.mat = model.matrix(~mutation)
    des.mat = cbind(des.mat, des.con[!isna,])
    computed = tryCatch({ 
      fit = lmFit(as.numeric(daf[match(gene, rownames(rna)), !isna]), des.mat)
      fit2 = eBayes(fit)
      tmp = topTable(fit2, coef='mutation', adjust.method='BH', number=nrow(daf), p.value=1.1, sort.by = 'none')  
      T},error = function(e) F, finally = {})
    if (computed){
      res = rbind(res, tmp)
      genes = c(genes, gene)
      transcriptome = c(transcriptome,  rownames(rna)[match(gene, rownames(rna))])
    }
  }
  res = data.frame(res, row.names = genes)
  res$FDR = p.adjust(res$P.Value, method = 'BH')
  rankd = order(-log10(res$FDR)-log10(res$P.Value), decreasing = T, na.last = T)
  res = res[rankd,]; genes = genes[rankd]; transcriptome = transcriptome[rankd] 
  valid = !is.na(res$FDR) & !is.na(res$P.Value)
  res = res[valid,]; genes = genes[valid]; transcriptome = transcriptome[valid] 
  colnames(res) = gsub('P.Value','P.value',colnames(res))
  if (any(valid)){
    xlsx::write.xlsx(cbind(Gene = genes, res[,-5], RNA = transcriptome), paste0('../data/RNA_regression_retro_cptac/Table.DNA.RNA.regression.linearLIMMA.RNAVsMut.',cancer,'.',mut.type,'.xlsx'), row.names = F)
  }
  # _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
} # analysis

