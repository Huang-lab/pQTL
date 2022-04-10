rm(list = ls(all.names = TRUE))
setwd('~/Library/CloudStorage/Box-Box/Huang_lab/manuscripts/pQTL/analysis')
library(limma)
library(readxl)

cancer = 'OV'
mut.types = c('missense', 'truncating', 'synonymous')
mut.type = 'missense'

for (mut.type in mut.types) {
  # pro = read.table(paste0('../../../Huang_lab_data/CPTAC2retrospective/', cancer, '/', cancer, '_PRO_formatted_normalized.txt'))
  # pro = read.table('../../../Huang_lab_data/CPTAC2retrospective/OV/OV_JHU_PRO_formatted_normalized.txt')
  pro = read.table('../../../Huang_lab_data/CPTAC2retrospective/OV/OV_PNNL_PRO_formatted_normalized.txt')
  cli = read.table(paste0('../../../Huang_lab_data/CPTAC2retrospective/', cancer, '/', cancer, '_clinical.txt'), quote = "", sep = '\t', header = T)
  mut = read.csv(paste0('../data/retro_cptac_mutation/', cancer, '_', mut.type, '.csv'), header = T)
  
  rownames(cli) <- cli$X
  rownames(mut) <- mut$X
  colnames(mut) <- sub("X","", colnames(mut))
  colnames(cli) <- sub("X","", colnames(cli))
  colnames(pro) <- sub("X","", colnames(pro))
  
  sam = intersect(colnames(pro), colnames(cli))
  sam = intersect(sam, colnames(mut))
  
  gen = intersect(rownames(mut), rownames(pro))
  
  mut = mut[gen, sam]
  pro = pro[gen, sam]
  cli = cli[, sam]
  
  des.con = c()
  # if (sum(!is.na(cli['DAYS_TO_BIRTH',])) > 1) {
  #   des.con = cbind(des.con, age = - as.numeric(cli['DAYS_TO_BIRTH',])/365)
  # }
  if (sum(!is.na(cli['birth_days_to',])) > 1) {
    des.con = cbind(des.con, age = - as.numeric(cli['birth_days_to',]) / 365)
  }
  # if (sum(!is.na(cli['GENDER',])) > 1) {
  #   des.con = cbind(des.con, sex = factor(as.character(cli['GENDER',])))
  # }
  # if (sum(!is.na(cli['ETHNICITY',])) > 1) {
  #   des.con = cbind(des.con, eth = factor(as.character(cli['ETHNICITY',])))
  # }
  if (sum(!is.na(cli['ethnicity',])) > 1) {
    des.con = cbind(des.con, eth = factor(as.character(cli['ethnicity',])))
  }

  daf = pro
  # regression
  res = c(); genes = c(); proteins = c() 
  gene = 'PIK3CA'
  for (gene in gen) {
    mutation = as.numeric(mut[gene,])
    isna = is.na(mutation) | rowSums(is.na(des.con))>0
    mutation = mutation[!isna]
    des.mat = model.matrix(~mutation)
    des.mat = cbind(des.mat, des.con[!isna,])
    computed = tryCatch({ 
      fit = lmFit(as.numeric(daf[match(gene, rownames(pro)), !isna]), des.mat)
      fit2 = eBayes(fit)
      tmp = topTable(fit2, coef='mutation', adjust.method='BH', number=nrow(daf), p.value=1.1, sort.by = 'none')  
      T},error = function(e) F, finally = {})
    if (computed){
      res = rbind(res, tmp)
      genes = c(genes, gene)
      proteins = c(proteins, rownames(pro)[match(gene, rownames(pro))])
    }
  }
  res = data.frame(res, row.names = genes)
  res$FDR = p.adjust(res$P.Value, method = 'BH')
  rankd = order(-log10(res$FDR)-log10(res$P.Value), decreasing = T, na.last = T)
  res = res[rankd,]; genes = genes[rankd]; proteins = proteins[rankd] 
  valid = !is.na(res$FDR) & !is.na(res$P.Value)
  res = res[valid,]; genes = genes[valid]; proteins = proteins[valid] 
  colnames(res) = gsub('P.Value','P.value',colnames(res))
  if (any(valid)){
    xlsx::write.xlsx(cbind(Gene = genes, res[,-5], Protein = proteins), paste0('../data/Pro_regression_retro_cptac/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.',cancer,'.',mut.type,'.xlsx'), row.names = F)
  }
}
