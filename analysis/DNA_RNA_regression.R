# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
# 
#   
#
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
rm(list = ls(all.names = TRUE))
setwd('~/Library/CloudStorage/Box-Box/Huang_lab/manuscripts/pQTL/analysis')
library(limma)
library(readxl)

# analysis = 'Poisson'
# analysis = 'LinearMut~RNA'
analysis = 'LinearRNA~Mut'

# pan.coh='BRCA' #***
pan.coh = c('BRCA', 'CRC', 'CCRCC', 'LUAD', 'OV', 'UCEC')
mut.types = c('missense', 'truncating', 'synonymous')

cancer = 'BRCA'
for (cancer in pan.coh){
  message('Processing ', cancer, ' ...')
  
  rna0 = read.table(paste0('../../../Huang_lab_data/PanCancerProteomicsData_HuangLab/',cancer,'.transcriptome.FPKM.formatted.txt.gz'))
  
  cli0 = read.table(paste0('../../../Huang_lab_data/PanCancerProteomicsData_HuangLab/',cancer,'.clinical.formatted.txt'), sep = '\t', header = T)
  
  mut.typ = 'missense'
  for (mut.typ in mut.types){
    rna = log2(1+rna0);
    cli = cli0; 
    mut = read.table(paste0('../../../Huang_lab_data/PanCancerProteomicsData_HuangLab/',cancer,'.WXS.SomaticVariant.',mut.typ,'.txt.gz'))
    # align 
    colnames(mut) = gsub('\\.Tumor.*$','',gsub('^X','',colnames(mut)))
    colnames(rna) = gsub('\\.Tumor.*$','',gsub('^X','',colnames(rna)))
    sam = intersect(cli$Participant.ID, colnames(mut))
    sam = intersect(sam, colnames(rna))
    
    gen = intersect(rownames(mut), rownames(rna))
    
    mut = mut[gen, sam]
    rna = rna[gen, sam]
    cli = cli[match(sam, cli$Participant.ID),]
    gen = intersect(rownames(mut)[which(rowSums(mut!=0, na.rm = T) > 0)],rownames(rna)[which(rowSums(rna!=0, na.rm = T) > 0)])
    mut = mut[gen, ]
    rna = rna[gen, ]
    
    
    if (analysis=='Poisson'){
      # _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
      # 
      #' #' Poission regression
      #' #' Dependent variable: the number mutations (missense/non-missense/truncating) in a given gene observed across the 
      #' #' patients. Predictors: the phosphorylation levels, gene expressions, protein expressions, and clinical/batch effects.
      # 
      # _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
      # Use binary mutation status instead of mutation counts
      mut = (mut>0)*1 
      res = c(); genes = c(); rnas = c() 
      gene = 'PIK3CA'
      for (gene in gen) {
        daf = data.frame(mut = as.numeric(mut[gene,]),
                         rna = as.numeric(rna[gene,]))
        formula = 'mut ~ rna'
        if ('Age.in.months.at.sampling' %in% colnames(cli)){
          if (length(levels(factor(cli$Age.in.months.at.sampling))) > 1) {
            daf$age = as.numeric(cli$Age.in.months.at.sampling)/12; formula = paste(formula, '+ age')}}
        if ('Gender' %in% colnames(cli)){
          if (length(levels(factor(cli$Gender))) > 1) {
            daf$sex = factor(as.character(cli$Gender)); formula = paste(formula, '+ sex')}}
        if ('Ethnicity' %in% colnames(cli)){
          if (length(levels(factor(cli$Ethnicity))) > 1) {
            daf$eth = factor(as.character(cli$Ethnicity)); formula = paste(formula, '+ eth')}}
        if ('TMT.plex' %in% colnames(cli)){
          if (length(levels(factor(cli$TMT.plex))) > 1) {
            daf$bat = factor(as.numeric(cli$TMT.plex)); formula = paste(formula, '+ bat')}}
        # summary(daf)
        computed = tryCatch({
          reg = summary(glm(formula, family="poisson", data=daf))$coefficients
          T}, error = function(e){e; F}, finally = {})
        if (computed & ('rna' %in% rownames(reg))) {
          res = rbind(res, reg['rna',])
          genes = c(genes, gene)
          rnas = c(rnas, gsub('\\:NP.*\\:','\\:',rownames(rna)[match(gene, rownames(rna))]) )
          # phosphosites = c(phosphosites, gsub('\\:NP.*\\:','\\:',rownames(pho)[match(gene, pho.gen)]) )
        }
      } # for gene
      res = data.frame(res, row.names = genes)
      colnames(res) = c('Beta', 'Std.Error', 'Z.value', 'P.value')
      res$FDR = p.adjust(res$P.value, method = 'BH')
      rankd = order(-log10(res$FDR)-log10(res$P.value), decreasing = T, na.last = T)
      res = res[rankd,]; genes = genes[rankd]; rnas = rnas[rankd] 
      valid = !is.na(res$FDR) & !is.na(res$P.value)
      res = res[valid,]; genes = genes[valid]; rnas = rnas[valid] 
      if (any(valid)){
        xlsx::write.xlsx(cbind(Gene = genes, RNA = rnas, res), paste0('../data/DNA_RNA_regression/Table.DNA.RNA.regression.Poisson.',cancer,'.',mut.typ,'.xlsx'), row.names = F)
      }
      # _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
      
      
    } else if (analysis=='LinearMut~RNA') {    
      # _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
      # 
      # LIMMA linear regression
      # 
      # _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
      # confounding covariates
      des.con = c()
      if ('Age.in.months.at.sampling' %in% colnames(cli)){
        if (length(levels(factor(cli$Age.in.months.at.sampling))) > 1) {
          des.con = cbind(des.con, age = as.numeric(cli$Age.in.months.at.sampling)/12)}}
      if ('Gender' %in% colnames(cli)){
        if (length(levels(factor(cli$Gender))) > 1) {
          des.con = cbind(des.con, sex = factor(as.character(cli$Gender)))}}
      if ('Ethnicity' %in% colnames(cli)){
        if (length(levels(factor(cli$Ethnicity))) > 1) {
          des.con = cbind(des.con, eth = factor(as.character(cli$Ethnicity)))}}
      if ('TMT.plex' %in% colnames(cli)){ 
        if (length(levels(factor(cli$TMT.plex))) > 1) {
          des.con = cbind(des.con, bat = factor(as.numeric(cli$TMT.plex)))}}
      daf = (mut>0)*1
      # regression
      res = c(); genes = c(); transcriptome = c() 
      gene = 'PIK3CA'
      for (gene in gen) {
        trans = as.numeric(rna[match(gene, rownames(rna)), ])
        isna = is.na(trans) | rowSums(is.na(des.con))>0
        trans = trans[!isna]
        des.mat = model.matrix(~trans)
        des.mat = cbind(des.mat, des.con[!isna,])
        
        computed = tryCatch({ 
          fit = lmFit(as.numeric(daf[gene,!isna]), des.mat)
          fit2 = eBayes(fit)
          tmp = topTable(fit2, coef='trans', adjust.method='BH', number=nrow(daf), p.value=1.1, sort.by = 'none')  
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
        xlsx::write.xlsx(cbind(Gene = genes, res[,-5], RNA = transcriptome), paste0('../data/DNA_RNA_regression/Table.DNA.RNA.regression.linearLIMMA.MutVsRNA.',cancer,'.',mut.typ,'.xlsx'), row.names = F)
      }
    } else if (analysis=='LinearRNA~Mut') { 
      mut = (mut>0)*1 
      # confounding covariates
      des.con = c()
      if ('Age.in.months.at.sampling' %in% colnames(cli)){
        if (length(levels(factor(cli$Age.in.months.at.sampling))) > 1) {
          des.con = cbind(des.con, age = as.numeric(cli$Age.in.months.at.sampling)/12)}}
      if ('Gender' %in% colnames(cli)){
        if (length(levels(factor(cli$Gender))) > 1) {
          des.con = cbind(des.con, sex = factor(as.character(cli$Gender)))}}
      if ('Ethnicity' %in% colnames(cli)){
        if (length(levels(factor(cli$Ethnicity))) > 1) {
          des.con = cbind(des.con, eth = factor(as.character(cli$Ethnicity)))}}
      if ('TMT.plex' %in% colnames(cli)){ 
        if (length(levels(factor(cli$TMT.plex))) > 1) {
          des.con = cbind(des.con, bat = factor(as.numeric(cli$TMT.plex)))}}
      daf = rna
      # regression
      res = c(); genes = c(); transcriptome = c()
      gene = 'PIK3CA'
      for (gene in gen) {
        mutation = as.numeric(mut[gene,])
        if (sum(mutation) < 3)
          next
        isna = is.na(mutation) | rowSums(is.na(des.con))>0
        mutation = mutation[!isna]
        des.mat = model.matrix(~mutation)
        des.mat = cbind(des.mat, des.con[!isna,])
        computed = tryCatch({ 
          fit = lmFit(as.numeric(daf[match(gene, rownames(rna)), !isna]), des.mat)
          # logFC.unsorted = fit$coefficients[,which(colnames(des.mat)=='mutation')]
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
        xlsx::write.xlsx(cbind(Gene = genes, res[,-5], RNA = transcriptome), paste0('../data/DNA_RNA_regression/Table.DNA.RNA.regression.linearLIMMA.RNAVsMut.',cancer,'.',mut.typ,'.xlsx'), row.names = F)
      }
      # _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
    } # analysis
    
  }
}

# put all data together 
cancer = 'BRCA'
mut.typ = 'missense' 
df = data.frame()
for (cancer in pan.coh) {
  for (mut.typ in mut.types) {
    dat <- read_excel(paste0('../data/DNA_RNA_regression/Table.DNA.RNA.regression.linearLIMMA.RNAVsMut.',cancer,'.',mut.typ,'.xlsx'))
    # dat = dat[dat$FDR <= 0.05,]
    dat$cancer = cancer
    dat$mutation = mut.typ
    df = rbind(df, dat)
    write.csv(df, '../data/DNA_RNA_regression/Table.DNA.RNA.regression.linearLIMMA.RNAVsMut.csv')
  }
}
