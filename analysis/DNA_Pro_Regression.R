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
# analysis = 'LinearMut~Pro'
analysis = 'LinearPro~Mut'

# pan.coh='BRCA' #***
pan.coh = c('BRCA', 'CRC', 'CCRCC', 'LUAD', 'OV', 'UCEC')
mut.types = c('missense', 'truncating', 'synonymous')

cancer = 'BRCA'
for (cancer in pan.coh){
  message('Processing ', cancer, ' ...')
  
  pro0 = read.table(paste0('../../../Huang_lab_data/PanCancerProteomicsData_HuangLab/',cancer,'.proteome.formatted.normalized.tumor.txt.gz'))
  cli0 = read.table(paste0('../../../Huang_lab_data/PanCancerProteomicsData_HuangLab/',cancer,'.clinical.formatted.txt'), sep = '\t', header = T)
  
  if (!any(!paste0('V',1:ncol(pro0)) %in% colnames(pro0))){
    colnames(pro0) = apply(pro0[1,], 1, as.character); pro0 = pro0[-1,]
  }
  
  mut.typ = 'missense' 
  for (mut.typ in mut.types){
    pro = pro0; cli = cli0
    mut = read.table(paste0('../../../Huang_lab_data/PanCancerProteomicsData_HuangLab/',cancer,'.WXS.SomaticVariant.',mut.typ,'.txt.gz'))
    # align 
    colnames(mut) = gsub('\\.Tumor.*$','',gsub('^X','',colnames(mut)))
    colnames(pro) = gsub('\\.Tumor.*$','',gsub('^X','',colnames(pro)))
    sam = intersect(intersect(cli$Participant.ID, colnames(mut)), colnames(pro))
    gen = intersect(rownames(mut), rownames(pro))
    mut = mut[gen, sam]
    pro = pro[gen, sam]
    cli = cli[match(sam, cli$Participant.ID),]
    gen = intersect(
      rownames(mut)[which(rowSums(mut!=0, na.rm = T) > 0)],
      rownames(pro)[which(rowSums(pro!=0, na.rm = T) > 0)])
    mut = mut[gen, ]
    pro = pro[gen, ]
    
    
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
      res = c(); genes = c(); proteins = c() 
      gene = 'TPCN1'
      for (gene in gen ) {
        daf = data.frame(mut = as.numeric(mut[gene,]),
                         pro = as.numeric(pro[gene,]))
        formula = 'mut ~ pro'
        # if (nrow(rna) > 0) {daf$rna = as.numeric(rna[gene,]); formula = paste(formula, '+ rna')}
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
        if (computed & ('pro' %in% rownames(reg))) {
          res = rbind(res, reg['pro',])
          genes = c(genes, gene)
          proteins = c(proteins, rownames(pro)[match(gene, rownames(pro))])
        }
      } # for gene
      res = data.frame(res, row.names = genes)
      colnames(res) = c('Beta', 'Std.Error', 'Z.value', 'P.value')
      res$FDR = p.adjust(res$P.value, method = 'BH')
      rankd = order(-log10(res$FDR)-log10(res$P.value), decreasing = T, na.last = T)
      res = res[rankd,]; genes = genes[rankd]; proteins = proteins[rankd] 
      valid = !is.na(res$FDR) & !is.na(res$P.value)
      res = res[valid,]; genes = genes[valid]; proteins = proteins[valid] 
      if (any(valid)){
        xlsx::write.xlsx(cbind(Gene = genes, Protein = proteins, res), paste0('../data/DNA_Pro_regression/Table.DNA.PRO.regression.Poisson.',cancer,'.',mut.typ,'.xlsx'), row.names = F)
      }
      # _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
      
      
    }  else if (analysis=='LinearMut~Pro') {    
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
      res = c(); genes = c(); proteins = c() 
      gene = 'PIK3CA'
      for (gene in gen) {
        prot = as.numeric(pro[match(gene, rownames(pro)), ])
        isna = is.na(prot) | rowSums(is.na(des.con))>0
        prot = prot[!isna]
        des.mat = model.matrix(~prot)
        des.mat = cbind(des.mat, des.con[!isna,])
        
        computed = tryCatch({ 
          fit = lmFit(as.numeric(daf[gene,!isna]), des.mat)
          # logFC.unsorted = fit$coefficients[,which(colnames(des.mat)=='phospho')]
          fit2 = eBayes(fit)
          tmp = topTable(fit2, coef='prot', adjust.method='BH', number=nrow(daf), p.value=1.1, sort.by = 'none')  
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
        xlsx::write.xlsx(cbind(Gene = genes, res[,-5], Protein = proteins), paste0('../data/DNA_Pro_regression/Table.DNA.PRO.regression.linearLIMMA.MutVsPro.',cancer,'.',mut.typ,'.xlsx'), row.names = F)
      }
    } else if (analysis=='LinearPro~Mut') { 
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
      daf = pro
      # regression
      res = c(); genes = c(); proteins = c() 
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
          fit = lmFit(as.numeric(daf[match(gene, rownames(pro)), !isna]), des.mat)
          # fit = lmFit(as.numeric(daf[match(gene, pho.gen), !isna]), des.mat)
          # logFC.unsorted = fit$coefficients[,which(colnames(des.mat)=='mutation')]
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
        xlsx::write.xlsx(cbind(Gene = genes, res[,-5], Protein = proteins), paste0('../data/DNA_Pro_regression/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.',cancer,'.',mut.typ,'.xlsx'), row.names = F)
      }
      # _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
    } # analysis
    
  }
}

# put all data together and select significant data
df = data.frame()
for (cancer in pan.coh) {
  for (mut.typ in mut.types) {
    dat <- read_excel(paste0('../data/DNA_Pro_regression/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.',cancer,'.',mut.typ,'.xlsx'))
    # dat = dat[dat$FDR <= 0.05,]
    dat$cancer = cancer
    dat$mutation = mut.typ
    df = rbind(df, dat)
    write.csv(df, '../data/DNA_Pro_regression/Table.DNA.PRO.regression.linearLIMMA.ProVsMut.csv')
  }
}


