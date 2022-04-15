# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
# 
#   
#
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
rm(list = ls(all.names = TRUE))
setwd('~/Library/CloudStorage/Box-Box/Huang_lab/manuscripts/pQTL/analysis')
library(readxl)
library(lmtest)

# pan.coh='BRCA' #***
pan.coh = c('BRCA', 'CRC', 'CCRCC', 'LUAD', 'OV', 'UCEC')
mut.types = c('missense', 'truncating', 'synonymous')
cancer = 'LUAD'
for (cancer in pan.coh){
  message('Processing ', cancer, ' ...')
  
  rna0 = read.table(paste0('../../../Huang_lab_data/PanCancerProteomicsData_HuangLab/',cancer,'.transcriptome.FPKM.formatted.txt.gz'))
  pro0 = read.table(paste0('../../../Huang_lab_data/PanCancerProteomicsData_HuangLab/',cancer,'.proteome.formatted.normalized.tumor.txt.gz'))
  
  if (!any(!paste0('V',1:ncol(pro0)) %in% colnames(pro0))){
    colnames(pro0) = apply(pro0[1,], 1, as.character); pro0 = pro0[-1,]
  }
  
  mut.typ = 'missense' 
  for (mut.typ in mut.types){
    rna = log2(1+rna0);
    pro = pro0; 
    mut = read.table(paste0('../../../Huang_lab_data/PanCancerProteomicsData_HuangLab/',cancer,'.WXS.SomaticVariant.',mut.typ,'.txt.gz'))
    # align 
    colnames(mut) = gsub('\\.Tumor.*$','',gsub('^X','',colnames(mut)))
    colnames(rna) = gsub('\\.Tumor.*$','',gsub('^X','',colnames(rna)))
    colnames(pro) = gsub('\\.Tumor.*$','',gsub('^X','',colnames(pro)))
    sam = intersect(colnames(pro), colnames(rna))
    sam = intersect(colnames(mut), sam)
    gen = intersect(rownames(mut), rownames(pro))
    gen = intersect(gen, rownames(rna))
    mut = mut[gen, sam]
    rna = rna[gen, sam]
    pro = pro[gen, sam]
    gen = intersect(
      rownames(mut)[which(rowSums(mut!=0, na.rm = T) > 0)],
      rownames(pro)[which(rowSums(pro!=0, na.rm = T) > 0)])
    gen = intersect(gen,rownames(rna)[which(rowSums(rna!=0, na.rm = T) > 0)])
    mut = mut[gen, ]
    rna = rna[gen, ]
    pro = pro[gen, ]
    
      mut = (mut>0)*1 
      daf = pro
      # regression
      res = c(); genes = c(); lrts = c()
      gene = 'MUC16'
      for (gene in gen) {
        trans = as.numeric(rna[match(gene, rownames(rna)), ])
        isna = is.na(trans) 
        trans = trans[!isna]
        des.mat = model.matrix(~trans)
        
        lm1 = lm((as.numeric(daf[match(gene, rownames(pro)), !isna]) ~ des.mat))
        
        mutation = as.numeric(mut[gene,])
        isna = is.na(mutation)
        mutation = mutation[!isna]
        des.mat = cbind(des.mat, mutation)
        
        if (sum(mutation) < 3)
          next
        
        lm2 = lm((as.numeric(daf[match(gene, rownames(pro)), !isna]) ~ des.mat))
        
        lrt = anova(lm1, lm2, test = "LRT")
        # lrt = lrtest(lm1, lm2)
        lrts = rbind(lrts, lrt$`Pr(>Chi)`[2])
        
        genes = c(genes, gene)
      }
      xlsx::write.xlsx(cbind(Gene = genes, lrt = lrts), paste0('../data/lrt/ProVsRNALrt',cancer,'.',mut.typ,'.xlsx'), row.names = F)
    
  }
}

lrt = c()
cancer = 'BRCA'
mut.typ = 'missense'
for (cancer in pan.coh) {
  for (mut.typ in mut.types){
    dat <- read_excel(paste0('../data/lrt/ProVsRNALrt',cancer,'.',mut.typ,'.xlsx'))
    colnames(dat)[2] = 'p.value'
    dat = dat[!is.na(dat$p.value),]
    dat$p.value = as.numeric(dat$p.value)
    dat_pvalue_0.1 = dat[dat$p.value < 0.05, ]
    dat_pvalue_0.1$cancer = cancer
    dat_pvalue_0.1$mutation = mut.typ
    dat_pvalue_0.1$FDR = p.adjust(dat_pvalue_0.1$p.value, method = 'BH')
    lrt = rbind(lrt, dat_pvalue_0.1)
    }
}
xlsx::write.xlsx(lrt, paste0('../data/lrt/ProVsRNALrt.xlsx'))
dat <- read_excel('../data/lrt/ProVsRNALrt.xlsx')
write.csv(dat, '../data/lrt/ProVsRNALrt.csv')

