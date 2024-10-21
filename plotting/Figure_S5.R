# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#     Compute the correlation between paired mRNA and protein expressions,
#       to generate Figure S5 and Table S7
#
#     Within cancer, analyze gene/protein expressions:
#       For each gene, calculate correlation coefficient b/w paired gene and protein expressions across the same sample set
#
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
rm(list = ls(all.names = TRUE))
setwd('~/Documents/workspace/pQTL/analysis/')

# Load data set ----
cancer.types = c('BRCA','CCRCC','CRC','LUAD','OV','UCEC')
rna = lapply(cancer.types, function(cancer){read.table(paste0('../../PanCancerProteomicsData_HuangLab/', cancer,'.transcriptome.FPKM.formatted.txt.gz'), sep = '\t')})
pro = lapply(cancer.types, function(cancer){read.table(paste0('../../PanCancerProteomicsData_HuangLab/', cancer,'.proteome.formatted.normalized.txt.gz'), sep = '\t', header = T)})
names(rna) = cancer.types
names(pro) = cancer.types
for (cancer in cancer.types) {colnames(rna[[cancer]]) = paste(cancer, gsub('^X','',colnames(rna[[cancer]])), sep = '.')}
for (cancer in cancer.types) {colnames(pro[[cancer]]) = paste(cancer, gsub('^X','',colnames(pro[[cancer]])), sep = '.')}

# Format ----
colnames(pro[['CRC']]) = sub('\\.Tissue\\..*$','',colnames(pro[['CRC']]))

# Compute ----
res = c()
cancer = 'LUAD' 
for (cancer in cancer.types){
  introw = intersect(rownames(rna[[cancer]]), rownames(pro[[cancer]]))
  intcol = intersect(colnames(rna[[cancer]]), colnames(pro[[cancer]]))
  gene = 'TP53'
  for (gene in introw){
    ct = tryCatch({cor.test(as.numeric(rna[[cancer]][gene,intcol]), as.numeric(pro[[cancer]][gene,intcol]), na.action('omit'))}, error=function(e){NA}, finally={})
    if (!any(is.na(ct))){
      res = rbind(res, data.frame(Cancer = cancer, Gene = gene, Coef = ct$estimate, P = ct$p.value, row.names = gene))
    }
  }
}

# Save ----
saveRDS(res, 'Genewise.correlations.Transcriptome.FPKM.formatted.Vs.Proteome.formatted.normalized.RDS')





# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#     Plot the mRNA-protein correlation coefficients, for Concordant/Discordant/Other genes in different cancer types
#     (Figure S5, Table S7)
#
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
rm(list = ls(all.names = TRUE))
setwd('~/Documents/workspace/pQTL/analysis/')
library(ggplot2)
library(ggrepel)

# Load concordant and discordant (spsQTL) gene sets ----
con = as.data.frame(readxl::read_xlsx('../doc/Supp Table/TableS3.concordant_e.pQTLs.xlsx', sheet = 1))
dis = as.data.frame(readxl::read_xlsx('../doc/Supp Table/TableS4.significant_spsQTLs.xlsx', sheet = 1))
con[paste(con$Gene, con$cancer) %in% intersect(paste(con$Gene, con$cancer), paste(dis$Gene, dis$cancer)),]
dis[paste(dis$Gene, dis$cancer) %in% intersect(paste(con$Gene, con$cancer), paste(dis$Gene, dis$cancer)),]


# Load correlation results ----
res = readRDS('Genewise.correlations.Transcriptome.FPKM.formatted.Vs.Proteome.formatted.normalized.RDS')
res$Category = 'Other'
res$Category[paste(res$Gene, res$Cancer) %in% paste(con$Gene, con$cancer)] = 'Concordant'
res$Category[paste(res$Gene, res$Cancer) %in% paste(dis$Gene, dis$cancer)] = 'Discordant'
col.val = c(setNames(c('blue', 'red'), setdiff(unique(res$Category), 'Other')), 'Other' = 'grey')

# Plot Figure S5 ----
pri = ggplot(res, aes(x = Category, y = Coef, color = Category)) +
  geom_violin(trim = F) +
  geom_jitter(aes(color = Category), width = 0.1, alpha = 0.3, size = .8) +
  geom_text_repel(aes(label = ifelse(Gene %in% c(con$Gene,dis$Gene) & Category != 'Other' & (abs(Coef)>.65 | !Cancer %in% c('CRC','UCEC')), Gene, NA)), min.segment.length = 0.05, max.iter = 10e6, box.padding = 0.05, max.overlaps = 20, max.time = 10) +
  facet_wrap(~ Cancer, ncol=3) +
  scale_color_manual(values = col.val) +
  theme_bw() +
  labs(x = "Category",
       y = "Coefficient",
       fill = "Category")
png(paste0('CorCoefs.ConcordantvsDiscordant.png'), w=8, h=5, units = 'in', res = 300); print(pri); dev.off()
# pdf(paste0('CorCoefs.ConcordantvsDiscordant.pdf'), w=8, h=5, useDingbats = F); print(pri); dev.off()

# Write Table S7 ----
xlsx::write.xlsx(res[!res$Category %in% 'Other',], 'CorCoefs.ConcordantvsDiscordant.xlsx', row.names = F)

# Report statistical significance of difference in corr_coef b/w groups ----
cancer = 'UCEC'
for (cancer in sort(unique(res$Cancer))){
  res.can = res[res$Cancer %in% cancer,]
  res.tmp = res.can; res.tmp$Category = sub('Discordant','Other',res.tmp$Category); rep.wil = wilcox.test(Coef ~ Category, data=res.tmp) # con vs all other
  message('Concordant genes vs all others [', cancer, ']: W = ', rep.wil$statistic, ', p = ', format.pval(rep.wil$p.value, eps = 1e-100, digits = 3))
  res.tmp = res.can; res.tmp$Category = sub('Concordant','Other',res.tmp$Category); rep.wil = wilcox.test(Coef ~ Category, data=res.tmp) # dis vs all other
  message('Discordant genes vs all others [', cancer, ']: W = ', rep.wil$statistic, ', p = ', format.pval(rep.wil$p.value, eps = 1e-100, digits = 3))
}

