# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#             Format WXS data downloaded from NIH GDC
#
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
rm(list = ls(all.names = TRUE))
setwd('~/Library/CloudStorage/Box-Box/Huang_lab/Huang_lab_data/CPTAC/WXS/')

workflow = 'SomaticVariant'
likelyFunctionalTypes = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Nonsense_Mutation","Splice_Site","Start_Codon_Del","Start_Codon_Del")
cancer.types = c('BRCA','CCRCC','CRC','LUAD','OV','UCEC')

cancer = 'UCEC'
for (cancer in cancer.types){
  # Load data
  map = read.table(paste(workflow, cancer, 'gdc_sample_sheet.tsv', sep = '/'), sep = '\t', header = T)
  map$Case.ID = as.factor(apply(as.matrix(map$Case.ID), 1, function(x){gsub('-','.',strsplit(as.character(x), split = ',')[[1]][1])}))
  table(map$Sample.Type)
  map.Sample.Type = array(NA,nrow(map))
  map.Sample.Type[which(regexpr('Primary Tumor',map$Sample.Type)>0)] = 'Tumor' #***
  any(!map.Sample.Type=='Tumor')
  
  sam.Hugo_Symbol = c()
  for (i in seq_len(nrow(map))) {
    if (file.exists(paste(workflow, cancer, 'data', map$File.ID[i], map$File.Name[i], sep = '/'))){
      loaded = F; trial = 1
      while (!loaded & trial < 5){
        loaded = tryCatch({mut = read.table(paste(workflow, cancer, 'data', map$File.ID[i], map$File.Name[i], sep = '/'), sep = '\t', header = T, quote = "", fill=F); T}
                          , error = function(e) F, finally = {})
        sam.Hugo_Symbol = c(sam.Hugo_Symbol, gsub( '\\..*$', '', mut$Hugo_Symbol))
        trial = trial + 1
      }
      if (!loaded) {message("File not loaded.")} else {
        message('[', trial, '] File missing: ', paste(workflow, cancer, 'data', map$File.ID[i], map$File.Name[i], sep = '/'))
      }
    } 
  }
  sam.Hugo_Symbol = unique(sam.Hugo_Symbol)
  
  tru = setNames(data.frame(matrix(0,length(sam.Hugo_Symbol),nrow(map)), row.names = sam.Hugo_Symbol), unique(gsub('-','.',paste(map$Case.ID, map.Sample.Type, sep='.'))))
  mis = tru; 
  syn = tru;
  
  
  # snv = table(alt$Gene.symbol, alt$Sample.ID)
  # Frame data
  i=3
  truHGVSps = c(); truHugos = c()
  misHGVSps = c(); misHugos = c()
  synHGVSps = c(); synHugos = c()
  for (i in seq_len(nrow(map))) {
    loaded = F; trial = 1
    if (file.exists(paste(workflow, cancer, 'data', map$File.ID[i], map$File.Name[i], sep = '/'))){
        while (!loaded & trial < 5){
          loaded = tryCatch({mut = read.table(paste(workflow, cancer, 'data', map$File.ID[i], map$File.Name[i], sep = '/'), sep = '\t', header = T, quote = "", fill=F); T}
                            , error = function(e) F, finally = {})
          trial = trial + 1
        } 
        mut$Hugo_Symbol = gsub( '\\..*$', '', mut$Hugo_Symbol)
        sam.idx = which(colnames(tru)==paste(map$Case.ID[i], map.Sample.Type[i], sep='.'))
        # Truncating mutations      
        idx = which(mut$Variant_Classification %in% likelyFunctionalTypes)
        for (j in idx){
          tru[mut$Hugo_Symbol[j], sam.idx] = mut$HGVSp_Short[j]
        }
        # Missense mutations
        idx = which(mut$Variant_Classification %in% 'Missense_Mutation')
        for (j in idx){
          mis[mut$Hugo_Symbol[j], sam.idx] = mut$HGVSp_Short[j]
        }
        # synonymous mutations
        idx = which(mut$Variant_Classification %in% 'Silent')
        for (j in idx){
          syn[mut$Hugo_Symbol[j], sam.idx] = mut$HGVSp_Short[j]
        }
        loaded = T
        if (!loaded) {message("File not loaded.")} else {
          trial = trial + 1
          message('[', trial, '] File missing: ', paste(workflow, cancer, 'data', map$File.ID[i], map$File.Name[i], sep = '/'))
        }
      } 
  }
 
  write.csv(tru, paste0('../../../manuscripts/pQTL/data/HGVSp/',cancer,'TruncatingHGVSp.csv'))
  write.csv(mis, paste0('../../../manuscripts/pQTL/data/HGVSp/',cancer,'MissenseHGVSp.csv'))
  write.csv(syn, paste0('../../../manuscripts/pQTL/data/HGVSp/',cancer,'SynonymousHGVSp.csv'))
}


