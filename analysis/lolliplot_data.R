rm(list = ls(all.names = TRUE))
setwd('~/Library/CloudStorage/Box-Box/Huang_lab/manuscripts/pQTL/analysis')

path = '../../../Huang_lab_data/CPTAC/WXS/SomaticVariant'
cancer.types = c('BRCA','CCRCC','CRC','LUAD','OV','UCEC')

cancer = 'CRC'
for (cancer in cancer.types){
  data = c()
  
  map = read.table(paste(path, cancer, 'gdc_sample_sheet.tsv', sep = '/'), sep = '\t', header = T)
  map$Case.ID = as.factor(apply(as.matrix(map$Case.ID), 1, function(x){gsub('-','.',strsplit(as.character(x), split = ',')[[1]][1])}))
  table(map$Sample.Type)
  map.Sample.Type = array(NA,nrow(map))
  map.Sample.Type[which(regexpr('Primary Tumor',map$Sample.Type)>0)] = 'Tumor' #***
  any(!map.Sample.Type=='Tumor')
  
  sam.Hugo_Symbol = c()
  for (i in seq_len(nrow(map))) {
    if (file.exists(paste(path, cancer, 'data', map$File.ID[i], map$File.Name[i], sep = '/'))){
      loaded = F; trial = 1
      while (!loaded & trial < 5){
        loaded = tryCatch({mut = read.table(paste(path, cancer, 'data', map$File.ID[i], map$File.Name[i], sep = '/'), sep = '\t', header = T, quote = "", fill=F); T}
                          , error = function(e) F, finally = {})
        sam.Hugo_Symbol = c(sam.Hugo_Symbol, gsub( '\\..*$', '', mut$Hugo_Symbol))
        trial = trial + 1
      } 
      if (!loaded) {message("File not loaded.")}
    } else {
      message('[', trial, '] File missing: ', paste(path, cancer, 'data', map$File.ID[i], map$File.Name[i], sep = '/'))
    }
  }
  sam.Hugo_Symbol = unique(sam.Hugo_Symbol)
  
  tru = setNames(data.frame(matrix(0,length(sam.Hugo_Symbol),nrow(map)), row.names = sam.Hugo_Symbol), gsub('-','.',paste(map$Case.ID, map.Sample.Type, sep='.')))
  mis = tru; nmi = tru
  
  i=55
  for (i in seq_len(nrow(map))) {
    loaded = F; trial = 1
    while (!loaded & trial < 5){
      if (file.exists(paste(path, cancer, 'data', map$File.ID[i], map$File.Name[i], sep = '/'))){
        loaded = F; trial = 1
        while (!loaded & trial < 5){
          loaded = tryCatch({mut = read.table(paste(path, cancer, 'data', map$File.ID[i], map$File.Name[i], sep = '/'), sep = '\t', header = T, quote = "", fill=F); T}
                            , error = function(e) F, finally = {})
          trial = trial + 1
        }
        gene = gsub('\\..*$', '', mut$Hugo_Symbol)
        refseq = sub("\\..*", "", mut$RefSeq)
        chromosome = mut$Chromosome
        start = mut$Start_Position
        aachange = mut$Amino_acids
        sample = mut$Tumor_Sample_Barcode
        class = gsub('Missense_Mutation', 'missense', mut$Variant_Classification)
        class = gsub('Frame_Shift_Del', 'frameshift', class)
        class = gsub('Frame_Shift_Ins', 'frameshift', class)
        class = gsub("5'UTR", 'utr_5', class)
        class = gsub("3'UTR", 'utr_3', class)
        class = gsub('In_Frame_Del', 'proteinDel', class)
        class = gsub('Nonsense_Mutation', 'nonsense', class)
        class = gsub('In_Frame_Ins', 'proteinIns', class)
        dat = cbind(gene, refseq, chromosome, start, aachange, sample, class)
        data = rbind(data, dat)
      } else {
        trial = trial + 1
        message('[', trial, '] File missing: ', paste(workflow, cancer, 'data', map$File.ID[i], map$File.Name[i], sep = '/'))
      }
    }
  }
        write.table(data, paste0("../data/lolliplot/", cancer, "_lolliplot_data.txt"), row.names=FALSE,sep="\t", quote = FALSE)
}

for (cancer in cancer.types) {
  dat = read.table(paste0("../data/lolliplot/", cancer, "_lolliplot_data.txt"), fill = TRUE, header = TRUE, sep="\t")
  df = paste(dat$chromosome, dat$start, dat$start)
  write.table(df, paste0('../data/lolliplot/', cancer, '_chr.txt'), row.names=FALSE, col.names = F, sep="\t", quote = FALSE)
}

cancer = 'UCEC'
dat = read.table(paste0("../data/lolliplot/", cancer, "_lolliplot_data.txt"), fill = TRUE, header = TRUE, sep="\t")
dat = dat[!(dat$start == 149583933),]
bed = read.table(paste0("../data/lolliplot/", cancer, ".bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
dat$start = bed$V2
write.table(dat, paste0('../data/lolliplot/', cancer, '_lolliplot_data.txt'), row.names=FALSE,sep="\t", quote = FALSE)

