library(data.table)
library(tidyverse)
setwd('../code/')

#Read gw significant snps
loci <- fread("../Output/FUMA_cQTL/snps.txt")
loci <- loci[!is.na(loci$gwasP),]
#Only loci with MAF > 0.05 were included in 500FG
loci <- loci[loci$MAF > 0.05,]
  

f500 <- list.files("../500FG_cQTL/",pattern = ".tsv")%>%
  .[grepl("IL1b|IL.10|IL6|IL.1Ra", .)]

results <- lapply(f500, function(name){
	print(name)
  path <- paste0("../500FG_cQTL/",name)
  #Read file
  gwas <- fread(path)
  #Adjust colnames
  if(ncol(gwas) == 11)  colnames(gwas) <- c("SNP", "EffectAllele", "beta", "p.value", "CHROM", "POS", "gene", "t.stat", "FDR", "CHROM2","AlternativeAllele")
  gwas <- dplyr::rename(gwas, "chr" = CHROM, "pos" = POS,"p.gwas" = p.value)
  loci_gwas <- left_join(loci, gwas, by = c("chr","pos"))
  #Get minimum p-value per locus
  loci_gwas <- loci_gwas[,c("GenomicLocus","rsID","p.gwas")]
  loci_gwas <- loci_gwas %>% group_by(GenomicLocus) %>% summarise(p.gwas = ifelse(all(is.na(p.gwas)), NA, min(p.gwas, na.rm = T)),
                                                                  SNP = ifelse(all(is.na(p.gwas)), NA, rsID[which.min(p.gwas)]))
  gwas.p <- loci_gwas$p.gwas
  names(gwas.p) <- loci_gwas$GenomicLocus
  return(gwas.p)
})
names(results) <- gsub("_hg\\d+.tsv","", f500)
results <- as.data.frame(results)
fwrite(results,"../Output/GWASCatalogue/loci_gwas_500FG.tsv", sep = "\t", row.names = T)


results <- fread("../Output/GWASCatalogue/loci_gwas_500FG.tsv")%>%
  as.data.frame()%>%column_to_rownames("V1")

#Replication rate at different p-value thresholds

apply(results, 1, function(x) sum(na.omit(x)<1e-5))

#Subset only 24h measurements and cytokines that are present

results <- as.data.frame(results)%>%
  .[,grepl("24h", colnames(.))]%>%
  .[,grepl("IL", colnames(.))]%>%
  .[,!grepl("IL.8", colnames(.))]%>%
  .[,grepl("PBMC", colnames(.))]


#All pvalues below 0.05, convert to NA. This is already like that in all cytokines except IL10
results <- apply(results, 2, function(x) ifelse(x>0.05, NA, x))

#Fill in missing loci
if(!all(unique(loci$GenomicLocus) %in% rownames(results))) stop('Some loci do not have a nominal hit in 500FG, add them to matrix')

#Heatmap
mat <- -log10(as.matrix(results))

#Discard traits without significant assoc
mat = mat[,apply(mat, 2, function(x) any(!is.na(x) & x > -log10(0.05)))]
thr = 0.05/nrow(mat)/ncol(mat)
thr = 1e-4
sig <- ifelse(mat > -log10(5e-8), "***",
              ifelse(mat > -log10(thr), "**",
                     ifelse(mat > -log10(0.05), "*", "")))

sig[is.na(sig)] <- ""
colnames(mat) <- gsub('_24h.tsv', '', colnames(mat))
pdf("../code_figures/Figures/Suppl.loci_gwas_500FG_heatmap.pdf", width =12, height = 8)
ComplexHeatmap::Heatmap(mat,rect_gp = grid::gpar(col = "black", lwd = .5),
                        cluster_rows = F, show_row_dend = F, show_column_dend = F, cluster_columns = F,
                        row_names_side = "left",column_names_side = "top",show_row_names = T,
                        col = circlize::colorRamp2(c(  0, #-log10(0.05)-0.0001 ,
                                                       -log10(0.05),-log10(0.05/(nrow(loci)*length(f500))), -log10(5e-8)),
                                                   c( "white","#F39B7FFF","#E64B35FF", "darkred")),
                        na_col = "white",
                        cell_fun = function(j,i,x,y,w,h,col)grid::grid.text(sig[i,j],x,y),
                        name = "-log10 p-value")
dev.off()
