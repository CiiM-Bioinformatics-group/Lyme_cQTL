library(data.table)
library(tidyverse)

#Read gw significant snps
loci <- fread("../Output/FUMA_cQTL/snps.txt")
loci <- loci[!is.na(loci$gwasP),]

fgwas <- list.files("../Output/GWASCatalogue/",pattern = ".gz")

#remove IL10 IL6 and DRAM1
fgwas <- fgwas[!grepl("IL10|IL6|DRAM1", fgwas)]


results <- lapply(fgwas, function(name){
  path <- paste0("../Output/GWASCatalogue/",name)
  #Read file
  gwas <- fread(path)
  #Quick fix for finngen
  if(grepl("finngen", name)){
    colnames(gwas) <- c("chromosome", "base_pair_location", "ref", "alt", "rsids", "nearest_genes", "p_value", "mlogp", "beta",
                                                 "sebeta", "af_alt", "af_alt_cases", "af_alt_controls")
    name <- paste0(name, "_", "hg38")}
  #Adjust colnames
  gwas <- dplyr::rename(gwas, "chr" = chromosome, "pos" = base_pair_location,"p.gwas" = p_value)
  #Liftover function
  gwas_liftover <- function(gwas){
    #Liftover positions to hg38 (https://master.bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html)
    library(rtracklayer)
    library(GenomicRanges)
    path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
    ch = import.chain(path)
    gwas = makeGRangesFromDataFrame(gwas, keep.extra.columns = T, seqnames.field = "chr", start.field = "pos",
                                    end.field = "pos")
    seqlevelsStyle(gwas) = "UCSC"
    gwas = liftOver(gwas, ch)
    gwas = unlist(gwas)
    genome(gwas) = "hg19"
    gwas <- data.frame(gwas)
    gwas <- gwas[gwas$seqnames != "chrX" & gwas$seqnames != "chrY",]
    gwas$chr <- as.integer(gsub("chr","",gwas$seqnames))
    gwas$pos <- gwas$start
    return(gwas)
  }
  #Liftover those from hg38
  if(grepl("hg38",name)) gwas <- gwas_liftover(gwas)
  #Extract p-value per loci
  gwas$chr <- as.integer(gwas$chr)

  #sapply(1:nrow(loci), function(n){
  #  snp = paste0(loci$chr[n], ":", loci$pos[n])
  #  matched = SNPsnap[Input_SNP == snp]
  #  matched.p <- sapply(matched, function(x){min(gwas$p.gwas[paste0(gwas$chr,":",gwas$pos) == x])})
  #})

  loci_gwas <- left_join(loci, gwas, by = c("chr","pos"))
  #Get only variant with the highest MAF in loci that have more than 1 allele
  #if(any(duplicated(loci_gwas$GenomicLocus))){
    #Which locus has duplicated
  #  dup <- which(loci_gwas$GenomicLocus == loci_gwas$GenomicLocus[duplicated(loci_gwas$GenomicLocus)])
    #Which is the allele with the lower MAF, to discard
  #  min <- which.min(loci_gwas$effect_allele_frequency[dup])
  #  dis <- seq(1,nrow(loci_gwas))[dup][min]
    #Remove that one
  #  loci_gwas <- loci_gwas[-dis,]
  #}

  #Get minimum p-value per locus
  loci_gwas <- loci_gwas[,c("GenomicLocus","rsID","p.gwas")]
  loci_gwas <- loci_gwas %>% group_by(GenomicLocus) %>% summarise(p.gwas = ifelse(all(is.na(p.gwas)), NA, min(p.gwas, na.rm = T)),
                                                                  SNP = ifelse(all(is.na(p.gwas)), NA, rsID[which.min(p.gwas)]))
  gwas.p <- loci_gwas$p.gwas
  names(gwas.p) <- loci_gwas$GenomicLocus
  return(gwas.p)
})
names(results) <- gsub("_hg\\d+.tsv.gz","", fgwas)
results <- as.data.frame(results)
fwrite(results,"../Output/GWASCatalogue/loci_gwas.tsv", sep = "\t", row.names = T)


results <- fread("../Output/GWASCatalogue/loci_gwas.tsv")%>%
  as.data.frame()%>%column_to_rownames("V1")

mat <- -log10(as.matrix(results))
rownames(mat) <- str_split(rownames(mat), pattern = " ", simplify = T)[,1]
#Some slight changes
mat <- mat[,!grepl("BMI|finngen", colnames(mat))]
colnames(mat) <- gsub("_", " ", colnames(mat))
#Discard traits without significant assoc
mat = mat[,apply(mat, 2, function(x) any(!is.na(x) & x > -log10(0.05)))]
thr = 0.05/nrow(mat)/ncol(mat)
thr = 1e-4
sig <- ifelse(mat > -log10(5e-8), "***",
              ifelse(mat > -log10(thr), "**",
                     ifelse(mat > -log10(0.05), "*", "")))



sig[is.na(sig)] <- ""
pdf("../code_figures/Figures/loci_gwas_heatmap.pdf", width =6, height = 12)
ComplexHeatmap::Heatmap(mat,rect_gp = grid::gpar(col = "black", lwd = .5),
                        cluster_rows = F, show_row_dend = F, show_column_dend = F, cluster_columns = F,
                        row_names_side = "left",column_names_side = "top",show_row_names = T,
                        col = circlize::colorRamp2(c(  0, #-log10(0.05)-0.0001 ,
                                                       -log10(0.05),-log10(0.05/(nrow(loci)*length(fgwas))), -log10(5e-8)),
                                         c( "white","#F39B7FFF","#E64B35FF", "darkred")),
                        na_col = "white",
                        cell_fun = function(j,i,x,y,w,h,col)grid::grid.text(sig[i,j],x,y),
                        name = "-log10 p-value")
dev.off()
