library(coloc)
library(data.table)
library(tidyverse)
library(LDlinkR)
library(susieR)
library(TwoSampleMR)
library(parallel)
source("../resources/Functions.R")


#The top SNP is rs4110197
#Boxplot of different things based on the snp
cQTL <- fread(cmd = " egrep 'LP|11:61997649:A:G%rs4110197' ../data/lyme-merged_dosage.tsv")

#Load whole pheno
cytokines <- read.csv2("../Raw/Phenotype/LymeProspect_cytokines_clinical_traits.csv")%>%
  .[,!grepl("rs|diss|combined",colnames(.))]%>%
  .[,-1]%>%
  mutate(sex = factor(sex,levels = c("F","M"), labels = c(0,1)))%>%
  column_to_rownames("Lab_ID")%>%
  mutate_all(as.numeric)%>%
  t()%>%  as.data.frame()


cytokines["t0_wba_il6_il10_bbmix_moi10",] <- as.numeric(cytokines["t0_wba_il6_bbmix_moi10",])/as.numeric(cytokines["t0_wba_il10_bbmix_moi10",])
cytokines["t6_wba_il6_il10_bbmix_moi10",] <- as.numeric(cytokines["t6_wba_il6_bbmix_moi10",])/as.numeric(cytokines["t6_wba_il10_bbmix_moi10",])


cytokines <- rownames_to_column(cytokines, "V1")
#Figure 4.A cis IL6 in Borrelia mix moi30 whole blood

cQTL_boxplot <- function(Cyt, write = T, col = 1){
  #Read genotype and cytokine
  Cytokine_measured<- cytokines[grepl(Cyt, cytokines$V1),]
  if(nrow(cQTL) > 1) stop("Too many snps")
  if(nrow(Cytokine_measured) > 1) stop("Too many snps")

  #Match column names
  Cytokine_measured <- Cytokine_measured[,colnames(cQTL)]

  #Bind
  df <- rbind(cQTL, Cytokine_measured)%>%
    remove_rownames()%>%
    column_to_rownames("V1")%>%
    t()%>%
    as.data.frame()

  split = str_split(colnames(df)[1], "[^[:alnum:]]+", simplify = T)
  ref = split[,3]
  alt = split[,4]
  snp = split[,5]

  df$snp <- ifelse(df[,1] <  0.5, paste0(ref,ref), ifelse(df[,1] < 1.5, paste0(ref,alt), paste0(alt,alt)))
  df$snp <- factor(df$snp, levels = c(paste0(ref,ref),paste0(ref,alt),paste0(alt,alt)))
  colnames(df)[2] <-"Cytokine"
  if(!Cyt %in% c("C6_t0", "C6_t6")) df$Cytokine <- log2(df$Cytokine)
  require(ggpubr)
  p <- ggboxplot(df, x = "snp",y = "Cytokine", fill  = ggsci::pal_npg()(3)[col], add = c("point"), add.params = list(position = position_jitter(w = 0.05)))+
    ylab(paste(Cyt))+xlab("rs4110197")+
    theme_bw()+theme(panel.grid = element_blank())+
    stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)), label = 'p.signif')
  if(write) ggsave(paste0("Figures/FinnLymeGWAS",Cyt,".pdf"), plot = p, width = 5, height = 5)
  p

}


cQTL_boxplot("C6_t0", col = 1)
cQTL_boxplot("C6_t6", col = 1)
cQTL_boxplot("t0_wba_il6_bbmix_moi10", col = 3)
cQTL_boxplot("t6_wba_il6_bbmix_moi10", col = 3)
cQTL_boxplot("t0_wba_il10_bbmix_moi10", col = 3)
cQTL_boxplot("t6_wba_il10_bbmix_moi10", col = 3)

cQTL_boxplot("t0_wba_il6_il10_bbmix_moi10", col = 2)
cQTL_boxplot("t6_wba_il6_il10_bbmix_moi10", col = 2)


#Heatmap of the significance of that SNP and the missense SNP on our data
cQTL.sumstats <- fread(cmd = " egrep 'LP|11:61997649:A:G%rs4110197|11:62010863:C:T%rs2232950' ../Output/Cytokines_goodqual.txt")
colnames(cQTL.sumstats) <- c("SNP","gene", "beta","t","p-value")
#Change this a bit
ind_lognorm <- cQTL.sumstats%>%
  mutate("log.pval" = sign(beta)*-log10(`p-value`))%>% #Reflecting cQTL direction
  dplyr::select(c("SNP","gene","log.pval"))%>%
  pivot_wider(names_from = gene, values_from = log.pval, values_fill = 0)%>%
  column_to_rownames("SNP")

ca = ComplexHeatmap::HeatmapAnnotation(Time = str_split(colnames(ind_lognorm), "_", simplify = T)[,1],
                                   `Cell system`= toupper(str_split(colnames(ind_lognorm), "_", simplify = T)[,2]),
                                   Cytokine = gsub("il","IL",str_split(colnames(ind_lognorm), "_", simplify = T)[,3]),
                                   Stimulation = str_split(colnames(ind_lognorm), "_", simplify = T)[,4],
                                   col = anno_col)

pdf("../Output/Finngen/FinnLymeGWAS_Heatmap.pdf", width = 15, height = 7)
ComplexHeatmap::Heatmap(as.matrix(ind_lognorm),top_annotation = ca,
                        col = circlize::colorRamp2(c(log10(1e-4),log10(0.05), log10(0.05)+0.0001, 0,  -log10(0.05)-0.0001 ,-log10(0.05), -log10(1e-4)),
                                                  c("blue","lightblue","white", "white","white","orange", "red")))
dev.off()




#Fine map this locuss to find the causal variants per trait
#My finemapping function
genesusie <-function(gene, chr, min, max){
  #gene <- sapply(indloci$gene, function(x) str_split(x, ";"))%>%unlist()%>%unique()%>% gsub("-","_",.)
  #First, I apply it for one gene to get the correlation matrix
  gene.1 <- gene[1]
  #Read the summary statistics for a given gene in locus 8
  locus.1 <- fread(cmd = paste0("awk -F ':'  '{if($1 ==",chr,
                                " && $2 > ",min,
                                " && $2 < ",max,
                                ")print}' ../Output/Cytokines_goodqual_pergene/",gene.1))
  #Extracting the dosage and correlation of locus8
  dosage <- fread(cmd = paste0("egrep 'LP|",paste(locus.1$V1, collapse = "|"),
                               "' ../data/lyme-merged_dosage.tsv"))
  dosage <- as.data.frame(dosage)%>%
    column_to_rownames("V1")
  R <- cor(t(dosage))


  #Then, apply through all genes
  lapply(gene, function(gene){
  print(gene)
  #Read the summary statistics for a given gene in locus
  locus <- fread(cmd = paste0("awk -F ':'  '{if($1 ==",chr,
                              " && $2 > ",min,
                              " && $2 < ",max,
                              ")print}' ../Output/Cytokines_goodqual_pergene/",gene))


  #Same roworder in both
  locus <- as.data.frame(locus)%>%
    column_to_rownames("V1")
  locus <- locus[rownames(R),]

  #Caluclate z from P-value and beta
  locus$z <- qnorm(locus$V5/2, lower.tail = F)*sign(locus$V3)
  #Plotz
  p1 <- susie_plot(locus$z, y = "z", b=locus$V3)

  #Fine mapping
  fitted_rss <- susie_rss(locus$z, R, L = 50)

  #Transform into df

  Susiedf <- data.frame("PIP" = fitted_rss$pip,  "CS" = NA)%>%
    rownames_to_column("SNP")%>%
    separate(SNP, into =c("CHR","BP","REF","ALT","rsID"),remove = F)

  for (x in names(fitted_rss$sets$cs)){
    Susiedf$CS[fitted_rss$sets$cs[[x]]] <- x
  }

  p2 <- ggplot(Susiedf, aes(x = as.numeric(BP), y = PIP, color = CS))+
    geom_point()+
    geom_text(data = Susiedf[which.max(Susiedf$PIP),], aes(label = rsID, y = PIP+0.05) )+
    theme_bw()+
    ggtitle(gene)
  Susiedf$gene <- gene
  return(list("df" = Susiedf, "p1" = p1, "p2" = p2))
  })

}


locus.fine <- genesusie(unique(cQTL.sumstats$gene[cQTL.sumstats$`p-value` < 0.05]), chr = "11", min =61997649-5e5, max= 61997649+5e5)

locus.fine.bind <- lapply(locus.fine, function(x) x$df)%>%
  bind_rows()

#SNPS that are fine-mapped in any gene
pip.snp <- unique(locus.fine.bind$rsID[locus.fine.bind$PIP >0.03])
#Only IL-6 and IL-10 results
top.expanded <- locus.fine.bind %>% filter(#grepl("t._wba_il[6|10]_bbmix", gene) &
  rsID %in% pip.snp)%>%
  pivot_wider(id_cols = c(gene), names_from = SNP, values_from = PIP)%>%
  column_to_rownames("gene")%>%
  as.matrix()


#Heatmap annotation
ra = ComplexHeatmap::rowAnnotation(Time = str_split(rownames(top.expanded), "_", simplify = T)[,1],
                                   `Cell system`= str_split(rownames(top.expanded), "_", simplify = T)[,2],
                                   Cytokine = str_split(rownames(top.expanded), "_", simplify = T)[,3],
                                   Stimulation = str_split(rownames(top.expanded), "_", simplify = T)[,4],
                                   col = anno_col)
#Change colnames
colnames(top.expanded) <-str_split(colnames(top.expanded), pattern = "%", simplify = T)[,2]

library(circlize)
pdf(paste0("../Output/Finngen/FinnLymeGWAS_finemapping_heatmap.pdf"), width = 6, height = 10)
print(ComplexHeatmap::Heatmap(top.expanded, show_row_dend = F, show_column_dend = F,
                              col =colorRamp2(c(0,0.2),c("white","red")),
                              right_annotation = ra, show_row_names = F))
dev.off()




#Let's see how our cQTL can help us prioritize the finngen results_coloc

#1. Heatmap


#Read gw significant snps
loci <- fread("../Output/FUMA_cQTL/snps.txt")
loci <- loci[!is.na(loci$gwasP),]


fgwas <- list.files("../Output/GWASCatalogue/",pattern = "finngen")

library(rtracklayer)
library(GenomicRanges)

#Liftover function
gwas_liftover <- function(gwas){
  #Liftover positions to hg18 (https://master.bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html)
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


#Get min p-value per locus in the finngen gwas
results <- lapply(fgwas, function(name){
  path <- paste0("../Output/GWASCatalogue/",name)
  #Read file
  gwas <- fread(path)
  #Adjust colnames
  colnames(gwas) <- c("chr", "pos","ref","alt","rsids","nearest_genes","p.gwas","mlogp","beta","sebeta", "af_alt", "af_alt_cases", "af_alt_controls")

  #Liftover those from hg38
  gwas <- gwas_liftover(gwas)
  #Extract p-value per loci
  gwas$chr <- as.integer(gwas$chr)

  loci_gwas <- left_join(loci, gwas, by = c("chr","pos"))

  #Get minimum p-value per locus
  loci_gwas <- loci_gwas[,c("GenomicLocus","rsID","p.gwas")]
  loci_gwas <- loci_gwas %>% group_by(GenomicLocus) %>% summarise(p.gwas = ifelse(all(is.na(p.gwas)), NA, min(p.gwas, na.rm = T)),
                                                                  SNP = ifelse(all(is.na(p.gwas)), NA, rsID[which.min(p.gwas)]))
  gwas.p <- loci_gwas$p.gwas
  names(gwas.p) <- loci_gwas$GenomicLocus
  return(gwas.p)
})
names(results) <- gsub(".gz","", fgwas) %>% gsub("finngen_R7_AB1_", "", .)
results <- as.data.frame(results)


mat <- -log10(as.matrix(results))
rownames(mat) <- str_split(rownames(mat), pattern = " ", simplify = T)[,1]
#Discard traits without significant assoc
mat = mat[,apply(mat, 2, function(x) any(!is.na(x) & x > -log10(0.05)))]
thr = 0.05/nrow(mat)
sig <- ifelse(mat > -log10(5e-8), "***",
              ifelse(mat > -log10(thr), "**",
                     ifelse(mat > -log10(0.05), "*", "")))

sig[is.na(sig)] <- ""
pdf("../Output/Finngen/loci_finngen_gwas_heatmap.pdf", width =6, height = 12)
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



#2. Coloc and MR
indloci <- fread("../Output/FUMA_cQTL/Indloci.annotated.txt")
#Function to remove, from a dataframe containing SNPid the duplicated SNPid with the highest p-value
remove_dupsnp <- function(df){
  dup <- df$SNPid[duplicated(df$SNPid)]
  discard <- sapply(dup, function(snp)df %>% arrange(-P) %>% filter(SNPid ==snp) %>% pull(SNP) %>%.[1])
  df <- df[!df$SNP %in% discard,]
  return(df)
}

gwas_name <- "OTHER_SPIROCHAETAL"
gw <- fread("../Output/Cytokines_goodqual_5e-8.txt")
gene.it <- unique(gw$gene) %>% .[grepl("bbmix", .)]
#	      	      gene.it <- c("t0_pbmc_il10_bbmix_10_5", "t6_pbmc_il6_bbmix_10_5")
#Read gwas
gwas <- fread(paste0("../Output/GWASCatalogue/finngen_R7_AB1_",gwas_name,".gz"))
#Adjust colnames
colnames(gwas) <- c("chr", "pos","other_allele","effect_allele","rsids","nearest_genes","p.gwas",
                    "mlogp","beta","sebeta", "effect_allele_frequency", "af_alt_cases", "af_alt_controls")
#liftOver
gwas <- gwas_liftover(gwas)


#Read QTL MAF: Read only MAF > 0.05
MAF <- fread(cmd = paste0("zgrep -v 'ALT_Frq' ../PostImp/*.info.gz | awk '{if($5 > 0.05) print}' | cut -f 2- -d ':' "))
colnames(MAF) <- c("SNPid","Ref","Alt","Alt_frq","MAF","AvgCall","Rsq","Genotyped","-","--","---","----","-----")
#Delete MAF below 0.05
MAF <- filter(MAF, MAF > 0.05)
res <- mclapply(gene.it, function(gene){
  print(gene)

  locus <- fread(cmd = paste0("cat ../Output/Cytokines_goodqual_pergene/",gene," | tr ':' '\t' | tr '%' '\t'"))
  colnames(locus) <- c("CHR", "BP","REF","ALT","rsID","gene","beta","t","pval")

  #Three different p-value thresholds

  lapply(c(1e-4, 1e-5, 1e-6, 1e-7), function(pthr){
    #SNP info only from interesting snps
    locus <- locus %>% filter(pval < pthr )
    #No sig snps, empty df
    if(nrow(locus) == 0) return(data.frame("id.exposure"=NA,
                                           "id.outcome" = NA, "outcome" = gwas_name, "exposure" = gene,  "method" = NA,
                                           "nsnp" = NA, "b" = NA, "se" = NA, "pval" = NA))


    gwas <- gwas %>% filter(rsids %in% unique(locus$rsID))
    gwas <- gwas %>% filter(nchar(effect_allele)  ==1 & nchar(other_allele) == 1)
    gwas$SNPid <- paste(gwas$chr, gwas$pos, sep = ":")



    #Coloc dataset for the gwas
    gwas <- dplyr::rename(gwas, "SNP" = rsids, "P" = p.gwas)
    gwas <- remove_dupsnp(gwas)
    gwas <- gwas[!is.na(gwas$effect_allele_frequency),]
    gwas$MAF <- ifelse(gwas$effect_allele_frequency > 0.5, 1 - gwas$effect_allele_frequency, gwas$effect_allele_frequency)


    #coloc dataset for the locus
    locus$SNPid <- paste(locus$CHR, locus$BP, locus$REF, locus$ALT, sep = ":")
    print("Locus:")
    print(head(locus))
    locus_MAF <- inner_join(locus, MAF, by = c("SNPid"))
    locus_MAF <- remove_dupsnp(locus_MAF)

    #locus_MAF <- locus_MAF %>% tidyr::separate(SNP, into = c("CHR","BP","REF","ALT","rsID"), remove = F)

    #Mendelian randomisation
    #Exposure dataframe

    locus_MAF$se <- abs(locus_MAF$beta/qnorm(locus_MAF$pval/2, lower.tail = F))
    print("Locus MAF:")
    print(head(locus_MAF))
    #locus_MAF <- locus_MAF[locus_MAF$P < 5e-8,]
    exposure <- format_data(locus_MAF, type = "exposure", snp_col = "rsID", beta_col = "beta", eaf_col = "Alt_frq",
                            effect_allele_col = "ALT", other_allele_col = "REF", pval_col = "pval", se_col = "se",
                            chr_col = "CHR", pos_col = "BP")

    #Perform ld clumping
    exposure <- clump_data(exposure)

    #Outcome dataframe
    outcome <- format_data(gwas, type = "outcome", eaf_col = "effect_allele_frequency", se_col = "sebeta",
                           pval_col = "P", chr_col = "chr", pos_col = "pos")
    #Harmonise
    dat <- harmonise_data(exposure_dat = exposure, outcome_dat = outcome)

    #Run MR
    mr <- mr(dat, method_list=c("mr_egger_regression", "mr_ivw"))
    mr$exposure <- gene
    mr$outcome <- "OTHER_SPIROCHAETAL"
    mr$pvalue.thr <- pthr

    theme_box <- theme_bw()+theme(panel.grid = element_blank())
    mrplot <- mr_scatter_plot(mr, dat)[[1]]+theme_box

    if(nrow(mr %>% filter(pval < 0.05 & nsnp > 2)) > 0) {

      res_single <- mr_singlesnp(dat)
      p4 <- mr_funnel_plot(res_single)
      p2 <- mr_forest_plot(res_single)
      ggsave(paste0("Figures/MR_multisnp_forestplot_",gene, "_", pthr,"_finngen.pdf"), p2[[1]]+theme_minimal()+theme(legend.position = 0))
      ggsave(paste0("../Output/Finngen/MR_multisnp_funnelplot_",gene, "_", pthr,"_finngen.pdf"), p4[[1]]+theme_minimal())
      ggsave(paste0("Figures/MR_multisnp_sedotplot_",gene, "_", pthr,"_finngen.pdf"), mrplot+theme(legend.position = 0))

    }
    return(mr)
  })
}, mc.cores = length(gene.it)) %>%
  bind_rows()

fwrite(res, "../Output/Finngen/MR_multisnp_finngen.txt", sep = "\t")
