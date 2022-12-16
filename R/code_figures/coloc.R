library(coloc)
library(data.table)
library(tidyverse)
library(LDlinkR)
library(susieR)
library(TwoSampleMR)
source("../resources/Functions.R")

#Set colors for the diseases
MR.colors <- ggsci::pal_locuszoom()(3)
names(MR.colors) <- c("Age_macular_degeneration", "Allergic_disease", "MS")
MR.colors[c("IBD", "UC","IL10","IL6", "DRAM1","BMI")] <- ggsci::pal_npg()(3)[3]

indloci <- fread("../Output/FUMA_cQTL/Indloci.annotated.txt")

#Function to remove, from a dataframe containing SNPid the duplicated SNPid with the highest p-value
remove_dupsnp <- function(df){
  dup <- df$SNPid[duplicated(df$SNPid)]
  discard <- sapply(dup, function(snp)df %>% arrange(-P) %>% filter(SNPid ==snp) %>% pull(SNP) %>%.[1])
  df <- df[!df$SNP %in% discard,]
  return(df)
}
#Function to plot locus comparer
locuscomparer2 <- function(df1, df2, c1, p1 , c2 , p2,
                           chr1, chr2, pos1, pos2,
                           name1, name2){

  df1$SNP <- df1 %>% pull(c1)
  df1$P1 <- df1 %>% pull(p1)
  df1$chr1 <- df1 %>% pull(chr1)
  df1$pos1 <- df1 %>% pull(pos1)
  df2$SNP <- df2 %>% pull(c2)
  df2$P2 <- df2 %>% pull(p2)
  df2$chr2 <- df2 %>% pull(chr2)
  df2$pos2 <- df2 %>% pull(pos2)

  p1 <- inner_join(df1, df2, by = "SNP")%>%
    ggplot(aes(x=-log10(P1), y= -log10(P2)))+
    geom_point(color = "darkgreen")+
    theme_bw()+xlab(paste0("-log10(",name1,")"))+ylab(paste0("-log10(",name2,")"))

  p1 <- ggplot(df1, aes(x=as.numeric(pos1), y=-log10(P1)))+
    geom_point(color = "blue")+
    theme_bw()+
    xlab("")+ylab("")+ggtitle(name1)+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p2 <- ggplot(df2, aes(x=as.numeric(pos2), y=-log10(P2)))+
    geom_point(color = "red")+
    theme_bw()+
    xlab("Position")+ylab("")+ggtitle(name2)+
    xlim(layer_scales(p2)$x$range$range)+
    theme(axis.text.x = element_text(angle = 90))

  library(patchwork)

  p1+p2/p3

}



comparisons <- list(c("2","Age_macular_degeneration",1295,456348,"t0_pbmc_il1ra_cand"), #Accession number for GWASCatalogue GCST90043776
                    c("3","IBD",12882,34652,"t6_wba_il10_bbmix_moi30"), #GCST003043
                    c("3","UC",6968,27432,"t6_wba_il10_bbmix_moi30"), #GCST003045
                    c("7","Allergic_disease",35890,484598,"t6_wba_il10_bbmix_moi10"), #GCST90038661
                    c("7","BMI","quant",315.347,"t0_wba_il6_bbmix_moi3"), #GCST006368
                    c("15", "MS",1683,325757,"t6_pbmc_il6_bbmix_10_5"), #GCST90014448

                    c("3","IL10","quant",10000,"t6_wba_il10_bbmix_moi30"),  #N = Approximate number from eqtlgen meta-analysis
                    c("16","IL6","quant",10000,"t0_wba_il6_bbmix_moi30"))  #N = Approximate number from eqtlgen meta-analysis
                    #c("28","DRAM1","quant",10000,"t6_wba_il10_lps"))  #N = Approximate number from eqtlgen meta-analysis




results_coloc = lapply(comparisons, function(compar){
  if(length(compar) < 5) stop("Lacking arguments")
  locus_name <- compar[1]
  gwas_name <- compar[2]
  gwas_cases <- compar[3]
  gwas_n <- compar[4]
  gene <- compar[5]

  print(gwas_name)

  #Read gwas
  gwas <- fread(paste0("../Output/GWASCatalogue/",gwas_name,"_hg19.tsv.gz"))
  #Read locus
  ref <- indloci[indloci$GenomicLocus == locus_name,]
  chr = ref$chr
  min = ref$pos-5e5
  max = ref$pos+5e5

  locus <- fread(cmd = paste0("awk -F ':'  '{if($1 ==",chr,
                              " && $2 > ",min,
                              " && $2 < ",max,
                              ")print}' ../Output/Cytokines_goodqual_pergene/",gene))

  colnames(locus) <- c("SNP","gene","beta","t","P")

  #SNP info only from interesting loci
  gwas <- gwas[chromosome == chr & base_pair_location > min & base_pair_location <  max,]
  gwas <- gwas[nchar(effect_allele)  ==1 & nchar(other_allele) == 1,]
  gwas$SNPid <- paste(gwas$chromosome, gwas$base_pair_location, sep = ":")

  #Removing the one duplicated snp
  gwas <- gwas[paste(gwas$variant_id, gwas$other_allele) != "rs6840818 G",]

  #Read QTL MAF

  MAF <- fread(cmd = paste0("zcat ../PostImp/chr", chr, ".info.gz | awk -F ':'  '{if($2 > ",min,
                            " && $2 < ",max,
                            ")print}'"))
  colnames(MAF) <- c("SNPid","Ref","Alt","Alt_frq","MAF","AvgCall","Rsq","Genotyped","-","--","---","----","-----")


  #Coloc dataset for the gwas
  gwas <- rename(gwas, "SNP" = variant_id, "P" = p_value)
  gwas <- remove_dupsnp(gwas)
  gwas <- gwas[!is.na(gwas$effect_allele_frequency),]
  gwas$MAF <- ifelse(gwas$effect_allele_frequency > 0.5, 1 - gwas$effect_allele_frequency, gwas$effect_allele_frequency)
  if(gwas_cases == "quant"){
    gwasD <- list("MAF" = gwas$MAF,
                  "snp" = gwas$SNPid,
                  "type" = "quant", "pvalues" = gwas$P,
                  "N" = as.numeric(gwas_n),
                  "beta" = gwas$beta) #Number from the GWASCatalogue accession
  }else{
    gwasD <- list("MAF" = gwas$effect_allele_frequency,
                  "snp" = gwas$SNPid,
                  "type" = "cc", "pvalues" = gwas$P,
                  "N" = as.numeric(gwas_n),
                  "s" = as.numeric(gwas_cases)/as.numeric(gwas_n),
                  "beta" = gwas$beta) #Number from the GWASCatalogue accession
  }


  check_dataset(gwasD)
  gwas.fine <- finemap.abf(gwasD)

  #coloc dataset for the locus
  locus$SNPid <- str_split(locus$SNP, pattern = "%", simplify = T)[,1]
  locus_MAF <- left_join(locus, MAF, by = c("SNPid"))
  locus_MAF <- locus_MAF[!grepl("G%rs6840818",locus_MAF$SNP),]
  locus_MAF$SNPid <- apply(str_split(locus_MAF$SNPid, pattern = ":", simplify = T)[,c(1,2)], 1, function(x) paste(x, collapse =":"))
  locus_MAF <- remove_dupsnp(locus_MAF)
  #Datasets for coloc
  QTL <- list("MAF" = locus_MAF$MAF,
              "snp" = locus_MAF$SNPid,
              "type" = "quant", "pvalues" = locus_MAF$P,
              "N" = 1063,
              "beta" = locus_MAF$beta)
  check_dataset(QTL)
  QTL.fine <- finemap.abf(QTL)

  #Run coloc
  my.res <- coloc.abf(QTL, gwasD)

  #Tried to run susi, no luck
  #Extracting the dosage and correlation of the locus
  #dosage <- fread(cmd = paste0("egrep 'LP|",paste(locus_MAF$SNP, collapse = "|"),
  #                             "' ../Samples/lyme-merged_dosage.tsv"))
  #dosage <- as.data.frame(dosage)%>%
  #  column_to_rownames("V1")
  #R <- cor(t(dosage))
  #Adding the matrix to the data
  #colnames(R) <- locus_MAF$SNPid
  #rownames(R) <- locus_MAF$SNPid
  #QTL$LD <- R
  #Remove the gwas SNPs that are not in the QTL SNPs
  #gwasD2 <- as.data.frame(gwasD) %>%
  #  filter(snp %in% colnames(R))%>%
  #  as.list()
  #gwasD2$LD <- R
  #check_dataset(QTL)
  #check_dataset(gwasD2)
  #Add the Z values based on the p-value and sign of beta
  #QTL$z <- qnorm(QTL$pvalues/2, lower.tail = F)*sign(QTL$beta)
  #gwasD2$z <- qnorm(gwasD2$pvalues/2, lower.tail = F)*sign(gwasD2$beta)


  #QTL.susi <- runsusie(QTL)
  #gwasD.susi <- runsusie(gwasD2)

  #my.res_susie <- coloc.susie(QTL, gwasD2)

  sensitivity(my.res,"H4 > 0.9")
  df <- subset(my.res$results,SNP.PP.H4>0.1)

  locus_MAF <- locus_MAF %>% tidyr::separate(SNP, into = c("CHR","BP","REF","ALT","rsID"), remove = F)


  #Plotting

  p1 <- ggplot(locus_MAF, aes(x=as.numeric(BP), y=-log10(P)))+
    geom_point(color = ggsci::pal_npg()(3)[2])+
    geom_text(data = locus_MAF[which.min(locus_MAF$P)], aes(x=as.numeric(BP),y=-log10(P), label=rsID), nudge_y = 0.2)+
    theme_bw()+
    ylab(paste0("Locus ",locus_name, "\n\n-log10(Pvalue)"))+xlab("")+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  p2 <- ggplot(gwas, aes(x=as.numeric(base_pair_location), y=-log10(P)))+
    geom_point(color =MR.colors[gwas_name])+
    geom_text(data = gwas[which.min(gwas$P)], aes(x=as.numeric(base_pair_location),y=-log10(P), label=SNP), nudge_y = 0.2)+
    theme_bw()+
    ylab(paste0(gsub("_"," ",gwas_name), "\n\n-log10(Pvalue)"))+xlab("")+
    theme(axis.text.x = element_text(angle = 90))+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  #Genes
  library(biomaRt)
  #Snp reference
  gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37) # we will need an additional mart for genes
  #extract genes
  out.bm.genes.region <- getBM(
    attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype'),
    filters = c('chromosome_name','start','end'),
    values = list(chr, min, max),
    mart = gene.ensembl)
  ## define plot range for x-axis
  plot.range <- c(min(min, out.bm.genes.region$start_position),
                  max(max, out.bm.genes.region$end_position))

  ## rank gene_biotype label
  out.bm.genes.region <- out.bm.genes.region %>% mutate(gene_biotype_fac = fct_relevel(as.factor(gene_biotype),
                                                                                       "protein_coding"), external_gene_name = fct_reorder2(external_gene_name,
                                                                                                                                            start_position, gene_biotype_fac, .desc = TRUE))%>%
    filter(gene_biotype == "protein_coding")
  #Plot genes
  p3 <- ggplot(data = out.bm.genes.region) +
    geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, colour = gene_biotype_fac, group = gene_biotype_fac)) +
    coord_flip() + ylab("") +
    ylim(plot.range) +
    geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, colour = gene_biotype_fac), fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5) +
    labs(caption = paste0("Data source: ", gene.ensembl@host, " + Data set: ", gene.ensembl@dataset), color = "Gene Biotype") +
    theme_bw()+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text.y = element_text(angle = 0),
          legend.position="bottom",
          panel.grid.major.y = element_blank()) +
    expand_limits(y=c(-1, 1))+
    scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out.bm.genes.region$gene_biotype_fac)-1)))
  library(patchwork)
  coloc.theme <-theme(plot.margin = margin(0, 0, 0, 0), panel.border = element_blank(), panel.grid = element_blank(),
                      axis.line = element_line(size = .2))
  p1 <- p1+xlim(plot.range)+coloc.theme
  p2 <- p2+xlim(plot.range)+coloc.theme
  p3 <- p3+coloc.theme+theme(axis.line.y  = element_blank())
  #Adding coloc results, if any
  if(nrow(df) >0){
    snp.tmp <- df$snp
    locus.tmp <- locus_MAF[locus_MAF$SNPid %in% snp.tmp,]
    locus.tmp$PP <- df$SNP.PP.H4[df$snp == locus.tmp$SNPid]
    p4 <- ggplot(locus.tmp, aes(x=as.numeric(BP), y=PP, size=PP))+
      geom_point(color = ggsci::pal_npg()(2)[1])+
      geom_text(aes(label=rsID), nudge_y = 0.1)+
      theme_bw()+ylim(c(0,1.1))+
      ylab(paste0("Locus ",locus_name," - ",gsub("_"," ",gwas_name), "\n\nPosterior Probability"))+xlab("")+
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
      theme(legend.position = 0)
    p4 <- p4+xlim(plot.range)+coloc.theme
    plot <- p1/p2/p4/p3+plot_layout(heights = c(2,2,1,1))
  }else{plot <- p1/p2/p3}
  path <- paste0("Figures/Locuscomparer_",gwas_name,"_",locus_name,".pdf")
  print(paste("Saving plot to",path))
  ggsave(plot = plot,filename = path, width = 5, height = 10)


  #Mendelian randomisation
  #Exposure dataframe
  locus_MAF$se <- abs(locus_MAF$beta/qnorm(locus_MAF$P/2, lower.tail = F))
  locus_MAF <- locus_MAF[locus_MAF$P < 5e-8,]
  exposure <- format_data(locus_MAF, type = "exposure", snp_col = "rsID", beta_col = "beta", eaf_col = "Alt_frq",
                          effect_allele_col = "ALT", other_allele_col = "REF", pval_col = "P", se_col = "se",
                          chr_col = "CHR", pos_col = "BP")

  #Perform ld clumping
  exposure <- clump_data(exposure)

  #Outcome dataframe
  outcome <- format_data(gwas, type = "outcome", eaf_col = "effect_allele_frequency", se_col = "standard_error",
                         pval_col = "P", chr_col = "chromosome", pos_col = "base_pair_location")
  #Harmonise
  dat <- harmonise_data(exposure_dat = exposure, outcome_dat = outcome)

  if(nrow(dat) == 0) return(list("Coloc" = my.res,"Coloc_PP.H4.0.1" =df,"Locuscomparer" = plot))

  #Run MR
  mr <- mr(dat)
  if(nrow(mr) == 0) return(list("Coloc" = my.res,"Coloc_PP.H4.0.1" =df,"Locuscomparer" = plot))

  res<-mr_singlesnp(dat)
  if(nrow(res) == 0) return(list("Coloc" = my.res,"Coloc_PP.H4.0.1" =df,"Locuscomparer" = plot))



  mrplot <- mr_scatter_plot(mr, dat)[[1]]+theme_bw()
  dat$id.exposure <- unique(locus_MAF$gene)
  dat$id.outcome <- gwas_name
  res$id.exposure <- unique(locus_MAF$gene)
  res$id.outcome <- gwas_name
  ggsave(paste0("../Output/MR/MR_",gwas_name,"_",locus_name,".pdf"), plot = mrplot)

  return(list("Coloc" = my.res,
              "Coloc_PP.H4.0.1" =df,
              "Locuscomparer" = plot,
              "MRplot" = mrplot,
              "MR" = res,
              "MR.dat" = dat))

})

saveRDS(results_coloc, "coloc.MR.rds")


results_coloc <- readRDS("coloc.MR.rds")


df <- lapply(results_coloc, function(x){x[["MR.dat"]]} )%>%
  bind_rows()

df2 <-  lapply(results_coloc, function(x){x[["MR"]]} )%>%
  bind_rows()
df <- left_join(df, df2, c("SNP",'id.exposure','id.outcome'))
fwrite(df, "../Output/MR/MR_diseases_table.txt", sep = '\t')

#Remove IBD
df <- subset(df, id.outcome %in% c('Age_macular_degeneration','Allergic_disease','MS','UC'))

ggplot(df,aes(y=beta.exposure, x = beta.outcome))+
  geom_pointrange(aes(ymin=beta.exposure-se.exposure,
                      ymax=beta.exposure+se.exposure, color = id.outcome))+
  geom_pointrange(aes(xmin=beta.outcome-se.outcome,
                      xmax=beta.outcome+se.outcome, color = id.outcome))+
  geom_vline(aes(xintercept=0), linetype = "dashed")+
  geom_hline(aes(yintercept=0), linetype = "dashed")+
  theme_bw()+theme(panel.grid = element_blank())+
  scale_color_manual(values = MR.colors)+
  xlab("SNP effect on outcome")+
  ylab("SNP effect on cQTL")
ggsave("Figures/MR_sedotplot.pdf", width = 7, height = 5)
