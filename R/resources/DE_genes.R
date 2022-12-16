library(data.table)
library(tidyverse)
library(VennDiagram)
source("~/bin/Functions.R")
library(RColorBrewer)
setwd('../code/')


#Read locitable

loci <- fread("../Output/FUMA_cQTL/Indloci.annotated.txt")

#Read DE genes

DE_24 <- fread("../Output/DE_genes/Borrelia burgdorferi_24h_DE_all_FDR005_fc2both.csv")
DE_4 <- fread("../Output/DE_genes/Borrelia burgdorferi_4h_DE_all_FDR005_fc2both.csv")

#Save to Table
DE_24 %>% select(gene.name, log2FoldChange)%>%arrange(-log2FoldChange)%>%fwrite('../Output/DE_genes/Borrelia burgdorferi_24h_DE_all_FDR005_fc2both_table.txt')
DE_4 %>% select(gene.name, log2FoldChange)%>%arrange(-log2FoldChange)%>%fwrite('../Output/DE_genes/Borrelia burgdorferi_4h_DE_all_FDR005_fc2both_table.txt')

#Which genes with known eqtl effect are misregulated
DE_ov24 <- c()
DE_ov4 <- c()
for(i in 1:nrow(loci)){
  g = unlist(str_split(loci$exp_gene[i], ";"))
  g1 = g[which(g %in% DE_24$gene.name)]
  g2 = g[which(g %in% DE_4$gene.name)]
  DE_ov24 <- c(DE_ov24, paste(g1, collapse = ";"))
  DE_ov4 <- c(DE_ov4, paste(g2, collapse = ";"))
}

loci$DE_stimulation_24 = DE_ov24
loci$DE_stimulation_4 = DE_ov4


loci %<>% 
  mutate(gene = strsplit(gene, ";"))%>%
  unnest(gene)%>%
  cyt_split("gene")%>%
  group_by(GenomicLocus)%>%
  summarise_all(function(x) paste(unique(x), collapse = ";"))
  

#Write again

fwrite(loci, "../Output/FUMA_cQTL/Indloci.annotated.txt", sep = "\t")


#

library(kableExtra)
options(knitr.kable.NA = "")

loci %>%
  select(GenomicLocus, nGWASSNPs, SNP, chr, pos, p, Cytokine, Stimulation, Time, Sample, exp_gene)%>%
  rename("Locus Nr." = GenomicLocus,
         "nSNPS" = nGWASSNPs)%>%
  mutate_all(function(x) gsub(";", ", ", x))%>%
  mutate(p = formatC(as.numeric(p), format = "e", digits = 2))%>%
  kable(align = "c", escape = T)%>%
  kable_styling()%>%
  column_spec(6, width = "10em")%>%
  add_header_above(c(" " = 2,"Lead SNP" = 4, " " = 5))%>%
  save_kable("../Output/FUMA_cQTL/Indloci.table.pdf")




## GSEA



#Only the loci that are associated to borrelia
indloci <- loci[grepl("bbmix", loci$gene),]%>%
  mutate("exp_gene" = gsub("NA;","", exp_gene)%>%gsub("-;","",.))

#genes that are related to borrelia
genes <- indloci%>%arrange(p)%>% .$exp_gene %>%
  gsub("NA;","", .)%>%gsub("-;",";",.)%>%
  na.omit()%>%paste(collapse = ";")%>%
  str_split(";")%>%.[[1]]

#My own annotation from Yang's DE list upon borrelia
DE_24 <- DE_24%>%.[.$padj < 0.05,]%>%
  mutate("group" = ifelse(log2FoldChange > 0, "DE_24_pos", "DE_24_neg"))%>%
  mutate("score" = sign(log2FoldChange)*-log10(padj))
DE_4 <- DE_4%>%.[.$padj < 0.05,]%>%
  mutate("group" = ifelse(log2FoldChange > 0, "DE_4_pos", "DE_4_neg"))%>%
  mutate("score" = sign(log2FoldChange)*-log10(padj))

#Bind both times to the same df
DE <- rbind(DE_24, DE_4)
DE <- DE[,c("group","gene.name")]


library(clusterProfiler)
library(enrichplot)

#Pull log2FC ordered genes from 24h DEGs
DE_24_all <- DE_24 %>% pull(abs(log2FoldChange), name = gene.name)%>%sort(decreasing = T)
#Remove Inf values if any (happened when calculating log)
DE_24_all = DE_24_all[!DE_24_all == Inf]
#Calculate enrichment of 24h DEGs in cQTLs
enr <- GSEA(DE_24_all, TERM2GENE = data.frame("group" = "cQTL", "gene" = genes))
#Save enrichment table
fwrite(as.data.frame(enr), "../Output/DE_genes/DEgenes_enr.txt", sep = "\t")
#Plot GSE
p <- gseaplot(enr, geneSetID = "cQTL", pvalue_table = T)
p
ggsave("../Output/DE_genes/DE_borrelia_24_enrichment.pdf")
