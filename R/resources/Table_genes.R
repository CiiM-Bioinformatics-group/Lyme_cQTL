library(data.table)
library(tidyverse)
library(VennDiagram)
source("~/bin/Functions.R")
library(RColorBrewer)
setwd('resources')

lognorm <- fread("../../Output/Cytokines_goodqual_5e-8.txt")
colnames(lognorm) <- c("SNP","gene","beta","t","P")

lognorm <- remove_badqual(lognorm, column = "gene")
lognorm <- separate(lognorm, SNP, into = c("chrpos","SNP"),"%")%>%
  unite("SNP_gene",c(SNP,gene), sep = "_", remove = F)%>% 
  mutate("notime" = gsub("t\\d_","",gene))%>%
  separate(gene, into = c("Time","Sample","Cyt","Stim"), sep ="_", remove = F)%>%
  separate(chrpos, into = c("chr","BP"), ":", remove = F)



##Read independents snps, r2 < 0.6

Indsnp = fread("../../Output/FUMA_cQTL/IndSigSNPs.txt")

##Add cytokines and stimulations in which they are mutated

Indsnp %<>% 
  rename("SNP" = rsID)%>%
  left_join(lognorm[,c("SNP","gene")], "SNP")%>%
  group_by(GenomicLocus, uniqID,SNP, chr,pos,p,nSNPs, nGWASSNPs)%>%
  summarise("gene" = paste(gene, collapse = ";"))

loci = Indsnp

library(phenoscanner)

ps <- phenoscanner(loci$SNP, catalogue = "GWAS")
pos <- ps$snps %>%
  select(c("snp","consequence","hgnc"))%>%
  rename("SNP" = snp)
eqtl <- phenoscanner(loci$SNP, catalogue = "eQTL")$results %>%
  select(c("snp","exp_gene"))%>%
  rename("SNP" = snp)
gwas <- ps$results %>%
  select(c("snp","trait"))%>%
  rename("SNP" = snp)

loci %<>%
  group_by(across(c(-exp_gene)))%>%
  summarise_at("exp_gene", function(x) paste(unique(x), collapse = ";"))%>%
  left_join(gwas, "SNP")%>%
  group_by(across(c(-trait)))%>%
  summarise_at("trait", function(x) paste(unique(x), collapse = ";"))


loci2 <- loci

cytoks = c()
stims = c()
for(i in 1:nrow(loci)){
  d = cyt_split(data.frame("g" = unlist(str_split(loci$gene[i],";"))), "g")
  cytoks = c(cytoks, paste(unique(d$Cytokine), collapse = ";"))
  stims = c(stims, paste(unique(d$Stim), collapse = ";"))
}

loci$Cytokines = cytoks
#loci$Stimulation.wli = loci$Stimulation
loci$Stimulation = stims

fwrite(loci, "../../Output/FUMA_cQTL/IndSigSNPs.annotated.txt",  sep = "\t")


# Squeeze into loci

splitpaste <- function(vector){
  vector =paste(vector, collapse = ";")
  vector = unlist(str_split(vector, ";"))
  return(paste(unique(vector), collapse = ";"))
}

loci.2.0 <- loci %>%
  group_by(GenomicLocus)%>%
  summarise(nGWASSNPs = sum(nGWASSNPs), nIndSNP = n(), 
            uniqID = uniqID[which.min(p)], SNP = SNP[which.min(p)], chr = chr[which.min(p)], pos = pos[which.min(p)], p = p[which.min(p)],
            gene = splitpaste(gene),Cytokines = splitpaste(Cytokines), Stimulation = splitpaste(Stimulation),
            #gene.wli = splitpaste(gene.wli), Stimulation.wli = splitpaste(Stimulation.wli),
            consequence = splitpaste(consequence),
            hgnc = splitpaste(hgnc), exp_gene = splitpaste(exp_gene),
            trait = splitpaste(trait)) 

fwrite(loci.2.0, "../../Output/FUMA_cQTL/Indloci.annotated.txt",  sep = "\t")



loci.2.0 <- fread("../../Output/FUMA_cQTL/Indloci.annotated.txt")

