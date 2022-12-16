library(data.table)
library(tidyverse)
library(VennDiagram)
source("~/bin/Functions.R")
library(RColorBrewer)
library(susieR)
setwd('resources')

#Read loci
indloci <- fread("../../Output/FUMA_cQTL/Indloci.annotated.txt")

#Function to remove, from a dataframe containing SNPid the duplicated SNPid with the highest p-value
remove_dupsnp <- function(df){
  dup <- df$SNPid[duplicated(df$SNPid)]
  discard <- sapply(dup, function(snp)df %>% arrange(-P) %>% filter(SNPid ==snp) %>% pull(SNP) %>%.[1])
  df <- df[!df$SNP %in% discard,]
  return(df)
}

parallel::mclapply(unique(indloci$GenomicLocus), function(locus){

#Extracting the information for the locus
ref <- indloci[indloci$GenomicLocus == locus,]
chr = ref$chr
min = ref$pos-5e5
max = ref$pos+5e5

#We need to apply it per affected measure.
#Option 1: Only cytokine stimulations with a gw snp
gene <- str_split(ref$gene, ";")[[1]]%>%
  gsub("-","_",.)
#Option 2: all genes
#gene <- sapply(indloci$gene, function(x) str_split(x, ";"))%>%unlist()%>%unique()%>% gsub("-","_",.)
#First, I apply it for one gene to get the correlation matrix
gene.1 <- gene[1]
#Read the summary statistics for a given gene in locus 8
locus.1 <- fread(cmd = paste0("awk -F ':'  '{if($1 ==",chr,
                            " && $2 > ",min,
                            " && $2 < ",max,
                            ")print}' ../../Output/Cytokines_goodqual_pergene/",gene.1))
#Extracting the dosage and correlation
dosage <- fread(cmd = paste0("egrep 'LP|",paste(locus.1$V1, collapse = "|"),
                             "' ../../data/lyme-merged_dosage.tsv"))
dosage <- as.data.frame(dosage)%>%
  column_to_rownames("V1")
R <- cor(t(dosage))

#Now, we loop through the genes to compute the fine mapping

genesusie <- lapply(gene, function(gene){
  print(gene)
  #Read the summary statistics for a given gene in locus 8
  locus <- fread(cmd = paste0("awk -F ':'  '{if($1 ==",chr,
                                " && $2 > ",min,
                                " && $2 < ",max,
                                ")print}' ../../Output/Cytokines_goodqual_pergene/",gene))
  
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
#Save susie results
fulldf = lapply(genesusie, function(x) x$df) %>% bind_rows() %>% filter(!is.na(CS))
fwrite(fulldf, paste0("../../Output/finemapping/finemapping_locus",locus,".tsv"), sep = "\t")

#Which conditions have bigger CS? What is the overlap between the CS variants?
ggplot(fulldf, aes(x=gene))+geom_bar()+theme_minimal()+
theme(axis.text.x = element_text(angle = 90))
ggsave(paste0("../../Output/finemapping/finemapping_CSsize_locus",locus,".pdf"))

#Compare CS size between time 0 and time 6
library(ggsignif)
fulldf%>%
  group_by(gene)%>%summarise("CSsize" = n())%>%
  cyt_split("gene")%>%
  ggplot(aes(x=Time, y=CSsize))+
  geom_jitter(fill = "black", alpha = .5, pch = 21, size = 1)+
  geom_violin(alpha = .5)+
  geom_boxplot(alpha = .5, width = .2)+
  geom_signif(comparisons = list(c("t0","t6")), map_signif_level = T)+
  theme_bw()+theme(panel.grid = element_blank())
ggsave(paste0("../../Output/finemapping/finemapping_CSsize_locus",locus,"_boxplot.pdf"))


df2 = fulldf %>% group_by(rsID)%>%summarise(group = n())%>%group_by(group)%>%mutate(value = n())%>%
  mutate(group = ifelse(group == 1, "1", ifelse(group < 5, "1-5", ifelse(group < 10, "5-10",ifelse(group < 14, "10-14","15")))))%>%
  mutate(group = factor(group, levels = c("1","1-5","5-10","10-14","15")))
ggplot(df2%>%group_by(group) %>% mutate(sum = sum(value)),
  aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=1) +
    geom_text(data = df2%>%group_by(group) %>% summarise(sum = sum(value)),aes(label = sum, y = sum),
     position = position_stack(vjust = 0.5))+
    coord_polar("y", start=0)+theme_minimal()
ggsave(paste0("../../Output/finemapping/finemapping_piechart_locus",locus,".pdf"))

#Top snp per cytokine stimulation
top = lapply(genesusie, function(x) x$df %>% arrange(-PIP) %>% slice(1)) %>% bind_rows()
top
#PIP per cytokine stimulation of all the snps that are top in any of them
top.expanded = lapply(genesusie, function(x) x$df %>% filter(SNP %in% top$SNP))%>%
  bind_rows()%>%
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
pdf(paste0("../../Output/finemapping/finemapping_locus",locus,"_heatmap.pdf"), width = 6, height = 10)
print(ComplexHeatmap::Heatmap(top.expanded, show_row_dend = F, show_column_dend = F,
                        col =colorRamp2(c(0,0.5,1),c("white","orange","red")),
                        right_annotation = ra, show_row_names = F))

dev.off()
}, mc.cores = 3)

#Plot change in sigificance levels (absolute T statistics) for that SNP between timepoints

#rs574 <- fread(cmd = "egrep 'rs5743618|gene' ../../Output/FUMA_cQTL/allsnps.pval.tsv")%>%
#  cyt_split("gene")%>%
#  filter(grepl("p3c|bbmix|bafzelii", gene))%>%
# mutate(absT = abs(`t-stat`))
#library(ggpubr)
#ggpaired(rs574, x = "Time", y = "t-stat", id = "notime",alpha = .5, width = .2)+
#  stat_compare_means(paired = T)+
#  geom_violin(alpha = .5)+
#  ylab("T-statistic")+xlab("Time")+
#  theme_bw()+theme(panel.grid = element_blank())
#ggsave("../../Output/finemapping/rs5743618_time_boxplot.pdf")
