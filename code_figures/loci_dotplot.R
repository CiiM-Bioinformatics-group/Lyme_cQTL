library(data.table)
library(tidyverse)
source("../resources/Functions.R")

#Read loci
loci <- fread("../Output/FUMA_cQTL/Indloci.annotated.txt")
loci$SNPid <- paste(loci$chr, loci$pos, sep = ":")

#Add MAF
MAF <- fread(cmd = paste0("zgrep -h -w ",paste(paste("-e", unique(loci$SNPid)), collapse = " "),
                          " ../PostImp/*.info.gz"))
colnames(MAF) <- c("SNPid","Ref","Alt","Alt_frq","MAF","AvgCall","Rsq","Genotyped","-","--","---","----","-----")

MAF <- MAF %>%
  separate(SNPid, c("chr","pos","Ref","Alt"))%>%
  mutate_at(c("chr","pos"), as.integer)

loci <- left_join(loci, MAF, by = c("chr","pos"))

#Add if they are study-wide significant or not
loci$SW <- loci$p < 1.55e-9
loci$Cytokine <- ifelse(grepl(";", loci$Cytokine), "More", loci$Cytokine)

nrow(unique(loci))
#One locus is duplicated, keep the good one

loci <- filter(loci, MAF > 0.01)

#Plot dotplot

ggplot(loci, aes(y=-GenomicLocus*100, x=1))+
  geom_point(aes( fill = Cytokine, size = 100*MAF, color = SW),
             pch=21, stroke =1)+
  scale_size(range = c(3,10))+
  scale_color_manual(values = c("white","black"))+
  scale_fill_manual(values = c(anno_col$Cytokine, "More" = "darkred"))+
  geom_text(aes(label = as.integer(MAF*100)), size = 3)+
  theme_void()+
  theme(legend.position = 0)
ggsave("Figures/loci.dotplot.row.pdf", height = 12)

ggplot(loci, aes(x=GenomicLocus*100, y=1))+
  geom_point(aes( fill = Cytokine, size = 100*MAF, color = SW),
             pch=21, stroke =1)+
  scale_size(range = c(3,10))+
  scale_color_manual(values = c("white","black"))+
  scale_fill_manual(values = c(anno_col$Cytokine, "More" = "darkred"))+
  geom_text(aes(label = as.integer(MAF*100)), size = 3)+
  theme_void()+
  theme(legend.position = 0)
ggsave("Figures/loci.dotplot.column.pdf", width  = 12)
