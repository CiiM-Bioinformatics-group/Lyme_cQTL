library(data.table)
library(tidyverse)
source("~/bin/Functions.R")
setwd("resources")

#gw <- fread("../Output/Cytokines_goodqual.txt")

#gw <- remove_badqual(gw, "gene")

#gw_man <- gw %>%
 # group_by(SNP) %>% summarise (`p-value` = min(`p-value`))%>%
 # separate(SNP, into=c("CHR","BP"), sep = ":", remove = F)%>%
 # dplyr::rename("P" = `p-value`) %>%
 # mutate("CHR" = as.integer(CHR), "BP" = as.integer(BP))

#fwrite(gw_man, "../Output/Cytokines_goodqual_manhattan.txt", sep = "\t")

manhattan <- fread("../../Output/Cytokines_goodqual_manhattan.txt")
Indloci <- fread("../../Output/FUMA_cQTL/GenomicRiskLoci.txt")
Indloci <- Indloci %>%
  mutate("min"= pos - 500000, "max" = pos + 500000)


match_loci <- function(CHR, BP, n) {
  print(paste0(format(n, digits = 2), " %"))
  df <- Indloci %>% filter(chr == CHR)
  if(nrow(df) == 0) return(0)
  df$pos <- as.numeric(df$pos)
  BP <- as.numeric(BP)
  if(all(BP < df$min | BP > df$max)) return(0)
  return(df$GenomicLocus[BP > df$min & BP < df$max])    
}

manhattan.2 <- manhattan %>%
  mutate("n" = 1:nrow(.))%>%
  rowwise()%>%
  mutate("loci"= match_loci(CHR, POS, n*100/nrow(.)))

fwrite(manhattan.2, "../../Output/Cytokines_goodqual_manhattan_loci.txt", sep = "\t")
