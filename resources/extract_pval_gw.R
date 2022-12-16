library("data.table")
library("tidyverse")
setwd('resources')
######ONLY IF FILES ARE NOT YET THERE############################################################
indsnp = fread("../../Output/FUMA_cQTL/snps.txt")

lb = fread(cmd = paste0("egrep -w 'SNP|",paste(indsnp$rsID, collapse = "|"),"' ../../Output/Cytokines_goodqual.txt"))
print('Sum stats read')
lb2 <- lb %>% separate("SNP", into = c("chrpos","rsID"), sep = "%", remove = F)                                                                                                #
lb2 <- lb2 %>% filter(rsID %in% indsnp$rsID)

fwrite(lb2, "../../Output/FUMA_cQTL/allsnps.pval.tsv", sep = "\t", row.names = T)
print('allsnps.pval.tsv written')


lb = fread(cmd = paste0("egrep -w 'LP|",paste(indsnp$rsID, collapse = "|"),"' ../../data/lyme-merged_dosage.tsv"))
print('dosage read')
lb2 <- lb %>% separate("V1", into = c("chrpos","SNP"), sep = "%")                                                                                                #
lb2 <- lb2 %>% filter(SNP %in% indsnp$rsID)                                                  #
lb2 <- lb2 %>% unite("V1",  c("chrpos","SNP"),sep = "%")%>%
  column_to_rownames(var= "V1")

fwrite(lb2, "../../Output/FUMA_cQTL/allsnps.dosage.tsv", sep = "\t", row.names = T)
print('allsnps.dosage.tsv written')

################################################################################################
