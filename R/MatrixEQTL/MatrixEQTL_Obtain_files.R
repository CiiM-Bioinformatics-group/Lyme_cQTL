library(tidyverse)
library(data.table)
library(reshape)
source('resources/Functions.R')



#Cytokine formatting
#Load raw phenotype data
pheno <- read.csv2("../Phenotype/LymeProspect_cytokines_clinical_traits.csv")%>%
  mutate(sex = factor(sex,levels = c("F","M"), labels = c(0,1)))%>%
  column_to_rownames("Lab_ID")%>% 
  mutate_all(as.numeric)%>%
  t()%>%  as.data.frame()


#Keep cytokine values
Cyt <- pheno[grepl("pbmc|wba", rownames(pheno)),] %>%
  apply(2, log2) %>% as.data.frame()

#Write all cytokines into one tsv
fwrite(Cyt, file = "../Phenotype/Cytokines.tsv", sep = "\t", row.names = T)

#Covariates formatting
Cov <- pheno[grepl("age|sex", rownames(pheno)),]

fwrite(Cov, file = "../Covariates/covariates.tsv", sep ="\t", row.names = T)

#C6 data formatting
C6 <- pheno[grepl("C6_", rownames(pheno)),]
C6 <- apply(C6,1,INT)%>%t()

fwrite(C6, file = "../Phenotype/C6_intnorm.tsv", sep ="\t", row.names = T)


##Clean everything, start over ################################
rm(list = ls(all.names = T))
gc()
#############################################################

source("resources/Functions.R")
###IMPORTANT: the exact same samples should be in the phenotype and genotype files

Cyt <- fread("../Phenotype/Cytokines.tsv")%>%column_to_rownames(var = colnames(.)[1])
C6 <- fread("../Phenotype/C6_intnorm.tsv")%>%column_to_rownames(var = colnames(.)[1])
#Remove badqual cytokines
Cyt <- remove_badqual(Cyt)

Cov <- fread("../Covariates/covariates.tsv") %>% column_to_rownames(var = colnames(.)[1])
if(grepl("LP", Cov[1,1])){
  colnames(Cov) <- Cov[1,]
  Cov <- Cov[-1,]
}


##Add institute and batch information
Cov["Institute",] <- str_split(colnames(Cov), pattern = "-", simplify = T)[,2]
Cov["Institute",] <- ifelse(Cov["Institute",] == "A", 1, 0)

#Add batch information
lb1 <- fread("lb1.dosage.tsv", nrows = 1)%>%
  colnames(.) %>% .[-1]


Cov["Batch",] <- ifelse(colnames(Cov) %in% lb1, 0, 1)


if(!all(colnames(score) == colnames(Cyt))) Cov <- Cov[,colnames(Cov)%in% colnames(Cyt)]


#Read genotype
lb<- fread("../data/lyme-merged_dosage.tsv")%>%column_to_rownames(var = colnames(.)[1])
print("Dosage Read!")
#Remove the extra stuff
colnames(lb) <- str_split(colnames(lb), pattern = "_", simplify = T)[,2]
print("Column names changed!")

##Remove samples in Cyt, Score and Cov that are not present in the genotype
Cyt <- Cyt[,colnames(Cyt) %in% colnames(lb)]
Cov <- Cov[,colnames(Cov) %in% colnames(lb)]
C6 <- C6[,colnames(C6) %in% colnames(lb)]

all(colnames(score) == colnames(Cyt) & colnames(Cyt) == colnames(Cov) & colnames(Cyt) == colnames(C6))


##Subset genotype data to only have samples that are present in cytokines

lb <- lb[,colnames(Cyt)]
C6 <- C6[,colnames(Cyt)]
Cov <- Cov[,colnames(Cyt)]
print("Dosage samples subsetted!")
print("Have all the files the same samples? T/F")
print(all(colnames(score) == colnames(Cyt) & colnames(Cyt) == colnames(lb) & colnames(Cov) == colnames(lb)))


##If everything is true, write the files adding rownames
##Write the formated files

fwrite(lb, "../data/lyme-merged_dosage.tsv", sep = "\t", quote = F, row.names = T, col.names = T)
fwrite(Cyt, "../Phenotype/Cytokines_goodqual.tsv", sep = "\t", quote = F, row.names = T, col.names = T)
fwrite(Cov, "../Covariates/covariates.tsv", sep = "\t", quote = F, row.names = T, col.names = T)
fwrite(C6, "../Phenotype/C6_intnorm.tsv", sep = "\t", quote = F, row.names = T, col.names = T)

