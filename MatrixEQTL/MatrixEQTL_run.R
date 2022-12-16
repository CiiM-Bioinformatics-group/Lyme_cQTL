library(data.table)
library(tidyverse)
require(MatrixEQTL)

#Prefix for the Phenotype file, stored in ../Phenotype
prefix = ""
#P-value threshold for the Output
p.val = 1
#Dosage filename
SNP_file_name = ""
#Covariate filename
covariates_fl = ""
##############################


pheno_fl <- paste0("../Phenotype/",prefix,".tsv")
output_fl <- paste0("../Output/",prefix,".txt")

##Check that files exist
if(!file.exists(SNP_file_name)) stop(paste("Geno file",SNP_file_name, "does not exist"))
if(!file.exists(pheno_fl)) stop(paste("Pheno file",pheno_fl, "does not exist"))

##Check that colnames are the same

if(any(colnames(fread(SNP_file_name, nrows = 1))[-1]!=colnames(fread(covariates_fl, nrows = 1))[-1])){
  stop("Covariates and genotype do not match")
}
if(any(colnames(fread(SNP_file_name, nrows = 1))[-1]!=colnames(fread(pheno_fl, nrows = 1))[-1])){
  stop("Phenotype and genotype do not match")
}

useModel <- modelLINEAR
errorCovariance <- numeric()

snps = SlicedData$new()
snps$fileOmitCharacters = "NA"
snps$fileDelimiter = "\t"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 200
snps$LoadFile(SNP_file_name)

pheno = SlicedData$new()
pheno$fileOmitCharacters = "NA"
pheno$fileDelimiter = "\t"
pheno$fileSkipRows = 1
pheno$fileSkipColumns = 1
pheno$fileSliceSize = 200
pheno$LoadFile(pheno_fl)

cvrt = SlicedData$new()
cvrt$fileOmitCharacters = "NA"
cvrt$fileDelimiter = "\t"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$fileSliceSize = 200
cvrt$LoadFile(covariates_fl)


quest_qtl = Matrix_eQTL_engine(
  snps = snps,
  gene = pheno,
  cvrt = cvrt,
  output_file_name = output_fl,
  pvOutputThreshold = p.val,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = T,
  min.pv.by.genesnp = F,
  noFDRsaveMemory = T,
  pvalue.hist = T
  )
