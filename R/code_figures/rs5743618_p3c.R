library(data.table)
library(tidyverse)
library(DESeq2)


#Read the genotype data for the snp of interest
geno <- fread(cmd = "egrep '#CHROM|4:38797027' ../../dosage.txt")
geno <- select(geno, -paste0("V",1:5))
geno <- as.data.frame(geno)

#Add the sample names
samples <- fread("../../sample.id.txt", header = F)
colnames(geno) <- gsub("^.*_0", "", samples$V1)

#Read the RNAseq data
RNA <- fread("../../transcriptome_matrix.txt")
RNA <- as.data.frame(RNA)
#remove tension
RNA <- RNA[,!colnames(RNA) %in%  outlier]


#Transfer samplenames
library(readxl)
Genotype_pheno <- read_excel("../../Genotype_pheno.xlsx")
Genotype_names <- Genotype_pheno$PatientID[Genotype_pheno$ProbenID == colnames(geno)]

#Get phenotype
Phenotype <- read_excel("../../Phenotype.xlsx")%>%
  select(Patient, Age, Gender)%>%
  unique()

#Subset genotype to those with patientIDs
geno <- geno[,Genotype_names != "NA"]
colnames(geno) <- Genotype_names[Genotype_names != "NA"]

#Add geno to phenotype
geno <- t(geno)%>%
  as.data.frame()%>%
  rownames_to_column("Patient")%>%
  dplyr::rename("geno" = V1)

Phenotype <- left_join(geno, Phenotype)
Phenotype <- Phenotype %>%
  column_to_rownames("Patient")

#Remove the mannlich so that DEseq does not complain
Phenotype$Gender <- gsub("ä","a", Phenotype$Gender)

#Convert geno in three factors
Phenotype$geno <- factor(ifelse(Phenotype$geno < 0.5, "CC", "CA/AA"), levels = c("CC","CA/AA"))
#Phenotype$geno <- factor(Phenotype$geno, levels = c('CC', 'CA', 'AA') )
#Start the analysis, P3C = 4, CPG = 2, polyIC = 3

  #Character ending to remove (.4)
  str <- "\\.4$"
  #RNA data for that stim
  P3C <- RNA[,c("genes", colnames(RNA)[grepl(str, colnames(RNA))])]%>%
    column_to_rownames("genes")
  colnames(P3C) <- gsub(str,"",colnames(P3C))
  
  #Select only patients present in both RNA and phenotype
  P3C <- P3C[,colnames(P3C) %in% rownames(Phenotype)]
  #Subset the phenotype to only present in both
  Pheno.P3C <- Phenotype[colnames(P3C),]
  P3C <- P3C[,colnames(P3C) %in% rownames(Pheno.P3C)]
  
  
  #Run DEseq
  P3C.DE <- DESeqDataSetFromMatrix(countData = P3C, colData = Pheno.P3C, design = ~Age+Gender+geno)%>%estimateSizeFactors()
  
  #Lets check this manually
  genes <- c("IL10","IL6","IL1B","IL1RN")
  P3C.cyt <- counts(P3C.DE, normalized = T)[RNA$gene_names %in% genes ,]
  #Change ENSG to gene symbol
  genes.convert <- RNA$gene_names
  names(genes.convert) <- RNA$genes
  rownames(P3C.cyt) <- genes.convert[rownames(P3C.cyt)]
  
  P3C.df <- as.data.frame(t(P3C.cyt))
  P3C.df <- cbind(P3C.df, Pheno.P3C)
  
  P3C.df <- pivot_longer(P3C.df, cols =genes,names_to = "gene", values_to = "rawcounts")
  
  library(ggpubr)
  p <- ggboxplot(P3C.df, x = "geno",y = "rawcounts", fill  = ggsci::pal_npg()(3)[1], add = c("point"), add.params = list(position = position_jitter(w = 0.05)))+
    facet_wrap(~gene, scales = "free") + theme_bw()+theme(panel.grid = element_blank(), strip.background = element_blank(),axis.line = element_line(colour = "black"),
							  panel.border = element_blank(),panel.background = element_blank())+
    stat_compare_means(comparisons = list(c(1,2)))+ylab('Normalized counts')+xlab("rs5743618")
  
  P3C.deseq <- DESeq(P3C.DE)
  res <- results(P3C.deseq)
  res <- res[order(res$pvalue),]
  res.df <- as.data.frame(res@listData)
  res.df$gene <- res@rownames
  #Add gene_symbol
  genes.convert <- RNA$gene_names
  names(genes.convert) <- RNA$genes
  res.df$gene_name <- genes.convert[res.df$gene]
  

ggsave("../code_figures/Figures/rs5743618_P3C_boxplot.pdf",p)


#Pathways enriched
library(clusterProfiler)
up <- bitr(res.df %>% filter(padj < 0.05 & log2FoldChange > 0) %>% pull(gene),
                 fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% pull(ENTREZID)
down <- bitr(res.df %>% filter(padj < 0.05 & log2FoldChange < 0) %>% pull(gene),
           fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% pull(ENTREZID)


enr.GO.up <- enrichGO(gene = up,
                      OrgDb = "org.Hs.eg.db",
                      keyType = "ENTREZID",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      qvalueCutoff = 0.05,
                      readable = TRUE)

enr.GO.down <- enrichGO(gene = down,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 0.05,
                        readable = TRUE)

enr.KEGG.up <- enrichKEGG(gene = up,
                          organism = 'hsa',
                          pAdjustMethod = "BH",
                          pvalueCutoff = 1,
                          qvalueCutoff = 0.05)

enr.KEGG.down <- enrichKEGG(gene = down,
                            organism = 'hsa',
                            pAdjustMethod = "BH",
                            pvalueCutoff = 1,
                            qvalueCutoff = 0.05)

res = list('GO Up' = enr.GO.up,
           'GO Down' = enr.GO.down,
           'KEGG Up' = enr.KEGG.up,
           'KEGG Down' = enr.KEGG.down)

merged <- merge_result(res)

# Construct the plot manually
df <- merged@compareClusterResult


df$old_cluster <- df$Cluster
df <- cbind(df,
            reshape2::colsplit(df$Cluster, pattern = ' ', names = c('enr', 'direction')))

# Mutate the df to the final df we use for plotting
df %<>%
  mutate(GeneRatio2 = sapply(df$GeneRatio, function(x) eval(parse(text=x)))) %>%
  mutate(Description = stringr::str_wrap(.$Description, width = 50)) %>%
  group_by(old_cluster) %>%
  dplyr::slice(1:10)

# order the df as we want and set the description factor last
# to ensure proper plotting order
df$enr <- factor(df$enr, levels = c('GO', 'KEGG'))
df$direction <- factor(df$direction, levels = c('Down', 'Up'))
df %<>% arrange(enr, direction)
df$Description <- str_to_sentence(df$Description) %>% factor(levels = unique(.))

# Hacky stuff
ratio <- nrow(df) * 0.5

print(ggplot(data = df) +
        geom_point(aes(x = direction, y = Description, size = GeneRatio2, color = p.adjust)) +
        DOSE::theme_dose() +
        theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45),
              strip.background = element_rect(fill = 'white', color = 'black'),
              axis.title.y = element_blank(),
              axis.title.x = element_blank(), aspect.ratio = ratio) +
        scale_colour_gradient(limits=c(0, 0.05), low="red", high = 'blue') +
        scale_x_discrete(drop=F) +
        labs(color = 'Adj. P', size = 'Gene ratio') +
        facet_wrap('enr'))


