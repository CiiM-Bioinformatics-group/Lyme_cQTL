library(data.table)
library(tidyverse)
library(rstatix)
library(ggpubr)
source("../resources/Functions.R")
library(magrittr)
library(parallel)

#Mirror locuszoom, read the details for the TLR locus, in order to make the locuszoom plot for each cytokine
indloci <- fread("../Output/FUMA_cQTL/Indloci.annotated.txt")
#Read locus
ref <- indloci[indloci$GenomicLocus == 8,]
chr = ref$chr
min = ref$pos-5e5
max = ref$pos+5e5

locus <- fread(cmd = paste0("awk -F ':'  '{if($1 ==",chr,
                            " && $2 > ",min,
                            " && $2 < ",max,
                            ")print}' ../Output/Cytokines_goodqual_0.05.txt"))

colnames(locus) <- c("SNP","gene","beta","t","P")
locus <- locus  %>% tidyr::separate(SNP, into = c("CHR","BP","REF","ALT","rsID"), remove = F)%>%cyt_split("gene")

#Only borrelia and p3c
locus <- locus[grepl("bbmix|p3c", locus$gene),]

#Split to get chr and so on
locus$logP <- -log10(locus$P)
locus$logP[locus$Time == "t6"] <- locus$logP[locus$Time == "t6"]*-1

#Plot per gene
coloc.theme <-theme(plot.margin = margin(0, 0, 0, 0), panel.border = element_blank(), panel.grid = element_blank(),
                    axis.line = element_line(size = .2), panel.background=element_blank())

t0 <- ggplot()+
 geom_point(data = locus[locus$Time == "t0",], aes(x=as.numeric(BP), y=logP, color = Cytokine))+
 geom_hline(data = data.frame(), aes(yintercept=c(-log10(5e-8), log10(5e-8)  )), linetype = "dashed")+
 ylim(c(-log10(0.05), max(abs(locus$logP))))+
 ylab(paste0("-log10(Pvalue)"))+xlab("")+
 scale_color_manual(values = anno_col$Cytokine)+
 theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 0)+coloc.theme

 t6 <- ggplot()+
  geom_point(data = locus[locus$Time == "t6",], aes(x=as.numeric(BP), y=logP, color = Cytokine))+
  geom_hline(data = data.frame(), aes(yintercept=c(-log10(5e-8), log10(5e-8)  )), linetype = "dashed")+
  ylim(c(-max(abs(locus$logP)), log10(0.05)))+
  ylab(paste0("-log10(Pvalue)"))+xlab("")+
  scale_color_manual(values = anno_col$Cytokine)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 0)+coloc.theme

#Plot Genes
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
  expand_limits(y=c(-1, 1))+
  theme_void()+theme(legend.position = 0)+
  scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out.bm.genes.region$gene_biotype_fac)-1)))

library(patchwork)
t0<-t0+xlim(plot.range)
t6<-t6+xlim(plot.range)

p <- t0/p3/t6+plot_layout(heights = c(4,1,4))

ggsave("../Output/TLR/TLR_locuszoom_mirror.pdf", p)
print("Mirror saved!")


#Mirror but only Fine-mapping results
susie <- fread("../Output/TLR/finemapping_locus8.tsv")%>%cyt_split("gene")
susie$PIP[susie$Time == "t6"] <- susie$PIP[susie$Time == "t6"]*-1

plot.range2 <- c(min(susie$BP)-1000, max(susie$BP)+1000)

susie0 <- ggplot()+
 geom_point(data = susie[susie$Time == "t0",] , aes(x=as.numeric(BP), y=PIP, color = Cytokine))+
 ylim(c(0, 1.05))+
 xlim(plot.range2)+
 geom_text(data = susie[susie$Time == "t0",] %>% filter(PIP == max(PIP)), aes(x=as.numeric(BP),label = rsID, y = PIP+0.05))+
 ylab(paste0("PIP"))+xlab("")+
 scale_color_manual(values = anno_col$Cytokine)+
 theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 0)+coloc.theme

 susie6 <- ggplot()+
  geom_point(data = susie[susie$Time == "t6",] , aes(x=as.numeric(BP), y=PIP, color = Cytokine))+
  ylim(c(-1.05, 0))+
  scale_x_continuous(position = "top", limits = plot.range2)+
  geom_text(data = susie[susie$Time == "t6",] %>% filter(PIP == min(PIP)), aes(x=as.numeric(BP),label = rsID, y = PIP-0.05))+
  ylab(paste0("PIP"))+xlab("")+
  scale_color_manual(values = anno_col$Cytokine)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 0)+coloc.theme

p3 <- p3+ylim(plot.range2)

psusie <- susie0/p3/susie6+plot_layout(heights = c(4,1,4))
ggsave("../Output/TLR/TLR_susie_mirror.pdf", psusie)
print("Susie saved!")

#Heatmap of independent significant snps and all snps only for the TLR Locus
library(ComplexHeatmap)
library(circlize)

ind_snp.1 <- fread("..//Output/FUMA_cQTL/snps.txt")

#Only independent or exonic SNPs will be shown
ind_snp.1 <- ind_snp.1 %>% filter(GenomicLocus == 8)%>%filter(func == "exonic" | rsID == IndSigSNP)

ind_snp <- paste(ind_snp.1$chr, ind_snp.1$rsID, sep = "_")

#Read all the pvals per snp
ind_lognorm <- fread("../Output/FUMA_cQTL/allsnps.pval.tsv")
#Read GWAS t align directions
#ind_lognorm_gwas <- fread(cmd = paste0("egrep 'CHR|",paste(unique(ind_lognorm$rsID), collapse = "|"), "' ../../Lyme_GWAS/GWAS/LymevsHealthy/LB1_LB2_mymetal_stderr.txt"))
#Risk allele per snp
#ind_lognorm_gwas$GWASRisk <- ifelse(ind_lognorm_gwas$z.meta > 0, ind_lognorm_gwas$A11, ind_lognorm_gwas$A21)
#Bind both
#ind_lognorm<- left_join(ind_lognorm,
#                       rename(ind_lognorm_gwas, "rsID" = SNP1) %>% dplyr::select(rsID,GWASRisk))%>%
#  separate(SNP, into = c("chr","pos","ref","alt","rsID"), remove = F)%>%
#  mutate(beta = ifelse(is.na(GWASRisk) ,beta, ifelse(ref == GWASRisk, -beta, beta )))

#Adjust df
ind_lognorm <- ind_lognorm%>%
  mutate("log.pval" = sign(beta)*-log10(`p-value`))%>% #Reflecting cQTL direction
  dplyr::select(c("SNP","gene","log.pval"))%>%
  separate(SNP, c("rsID"),sep = "_", remove = F) %>%
  pivot_wider(names_from = gene, values_from = log.pval, values_fill = 0)%>%
  separate(rsID, into = c("CHR:BP", "SNP"), "%")%>%
  separate(`CHR:BP`, into = c("chr"), ":")%>%
  unite(chr,SNP,col = "chr_SNP")%>%
  column_to_rownames(var = "chr_SNP")%>%
  t() %>%
  .[,order(as.numeric(str_split(colnames(.), "_", simplify = T)[,1]))]%>%
  .[order(gsub("t._","",rownames(.))),]%>%
  .[!grepl("C6", rownames(.)),]




ra = ComplexHeatmap::rowAnnotation(Time = str_split(rownames(ind_lognorm), "_", simplify = T)[,1],
                                   `Cell system`= toupper(str_split(rownames(ind_lognorm), "_", simplify = T)[,2]),
                                   Cytokine = gsub("il","IL",str_split(rownames(ind_lognorm), "_", simplify = T)[,3]),
                                   Stimulation = str_split(rownames(ind_lognorm), "_", simplify = T)[,4],
                                   col = anno_col)

#ca = HeatmapAnnotation(CHR = anno_text(str_split(colnames(ind_lognorm[,colnames(ind_lognorm) %in% ind_snp]), "_", simplify = T)[,1],
#                                       show_name = T, rot = 0, gp = gpar(fontsize = 8)))
#locia = HeatmapAnnotation(Loci = anno_text(ind_snp.1$GenomicLocus, show_name = T, rot =0, gp = gpar(fontsize = 8)))

pdf("../Output/TLR/TLR_Heatmap_indsig_snps.pdf", width = 7.5, height = 5.5)
ComplexHeatmap::Heatmap(ind_lognorm[,colnames(ind_lognorm) %in% ind_snp], cluster_columns = F, show_row_dend = F, left_annotation = ra,
                        show_row_names = F, cluster_rows = F,
                        col = colorRamp2(c(log10(5e-8),log10(0.05), log10(0.05)+0.0001, 0,  -log10(0.05)-0.0001 ,-log10(0.05), -log10(5e-8)),
                                         c("blue","lightblue","white", "white","white","orange", "red")),
                        #bottom_annotation = ca,show_column_names = F, top_annotation = locia,
                        name = "Signed\nLog10 pvalue")
dev.off()







#There are changes in the t statistics of the TLR locus between time 1 and 6. Lets check if beta is changing or the standard error or beta is changing
locus_split <- locus %>%
                  filter(P < 1e-5) %>% #Only suggestive
                  cyt_split("gene")%>%
                  mutate(se = beta/t)%>% #Add se
                  group_by(notime, SNP) %>% mutate(nt = n(),keep = (nt == 2) )%>%filter(keep)#Fix to get only snps significant in both timepoints

#Paired test of beta values per gene
for (i in c("se","beta","t")){
  print(ggpaired(locus_split,
           x = "Time", y = i, id = "SNP")+facet_wrap(~notime)+stat_compare_means(method = "wilcox",paired = T)+
    ggtitle(paste("Change in",i))+ylab(i))
}

#All right, t changes mainly because the effect, beta, is higher. The SNPs in the TLR locus have a higher effect after antibody treatment

#But, fine mapping is easier because the main snp is still very significant compared to the others. Is the change in beta of rs574 different?
locus_split_change <- locus_split %>%
  pivot_wider(id_cols = c("SNP", "notime"), names_from = "Time", values_from = "beta")%>%
  mutate(beta.change = t6 - t0)
ggplot(locus_split_change)+
  geom_density(aes(x=beta.change, fill = notime))+
  geom_vline(data = locus_split_change %>% filter (SNP == "4:38798648:C:A%rs5743618"), aes(xintercept=beta.change))+
  facet_wrap(~notime, scales = "free_y")+
  theme_minimal()

#Density plot with the time 0, time 6 and the snp in all of them
density_df_median <- locus_split_change%>%
  pivot_longer(cols = c(t0,t6,beta.change), names_to = "time", values_to = "beta")%>%
  group_by(SNP,time) %>% summarise(beta = median(beta))%>%
  mutate(time = factor(time, levels = c("t0","t6", "beta.change")))
ggplot(density_df_median)+
  geom_density(aes(x=beta))+
  geom_vline(data = density_df_median %>% filter (SNP == "4:38798648:C:A%rs5743618"), aes(xintercept=beta))+
  facet_wrap(~time, ncol = 1)+
  theme_minimal()+theme(panel.grid = element_blank())
ggsave("../Output/TLR/TLR_density_median_beta.pdf")


#The change in beta values of the top SNP is almost 0, compared to the other SNPs in the locus, which increase significantly
#look at il10_bbmix_pbmc, the change is 0, while the second part of the distribution changes more


#Next hypothesis, the snps in LD with the top one change more than the top snp, or the correlation in correlations changes

#Extracting the dosage and correlation of locus8
dosage <- fread(cmd = paste0("egrep 'LP|",paste(unique(locus$SNP), collapse = "|"),
                             "' ../data/lyme-merged_dosage.tsv"))
dosage <- as.data.frame(dosage)%>%
  column_to_rownames("V1")
R <- cor(t(dosage))
#Correlation for the first SNP
topR <- R["4:38798648:C:A%rs5743618",]%>%
  as.data.frame()%>%rownames_to_column("SNP")
colnames(topR)<- c("SNP","LD")
#Correlation of beta values
locus_split_change_R <- inner_join(topR, locus_split_change, by = "SNP")

density.margin <- lapply(unique(locus_split_change_R$notime), function(i){
 ggscatterhist(locus_split_change_R %>% mutate(LD = abs(LD))%>%
                        filter(notime == i),title = i, x="beta.change", y = "LD", margin.params = list(fill = "notime"))

})

ggplot(locus_split_change_R, aes(x=beta.change, y= abs(LD)))+
  geom_point()+theme_minimal()+
  geom_text(data  = locus_split_change_R %>%  filter (SNP == "4:38798648:C:A%rs5743618"), aes(label = "rs5743618"))+
    facet_wrap(~notime)+
    geom_smooth(method = "lm")


#The effect of the other SNPs on the main, missense snp
#Only those that are a bit linked to the top one: r2 >0.2
linked <- topR$SNP[topR$LD^2>.2]
#Do not do the r2 filter
linked <- topR$SNP

#No clumping and no r2 filter
locus.clump <- locus[locus$SNP %in% linked,]%>%
  rename("rsid" = rsID, "pval" = P)#%>%
  #ieugwasr::ld_clump(clump_r2 = 0.2)

dosage.clump <- dosage[unique(locus.clump$SNP),]

#LEt's add the covariates
covar <- fread("..//Covariates/covariates.tsv")%>%
  column_to_rownames("V1")

#Let's check the possible interactions
#For that, we need the cytokine values of a cytokine significantly affected by the TLR locus
#Let's add also the average of the four cytokines per measure

cyt.full <- fread("../Phenotype/Cytokines_goodqual.tsv")%>%
  column_to_rownames("V1")
cyt.full <- cyt.full[gsub("-", "_", unique(locus$gene)),]

#Create new measures, averaging the cytokine values for each condition
newcyts <- gsub("_il(6|1ra|10|1b)_", "_mean_", rownames(cyt.full))
newcyts.df <- lapply(unique(newcyts), function(x){
  tmp <- cyt.full[newcyts %in% x,]
  #Average
  tmp <- as.data.frame(t(apply(tmp, 2, function(y) mean(na.omit(y)))))
  rownames(tmp) <- x
  return(tmp)
})%>%bind_rows()
cyt.full <- rbind(cyt.full, newcyts.df)

print("computing interaction")
interaction <- mclapply(rownames(cyt.full), function(g){
  #Read the cytokine
  g <- gsub("-","_", g)
  print(g)
  cyt <- cyt.full[g,]
  all(colnames(cyt) == colnames(dosage.clump))

  #Let's map the interaction of both into one

  interaction.gene <- lapply(rownames(dosage.clump), function(x){
      print(x)
      #Check the LD of the SNP, if there is a negative correlation, flip it
      tmpR <- topR$LD[topR$SNP == x]
      if(tmpR < 0) test <-  2-as.numeric(dosage.clump[x,]) else test = as.numeric(dosage.clump[x,])
      tmp <- rbind(dosage["4:38798648:C:A%rs5743618",], test, cyt, covar)%>%
        t()%>%  as.data.frame()
      colnames(tmp) <- c("snp","test","cyt", rownames(covar))
      sum <- summary(lm(cyt ~ test*snp+age+sex+Batch+Institute, tmp))
      if(x == "4:38798648:C:A%rs5743618") return(data.frame("SNP" = x, "gene" = g, "beta" = NA, "P" = NA))
      return(data.frame("SNP" = x, "gene" = g, "beta" = sum$coefficients[8,1], "P" = sum$coefficients[8,4]))
      })%>%bind_rows()
  }, mc.cores = 30) %>% bind_rows()
fwrite(interaction, "../Output/TLR/TLR_interaction.txt", sep = "\t")

interaction <- fread("../Output/TLR/TLR_interaction.txt")%>%
  separate(SNP, into = c("CHR","BP","REF","ALT", "rsID"), remove = F)

#Locuszoom plot

locuszoom.plot <- function(df){
  ggplot(df, aes(x=as.numeric(BP), y=-log10(P)*sign(beta)))+
    geom_point(color = ggsci::pal_npg()(3)[2])+
    geom_vline(aes(xintercept = 38798648), linetype = "dashed")+
    geom_hline(aes(yintercept = c(-log10(0.05))), linetype = "dashed")+
    geom_hline(aes(yintercept = c(log10(0.05))), linetype = "dashed")+
    geom_text(data = df[df$rsID == df$rsID[which.min(df$P)]],
              aes(label=rsID))+
    theme_bw()+
    facet_wrap(~gene)+
    ylab("-log10(Pvalue)")+xlab("")+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
locuszoom.plot(interaction%>% filter(!grepl("mean", gene)))
locuszoom.plot(interaction%>% filter(grepl("mean", gene)))
#Locuszoom plot: only pbmc_bbmix_10_5
locuszoom.plot(interaction %>% filter(grepl("pbmc_il.+_bbmix_10_5", gene)))



#Only SNPs with a significant hit
nom.hit <- c(unique(interaction$rsID[interaction$P < 5e-4]), "rs5743618")
#Add the names
interaction.mat <- interaction %>%
  filter(rsID %in% nom.hit)%>%
  mutate(logP = -log10(P)*sign(beta))%>%
  select(SNP, gene, logP)%>%
  pivot_wider(names_from = SNP, values_from = logP)%>%
  column_to_rownames("gene")%>%
  as.matrix()

#Order interactions
#interaction <- interaction[order(str_split(colnames(interaction), pattern = ":", simplify = T)[,2])]

pdf("../Output/TLR/TLR_interaction_heatmap.pdf", width = 20, height = 12)
ComplexHeatmap::Heatmap(interaction.mat, cluster_columns = F ,
                        col = circlize::colorRamp2(c(log10(1e-4),log10(0.05), log10(0.05)+0.0001, 0,  -log10(0.05)-0.0001 ,-log10(0.05), -log10(1e-4)),
                                                      c("blue","lightblue","white", "white","white","orange", "red")))

dev.off()
#Clumped version
interaction.mat <- interaction %>%
  mutate(logP = -log10(P)*sign(beta))%>%
  select(SNP, gene, logP)%>%
  pivot_wider(names_from = SNP, values_from = logP)%>%
  column_to_rownames("gene")%>%
  as.matrix()

#Clumped version
#Clump SNPs
interaction.clump <- data.frame("rsid" = as.character(str_split(colnames(interaction.mat), pattern = "%", simplify = T)[,2]),
                                "pval" = apply(interaction.mat, 2, function(x) min(10^-abs(x))))%>%
  rownames_to_column("SNP")%>%
  na.omit()%>%
  ieugwasr::ld_clump(clump_r2 = 0.5)

interaction.mat.clump <- interaction.mat[,c(interaction.clump$SNP, "4:38798648:C:A%rs5743618")]

#Order interactions
#interaction <- interaction[order(str_split(colnames(interaction), pattern = ":", simplify = T)[,2])]

ComplexHeatmap::Heatmap(interaction.mat.clump, cluster_columns = F ,
                        col = circlize::colorRamp2(c(log10(1e-4),log10(0.05), log10(0.05)+0.0001, 0,  -log10(0.05)-0.0001 ,-log10(0.05), -log10(1e-4)),
                                                   c("blue","lightblue","white", "white","white","orange", "red")))

#4:38815502:A:C%rs4833103 and 4:38830350:A:G%rs5743810 are the more significant
#cyt <- fread(cmd = paste0("egrep 'LP|t0_pbmc_il6_bbmix_10_5' ../Phenotype/Cytokines_goodqual.tsv"))%>%
#  column_to_rownames("V1")
#int <- rbind(dosage["4:38798648:C:A%rs5743618",], as.numeric(dosage.clump["4:38815502:A:C%rs4833103",]), cyt)%>%
#  t()%>%  as.data.frame()
#colnames(int) <- c("snp","test","cyt")

#int$snp <- factor(ifelse(int$snp < 0.5, "CC", ifelse(int$snp < 1.5, "CA", "AA")), levels = c("CC","CA","AA"))
#int$test <- factor(ifelse(int$test < 0.5, "AA", ifelse(int$test < 1.5, "AC", "CC")), levels = c("AA","AC","CC"))

#ggboxplot(int, x="snp", y="cyt", color = "test", add = "jitter")+
#  theme_minimal()
