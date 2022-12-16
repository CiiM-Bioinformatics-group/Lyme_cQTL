library(data.table)
library(tidyverse)
source("../resources/Functions.R")

indloci <- fread("../Output/FUMA_cQTL/Indloci.annotated.txt")

locuszoom <- function(locus_name, gene, point_color){
  #Read locus
  ref <- indloci[indloci$GenomicLocus == locus_name,]
  chr = ref$chr
  min = ref$pos-5e5
  max = ref$pos+5e5

  locus <- fread(cmd = paste0("awk -F ':'  '{if($1 ==",chr,
                              " && $2 > ",min,
                              " && $2 < ",max,
                              ")print}' ../Output/Cytokines_goodqual_pergene/",gene, " | tr ':' '\t' | tr '%' '\t'"))

  colnames(locus) <- c("CHR", "BP", "REF", "ALT","rsID","gene","beta","t","P")

  #Plotting

  p1 <- ggplot(locus, aes(x=as.numeric(BP), y=-log10(P)))+
    geom_point(color = point_color)+
    geom_text(data = locus[which.min(locus$P)], aes(x=as.numeric(BP),y=-log10(P), label=rsID), nudge_y = 0.2)+
    theme_bw()+
    ylab(paste0("Locus ",locus_name, "\n\n-log10(Pvalue)"))+xlab("")+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


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
    geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position)) +
    coord_flip() + ylab("") +
    ylim(plot.range) +
    geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name), fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5) +
    #labs(caption = paste0("Data source: ", gene.ensembl@host, " + Data set: ", gene.ensembl@dataset), color = "Gene Biotype") +
    theme_bw()+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text.y = element_text(angle = 0),
          legend.position="bottom",
          panel.grid.major.y = element_blank()) +
    expand_limits(y=c(-1, 1))+
    scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out.bm.genes.region$gene_biotype_fac)-1)))
  library(patchwork)
  coloc.theme <-theme(plot.margin = margin(0, 0, 0, 0), panel.border = element_blank(), panel.grid = element_blank(),
                      axis.line = element_line(size = .2))
  p1 <- p1+xlim(plot.range)+coloc.theme
  p3 <- p3+coloc.theme+theme(axis.line.y  = element_blank())
  plot <- p1/p3+plot_layout(heights = c(3,1))
  plot
}
#TLR
locuszoom(7, "t6_pbmc_il6_bbmix_10_5", "darkgrey")
ggsave("Figures/Locuszoom_TLR.pdf", width = 7, height = 8)
#KL
locuszoom(26, "t6_wba_il6_bbmix_moi10", "black")
ggsave("Figures/Locuszoom_KL.pdf", width = 7, height = 8)
#EFR3A
locuszoom(19, "t0_pbmc_il1b_bbmix_10_5", "black")
ggsave("Figures/Locuszoom_EFR3A.pdf", width = 7, height = 8)
