library(data.table)
library(tidyverse)
library(rstatix)
library(ggpubr)
source("../resources/Functions.R")
library(magrittr)

#General parameters
#Colors
anno_col = list(Time = c("t0" = "lightgrey", "t6" = "grey"),
                `Cell system` = c("pbmc" = "#a0711c","wba" = "#933b27",
                                  "PBMC" = "#a0711c","WBA" = "#933b27"),
                Cytokine = c("il1ra" = "#ec9488","il1b" = "#eb5a46","il10" = "#b04632","il6" = "#f5d3ce",
                             "IL1ra" = "#ec9488","IL1b" = "#eb5a46","IL10" = "#b04632","IL6" = "#f5d3ce"),
                Stimulation = c("p3c" = "#fce8d2" , "lps" = "#fdc788", "cand" = "#ffab4a", "bbmix" = "#d29034",
                                "bafzelii" = "brown", "rpmi" = "white"))

#Read the set of all GW snps and their summary statistics
ind_snp <- fread("../Output/FUMA_cQTL/snps.txt")
#Remove those that are not in the original sum stats
ind_snp <- filter(ind_snp,!is.na(gwasP) )

#Get P-vals
lognorm <- fread("../Output/FUMA_cQTL/allsnps.pval.tsv")
ind_lognorm <- lognorm%>%
  filter(rsID %in% ind_snp$rsID)%>%
  cyt_split("gene")

#Add the genomic locus
ind_lognorm <- left_join(ind_lognorm, ind_snp[,c("rsID","GenomicLocus")])

#Add MAF
snpid <- unique(ind_lognorm$chrpos)
MAF <- fread(cmd = paste0("zgrep -h -w ",paste(paste("-e", snpid), collapse = " "),
                          " ../PostImp/*.info.gz"))
colnames(MAF) <- c("chrpos","Ref","Alt","Alt_frq","MAF","AvgCall","Rsq","Genotyped","-","--","---","----","-----")

ind_lognorm <- left_join(ind_lognorm, 
                         MAF[,c("chrpos","MAF")])

#Add number of patients used to calculated number of alleles
Cytokines <- fread("../Phenotype/Cytokines_goodqual.tsv")%>%
  as.data.frame()%>%
  column_to_rownames("V1")
NCytokines <- apply(Cytokines, 1, function(x)sum(!is.na(x))) 
names(NCytokines) <- gsub("_10_","-10-", names(NCytokines)) %>% gsub("_moi","-moi", .)
NCytokines <- data.frame("gene" = names(NCytokines), "Npatients" = NCytokines)
ind_lognorm <- left_join(ind_lognorm, NCytokines, "gene")

ind_lognorm$Nalleles <- ind_lognorm$MAF * ind_lognorm$Npatients


#A generic function to plot consistency between two groups of a condition

Consistency <- function(df, var, binary = F, g1, g2){
  #In case it is a data.table
  df = as.data.frame(df)
  #If we only want to compare two groups of the variable, set binary = T
  if(binary) df <- filter(df, grepl(paste(g1,g2, sep = "|"), get(var)))
  #Select the variables we can group by
  all.vars <- c("Time","Sample","Cytokine","Stimulation")
  use.vars <- all.vars[all.vars != var]
  #The levels in that variable
  levels <- unique(df[,var])
  #Rotate the df
  df <- df %>% 
    group_by(.dots = c(use.vars,"SNP"))%>%
    mutate(minNalleles = min(Nalleles))%>%#Transform number of alleles to minimum namber of alleles between timepoints
    pivot_wider(id_cols = c(all_of(use.vars),"rsID","SNP","minNalleles","GenomicLocus"),names_from = all_of(var),names_sep = ".", 
                values_from = c(beta,`p-value`,`t-stat`))
  #Only snps that are gw in one of the two
  df <- df %>%
    filter(get(paste0("p-value.", levels[1])) < 5e-8 | get(paste0("p-value.", levels[2])) < 5e-8)
  
  #The dotplot comparing both conditions
  if(var == "Cytokine") {
    df$Cytokine = as.factor(df$GenomicLocus)
    colors =c("7" = "skyblue","3" = "darkgreen","16" = "darkred","26" = "black","2" = "darkblue", "19" = "orange", "grey" )
  }else{
    colors = anno_col$Cytokine[grepl("il", names(anno_col$Cytokine))]
  }
  
  p1 <- ggplot(df)+
    geom_rect(data = data.frame(),aes(xmin = c(-Inf,0), xmax = c(0,Inf), ymin = c(Inf,0), ymax = c(0,-Inf), fill = ""), alpha = 0.1)+
    scale_fill_manual(values = "grey")+
    geom_vline(aes(xintercept = 0), linetype = "dashed")+geom_hline(aes(yintercept = 0), linetype = "dashed")+
    geom_point(aes(x=get(paste0("t-stat.", levels[1])), y= get(paste0("t-stat.", levels[2])), color = Cytokine,  size=minNalleles))+
    ylim(c(-12,12))+xlim(c(-12,12))+
    scale_color_manual(values = colors)+
    theme_bw()+coord_fixed()+theme(panel.grid = element_blank())+
    xlab(paste("T-stat", toupper(levels[1])))+
    ylab(paste("T-stat", toupper(levels[2])))+scale_size_continuous(range = c(0,4))

  p1
  
  #The dotplot comparing minimum T with number of alleles
  p2 <- df%>%
    rowwise()%>%
    mutate("t" = min(abs(get(paste0("t-stat.", levels[1]))),abs(get(paste0("t-stat.", levels[2])))))%>%
    ggplot(aes(y=t,x=minNalleles, size=minNalleles, color = Cytokine))+
    scale_color_manual(values = colors)+
    geom_point()+
    theme_bw()+theme(panel.grid = element_blank())+xlab('Minimum number of alleles')+ylab("Minimum T-statistic")+
    geom_smooth(aes(group = "a"), color ="darkgrey", linetype = "dashed")+scale_size_continuous(range = c(0,4))
  p2
  
  library(patchwork)
  p1/p2+theme(legend.position = 0)+
    plot_layout(heights = c(2,1))
}


Consistency(ind_lognorm, "Time")
ggsave("../code_figures/Figures/Consistency_time.pdf", width = 7, height = 6)
Consistency(ind_lognorm, "Cytokine", T, "il10", "il6")
ggsave("../code_figures/Figures/Consistency_IL6_IL10.pdf", width = 7, height = 6)
Consistency(ind_lognorm, "Stimulation", T, "bbmix-10-5", "p3c")
ggsave("../code_figures/Figures/Consistency_p3c_bbmix.pdf", width = 7, height = 6)


Consistency(ind_lognorm, "Stimulation", T, "bbmix-10-5", "bbmix-10-6")
Consistency(ind_lognorm, "Sample")
Consistency(ind_lognorm, "Cytokine", T, "il10", "il6")
Consistency(ind_lognorm, "Cytokine", T, "il1ra", "il6")
Consistency(ind_lognorm, "Cytokine", T, "il10", "il1b")






#Read independent significant snps
ind_snp.1 <- fread("../Output/FUMA_cQTL/IndSigSNPs.txt")
ind_lognorm <- fread("../Output/FUMA_cQTL/allsnps.pval.tsv")%>%
  filter(rsID %in% ind_snp.1$rsID)


ind_lognorm <- ind_lognorm %>% cyt_split("gene")


#Consistency heatmaps
#Consistency rates per timepoint.

#How I calculate it. Per cytokine stimulation and snp, I calculate the number of significant hits.
#2 indicated both timepoints are significant, 1 indicates only one and therefore it is not replicated.
#I subtract 1 converting it to 1 = replicated, 0 = not replicated
#Sum the number of 1s and divide it by the total
time_con <- lapply(c(0.05, 0.01,0.001,1e-4,1e-5,1e-6,1e-7,5e-8), function(x){
  d = ind_lognorm %>%
    filter(`p-value` < x)%>%
    group_by(notime, rsID)%>%
    summarise(n = n())%>%
    group_by(notime)%>%
    summarise(rep = sum(n-1)/(n()))
  colnames(d) <- c("notime", x)
  return(as.data.frame(d))
}) %>% reduce(full_join)%>%
  column_to_rownames("notime")


rownames(time_con) <- gsub("-", "_",rownames(time_con))
ra = ComplexHeatmap::rowAnnotation(`Cell system`= str_split(rownames(time_con), "_", simplify = T)[,1],
                                   Cytokine = str_split(rownames(time_con), "_", simplify = T)[,2],
                                   Stimulation = str_split(rownames(time_con), "_", simplify = T)[,3],
                                   col = anno_col)

library(circlize)
#pdf("../Plots/cQTL_post/Time_consistency_heatmap.pdf", width = 6, height = 10)
ComplexHeatmap::Heatmap(as.matrix(time_con), show_row_dend = F, show_column_dend = F, cluster_columns = F,
                        col =colorRamp2(c(0,1),c("white","darkgreen")), cluster_rows = T,
                        right_annotation = ra, show_row_names = F)

#dev.off()





#Barplot proportion of replicated across condition
#Considering the GW significant snps, how many can be replicated at p.value < 0.05 in the different conditions?
#Select the variables we can group by

lapply(c("Time","Sample","Cytokine","Stimulation"), function(var){
  #Code to group by the other variables
  all.vars <- c("Time","Sample","Cytokine","Stimulation")
  use.vars <- all.vars[all.vars != var]
  
  #Generate a data.frame indicating the number of groups in the condition in which a snp-cytokine pair is significant
  df <- lognorm%>%
    cyt_split("gene") %>%
    filter(`p-value` < 0.05)%>%
    group_by(.dots = c(use.vars,"rsID"))%>%
    summarise(n = n())%>%
    group_by(n)%>%
    summarise(sum = n())%>%
    mutate("Condition" = var)
})%>%bind_rows()%>%
  mutate(Condition = gsub("Sample", "Cell system", Condition) %>% factor(levels = c("Cell system", "Cytokine","Time","Stimulation")))%>%
  mutate(n = factor(n, levels = max(n):min(n)))%>%
  ggplot(aes(y=sum, x = Condition, fill = n))+geom_col(position = "fill")+coord_flip()+theme_minimal()+
  scale_fill_manual(values = gray.colors(9, rev = T))+ theme(panel.grid = element_blank())
#ggsave("../Plots/cQTL_post/Consistency_barplot_condition.pdf", width = 7, height  = 3)

  
  
