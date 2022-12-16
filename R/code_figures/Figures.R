library(data.table)
library(tidyverse)
library(rstatix)
library(ggpubr)
source("../resources/Functions.R")
library(magrittr)
setwd("../code_figures/")

#General parameters
#Colors
anno_col = list(Time = c("t0" = "lightgrey", "t6" = "grey"),
                `Cell system` = c("pbmc" = "#a0711c","wba" = "#933b27",
                                  "PBMC" = "#a0711c","WBA" = "#933b27"),
                Cytokine = c("il1ra" = "#ec9488","il1b" = "#eb5a46","il10" = "#b04632","il6" = "#f5d3ce",
                             "IL1ra" = "#ec9488","IL1b" = "#eb5a46","IL10" = "#b04632","IL6" = "#f5d3ce"),
                Stimulation = c("p3c" = "#fce8d2" , "lps" = "#fdc788", "cand" = "#ffab4a", "bbmix" = "#d29034",
                                "bafzelii" = "brown", "rpmi" = "white"))


#Loading and formating
#Read each cytokine file that was previously formatted for MatrixEqtl and bind them
Cyt1 <- as.data.frame(fread("../old/QTL/Phenotype/Cytokines_questionnaires_LB1.tsv"))

Cyt_ <- read.csv2("../Raw/Phenotype/LymeProspect_cytokines_clinical_traits.csv")%>%
  .[,-1]%>%
  mutate(sex = factor(sex,levels = c("F","M"), labels = c(0,1)))%>%
  column_to_rownames("Lab_ID")%>% 
  t()%>%  as.data.frame()

#Make sure only cytokines are present
Cyt <- Cyt_[grepl("pbmc|wba", rownames(Cyt_)),]

Cyt2.0 <- as.data.frame(fread("../Phenotype/Cytokines.tsv"))%>%column_to_rownames(colnames(.)[1])

#Transpose
Cyt_t <- as.data.frame(t(Cyt)) %>% mutate_all(as.numeric)


#Add cols cytokine, stimulation, Sample
Cyt_df <- Cyt_t %>% rownames_to_column("patient") %>% melt()%>%
  separate(variable, into = c("Time","Sample","Cytokine","Stimulation"), sep = "_", remove = F) %>%
  mutate("notime" = gsub("t\\d_","",variable))

#Add notime variable (Cytokine and stimulation), calculate normality and order by normality
Cyt_df %<>% group_by(notime) %>%
  mutate("norm" = shapiro.test(log2(value))$p.value, "notime" = as.factor(notime))%>%
  ungroup()%>%
  mutate(notime =  fct_reorder(notime,norm, min))

#Add logvalue and iqr for boxplot
Cyt_df <- Cyt_df  %>% mutate("logvalue" = log2(value))%>% group_by(variable) %>%
  mutate("iqr" = IQR(na.omit(value)), "iqrlog" = IQR(na.omit(logvalue)), "intvalue" = INT(value)) %>%
  ungroup()


#Figure 1 - Cytokine heatmap and PCA
# Figure 1.B - Heatmap of correlation
library(ComplexHeatmap)
library(circlize)
#First with time0
Cyt_t0 <- Cyt_t[,grepl("t0", colnames(Cyt_t))]
dist_cyt_t0 <- cor(as.matrix(Cyt_t0),use = "pairwise.complete.obs",method = "spearman")

#Then time6
Cyt_t6 <- Cyt_t[,grepl("t6", colnames(Cyt_t))]
dist_cyt_t6 <- cor(as.matrix(Cyt_t6),use = "pairwise.complete.obs",method = "spearman")

#Annotate measures
names_split <- str_split(colnames(Cyt_t0), "_", simplify = T)
row_ha <- rowAnnotation(Cytokine = gsub("il","IL",names_split[,3]),
                        #Time = names_split[,1],
                        `Cell system` = toupper(names_split[,2]),
                        Stimulation = names_split[,4],
                        col = anno_col)
#Heatmap t0
h1 =  Heatmap(dist_cyt_t0, row_title  = "", column_title = "Baseline",
              show_column_dend = F, show_column_names = F, show_row_names = F,
              show_row_dend = F,
              col = colorRamp2(c(1,0,-1), c("red", "white","blue")),
              left_annotation = row_ha, name = "Spearman\ncorrelation")
#Heatmap t6
h2 = Heatmap(dist_cyt_t6, row_title  ="", show_column_dend = F, column_title = "After treatment",
             show_column_names = F, show_row_names = F,
             show_row_dend = F,
             col = colorRamp2(c(1,0,-1), c("red", "white","blue")),  column_order = column_order(h1),
             show_heatmap_legend = F)
pdf("Figures/Figure1.Correlation_heatmap.pdf", width = 12, height = 8)
set.seed(222)
print(h1 + h2)
dev.off()


#Figure 1.C - PCA plot

#Imputation - needed for PCA
#Remove two stimulations that are exclusive
Cyt_t2 <- Cyt_t[!grepl("bafzelii|bbmix_10_3", colnames(Cyt_t))]

#Remove samples with too many missing
discard <- apply(Cyt_t2,1, function(x) sum(is.na(x))) %>% .[.>30] %>% names()
Cyt_t2 <- Cyt_t2[!rownames(Cyt_t2) %in% discard,]

##bnstruct will impute missing values based on the knn algorithm with k=20
library('bnstruct')
Cyt_t_imputed <- knn.impute(as.matrix(Cyt_t2),  k = 10)

#Use only normally distributed measures, which will be used for QTL
set.seed(123)
Cyt_PC <- remove_badqual(t(Cyt_t_imputed)) %>%
  apply(2, function(x) log2(x))%>%
  apply(1, function(x) scale(x))%>%t()%>%
  #apply(2, function(x) scale(x))%>%
  prcomp()

#Adjust df
PCdf <-Cyt_PC$x %>%
  as.data.frame()%>%
  mutate("label" = rownames(remove_badqual(t(Cyt_t_imputed))))%>%
  cyt_split("label")%>%
  rename("Cell system" = Sample)

#Get % of variance
sum <- summary(Cyt_PC)

#Plot PCA
PCplot <- PCdf%>%
  ggplot(aes(x=PC1,y=PC2,color =Cytokine,shape=`Cell system`))+geom_point()+
  scale_color_manual(values = anno_col$Cytokine)+
  theme_bw()+#xlim(c(-25,25))+
  xlab(paste0("PC1 (",sum[["importance"]][2,1]*100,"% of variance)"))+
  ylab(paste0("PC2 (",sum[["importance"]][2,2]*100,"% of variance)"))+
  theme(plot.margin = margin(), panel.grid = element_blank())

PCdens <- PCdf%>%
  ggplot()+geom_density(aes(x=PC1, fill = Cytokine),color=NA, alpha = .9)+
  theme_bw()+#xlim(c(-25,25))+
  scale_fill_manual(values = anno_col$Cytokine)+
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position = 0,
        axis.text.x = element_blank(), plot.margin = margin(), panel.grid = element_blank())+
  xlab("")

library(patchwork)
PCdens/PCplot+plot_layout(heights = c(1,3))

#Which PCs capture cytokines?
for(i in c(60,62,64,65)){
  print(colnames(PCdf)[i])
  print(cor(PCdf$PC4, as.numeric(as.factor(PCdf[,i]))))
}


ggpubr::ggscatterhist(data = PCdf, x = 'PC1', y = 'PC2',
                      color = 'Cytokine', fill = 'Cytokine',palette =  anno_col$Cytokine,
                      alpha = 1, margin.plot = 'density', legend = 'right', title = '',
                      margin.params = list(fill = "Cytokine", color = 'black', size = 0.2))


ggsave("Figures/Figure1.PCA_density.pdf", width = 8, height = 8)



#Figure 2 - Phenotypic and disease effects on cytokine production

#Formating specific for this analysis
Cyt_df.2 <- Cyt_df %>% ind_N(column = "variable")%>% ##Remove patient and measurements in which one of the timepoints is missing
  group_by(patient, notime) %>% #paired test will fail if you dont do that
  mutate("na" = any(is.na(value))) %>%
  filter(na == FALSE)

#Edit text for figures
edit.Stimulation <- function(x){
  x =  gsub("bbmix-","Borrelia mix ", x) %>% gsub(" 10-", " 10^", .) %>% gsub("rpmi", "RPMI", .) %>%
    gsub("p3c", "Pam3Cys", .) %>%  gsub("lps", "LPS", .) %>% gsub("bafzelii", "Borrelia afzelii", .) %>% gsub("cand", "C.albicans", .)
  x = as.factor(x)
}
edit.Cytokines <- function(x){
  x = gsub("il10","IL10", x) %>% gsub("il1b", "IL1b", .)%>%gsub("il6", "IL6", .)%>% gsub("il1ra", "IL1Ra", .)
}
edit.Sample <- function(x){
  x = gsub("pbmc","PBMCs",x)%>%gsub("wba", "Whole blood", .)
}
edit.Time <- function(x){
  x = gsub("t0", "Baseline", x)%>% gsub("t6", "After treatment", .)
  x = factor(x, levels = c("Baseline","After treatment"))
}


#Perform paired test, fdr adjusted
paired.test <- compare_means(value ~ Time, data = Cyt_df,
                             group.by = "notime", paired = T,
                             method = "wilcox.test", ref.group = "t0",
                             p.adjust.method = "fdr")


#Get median difference between both
diff = Cyt_df.2 %>%
  group_by(notime, patient) %>%
  pivot_wider(id_cols = c("patient","notime"),names_from = "Time", values_from = "value")%>%
  mutate("diff" = t6 - t0)%>%
  group_by(notime) %>% summarise(diff = median(diff))%>%
  filter(notime %in% paired.test$notime)


p <- list()

p$treatment <- paired.test %>%
  left_join(diff,"notime")%>% #Add diff
  mutate("P" = p.adj)%>%
  mutate("test" = gsub("bbmix_10_","bbmix-10-", notime)%>%gsub("bbmix_moi","bbmix-moi",.))%>%
  separate("test", into = c("Sample", "Cytokine", "Stimulation"), "_") %>%
  mutate("values" = -log10(as.numeric(P)),
         "signif" = ifelse(P>0.05, "",ifelse(P>0.01,"*",ifelse(P>0.001,"**","***"))))%>%
  mutate(Stimulation = edit.Stimulation(Stimulation), Cytokine = edit.Cytokines(Cytokine), Sample = edit.Sample(Sample))%>%
  ggplot(aes(x=Cytokine, y=Stimulation, fill=values*sign(diff)))+
  geom_tile(color=NA)+ggtitle("Treatment effect")+facet_grid(rows = vars(Sample),scales = "free", switch="both")+
  geom_text(aes(label = signif))+
  scale_color_manual(values = "black")+
  scale_fill_gradient2(limits = c(-4,4),oob=scales::squish,
                       high = scales::muted("red"),
                       mid = "white",
                       low = scales::muted("blue"),
                       midpoint = 0,
                       space = "Lab",
                       na.value = "grey50",
                       guide = "colourbar",
                       aesthetics = "fill",
                       name = "Signed\nlog10(p value)"
  )+ theme_bw()+
  theme(axis.ticks = element_blank(), strip.background = element_blank(), panel.grid = element_blank())


#Add other info to test

Cyt_t_age_sex <- Cyt_ %>% t() %>% as.data.frame() %>% rownames_to_column('patient') %>%
  mutate("institute" = ifelse(patient %like% "LP-A", "Amsterdam", "Radboud"), "batch" = ifelse(patient %in% colnames(Cyt1), "LB1", "LB2"))%>%
  mutate_at(c(colnames(.)[grepl("pbmc|wba", colnames(.))], "C6_t0", "C6_t6", "age"), as.numeric)%>%
  mutate('diss_or_EM' = factor(diss_or_EM, levels = c("EM","DISS")))

Cyt_t_age_sex <- Cyt_t_age_sex %>%
  left_join(fread('../Raw/Phenotype/LymeProspect_GWAS_BloodbeforeAB.csv')%>% as.data.frame()%>%dplyr::rename('patient'=Lab_ID))


for(trait in c("age","sex", "C6_t0", "C6_t6", "sex","diss_or_EM", 'EMdiam','y_combined_s2_f2_6mt', 'BloodbeforeABstart_atBL')){
  plist <- sapply(Cyt_t_age_sex[,grepl("pbmc|wba", colnames(Cyt_t_age_sex))], function(y){
    
    x <- Cyt_t_age_sex[,trait]
    
    if(is.numeric(x)){
      res <- cor.test(y, x, use = "pairwise.complete.obs", method = "spearman")
      value <- cor(y, x, use = "pairwise.complete.obs", method = "spearman")
    }else if(length(unique(na.omit(x))) == 2){
      df <- na.omit(data.frame("x" = x, "y" = y))
      res <- wilcox.test(y ~ x, df)
      value <- log(mean(df$y[df$x == unique(x)[2]])/mean(df$y[df$x == unique(x)[1]]))
    }else{
      df <- na.omit(data.frame("x" = x, "y" = y))
      res <- kruskal.test(y ~ x, df)
      value <- lm(log2(y) ~ as.numeric(as.factor(x)), df)$coefficients[2]
    }
    
    a <- res
    return(data.frame("sign" = sign(value),
                      "praw" = a$p.value))
    
  }) %>% t() %>% as.data.frame()
  
  #Apply fdr correction
  plist$padj = p.adjust(plist$praw, method = 'fdr')
  #Add significance
  plist <- mutate(plist, "p" = ifelse(padj>0.05, "",ifelse(padj>0.01,"*",ifelse(padj>0.001,"**","***"))))
  #Add value
  plist$values <- unlist(plist$sign)*(-log10(plist$padj))
  
  #plot
  p[[trait]] <- plist %>% rownames_to_column(var = "test") %>%
    mutate("test" = gsub(".rho", "", test)%>%gsub("bbmix_10_","bbmix-10-", .)
           %>%gsub("bbmix_moi","bbmix-moi",.))%>%
    separate("test", into = c("Time", "Sample", "Cytokine", "Stimulation"), "_") %>%
    mutate("values" = as.numeric(values))%>%
    
    mutate(Stimulation = edit.Stimulation(Stimulation), Cytokine = edit.Cytokines(Cytokine),
           Sample = edit.Sample(Sample), Time = edit.Time(Time))%>%
    ggplot(aes(x=Cytokine, y=Stimulation, fill=values))+
    geom_tile(color = NA)+ggtitle(str_to_sentence(trait))+facet_grid(rows = vars(Sample),cols = vars(Time), scales = "free", switch="both")+
    geom_text(aes(label = p))+
    scale_color_manual(values = "black")+
    scale_fill_gradient2(limits = c(-4,4),oob=scales::squish,
                         high = scales::muted("red"),
                         mid = "white",
                         low = scales::muted("blue"),
                         midpoint = 0,
                         space = "Lab",
                         na.value = "grey50",
                         guide = "colourbar",
                         aesthetics = "fill",
                         name = "Signed\nlog10(p value)"
    )+ theme_bw()+
    theme(axis.ticks = element_blank(), strip.background = element_blank(), panel.grid = element_blank())
  #ggsave(plot = p[[trait]],filename =  paste0("../Plots/cQTL_descriptive/",trait,"_cyt_cor.pdf"), height = 7, width = 6)
}

#Correlate C6 levels with each timepoint
p$C6_t0$data <- p$C6_t0$data %>% filter(Time == "Baseline")
p$C6_t6$data <- p$C6_t6$data %>% filter(Time == "After treatment")

#Figure 2.A Phenotype effect
library(patchwork)
p$age+
  theme(legend.position = 0)+
  p$sex+theme(axis.text.y = element_blank(), strip.text.y = element_blank())+ylab("")
ggsave("Figures/Figure2.General_phenotype.pdf", height = 7, width =10)

p$treatment+theme(legend.position = 0)+
  p$diss_or_EM+theme(axis.text.y = element_blank(), strip.text.y = element_blank(), legend.position = 0)+ylab("")+ggtitle("Disseminated vs EM")+
  p$C6_t0+theme(axis.text.y = element_blank(), strip.text.y = element_blank(), legend.position = 0)+ylab("")+ggtitle("C6 levels at Baseline")+
  p$C6_t6+theme(axis.text.y = element_blank(), strip.text.y = element_blank())+ylab("")+ggtitle("C6 levels after treatment")+
  plot_layout(widths = c(1,2,1,1))
ggsave("Figures/Figure2.Disease_phenotype.pdf", height = 7, width =12.5)


#Figure 2.c Boxplots/dotplots for top

age_dp <- ggplot(Cyt_t_age_sex, aes(x=age, y= log2(t6_pbmc_il1ra_bbmix_10_4)))+
  geom_point(fill = "black", alpha = .5, pch = 21, size = 1)+
  theme_bw()+geom_smooth(color = "black")+theme(panel.grid = element_blank())
sex_bp <- ggplot(Cyt_t_age_sex%>%
                   filter(!is.na(sex))%>%
                   mutate(sex = ifelse(sex == 0 ,"Female","Male")),
                 aes(x=sex, y= log2(t0_pbmc_il1ra_p3c)))+
  geom_jitter(fill = "black", alpha = .5, pch = 21, size = 1)+
  geom_violin(alpha = .5)+
  geom_boxplot(alpha = .5, width = .2)+
  theme_bw()+theme(panel.grid = element_blank())

treatment_bp <- ggplot(Cyt_df.2%>%
                         filter(notime == "wba_il10_bbmix_moi30"), aes(x=Time, y=logvalue))+
  geom_jitter(fill = "black", alpha = .5, pch = 21, size = 1)+
  geom_violin(alpha = .5)+
  geom_boxplot(alpha = .5, width = .2)+
  theme_bw()+theme(panel.grid = element_blank())

diss_or_EM_bp <- ggplot(Cyt_t_age_sex%>%
                          filter(!is.na(diss_or_EM)), aes(x=diss_or_EM, y= log2(t0_pbmc_il1ra_p3c)))+
  geom_jitter(fill = "black", alpha = .5, pch = 21, size = 1)+
  geom_violin(alpha = .5)+
  geom_boxplot(alpha = .5, width = .2)+
  theme_bw()+theme(panel.grid = element_blank())

C6_t0_dp <- ggplot(Cyt_t_age_sex, aes(x=C6_t0, y= log2(t0_wba_il10_bbmix_moi10)))+
  geom_point(fill = "black", alpha = .5, pch = 21, size = 1)+
  theme_bw()+geom_smooth(color = "black")+theme(panel.grid = element_blank())
C6_t6_dp <- ggplot(Cyt_t_age_sex, aes(x=C6_t6, y= log2(t6_wba_il10_bbmix_moi10)))+
  geom_point(fill = "black", alpha = .5, pch = 21, size = 1)+
  theme_bw()+geom_smooth(color = "black")+theme(panel.grid = element_blank())

age_dp+sex_bp
ggsave("Figures/Figure2.General_phenotype_plots.pdf", height = 7, width = 10)
treatment_bp+diss_or_EM_bp+C6_t0_dp+C6_t6_dp+plot_layout(nrow = 1)
ggsave("Figures/Figure2.Disease_phenotype_plots.pdf", height = 7, width = 20)

#In disseminated there are more men than in EM, let's plot the same stratified
diss_or_EM_bp_sex <- Cyt_t_age_sex%>%
  mutate(sex = ifelse(sex == 0 ,"Female","Male"))%>%
  mutate(Cyt =as.numeric(log2(t0_pbmc_il1ra_p3c)))%>%
  select(diss_or_EM, Cyt, sex)%>%na.omit()%>%
  mutate(xvar = factor(paste(sex, diss_or_EM), levels = c('Female EM', 'Female DISS', 'Male EM', 'Male DISS')))%>%
  ggboxplot(x='xvar', y= 'Cyt', add = 'jitter', outlier.shape = NA)+
  stat_compare_means(comparisons = list(c(1,2),c(3,4),c(1,3),c(2,4)), method.args = list(exact=F))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab('log2 t0_pbmc_il1ra_p3c')
ggsave("Figures/Figure2.Diss_EM_sex_bp.pdf", height = 6.5, width = 4.5)



#Figure. Antibiotics effect, heatmap and boxplot

p$BloodbeforeABstart_atBL
ggsave("Figures/Figure2.Antibody_start.pdf", height = 7, width =6)

Antibody_start_dp <- ggplot(Cyt_t_age_sex, aes(x=BloodbeforeABstart_atBL, y= log2(t0_pbmc_il1ra_bbmix_10_5)))+
  geom_jitter(fill = "black", alpha = .5, pch = 21, size = 1)+
  geom_violin(alpha = .5)+
  geom_boxplot(alpha = .5, width = .2)+
  theme_bw()+theme(panel.grid = element_blank())
Antibody_start_dp
ggsave("Figures/Figure2.Antibiotics_start_bp.pdf", height = 7, width = 5)

#Are disseminated, sex and antibiotics indepdendent
#Men are more disseminated. 2.5% of women, 5.9% of men have disseminated
chisq.test(table(Cyt_t_age_sex$sex, Cyt_t_age_sex$diss_or_EM)) 
#Men are less likely to have started antibiotics, 50% of no-antibiotics are men, 40% of antibiotics are men
chisq.test(table(Cyt_t_age_sex$sex, Cyt_t_age_sex$BloodbeforeABstart_atBL))
#Disseminated are more likely to not have started antibiotics (?). 12% of EM have not started antibiotics, 34% of DISS have not started
chisq.test(table(Cyt_t_age_sex$diss_or_EM, Cyt_t_age_sex$BloodbeforeABstart_atBL))

#In Non-Antibiotics start there are more men than and DISS, let's plot the same stratified
Antibiotics_start_bp_sex <- Cyt_t_age_sex%>%
  mutate(sex = ifelse(sex == 0 ,"Female","Male"))%>%
  mutate(Cyt =as.numeric(log2(t0_pbmc_il1ra_p3c)))%>%
  select(BloodbeforeABstart_atBL, Cyt, sex,diss_or_EM)%>%na.omit()%>%
  mutate(xvar = factor(paste(sex, BloodbeforeABstart_atBL, sep = '\n'), levels = c('Female\nAntibiotics', 'Female\nNot_started_antibiotics', 'Male\nAntibiotics', 'Male\nNot_started_antibiotics')))%>%
  ggboxplot(x='xvar', y= 'Cyt', add = 'jitter', outlier.shape = NA)+
  stat_compare_means(comparisons = list(c(1,2),c(3,4),c(1,3),c(2,4)), method.args = list(exact=F))+
  theme_bw()+theme(panel.grid = element_blank())+
  facet_wrap(~diss_or_EM)+
  ylab('log2 t0_pbmc_il1ra_p3c')
Antibiotics_start_bp_sex
ggsave("Figures/Figure2.Antibiotic_start_sex_diss_EM_bp.pdf", height = 8, width = 10)


#Figure 3. Mannhattan plot and Heatmap

#Read results
lognorm <- fread("../Output/Cytokines_goodqual_5e-8.txt")

colnames(lognorm) <- c("SNP","gene","beta","t","P")

lognorm <- remove_badqual(lognorm, column = "gene")

lognorm <- separate(lognorm, SNP, into = c("chrpos","SNP"),"%")%>%
  unite("SNP_gene",c(SNP,gene), sep = "_", remove = F)%>%
  mutate("notime" = gsub("t\\d_","",gene))%>%
  separate(gene, into = c("Time","Sample","Cyt","Stim"), sep ="_", remove = F)%>%
  separate(chrpos, into = c("chr","BP"), ":", remove = F)

#Read dataframe for Manhattan
manhattan <- fread("../Output/Cytokines_goodqual_manhattan_loci.txt")

#Merge the chr:pos:ref:alt, to allow for identifying "." SNPs
manhattan <- unite(manhattan,"SNP", c("CHR","POS","REF","ALT", "SNP"), sep = ":", remove = F)

#Read independent loci
loc <- fread("../Output/FUMA_cQTL/Indloci.annotated.txt")

manhattan <- manhattan %>%
  left_join(loc%>%
              select(GenomicLocus, Cytokine)%>%rename("loci" = "GenomicLocus"))

source("~/bin/ggmanhattan.R")

ggmanhattan(manhattan, bp = "POS", significance = c(5e-8,1.55e-9), scale_color = scale_color_manual(values = rep(c("lightgrey","darkgrey"), 11)),
            highlight = list(manhattan %>%filter(Cytokine == "il6")%>% pull(SNP),
                             manhattan %>%filter(Cytokine == "il10")%>% pull(SNP),
                             manhattan %>%filter(Cytokine == "il1ra")%>% pull(SNP),
                             manhattan %>%filter(Cytokine == "il1b")%>% pull(SNP),
                             manhattan %>%filter(grepl(";",Cytokine))%>% pull(SNP)),
            highlight_col = list("#f5d3ce","#b04632","#ec9488","#eb5a46","darkred"),
            lead_snp = manhattan %>%
              filter(grepl(paste(loc$SNP, collapse = "|"), SNP) & loci != 0) %>% pull(SNP),
            annotate_snp = T
)
ggsave("Figures/Figure3.manhattan.png", width = 12, height = 10)
ggsave("Figures/Figure3.manhattan.pdf", width = 12, height = 10)



#Heatmap
library(ComplexHeatmap)
library(circlize)

ind_snp.1 <- fread("../Output/FUMA_cQTL/GenomicRiskLoci.txt")

ind_snp <- paste(ind_snp.1$chr, ind_snp.1$rsID, sep = "_")

#Read all the pvals per snp
ind_lognorm <- fread("../Output/FUMA_cQTL/allsnps.pval.tsv")
#Read GWAS t align directions

#Option 1: My GWAS
#ind_lognorm_gwas <- fread(cmd = paste0("egrep 'CHR|",paste(unique(ind_lognorm$rsID), collapse = "|"), "' ../GWAS/LymevsHealthy/LB1_LB2_mymetal_stderr.txt"))
#Risk allele per snp
#ind_lognorm_gwas$GWASRisk <- ifelse(ind_lognorm_gwas$z.meta > 0, ind_lognorm_gwas$A11, ind_lognorm_gwas$A21)

#Option 2: Finngen spirochaetal
ind_lognorm_gwas <- fread(cmd = paste0("zgrep -w -e chrom ",paste(paste("-e", unique(ind_lognorm$rsID)), collapse = " "), " ../Output/GWASCatalogue/finngen_R7_AB1_OTHER_SPIROCHAETAL.gz"))
#Risk allele per snp
ind_lognorm_gwas$GWASRisk <- ifelse(ind_lognorm_gwas$beta > 0, ind_lognorm_gwas$alt, ind_lognorm_gwas$ref)

#Bind both
ind_lognorm<- left_join(ind_lognorm,
                        rename(ind_lognorm_gwas, "rsID" = rsids) %>% select(rsID,GWASRisk))%>%
  separate(SNP, into = c("chr","pos","ref","alt","rsID"), remove = F)%>%
  mutate(beta = ifelse(is.na(GWASRisk) ,beta, ifelse(ref == GWASRisk, -beta, beta )))


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

ca = HeatmapAnnotation(CHR = anno_text(str_split(colnames(ind_lognorm[,colnames(ind_lognorm) %in% ind_snp]), "_", simplify = T)[,1],
                                       show_name = T, rot = 0, gp = gpar(fontsize = 8)))
locia = HeatmapAnnotation(Loci = anno_text(ind_snp.1$GenomicLocus, show_name = T, rot =0, gp = gpar(fontsize = 8)))

pdf("Figures/Figure3.Heatmap_indloci.pdf", width = 12, height = 10)
ComplexHeatmap::Heatmap(ind_lognorm[,colnames(ind_lognorm) %in% ind_snp], cluster_columns = F, show_row_dend = F, left_annotation = ra,
                        show_row_names = F, cluster_rows = F,
                        col = colorRamp2(c(log10(5e-8),log10(0.05), log10(0.05)+0.0001, 0,  -log10(0.05)-0.0001 ,-log10(0.05), -log10(5e-8)),
                                         c("blue","lightblue","white", "white","white","orange", "red")),
                        bottom_annotation = ca,show_column_names = F, top_annotation = locia,
                        name = "Signed\nLog10 pvalue")
dev.off()




#Table 1.

loci <- fread("../Output/FUMA_cQTL/Indloci.annotated.txt")

library(kableExtra)
options(knitr.kable.NA = "")

table1.loci <- loci %>%
  mutate(exp_gene = gsub("NA;","", exp_gene)%>% gsub(";NA","",.)) %>%#Remove the NAs in exp_gene
  rowwise()%>%
  mutate("DE" = paste(DE_stimulation_4, DE_stimulation_24, sep =";") %>% str_split(";") %>% unlist() %>% unique()%>%paste(collapse = ";")) %>%
  select(GenomicLocus, nGWASSNPs, SNP, chr, pos, p, Cytokine, Stimulation, Time, Sample, exp_gene, DE)%>%
  rename("Locus Nr." = GenomicLocus,
         "nSNPS" = nGWASSNPs,
         "Cell System" = Sample,
         "eQTL" = exp_gene)%>%
  mutate("Cytokine" = gsub("il","IL",Cytokine), "Time" = gsub("t0","Baseline",Time) %>% gsub("t6","After treatment",.))%>%
  mutate_all(function(x) gsub(";", ", ", x))%>%
  mutate(p = formatC(as.numeric(p), format = "e", digits = 2))

#One big table
table1.loci %>% kable(align = "c", escape = T)%>%
  kable_styling()%>%
  column_spec(6, width = "10em")%>%
  add_header_above(c(" " = 2,"Lead SNP" = 4, " " = 6))%>%
  save_kable("Figures/Table1.Indloci.table.pdf")

#One table for SW loci
table1.loci %>% filter(as.numeric(p) < 1.55e-9)%>%
  kable(align = "c", escape = T)%>%
  kable_styling()%>%
  column_spec(6, width = "10em")%>%
  add_header_above(c(" " = 2,"Lead SNP" = 4, " " = 6))%>%
  save_kable("Figures/Table1.Indloci.SW.pdf")

#One table for GW loci
table1.loci %>% filter(as.numeric(p) > 1.55e-9)%>%
  kable(align = "c", escape = T)%>%
  kable_styling()%>%
  column_spec(6, width = "10em")%>%
  add_header_above(c(" " = 2,"Lead SNP" = 4, " " = 6))%>%
  save_kable("Figures/Table1.Indloci.GW.pdf")

#One table per cytokine
#table1.loci %>%
#  group_by(Cytokine)%>%
#  group_walk(~{.x%>%
#      kable(align = "c", escape = T)%>%
#      kable_styling()%>%
#      column_spec(6, width = "10em")%>%
#      add_header_above(c(" " = 2,"Lead SNP" = 4, " " = 5))%>%
#      save_kable(paste0("Figures/Table1/",unique(.y$Cytokine),".Indloci.table.pdf"))
#  })


#Figure 4. cQTL examples
dosages <- fread("../Output/FUMA_cQTL/allsnps.dosage.tsv")
cytokines <- as.data.frame(fread("../Phenotype/Cytokines_goodqual.tsv"))
#Figure 4.A cis IL6 in Borrelia mix moi30 whole blood

boxplot_theme <- function(){
	theme_bw()+theme(panel.grid = element_blank(), strip.background = element_blank(),axis.line = element_line(colour = "black"),
                                                          panel.border = element_blank(),panel.background = element_blank())
	}

cQTL_boxplot <- function(rsID, Cyt, write = T){
  #Read genotype and cytokine
  cQTL <- dosages[grepl(paste0(rsID,"$"), dosages$V1),]
  Cytokine_measured<- cytokines[grepl(Cyt, cytokines$V1),]
  if(nrow(cQTL) > 1) stop("Too many snps")
  if(nrow(Cytokine_measured) > 1) stop("Too many snps")
  
  #Match column names
  Cytokine_measured <- Cytokine_measured[,colnames(cQTL)]
  
  #Bind
  df <- rbind(cQTL, Cytokine_measured)%>%
    remove_rownames()%>%
    column_to_rownames("V1")%>%
    t()%>%
    as.data.frame()
  
  split = str_split(colnames(df)[1], "[^[:alnum:]]+", simplify = T)
  ref = split[,3]
  alt = split[,4]
  snp = split[,5]
  
  df$snp <- ifelse(df[,1] <  0.5, paste0(ref,ref), ifelse(df[,1] < 1.5, paste0(ref,alt), paste0(alt,alt)))
  df$snp <- factor(df$snp, levels = c(paste0(ref,ref),paste0(ref,alt),paste0(alt,alt)))
  colnames(df)[2] <-"Cytokine"
  df$Cytokine <- log2(df$Cytokine)
  require(ggpubr)
 p <- ggboxplot(df, x = "snp",y = "Cytokine", fill  = ggsci::pal_npg()(3)[1], 
                add = c("point"), add.params = list(position = position_jitter(w = 0.05)))+
   ylab(paste("log2",Cyt))+xlab(rsID)+
   boxplot_theme()+
    stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)), label = 'p.signif')
  if(write) ggsave(paste0("Figures/Figure5.Boxplot",Cyt,rsID,".pdf"), plot = p, width = 6, height = 6) else p
}
cQTL_boxplot("rs35345753", "t0_wba_il6_bbmix_moi30")
cQTL_boxplot("rs1518111", "t6_wba_il10_bbmix_moi30")
cQTL_boxplot("rs488380", "t0_pbmc_il1ra_cand")
cQTL_boxplot("rs5743618", "t0_pbmc_il1b_bbmix_10_5")

#the snp that is not replicated in borrelia afzelii
#Selecting only patients that are measured in borrelia afzelii
keep = colnames(cytokines)[!is.na(cytokines[cytokines$V1 == "t6_pbmc_il10_bafzelii",])]
dosages = as.data.frame(dosages)[,colnames(dosages) %in% keep]
cQTL_boxplot("rs72653590", "t6_pbmc_il10_bafzelii")
cQTL_boxplot("rs72653590", "t0_pbmc_il10_bafzelii")
cQTL_boxplot("rs72653590", "t6_pbmc_il10_bbmix_10_5")
cQTL_boxplot("rs72653590", "t0_pbmc_il10_bbmix_10_5")



#Figures to compare il1ra vs p3c in locus 8, also for il6
cQTL_boxplot("rs5743618", "t0_pbmc_il1ra_p3c")
cQTL_boxplot("rs5743618", "t0_pbmc_il1ra_bbmix_10_5")
cQTL_boxplot("rs5743618", "t0_pbmc_il6_p3c")
cQTL_boxplot("rs5743618", "t0_pbmc_il6_bbmix_10_5")

#Also time 6
cQTL_boxplot("rs5743618", "t6_pbmc_il1ra_p3c")
cQTL_boxplot("rs5743618", "t6_pbmc_il1ra_bbmix_10_5")
cQTL_boxplot("rs5743618", "t6_pbmc_il6_p3c")
cQTL_boxplot("rs5743618", "t6_pbmc_il6_bbmix_10_5")

cQTL_boxplot("rs5743618", "t0_pbmc_il1ra_bbmix_10_6")
cQTL_boxplot("rs5743618", "t0_pbmc_il1ra_p3c")
