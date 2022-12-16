#!/usr/bin/env Rscript
if(!require("plinkQC")){
  install.packages("plinkQC",repos = "http://cran.us.r-project.org")
}
library("plinkQC")
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

##Change sample name
name <- args[1]
snpQC <- as.logical(args[2])

###
package.dir <- paste0(find.package('plinkQC'),"/","extdata")
indir <- 'data'
qcdir <- 'QC'
path2plink <- "plink"

fail_sample <- perIndividualQC(indir = indir, qcdir = qcdir, name=name,
                               refSamplesFile=paste(package.dir, "/Genomes1000_ID2Pop.txt",
                                                    sep=""),
                               refColorsFile=paste(package.dir, "/Genomes1000_PopColors.txt",
                                                   sep=""),
                               prefixMergedDataset=paste0(name,"_1000g"),
                               path2plink=path2plink, do.run_check_ancestry = FALSE,
                               interactive=TRUE, verbose=TRUE,
                               europeanTh = 6, imissTh = 0.03,highIBDTh = 0.185,
                               hetTh = 3)


ggsave(paste0("Plots/",name,".sampleQC.pdf"),fail_sample$p_sampleQC, width = 10, height = 12)
overview_individuals <- overviewPerIndividualQC(fail_sample)
ggsave(paste0("Plots/",name,".overview_samples.pdf"),overview_individuals$p_overview)

##If snpQc is required also, when analyzing the Healthy cohort, samples and snps
# are cleaned after snpQC. If it is not required, samples are cleaned before

if(!snpQC){

##In order to perform quality control only on samples
samplesID <- cleanData(indir=indir, qcdir=qcdir, name=name, path2plink=path2plink,
                       verbose=TRUE, showPlinkOutput=T, filterSNPMissingness = F,
                       filterHWE = F, filterMAF = F)
}

if(snpQC){

fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name,
                            path2plink=path2plink,
                            verbose=TRUE, interactive=TRUE,
                            showPlinkOutput=FALSE,
                            mafTh = 0.01, hweTh = 1e-4)

ggsave(paste0("Plots/",name,".snpQC.pdf"),fail_markers$p_markerQC, width = 10, height = 12)

overviewPerMarkerQC <- function (results_perMarkerQC, interactive = FALSE)
{
  if (length(perMarkerQC) == 2 && !all(names(results_perMarkerQC) ==
                                       c("fail_list", "p_markerQC"))) {
    stop("results_perMarkerQC not direct output of perMarkerQC")
  }
  fail_list <- results_perMarkerQC$fail_list
  unique_markers_fail <- unique(unlist(fail_list))
  fail_counts <- UpSetR::fromList(fail_list)
  p <- UpSetR::upset(fail_counts, order.by = "freq", empty.intersections = "on",
                     text.scale = 1.2, mainbar.y.label = "Markers failing multiple QC checks",
                     sets.x.label = "Marker fails per QC check", main.bar.color = "#1b9e77",
                     matrix.color = "#1b9e77", sets.bar.color = "#d95f02")
  p_overview <- cowplot::plot_grid(NULL, p$Main_bar, p$Sizes,
                                   p$Matrix, nrow = 2, align = "v", rel_heights = c(3, 1),
                                   rel_widths = c(2, 3))
  if (interactive) {
    print(p_overview)
  }
  nr_fail_markers <- length(unique_markers_fail)
  return(list(nr_fail_samples = nr_fail_markers, fail_QC = fail_counts,
              p_overview = p_overview))
}

overview_marker <- overviewPerMarkerQC(fail_markers)
overview_marker$p_overview

ggsave(paste0("Plots/",name,".overview_snps.pdf"),overview_marker$p_overview)

exclude_ancestry <- check_ancestry(indir=indir, qcdir=qcdir, name=name,refSamplesFile=paste(package.dir, "/Genomes1000_ID2Pop.txt",
                                                                                            sep=""),
                                   refColorsFile=paste(package.dir, "/Genomes1000_PopColors.txt",
                                                       sep=""),
                                   prefixMergedDataset=paste0(name,"_1000g"),
                                   path2plink=path2plink, run.check_ancestry = FALSE,
                                   interactive=TRUE,
                                   europeanTh = 6)
ggsave(paste0("Plots/",name,".PCA.ancestry.pdf"), exclude_ancestry$p_ancestry, width = 15, height = 12)



Ids  <- cleanData(indir=indir, qcdir=qcdir, name=name, path2plink=path2plink,
                  verbose=TRUE, showPlinkOutput=FALSE, mafTh = 0.01, hweTh = 1e-4)

}
