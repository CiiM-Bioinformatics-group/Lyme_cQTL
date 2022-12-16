PCA_hm <- function(df, df2, colname){
  ##df = df on which PCA is performed
  ##df2 = df containing data for plotting and heatmap. Rownames of PCA are stored in column "colname" of this df
  PC <- prcomp(df, scale. = T, center = T)
  PC_df <- PC %>% .$x %>% as.data.frame() %>% rownames_to_column(var =colname)

  PC_df <- PC_df %>% select(1:21) %>%
    left_join(df2, colname)


  comp <- colnames(PC_df)[grepl("PC", colnames(PC_df))]
  vars <- colnames(PC_df)[!grepl("PC", colnames(PC_df))]

  pmat <- matrix(nrow = length(vars), ncol = length(comp),
                 dimnames = list(vars, comp))

  for (pc in comp){
    for (var in vars){
      cat(pc,"\n", var, "\n")
      if(is.numeric(PC_df[,var])){
        res <- cor.test(PC_df[,pc], PC_df[,var], use = "pairwise.complete.obs", method = "spearman")
      }else if(length(unique(na.omit(PC_df[,var]))) == 2){
        df <- na.omit(PC_df[,c(var,pc)])
        colnames(df) <- c("x","y")
        res <- wilcox.test(y ~ x, df)
      }else{
        df <- na.omit(PC_df[,c(var,pc)])
        colnames(df) <- c("x","y")
        res <- kruskal.test(y ~ x, df)
      }
      pmat[var,pc] = -log10(res$p.value)
    }
  }
  return(pmat)
}

PCA_wrap <- function(df){
  ##df = df on which PCA is performed
  PC <- prcomp(df, scale. = T, center = T)
  PC_df <- PC %>% .$x %>% as.data.frame() %>% rownames_to_column(var ="Individual")
  return(PC_df)
}

cyt_split <- function(df, col){
  df <- mutate_at(df, col, function(col) gsub("_10_","-10-", col) %>% gsub("_moi","-moi", .))
  df <-tidyr::separate(df, col, into = c("Time","Sample","Cytokine","Stimulation"), sep = "_", remove = F)
  df <- tidyr::unite(df, "Cyt_stim", c("Cytokine","Stimulation"), sep = "_", remove = F)
  df <- tidyr::unite(df, "notime", c("Sample","Cytokine","Stimulation"), sep = "_", remove = F)
  return(df)
}

remove_badqual <- function(df, column = F,  dataset = "notnorm") {
  ##Column can be, either the column name where the cytokine values are stored
  ## or F if cytokines are in rows
  df <- as.data.frame(df)
  if(dataset == "badqual") badqual <- fread("badqual.txt")%>% as.data.frame() %>%
      .[,1] #read badqual file
  if(dataset == "notnorm") badqual <- fread("notnorm.txt") %>% as.data.frame() %>%
      .[,1]#read badqual file
  badqual <- c(paste0("t0_", badqual), paste0("t6_", badqual))
  #Remove questionnaires
  if(column == FALSE){
    gene <- row.names(df)
  }else{
    gene <- df[,column]
  }
  df <- subset(df, grepl("pbmc|wba", gene) & !gene %in% badqual  )
  return(df)

}

ind_N <- function(df, column = F, value = "value", id = "patient", thr = 0.5){
  df. <- df
  if(column != F){
    df <- df[,c(id, column, value)]
    df <- pivot_wider(df, names_from = all_of(column), values_from = value)%>%
      column_to_rownames(var = id)
  }
  mat <- cor(df, method = "spearman", use = "pairwise.complete.obs")
  mat[is.na(mat)] <- 0
  cl <- hclust(as.dist(1- abs(mat)))
  gr <- cutree(cl, h =1- thr)
  df.$N <- length(unique(gr))
  return(df.)
}

INT <- function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
hardINT <- function(x) qnorm((rank(x,na.last="keep",ties.method = "random")-0.5)/sum(!is.na(x)))


metal <- function(df1, df2, method = "Inverse-variance",
                  P = "P", E = "E", N = "N", A1 = "A1", A2 = "A2", SE = "SE",
                  CHR= "CHR", BP = "BP", SNP = "SNP"){

  #Method can be one of: "Sample-size" or "Inverse-variance"
  ##Weighted z-score meta-analysis based on the approach taken by Metal meta-analysis.
  # Reference: http://csg.sph.umich.edu/abecasis/publications/pdf/Bioinformatics.vol.26-pp.2190.pdf

  df = dplyr::inner_join(df1, df2, by = c("CHR", "BP"), suffix = c("1","2"))
  #If reference and alternative are flipped, change direction
  df$E2 <- ifelse(df$A11 == df$A22 & df$A12 == df$A21, -df$E2,df$E2)

  if(method == "Sample-size"){
    P1 = df[,P1] #P-value from the first dataset
    P2 = df[,P2] #P-value from the second dataset

    E1 = sign(df[,E1]) #Effect direction from the first dataset
    E2 = sign(df[,E2]) #Effect direction from the second dataset

    W1 = sqrt(df[,N1])/sqrt(df[,N1] + df[,N2]) #Weight for the first dataset
    W2 = sqrt(df[,N2])/sqrt(df[,N1] + df[,N2]) #Weight for the second dataset

    df$z.meta = qnorm(P1/2, lower.tail = F)*E1*W1+qnorm(P2/2, lower.tail = F)*E2*W2 #Weighted Z-score
  }else if(method == "Inverse-variance"){
    SE1 = df[,SE1]
    SE2 = df[,SE2]

    W1 = 1/(SE1^2)
    W2 = 1/(SE2^2)

    SE = sqrt(1/(W1+W2))

    E1 = df[,E1]
    E2 = df[,E2]

    E = (E1*W1+E2*W2)/(W1+W2)

    df$z.meta = E/SE
  }

  df$p.meta = 2*pnorm(abs(df$z.meta), lower.tail = F) #P-value from the weighted Z-score


  return(df)
}


minmax <- function(x){
  if(is.na(x)){
    return(NA)
  } else {
    return((x - min(x, na.rm = T)) /(max(x, na.rm = T)-min(x, na.rm = T)))
  }}
metal_gwas <- function(df, P1 = "P.lb1", P2 = "P.lb2", E1 = "OR.lb1", E2 = "OR.lb2", N1 = "n1", N2 = "n2"){

  ##Weighted z-score meta-analysis based on the approach taken by Metal meta-analysis.
  # Reference: http://csg.sph.umich.edu/abecasis/publications/pdf/Bioinformatics.vol.26-pp.2190.pdf

  df = as.data.frame(df) #Transform the dataframe in case it is a data.table object

  P1 = df[,P1] #P-value from the first dataset
  P2 = df[,P2] #P-value from the second dataset

  E1 = sign(log(df[,E1])) #Effect direction from the first dataset
  E2 = sign(log(df[,E2])) #Effect direction from the second dataset

  W1 = sqrt(df[,N1])/sqrt(df[,N1] + df[,N2]) #Weight for the first dataset
  W2 = sqrt(df[,N2])/sqrt(df[,N1] + df[,N2]) #Weight for the second dataset

  df$z.meta = qnorm(P1/2, lower.tail = F)*E1*W1+qnorm(P2/2, lower.tail = F)*E2*W2 #Weighted Z-score

  df$p.meta = 2*pnorm(abs(df$z.meta), lower.tail = F) #P-value from the weighted Z-score


  return(df)
}

proper_QTL <- function(genotype_fl, prefix, covariates_fl =  "../Covariates/covariates", p.val = 0.05, batch2 = "alone"){
  require(MatrixEQTL)

  useModel <- modelLINEAR

  if(genotype_fl == "LB1") SNP_file_name = "../Sample/lyme-batch1_500FG_imputed.dosage.tsv"
  if(genotype_fl == "LB2") {
    if(batch2 == "alone"){
      SNP_file_name = "../Sample/lyme-batch2_imputed_alone.dosage.tsv"
    }
    if(batch2 == "300BCG"){
      SNP_file_name = "../Sample/lyme-batch2_300BCG_imputed.dosage.tsv"
    }}


  covariates_fl = paste0(covariates_fl, "_", genotype_fl, ".tsv")
  pheno_fl <- paste0("../Phenotype/",prefix, "_", genotype_fl,".tsv")
  output_fl <- paste0("../Output/",prefix, "_", genotype_fl,".txt")

  if(!file.exists(SNP_file_name)) stop(paste("Geno file",SNP_file_name, "does not exist"))
  if(!file.exists(pheno_fl)) stop(paste("Pheno file",pheno_fl, "does not exist"))

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
    pvalue.hist = "qqplot",
    noFDRsaveMemory = F
  )
  saveRDS(quest_qtl, paste0("../Output/",prefix, "_", genotype_fl, ".RDS"))
}


meta_analysis <- function(file, pval = NULL){
  if(!is.null(pval)) pval = paste0("_",pval)
  #Read pheno file to count N in each dataset
  N1 = fread(paste0("../Phenotype/",file, "_LB1.tsv")) %>% column_to_rownames(var = colnames(.)[1])
  N2 = fread(paste0("../Phenotype/",file, "_LB2.tsv")) %>% column_to_rownames(var = colnames(.)[1])

  #Read the Output file
  df1.1 = fread(paste0("../Output/",file,pval, "_LB1.txt"))
  df2.1 = fread(paste0("../Output/",file,pval, "_LB2.txt"))

  #Combine output files
  df = inner_join(df1.1, df2.1, by = c("SNP", "gene"), suffix = c(".lb1",".lb2"))

  #Read N for each trait
  n1_gene <- sapply(unique(df$gene), function(x) length(na.omit(t(N1[x,]))))
  n2_gene <- sapply(unique(df$gene), function(x) length(na.omit(t(N2[x,]))))

  #Combine and meta-analyse each snp trait pair
  df <- df %>% mutate("n1" = n1_gene[gene], "n2" = n2_gene[gene]) %>%
    metal()

  return(df)
}

#Colors for heatmap of cytokines and stimulations
anno_col = list(Time = c("t0" = "lightgrey", "t6" = "grey"),
                `Cell system` = c("pbmc" = "#a0711c","wba" = "#933b27",
                                  "PBMC" = "#a0711c","WBA" = "#933b27"),
                Cytokine = c("il1ra" = "#ec9488","il1b" = "#eb5a46","il10" = "#b04632","il6" = "#f5d3ce",
                             "IL1ra" = "#ec9488","IL1b" = "#eb5a46","IL10" = "#b04632","IL6" = "#f5d3ce"),
                Stimulation = c("p3c" = "#fce8d2" , "lps" = "#fdc788", "cand" = "#ffab4a", "bbmix" = "#d29034",
                                "bafzelii" = "brown", "rpmi" = "white"))

##########################
#Code copied from https://github.com/cran/GenABEL/blob/master/R/rntransform.R

"ztransform" <- 
  function(formula,data,family=gaussian) {
    if (missing(data)) {
      if(is(formula,"formula")) 
        data <- environment(formula)
      else  
        data <- environment()
      #		wasdata <- 0
    } else {
      if (is(data,"gwaa.data")) {
        data <- data@phdata
      } 
      else if (!is(data,"data.frame")) {
        stop("data argument should be of gwaa.data or data.frame class")
      }
      #		attach(data,pos=2,warn.conflicts=FALSE)
      #		wasdata <- 1
    }
    
    if (is.character(family)) 
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
    
    if ( is(try(formula,silent=TRUE),"try-error") ) { 
      formula <- data[[as(match.call()[["formula"]],"character")]] 
    }
    
    if (is(formula,"formula")) {
      #		mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
      mf <- model.frame(formula,data,na.action=na.pass,drop.unused.levels=TRUE)
      mids <- complete.cases(mf)
      mf <- mf[mids,]
      y <- model.response(mf)
      desmat <- model.matrix(formula,mf)
      lmf <- glm.fit(desmat,y,family=family)
      #		if (wasdata) 
      #			mids <- rownames(data) %in% rownames(mf)
      #		else 
      resid <- lmf$resid
      #		print(formula)
    } else if (is(formula,"numeric") || is(formula,"integer") || is(formula,"double")) {
      y <- formula
      mids <- (!is.na(y))
      y <- y[mids]
      resid <- y
      if (length(unique(resid))==1) stop("trait is monomorphic")
      if (length(unique(resid))==2) stop("trait is binary")
    } else {
      stop("formula argument must be a formula or one of (numeric, integer, double)")
    }
    y <- (resid-mean(resid))/sd(resid)
    #	if (wasdata==1) detach(data)
    tmeas <- as.logical(mids)
    out <- rep(NA,length(mids))
    out[tmeas] <- y
    out
  }

"rntransform" <-
  function(formula,data,family=gaussian) {
    
    if ( is(try(formula,silent=TRUE),"try-error") ) { 
      if ( is(data,"gwaa.data") ) data1 <- phdata(data)
      else if ( is(data,"data.frame") ) data1 <- data
      else stop("'data' must have 'gwaa.data' or 'data.frame' class")
      formula <- data1[[as(match.call()[["formula"]],"character")]] 
    }
    
    var <- ztransform(formula,data,family)
    out <- rank(var) - 0.5
    out[is.na(var)] <- NA
    mP <- .5/max(out,na.rm=T)
    out <- out/(max(out,na.rm=T)+.5)
    out <- qnorm(out)
    out
  }




