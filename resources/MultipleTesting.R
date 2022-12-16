MeffLi <- function(mat){
  #mat: a matrix having the variables (traits, cytokines, metabolites...) as columns
  #Formula 5 taken from J Li and L Ji Heredity 2005
  
  #Calculate correlation matrix
  corr.matrix <- cor(mat, method = "spearman")
  #Eigenvalues for correlation matrix
  evals<-eigen(t(corr.matrix),symmetric=T)$values
  #Formula from the paper to calculate the number of effective tests (Meff)
  f <- function(x){
    ifelse(x < 0, 0, ifelse(x >= 1, 1, 0)) + (x - floor(x))
  }
  Meff <- sum(f(evals))
  #Alpha value calculation
  alpha <- 1 - (1-0.05)^(1/(Meff))
  alpha.gw <- alpha/1e6
  
  return(list("Meff" = Meff, "alpha" = alpha, "alpha.gw" = alpha.gw))
}
