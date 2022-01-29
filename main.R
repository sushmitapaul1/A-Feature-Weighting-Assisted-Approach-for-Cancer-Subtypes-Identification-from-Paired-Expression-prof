## Calculating feature weights based on proposed approach.
## "robust_weight.m" contains the MATLAB codes to generate weights based on robust regression.
## The MATLAB function is Called from R for easy working. 
## First feature weights are calculated for each of the expression profiles. 
## All the expression values should be greater than or equal to 0. Features containing
## more than 60% of zero values should be removed.

matlab.lines <- c(
  "x = csvread('mirna.csv')",
  "weight = robust_weight(x)",
  "csvwrite('weight_mirna.csv', weight)",
  "y = csvread('mrna.csv')",
  "weight_1 = robust_weight(y)",
  "csvwrite('weight_mrna.csv', weight_1)")
writeLines(matlab.lines, con ="robust_weights.m")
system("matlab -nodisplay -r \"run('robust_weights.m'); exit\"")

weight_mirna <- read.csv('weight_mirna.csv',header = F)
weight_mirna <- weight_mirna$V1

weight_mrna <- read.csv('weight_mrna.csv',header = F)
weight_mrna <- weight_mrna$V1

###########################################################################################
## The obtained normalized weights are further used for cancer subtypes identification as:
library(SNFtool)

HIST_SUB <- function(X1,X2, k=number_of_cluster)
  ##X1 and X2 are miRNA and mRNA expression matrices (Row is feature, column is sample)
{
  cat("\n Weighted Sample Similarity Network Generation...\n")
  Dist1 <- distanceWeighted(t(X1), weight_mirna)
  W1 <- affinityMatrix(Dist1, 20, 0.8)
  Dist2 <- distanceWeighted(t(X2), weight_mrna)
  W2 <- affinityMatrix(Dist2, 20, 0.8)
  cat("\n Sample Similarity Network Fusion and Sample Clustering... \n") 
  W <- SNF(list(W1,W2), 20, 20)
  group <- spectralClustering(W,k)
  result <- list(W=W,group=group)

}

distanceWeighted <- function(X,weight)
  ##X is the expression Matrix(Row is feature, column is sample) 
  ## weight is the normalized feature weight obtained by robust regression
{
  X_row <- nrow(X)      
  weight_diag <- diag(weight)
  X2 <- (X^2)%*%weight_diag 
  sumsqX <- rowSums(X2)
  X1 <- X%*%weight_diag
  XY <- X1 %*% t(X)
  XX <- matrix(rep(sumsqX, times = X_row), X_row, X_row)
  res <- XX+t(XX)-2*XY
  res[res < 0] <- 0
  diag(res) <- 0
  res <- sqrt(res)
  return(res)
}

######## Execution ##########################################################
miRNA_exp <- as.matrix(read.csv("mirna.csv", sep=",",header=F))
mRNA_exp <- as.matrix(read.csv("mrna.csv", sep=",",header=F))

miRNA_Names <- read.csv(file="miRNA_names.csv",check.names = F,header = T)
mRNA_Names <- read.csv(file="mRNA_names.csv",check.names = F,header = T)
samples <- read.csv(file="samples.csv",check.names = F,header = T)

colnames(miRNA_exp) = colnames(mRNA_exp) <- samples$x
rownames(miRNA_exp) <- miRNA_Names$x
rownames(mRNA_exp) <- mRNA_Names$x

clusters <- HIST_SUB(miRNA_exp,mRNA_exp,2)
cat("\n Clusters: \n",clusters$group)
write.csv(cbind(samples$x,clusters$group), file="cancer_subtypes.csv")

