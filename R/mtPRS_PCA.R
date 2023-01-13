#' Construct mtPRS by performing PCA on genetic correlation matrix among traits
#'
#' Combine multiple individual single-trait PRSs (stPRSs) with weights calculated by performing principal component analysis (PCA) on genetic correlation matrix among traits
#' @param dat a list of inputs including disease GWAS summary statistics of K traits from base cohort, individual-level data from base cohort, genetic correlation matrix, and true effect sizes
#' @param pcut p-value cutoff for C+T method to construct individual stPRS
#' @param varcut choose top principal components (PCs) until explaining varcut-percent variance
#' @param K number of traits
#' @param phenotype the type of phenotype in target cohort, either "dis" or "pgx"
#' @details mtPRS_PCA needs disease GWAS summary statistics from K traits, and genetic correlation matrix among K traits
#' @return a list of weights, individual stPRSs, and mtPRS-PCA
#' @references Zhai, S., Guo, B., Wu, B., Mehrotra, D.V., and Shen, J., 2023. Integrating multiple traits for improving polygenic risk prediction in disease and pharmacogenomics GWAS (submitted).
#' @author Song Zhai
#' @export
#' @examples
#' \donttest{
#' dat <- generate_pgx_data(structure = "clustered", sparseness = "more",
#' rho_DT = c(0.5,0.5,0.5,0.5), rho_T = 0.5, rho_E = 0.3, rho_C = 0.1,
#' K=4, m=1000, pcausal=0.1, blocksize=100,
#' gamma=1, samplesize=700, h2_base=0.3, h2_target=0.3)
#' re <- mtPRS_PCA(dat, pcut = 0.05, varcut = 0.8, K = 4, phenotype = "pgx")
#' hist(re$mtPRS)
#' }
#'
mtPRS_PCA <- function(dat, pcut=0.05, varcut=0.8, K=4, phenotype="pgx"){
  base <- dat$base # base is a list of disease GWAS summary statistics of K traits
  target <- dat$target # target is a data frame
  if(phenotype == "pgx"){G.target <- target[,-c(1:2)] %>% as.matrix()}
  if(phenotype == "dis"){G.target <- target[,-1] %>% as.matrix()}
  corr <- dat$corr # corr is the genetic correlation matrix among K traits, either calculated from LD score regression or obtained from public database
  X <- calculate_PRS(K, base, pcut, G.target)
  w <- calculate_weight(K, corr, varcut)
  PRS0 <- X %*% w %>% as.vector()
  re <- list(w=w, stPRS=X, mtPRS=PRS0)
  re
}

calculate_PRS <- function(K, base, pcut, G.target){
  X <- double()
  for(k in 1:K){
    ss <- base[[k]]
    index <- which(ss$p > pcut)
    beta <- ss$beta; beta[index] <- 0
    prs <- G.target %*% beta %>% as.vector()
    X <- cbind(X, prs)
  }
  for (k in 1:K) {
    if(sum(X[,k] != 0) > 0){X[,k] <- scale(X[,k])}
  }
  colnames(X) <- names(base)
  return(X)
}

calculate_weight <- function(K, corr, varcut){
  D <- nearPD(corr[1:K,1:K], corr = TRUE); mat1 <- D$mat %>% as.matrix()
  x = princomp(covmat=mat1, cor=TRUE)
  pc <- x$sdev^2/sum(x$sdev^2)
  pcc <- 0; ct <- 0
  while(pcc < varcut){
    ct <- ct + 1
    pcc <- pcc + pc[ct]
  }
  e1 <- eigen(mat1)
  w <- apply(e1$vectors[,1:ct], 1, sum)
  w[which(is.na(w))] <- 0
  w
}
