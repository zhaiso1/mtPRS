#' Multi-trait PRS association test
#'
#' Aggregate p-values from mtPRS-PCA, mtPRS-ML (Machine Learning), and stPRSs together with Cauchy Combination Test (CCT) for the robust multi-trait PRS association test
#' @param dat a list of inputs including disease GWAS summary statistics of K traits from base cohort, individual-level data from base cohort, genetic correlation matrix, and true effect sizes
#' @param pcut p-value cutoff for C+T method to construct individual stPRS
#' @param varcut choose top principal components (PCs) until explaining varcut-percent variance
#' @param K number of traits
#' @param phenotype the type of phenotype in target cohort, either "dis" or "pgx"
#' @details mtPRS_O needs disease GWAS summary statistics from K traits, and genetic correlation matrix among K traits
#' @return if the phenotype is from disease GWAS, then return main p-value; if the phenotype is from PGx GWAS (with both T and C arms), then return main p-value from 2df test, interaction p-value from 2df test, main p-value in T arm only, and main p-value in C arm only
#' @references Zhai, S., Guo, B., Wu, B., Mehrotra, D.V., and Shen, J., 2023. Integrating multiple traits for improving polygenic risk prediction in disease and pharmacogenomics GWAS (submitted).
#' @author Song Zhai
#' @export
#' @examples
#' \donttest{
#' dat <- generate_pgx_data(structure = "clustered", sparseness = "more",
#' rho_DT = c(0.5,0.5,0.5,0.5), rho_T = 0.5, rho_E = 0.3, rho_C = 0.1,
#' K=4, m=1000, pcausal=0.1, blocksize=100,
#' gamma=1, samplesize=700, h2_base=0.3, h2_target=0.3)
#' re <- mtPRS_O(dat, pcut = 0.05, varcut = 0.8, K = 4, phenotype = "pgx")
#' print(re)
#' }
#'
mtPRS_O <- function(dat, pcut, varcut, K, phenotype){
  #devtools::install_github("yaowuliu/ACAT")
  #library(ACAT)

  re <- mtPRS_PCA(dat, pcut, varcut, K, phenotype)
  mtPRSxPCA <- re$mtPRS
  stPRS <- re$stPRS
  mtPRSxML <- mtPRS_ML(dat, pcut, K, phenotype)
  mtPRSxMLxprog <- mtPRSxML$PRSprog
  mtPRSxMLxpred <- mtPRSxML$PRSpred
  if(phenotype == "dis"){
    target <- dat$target; Y <- target[,1]

    p0 <- double()
    for (k in 1:K) {fit0 <- summary(lm(Y ~ stPRS[,k])); p0 <- c(p0, fit0$coefficients[2,4])}
    fit1 <- summary(lm(Y ~ mtPRSxPCA)); p1 <- fit1$coefficients[2,4]
    fit2 <- summary(lm(Y ~ mtPRSxMLxprog)); p2 <- fit2$coefficients[2,4]

    pp <- c(p0,p1,p2)
    re <- ACAT(pp)
  }
  if(phenotype == "pgx"){
    target <- dat$target; Y <- target[,1]; Tr <- target[,2]

    p01 <- p02 <- p03 <- p04 <- double()
    for (k in 1:K) {
      fit0 <- summary(lm(Y ~ Tr*stPRS[,k]))
      p01 <- c(p01, fit0$coefficients[3,4])
      p02 <- c(p02, fit0$coefficients[4,4])

      index <- which(Tr==1)
      fit0 <- summary(lm(Y[index] ~ stPRS[index,k])); p03 <- c(p03, fit0$coefficients[2,4])
      index <- which(Tr==0)
      fit0 <- summary(lm(Y[index] ~ stPRS[index,k])); p04 <- c(p04, fit0$coefficients[2,4])
    }

    fit1 <- summary(lm(Y ~ Tr*mtPRSxPCA))
    p11 <- fit1$coefficients[3,4]
    p12 <- fit1$coefficients[4,4]
    index <- which(Tr==1)
    fit1 <- summary(lm(Y[index] ~ mtPRSxPCA[index])); p13 <- fit1$coefficients[2,4]
    index <- which(Tr==0)
    fit1 <- summary(lm(Y[index] ~ mtPRSxPCA[index])); p14 <- fit1$coefficients[2,4]

    fit2 <- summary(lm(Y ~ Tr + mtPRSxMLxprog + Tr:mtPRSxMLxpred))
    p21 <- fit2$coefficients[3,4]
    p22 <- fit2$coefficients[4,4]
    mtPRSxML0 <- mtPRSxMLxprog + mtPRSxMLxpred
    index <- which(Tr==1)
    fit2 <- summary(lm(Y[index] ~ mtPRSxML0[index])); p23 <- fit2$coefficients[2,4]
    index <- which(Tr==0)
    fit2 <- summary(lm(Y[index] ~ mtPRSxML0[index])); p24 <- fit2$coefficients[2,4]

    p_2df_prog <- c(p01,p11,p21)
    p_2df_pred <- c(p02,p12,p22)
    p_T_prog <- c(p03,p13,p23)
    p_C_prog <- c(p04,p14,p24)

    re <- c(ACAT(p_2df_prog),ACAT(p_2df_pred),ACAT(p_T_prog),ACAT(p_C_prog))
    names(re) <- c("p_2df_prog","p_2df_pred","p_T_prog","p_C_prog")
  }
  re
}

mtPRS_ML <- function(dat, pcut = 0.05, K = 4, phenotype){
  base <- dat$base
  target <- dat$target
  if(phenotype == "pgx"){
    G.target <- target[,-c(1:2)] %>% as.matrix()
    Y <- target[,1]; Tr <- target[,2]
    X <- calculate_PRS(K, base, pcut, G.target)
    T_X_TX <- cbind(Tr,X,Tr*X)

    cv_model <- cv.glmnet(T_X_TX, Y, alpha = 0.5)
    best_lambda <- cv_model$lambda.min
    best_model <- glmnet(T_X_TX, Y, alpha = 0.5, lambda = best_lambda, intercept = TRUE)

    betax <- best_model$beta; w.prog <- betax[2:(K+1)]; w.pred <- betax[(K+2):(2*K+1)]
    w.prog[which(is.na(w.prog))] <- 0; w.pred[which(is.na(w.pred))] <- 0

    PRS.prog <- X %*% w.prog %>% as.vector(); PRS.pred <- X %*% w.pred %>% as.vector()
    re <- list(wprog=w.prog,wpred=w.pred,PRSprog=PRS.prog,PRSpred=PRS.pred)
  }
  if(phenotype == "dis"){
    G.target <- target[,-1] %>% as.matrix()
    Y <- target[,1]
    X <- calculate_PRS(K, base, pcut, G.target)

    cv_model <- cv.glmnet(X, Y, alpha = 0.5)
    best_lambda <- cv_model$lambda.min
    best_model <- glmnet(X, Y, alpha = 0.5, lambda = best_lambda, intercept = TRUE)

    w <- best_model$beta; w <- w[1:4]; w[which(is.na(w))] <- 0

    PRS <- X %*% w %>% as.vector()

    re <- list(w=w,PRSprog=PRS,PRSpred=PRS)
  }
  re
}

