#' Conditional score test for individual variants in a given variant-set
#'
#' The \code{Indiv_Score_Test_Region_cond} function takes in genotype,
#' the genotype of variants to be adjusted for in conditional analysis, and
#' the object from fitting the null model to analyze the conditional associations between
#' a quantitative/dichotomous phenotype and all individual variants in
#' a given variant-set by using score test, adjusting for a given list of variants.
#' @param genotype an n*p genotype matrix (dosage matrix) of the target sequence,
#' where n is the sample size and p is the number of genetic variants.
#' @param genotype_adj an n*p_adj genotype matrix (dosage matrix) of the target
#' sequence, where n is the sample size and p_adj is the number of genetic variants to be
#' adjusted for in conditional analysis (or a vector of a single variant with length n
#' if p_adj is 1).
#' @param obj_nullmodel an object from fitting the null model, which is the
#' output from either \code{\link{fit_null_glm}} function for unrelated samples or
#' \code{\link{fit_null_glmmkin}} function for related samples. Note that \code{\link{fit_null_glmmkin}}
#' is a wrapper of the \code{\link{glmmkin}} function from the \code{\link{GMMAT}} package.
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @param rv_num_cutoff_max the cutoff of maximum number of variants of analyzing
#' a given variant-set (default = 1e+09).
#' @param method_cond a character value indicating the method for conditional analysis.
#' \code{optimal} refers to regressing residuals from the null model on \code{genotype_adj}
#' as well as all covariates used in fitting the null model (fully adjusted) and taking the residuals;
#' \code{naive} refers to regressing residuals from the null model on \code{genotype_adj}
#' and taking the residuals (default = \code{optimal}).
#' @return A data frame with p rows corresponding to the p genetic variants in the given variant-set
#' and three columns: \code{Score_cond} (the conditional score test statistic adjusting for variants
#' in \code{genotype_adj}), \code{SE_cond} (the standard error associated with the
#' conditional score test statistic), and \code{pvalue_cond} (the conditional score test p-value).
#' If a variant in the given variant-set has minor allele frequency = 0 or
#' greater than \code{rare_maf_cutoff}, the corresponding row will be \code{NA}. If a variant in
#' the given variant-set has standard error equal to 0, the p-value will be set as 1.
#' @references Chen, H., et al. (2016). Control for population structure and relatedness for binary traits
#' in genetic association studies via logistic mixed models. \emph{The American Journal of Human Genetics}, \emph{98}(4), 653-666.
#' (\href{https://doi.org/10.1016/j.ajhg.2016.02.012}{pub})
#' @references Sofer, T., et al. (2019). A fully adjusted two-stage procedure for rank-normalization
#' in genetic association studies. \emph{Genetic Epidemiology}, \emph{43}(3), 263-275.
#' (\href{https://doi.org/10.1002/gepi.22188}{pub})
#' @export

Indiv_Score_Test_Region_cond <- function(genotype,genotype_adj,obj_nullmodel,
                                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                         rv_num_cutoff_max=1e9,
                                         method_cond=c("optimal","naive")){

  method_cond <- match.arg(method_cond) # evaluate choices
  if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(dim(genotype)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }

  results <- data.frame(Score_cond = rep(NA, dim(genotype)[2]),
                        SE_cond = rep(NA, dim(genotype)[2]),
                        pvalue_cond = rep(NA, dim(genotype)[2]))

  if(inherits(genotype, "sparseMatrix")){
    genotype <- as.matrix(genotype)
  }

  if(inherits(genotype_adj, "numeric")){
    genotype_adj <- matrix(genotype_adj, ncol=1)
  }

  if(dim(genotype)[1] != dim(genotype_adj)[1]){
    stop(paste0("Dimensions don't match for genotype and genotype_adj!"))
  }
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[,RV_label]

  rm(genotype,MAF)
  gc()

  if(sum(RV_label) >= rv_num_cutoff_max){
    stop(paste0("Number of rare variant in the set is more than ",rv_num_cutoff_max,"!"))
  }

  if(sum(RV_label) >= rv_num_cutoff){
    G <- as(Geno_rare,"dgCMatrix")
    rm(Geno_rare)
    gc()

    if(obj_nullmodel$relatedness){
      if(!obj_nullmodel$sparse_kins){
        P <- obj_nullmodel$P

        residuals.phenotype <- obj_nullmodel$scaled.residuals
        if(method_cond == "optimal"){
          residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj+obj_nullmodel$X-1)
        }else{
          residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj)
        }
        residuals.phenotype <- residuals.phenotype.fit$residuals
        X_adj <- model.matrix(residuals.phenotype.fit)
        PX_adj <- P%*%X_adj
        P_cond <- P - X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj) -
          PX_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj) +
          X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj)%*%X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj)
        rm(P)
        gc()

        results[RV_label,] <- do.call(cbind,Indiv_Score_Test_SMMAT(G,P_cond,residuals.phenotype))
      }else{
        Sigma_i <- obj_nullmodel$Sigma_i
        Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
        cov <- obj_nullmodel$cov

        residuals.phenotype <- obj_nullmodel$scaled.residuals
        if(method_cond == "optimal"){
          residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj+obj_nullmodel$X-1)
        }else{
          residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj)
        }
        residuals.phenotype <- residuals.phenotype.fit$residuals
        X_adj <- model.matrix(residuals.phenotype.fit)

        results[RV_label,] <- do.call(cbind,Indiv_Score_Test_SMMAT_sparse_cond(G,Sigma_i,Sigma_iX,cov,X_adj,residuals.phenotype))
      }
    }else{
      X <- model.matrix(obj_nullmodel)
      working <- obj_nullmodel$weights
      sigma <- sqrt(summary(obj_nullmodel)$dispersion)
      if(obj_nullmodel$family[1] == "binomial"){
        P <- diag(working) - X%*%solve(t(X)%*%diag(working)%*%X)%*%t(X)
      }else if(obj_nullmodel$family[1] == "gaussian"){
        P <- diag(length(working)) - X%*%solve(t(X)%*%X)%*%t(X)
      }

      residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values
      if(method_cond == "optimal"){
        residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj+model.matrix(obj_nullmodel)-1)
      }else{
        residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj)
      }
      residuals.phenotype <- residuals.phenotype.fit$residuals
      X_adj <- model.matrix(residuals.phenotype.fit)
      PX_adj <- P%*%X_adj
      P_cond <- P - X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj) -
        PX_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj) +
        X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj)%*%X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj)
      rm(P)
      gc()

      results[RV_label,] <- do.call(cbind,Indiv_Score_Test_SMMAT(G,P_cond,residuals.phenotype))
    }

    return(results)
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }

}

