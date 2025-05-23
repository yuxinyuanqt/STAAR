#' STAAR-SPA procedure using omnibus test (dgCMatrix version)
#'
#' The \code{STAAR_Binary_SPA_sp} function is a sparse matrix (dgcMatrix) version of \code{STAAR_Binary_SPA}.
#' It takes in a sparse genotype matrix, the object from fitting the null model, and functional annotation data
#' to analyze the association between a imbalanced case-control phenotype and a variant-set by using STAAR-SPA procedure.
#' For each variant-set, the STAAR-B p-value is a p-value from an omnibus test
#' that aggregated Burden(1,25) and Burden(1,1) together with p-values of each test weighted by each annotation
#' using Cauchy method.
#' 
#' @param genotype_sp a sparse genotype matrix (dgCMatrix) of dimension n*p,
#' where n is the sample size and p is the number of genetic variants.
#' The matrix should be extracted using \code{STAARpipeline::Genotype_sp_extraction},
#' and it has been flipped: only the minor allele is stored (coded as 1), while the major allele is not stored (coded as 0).
#' @param MAF a numeric vector of minor allele frequencies for the variants in \code{genotype_sp}.
#' It can be computed using \code{Matrix::colMeans(genotype_sp,na.rm = TRUE)/2}.
#' @param obj_nullmodel an object from fitting the null model, which is the
#' output from either \code{\link{fit_null_glm}} function for unrelated samples or
#' \code{\link{fit_null_glmmkin}} function for related samples. Note that \code{\link{fit_null_glmmkin}}
#' is a wrapper of the \code{\link{glmmkin}} function from the \code{\link{GMMAT}} package.
#' @param annotation_phred a data frame or matrix of functional annotation data
#' of dimension p*q (or a vector of a single annotation score with length p).
#' Continuous scores should be given in PHRED score scale, where the PHRED score
#' of j-th variant is defined to be -10*log10(rank(-score_j)/total) across the genome. (Binary)
#' categorical scores should be taking values 0 or 1, where 1 is functional and 0 is
#' non-functional. If not provided, STAAR will perform the
#' Burden(1,25) and Burden(1,1) tests (default = NULL).
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @param rv_num_cutoff_max the cutoff of maximum number of variants of analyzing
#' a given variant-set (default = 1e+09).
#' @param tol a positive number specifying tolerance, the difference threshold for parameter
#' estimates in saddlepoint apporximation algorithm below which iterations should be stopped (default = ".Machine$double.eps^0.25").
#' @param max_iter a positive integers pecifying the maximum number of iterations for applying the saddlepoint approximation algorithm (default = "1000").
#' @param SPA_p_filter logical: are only the variants with a normal approximation based p-value smaller than a pre-specified threshold use the SPA method to recalculate the p-value, only used for imbalanced case-control setting (default = FALSE).
#' @param p_filter_cutoff threshold for the p-value recalculation using the SPA method, only used for imbalanced case-control setting (default = 0.05)
#' @return A list with the following members:
#' @return \code{num_variant}: the number of variants with minor allele frequency > 0 and less than
#' \code{rare_maf_cutoff} in the given variant-set that are used for performing the
#' variant-set using STAAR.
#' @return \code{cMAC}: the cumulative minor allele count of variants with
#' minor allele frequency > 0 and less than \code{rare_maf_cutoff} in the given variant-set.
#' @return \code{RV_label}: the boolean vector indicating whether each variant in the given
#' variant-set has minor allele frequency > 0 and less than \code{rare_maf_cutoff}.
#' @return \code{results_STAAR_B}: the STAAR-B p-value that aggregated Burden(1,25) and Burden(1,1) together
#' with p-values of each test weighted by each annotation using Cauchy method.
#' @return \code{results_STAAR_B_1_25}: a vector of STAAR-B(1,25) p-values,
#' including Burden(1,25) p-value weighted by MAF, the Burden(1,25)
#' p-values weighted by each annotation, and a STAAR-B(1,25)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_STAAR_B_1_1}: a vector of STAAR-B(1,1) p-values,
#' including Burden(1,1) p-value weighted by MAF, the Burden(1,1)
#' p-values weighted by each annotation, and a STAAR-B(1,1)
#' p-value by aggregating these p-values using Cauchy method.
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @references Liu, Y., et al. (2019). Acat: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(3), 410-421.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.01.002}{pub})
#' @references Li, Z., Li, X., et al. (2020). Dynamic scan procedure for
#' detecting rare-variant association regions in whole-genome sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(5), 802-814.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.03.002}{pub})
#' @export

STAAR_Binary_SPA_sp <- function(genotype_sp,MAF=NULL,obj_nullmodel,annotation_phred=NULL,
                                rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                rv_num_cutoff_max=1e9,
                                tol=.Machine$double.eps^0.25,max_iter=1000,
                                SPA_p_filter=FALSE,p_filter_cutoff=0.05){
  
  if(!inherits(genotype_sp, "matrix") && !inherits(genotype_sp, "Matrix")){
    stop("genotype is not a matrix!")
  }
  
  if(!inherits(genotype_sp, "sparseMatrix")){
    genotype_sp <- as(genotype_sp,"dgCMatrix")
  }
  
  if(dim(genotype_sp)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }
  
  if(is.null(MAF)){
    MAF <- Matrix::colMeans(genotype_sp,na.rm = TRUE)/2
  }
  
  annotation_phred <- as.data.frame(annotation_phred)
  if(dim(annotation_phred)[1] != 0 & (dim(genotype_sp)[2] != dim(annotation_phred)[1] | dim(genotype_sp)[2] != length(MAF))){
    stop(paste0("Dimensions don't match for genotype, MAF and annotation!"))
  }
  
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0)&(!is.na(MAF)))
  G <- genotype_sp[,RV_label]
  rm(genotype_sp)
  gc()
  
  annotation_phred <- annotation_phred[RV_label,,drop=FALSE]
  
  if(sum(RV_label) >= rv_num_cutoff_max){
    stop(paste0("Number of rare variant in the set is more than ",rv_num_cutoff_max,"!"))
  }
  
  if(sum(RV_label) >= rv_num_cutoff){
    MAF <- MAF[RV_label]
    if(any(is.na(G@x)))
    {
      G <- na.replace.sp(G,m=2*MAF)
    }
    
    annotation_rank <- 1 - 10^(-annotation_phred/10)
    
    ## beta(1,25)
    w_1 <- dbeta(MAF,1,25)
    ## beta(1,1)
    w_2 <- dbeta(MAF,1,1)
    if(dim(annotation_phred)[2] == 0){
      ## Burden, SKAT, ACAT-V
      w_B <- as.matrix(cbind(w_1,w_2))
    }else{
      ## Burden
      w_B_1 <- annotation_rank*w_1
      w_B_1 <- cbind(w_1,w_B_1)
      w_B_2 <- annotation_rank*w_2
      w_B_2 <- cbind(w_2,w_B_2)
      w_B <- cbind(w_B_1,w_B_2)
      w_B <- as.matrix(w_B)
    }
    
    ## use SPA for all p-value calculation
    if(!SPA_p_filter){
      if(obj_nullmodel$relatedness){
        if(!obj_nullmodel$sparse_kins){
          
          residuals.phenotype <- obj_nullmodel$scaled.residuals
          muhat <- obj_nullmodel$fitted.values
          
          pvalues <- STAAR_B_Binary_SPA_sp(G=G,XW=obj_nullmodel$XW,XXWX_inv=obj_nullmodel$XXWX_inv,
                                           residuals=residuals.phenotype,
                                           muhat=muhat,weights_B=w_B,
                                           tol=tol,max_iter=max_iter)
        }else{
          
          residuals.phenotype <- obj_nullmodel$scaled.residuals
          muhat <- obj_nullmodel$fitted.values
          
          pvalues <- STAAR_B_Binary_SPA_sp(G=G,XW=as.matrix(obj_nullmodel$XSigma_i),XXWX_inv=as.matrix(obj_nullmodel$XXSigma_iX_inv),
                                           residuals=residuals.phenotype,
                                           muhat=muhat,weights_B=w_B,
                                           tol=tol,max_iter=max_iter)
        }
      }else{
        
        residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values
        muhat <- obj_nullmodel$fitted.values
        
        pvalues <- STAAR_B_Binary_SPA_sp(G=G,XW=obj_nullmodel$XW,XXWX_inv=obj_nullmodel$XXWX_inv,
                                         residuals=residuals.phenotype,
                                         muhat=muhat,weights_B=w_B,
                                         tol=tol,max_iter=max_iter)
      }
    }else{
      if(obj_nullmodel$relatedness){
        if(!obj_nullmodel$sparse_kins){
          
          residuals.phenotype <- obj_nullmodel$scaled.residuals
          muhat <- obj_nullmodel$fitted.values
          
          P <- obj_nullmodel$P
          
          pvalues <- STAAR_B_SPA_SMMAT_sp(G=G,XW=obj_nullmodel$XW,XXWX_inv=obj_nullmodel$XXWX_inv,
                                          residuals=residuals.phenotype,
                                          muhat=muhat,weights_B=w_B,
                                          tol=tol,max_iter=max_iter,
                                          p_filter_cutoff=p_filter_cutoff,
                                          P=P)
        }else{
          
          residuals.phenotype <- obj_nullmodel$scaled.residuals
          muhat <- obj_nullmodel$fitted.values
          
          Sigma_i <- obj_nullmodel$Sigma_i
          Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
          cov <- obj_nullmodel$cov
          
          pvalues <- STAAR_B_SPA_SMMAT_sparse_sp(G=G,XW=as.matrix(obj_nullmodel$XSigma_i),XXWX_inv=as.matrix(obj_nullmodel$XXSigma_iX_inv),
                                                 residuals=residuals.phenotype,
                                                 muhat=muhat,weights_B=w_B,
                                                 tol=tol,max_iter=max_iter,
                                                 p_filter_cutoff=p_filter_cutoff,
                                                 Sigma_i=Sigma_i,Sigma_iX=Sigma_iX,cov=cov)
        }
      }else{
        
        residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values
        muhat <- obj_nullmodel$fitted.values
        
        X <- obj_nullmodel$X
        working <- obj_nullmodel$weights
        sigma <- sqrt(summary(obj_nullmodel)$dispersion)
        if(obj_nullmodel$family[1] == "binomial"){
          fam <- 1
        }else if(obj_nullmodel$family[1] == "gaussian"){
          fam <- 0
        }
        
        pvalues <- STAAR_B_SPA_sp(G=G,XW=obj_nullmodel$XW,XXWX_inv=obj_nullmodel$XXWX_inv,
                                  residuals=residuals.phenotype,
                                  muhat=muhat,weights_B=w_B,
                                  tol=tol,max_iter=max_iter,
                                  p_filter_cutoff=p_filter_cutoff,
                                  X=X,working=working,sigma=sigma,fam=fam)
      }
    }
    
    num_variant <- sum(RV_label) #dim(G)[2]
    cMAC <- sum(G@x)
    num_annotation <- dim(annotation_phred)[2]+1
    
    ## STAAR-B
    if(sum(is.na(pvalues))>0){
      ## all NAs
      if(sum(is.na(pvalues))==length(pvalues)){
        results_STAAR_B <- 1
      }else{
        ## not all NAs
        pvalues_sub <- pvalues[!is.na(pvalues)]
        if(sum(pvalues_sub[pvalues_sub<1])>0){
          ## not all ones
          results_STAAR_B <- CCT(pvalues_sub[pvalues_sub<1])
        }else{
          results_STAAR_B <- 1
        }
      }
    }else{
      if(sum(pvalues[pvalues<1])>0){
        results_STAAR_B <- CCT(pvalues[pvalues<1])
      }else{
        results_STAAR_B <- 1
      }
    }
    
    ## STAAR-B-1-25
    pvalues_sub <- pvalues[1:num_annotation]
    if(sum(is.na(pvalues_sub))>0){
      if(sum(is.na(pvalues_sub))==length(pvalues_sub)){
        pvalues_STAAR_B_1_25 <- 1
      }else{
        ## not all NAs
        pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
        if(sum(pvalues_sub[pvalues_sub<1])>0){
          ## not all ones
          pvalues_STAAR_B_1_25 <- CCT(pvalues_sub[pvalues_sub<1])
        }else{
          pvalues_STAAR_B_1_25 <- 1
        }
      }
    }else{
      if(sum(pvalues_sub[pvalues_sub<1])>0){
        pvalues_STAAR_B_1_25 <- CCT(pvalues_sub[pvalues_sub<1])
      }else{
        pvalues_STAAR_B_1_25 <- 1
      }
    }
    results_STAAR_B_1_25 <- c(pvalues[1:num_annotation],pvalues_STAAR_B_1_25)
    results_STAAR_B_1_25 <- data.frame(t(results_STAAR_B_1_25))
    
    ## STAAR-B-1-1
    pvalues_sub <- pvalues[(num_annotation+1):(2*num_annotation)]
    if(sum(is.na(pvalues_sub))>0){
      if(sum(is.na(pvalues_sub))==length(pvalues_sub)){
        pvalues_STAAR_B_1_1 <- 1
      }else{
        ## not all NAs
        pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
        if(sum(pvalues_sub[pvalues_sub<1])>0){
          ## not all ones
          pvalues_STAAR_B_1_1 <- CCT(pvalues_sub[pvalues_sub<1])
        }else{
          pvalues_STAAR_B_1_1 <- 1
        }
      }
    }else{
      if(sum(pvalues_sub[pvalues_sub<1])>0){
        pvalues_STAAR_B_1_1 <- CCT(pvalues_sub[pvalues_sub<1])
      }else{
        pvalues_STAAR_B_1_1 <- 1
      }
    }
    results_STAAR_B_1_1 <- c(pvalues[(num_annotation+1):(2*num_annotation)],pvalues_STAAR_B_1_1)
    results_STAAR_B_1_1 <- data.frame(t(results_STAAR_B_1_1))
    
    if(dim(annotation_phred)[2] == 0){
      colnames(results_STAAR_B_1_25) <- c("Burden(1,25)","STAAR-B(1,25)")
      colnames(results_STAAR_B_1_1) <- c("Burden(1,1)","STAAR-B(1,1)")
    }else{
      colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                          paste0("Burden(1,25)-",colnames(annotation_phred)),
                                          "STAAR-B(1,25)")
      colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                         paste0("Burden(1,1)-",colnames(annotation_phred)),
                                         "STAAR-B(1,1)")
    }
    
    return(list(num_variant = num_variant,
                cMAC = cMAC,
                RV_label = RV_label,
                results_STAAR_B = results_STAAR_B,
                results_STAAR_B_1_25 = results_STAAR_B_1_25,
                results_STAAR_B_1_1 = results_STAAR_B_1_1))
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }
  
}