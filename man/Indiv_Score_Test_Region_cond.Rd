% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Indiv_Score_Test_Region_cond.R
\name{Indiv_Score_Test_Region_cond}
\alias{Indiv_Score_Test_Region_cond}
\title{Conditional score test for individual variants in a given variant-set}
\usage{
Indiv_Score_Test_Region_cond(
  genotype,
  genotype_adj,
  obj_nullmodel,
  rare_maf_cutoff = 0.01,
  rv_num_cutoff = 2,
  rv_num_cutoff_max = 1e+09,
  method_cond = c("optimal", "naive")
)
}
\arguments{
\item{genotype}{an n*p genotype matrix (dosage matrix) of the target sequence,
where n is the sample size and p is the number of genetic variants.}

\item{genotype_adj}{an n*p_adj genotype matrix (dosage matrix) of the target
sequence, where n is the sample size and p_adj is the number of genetic variants to be
adjusted for in conditional analysis (or a vector of a single variant with length n
if p_adj is 1).}

\item{obj_nullmodel}{an object from fitting the null model, which is the
output from either \code{\link{fit_null_glm}} function for unrelated samples or
\code{\link{fit_null_glmmkin}} function for related samples. Note that \code{\link{fit_null_glmmkin}}
is a wrapper of the \code{\link{glmmkin}} function from the \code{\link{GMMAT}} package.}

\item{rare_maf_cutoff}{the cutoff of maximum minor allele frequency in
defining rare variants (default = 0.01).}

\item{rv_num_cutoff}{the cutoff of minimum number of variants of analyzing
a given variant-set (default = 2).}

\item{rv_num_cutoff_max}{the cutoff of maximum number of variants of analyzing
a given variant-set (default = 1e+09).}

\item{method_cond}{a character value indicating the method for conditional analysis.
\code{optimal} refers to regressing residuals from the null model on \code{genotype_adj}
as well as all covariates used in fitting the null model (fully adjusted) and taking the residuals;
\code{naive} refers to regressing residuals from the null model on \code{genotype_adj}
and taking the residuals (default = \code{optimal}).}
}
\value{
A data frame with p rows corresponding to the p genetic variants in the given variant-set
and three columns: \code{Score_cond} (the conditional score test statistic adjusting for variants
in \code{genotype_adj}), \code{SE_cond} (the standard error associated with the
conditional score test statistic), and \code{pvalue_cond} (the conditional score test p-value).
If a variant in the given variant-set has minor allele frequency = 0 or
greater than \code{rare_maf_cutoff}, the corresponding row will be \code{NA}. If a variant in
the given variant-set has standard error equal to 0, the p-value will be set as 1.
}
\description{
The \code{Indiv_Score_Test_Region_cond} function takes in genotype,
the genotype of variants to be adjusted for in conditional analysis, and
the object from fitting the null model to analyze the conditional associations between
a quantitative/dichotomous phenotype and all individual variants in
a given variant-set by using score test, adjusting for a given list of variants.
}
\references{
Chen, H., et al. (2016). Control for population structure and relatedness for binary traits
in genetic association studies via logistic mixed models. \emph{The American Journal of Human Genetics}, \emph{98}(4), 653-666.
(\href{https://doi.org/10.1016/j.ajhg.2016.02.012}{pub})

Sofer, T., et al. (2019). A fully adjusted two-stage procedure for rank-normalization
in genetic association studies. \emph{Genetic Epidemiology}, \emph{43}(3), 263-275.
(\href{https://doi.org/10.1002/gepi.22188}{pub})
}
