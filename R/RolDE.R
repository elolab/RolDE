#' Robust longitudinal Differential Expression
#'
#' Detects longitudinal differential expression between two conditions (or groups) in time point
#' aligned data or in data with non-aligned time points. A rank product from the results of three independent modules,
#' RegROTS, DiffROTS and PolyReg, is determined to indicate the strength of differential expression of features
#' between the conditions / groups. \code{RolDE} tolerates a fair amount of missing values and is especially suitable
#' for noisy proteomics data.
#'
#' @param data the preprocessed normalized data as as a numerical matrix or as a SummarizedExperiment instance. Features (rows) and variables (columns)
#' of the data must have unique identifiers. If \code{data} is a SummarizedExperiment object, the design matrix must be included in the
#' \code{colData} argument of the \code{data} object.
#' @param des_matrix the design matrix for the \code{data}. Rows correspond to columns of
#' the \code{data}. Must contain four character columns (see included example design matrices \code{des_matrix1} or \code{des_matrix3}).
#' First column should contain sample (column) names of the \code{data}. Second column should indicate condition status (for each sample),
#' for which longitudinal differential expression is to be examined. Third column should indicate time point (time point aligned data)
#' or time value (non-aligned time point data) for each sample. Fourth column should provide the replicate (individual)
#' information for each sample as a numerical value. Each replicate or indivdual should have a distinct number. If \code{data} is a
#' SummarizedExperiment object, the design matrix must be included in the \code{colData} argument of the \code{data} object.
#' @param aligned logical; are the time points in different conditions and replicates (individuals) in the \code{data} aligned (fixed)? In aligned time point data,
#' the time points should be the same for each replicate (individual).
#' @param n_cores a positive integer. The number of threads used for parallel computing. If set to 1 (the default), no parallel computing is used.
#' @param model_type a string indicating the type of regression to be used for the PolyReg module and
#' the maximum level for which random effects should be allowed in the case of mixed models. Default "auto" for automatic selection.
#' @param sigValSampN a positive integer indicating the number of permutations for significance value
#' calculations. The overall used number will be \code{sigValSampN} * the number of rows in the \code{data}. Or set \code{sigValSampN} to
#' 0 to turn significance calculations off. Should be > 100000.
#' @param sig_adj_meth The multiple test hypothesis correction method for the estimated significance values. Only
#' relevant if \code{sigValSampN} is not 0.

#' @return \code{RolDE} returns a list with the following components: \cr
#' \cr
#' \code{RolDE_Results} a dataframe with the \code{RolDE} main results. Contains the \code{RolDE} rank product, the estimated
#' significance values (if \code{sigValSampN} is not set to 0) and the multiple hypothesis adjusted estimated significance values.\cr
#' \cr
#' \code{RegROTS_Results} a data frame of results for the RegROTS module. RegROTS internal rank products.\cr
#' \cr
#' \code{RegROTS_P_Values} a data frame of significance values for all the RegROTS runs.\cr
#' \cr
#' \code{DiffROTS_Results} a data frame of results for the DiffROTS module. DiffROTS internal rank products.\cr
#' \cr
#' \code{DiffROTS_P_Values} a data frame of the significance values for all the DiffROTS runs.\cr
#' \cr
#' \code{PolyReg_Results} a data frame of results for the PolyReg module. The representative (minimum) condition - related significance values.\cr
#' \cr
#' \code{PolyReg_P_Values} a data frame of all the condition - related significance values for the PolyReg module.\cr
#' \cr
#' \code{ROTS_Runs} a list containing the samples in the different ROTS runs for the RegROTS and DiffROTS (time point aligned data) modules.\cr
#' \cr
#' \code{Method_Degrees} a list containing the used degrees for the RegROTS and the PolyReg (and DiffROTS in non-aligned time point data) modules.\cr
#' \cr
#' \code{Input} a list of all the used inputs for \code{RolDE}.\cr
#'
#' @details \code{RolDE}, is a composite method, consisting of three independent
#' modules with different approaches to detecting longitudinal differential expression.
#' The combination of these diverse modules allows RolDE to robustly detect
#' varying differences in longitudinal trends and expression levels in
#' diverse data types and experimental settings.
#'
#' The *RegROTS* module merges the power of regression modelling with the power
#' of the established differential expression method Reproducibility Optimized Test
#' Statistic (ROTS) (Elo et al., Suomi et al.). A polynomial regression
#' model of protein expression over time is fitted separately for each replicate
#' (individual) in each condition. Differential expression between two replicates
#' (individuals) in different conditions is examined by comparing the coefficients
#' of the replicate-specific regression models. If all coefficient differences
#' are zero, no longitudinal differential expression between the two replicates (individuals)
#' in different conditions exist. For a thorough exploration of differential
#' expression between the conditions, all possible combinations of replicates (individuals)
#' in different conditions are examined.
#'
#' In the *DiffROTS* module the expression of replicates (individuals) in different
#' conditions are directly compared at all time points. Again, if the expression
#' level differences at all time points are zero, no differential expression
#' between the examined replicates (individuals) in different conditions exist. Similarly to
#' the RegROTS module, differential expression is examined between all possible
#' combinations of replicates (individuals) in the different conditions. In non-aligned
#' time point data, the expression level differences between the conditions
#' is examined when accounting for time-associated trends of varying complexity
#' in the data. More specifically, the expression level differences between
#' the conditions are examined when adjusting for increasingly complex
#' time-related expression trends of polynomial degrees d=0,1,.,d where d
#' is the maximum degree for the polynomial and the same degree as is used
#' for the PolyReg module.
#'
#' In the *PolyReg* module, polynomial regression modelling is used to
#' detect longitudinal differential expression. Condition is included as
#' a categorical factor within the models and by investigating the
#' condition related intercept and the polynomial termns at different
#' levels of the condition factor, average differences in expression
#' levels as well as differences in longitudinal expression patterns
#' between the conditions can be examined.
#'
#' Finally, to conclusively detect any differential expression,
#' the detections from the different modules are combined using the rank product.
#' For more details about the method, see the original \code{RolDE} publication (Valikangas et al.).
#'
#' By bare minimum, the user should provide \code{RolDE} the data in a normalized numerical matrix,
#' adjusted for confounding effects if needed, together with a suitable design matrix for the data. If
#' the time points in the data are non-aligned, the user should set the parameter \code{aligned} to \code{FALSE}.
#' Other parameter values \code{RolDE} determines automatically by default. The default values should be suitable
#' for a typical longitudinal differential expression analysis but the user is given control of many of the
#' parameters for \code{RolDE}.
#'
#' By default, \code{RolDE} assumes aligned time points in the data. If the time points
#' in the data are non-aligned, the user should set the parameter \code{aligned} to \code{FALSE}.
#'
#' Parallel processing can be enabled by setting the parameter \code{n_cores} as larger than the default 1 (highly recommended). With
#' parallel processing using multiple threads, the run time for \code{RolDE} can be significantly decreased. The parameter
#' \code{n_cores} controls the number of threads available for parallel processing.
#'
#' By default, \code{RolDE} uses fixed effects only regression with a common intercept and slope for the replicates (individuals) when time points
#' in the data are aligned and mixed effects models with a random effect for the individual baseline (intercept) if the time points are non aligned
#' for the PolyReg and the DiffROTS (only in data with non aligned time points) modules. This behaviour is controlled with the parameter \code{model_type}
#' and the default behaviour is induced when \code{model_type} is allowed to be "auto". However, the user can choose to use mixed effects regression modelling
#' when appropriate by setting the parameter \code{model_type} as "mixed0" for random effects for the individual baseline and
#' setting \code{model_type} as "mixed1" for an individual baseline and slope. Fixed effects only models can be chosen to be used by setting
#' as "fixed". Valid inputs for \code{model_type} are "auto" (the default), "mixed0", "mixed1" and "fixed".
#'
#' If the interest is only in ordering the features based on the strength of longitudinal differential expression between the conditions, \code{sigValSampN} can
#' be set to 0 to disable significance value estimation and to reduce the computational time used by \code{RolDE}. Otherwise, Parameter \code{sigValSampN}
#' indicates how many permutations should be performed when estimating the significance values. A larger value will lead to more accurate estimates but increases the
#' required computational time. The total number of permutataions for the significance value estimation will be approximately \code{sigValSampN}.
#' The default value used by \code{RolDE} is 500 000. The realized value of permutations might be sightly different, depending on the number of features
#' in the data. Using parallel processing greatly decreases the time needed for the significance value calculations. The estimated significance
#' values can be adjusted by any method in the \code{p.adjust} method in the \code{stats} package. Alternatively, q-values as defined by Storey et al. in the
#' Bioconductor package \code{qvalue} can be used. Valid values for \code{sig_adj_meth} are then:
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr","none", "qvalue". The default value is "fdr".
#' For more details about \code{RolDE}, see the original \code{RolDE} publication (Valikangas et al.)
#'
#' Please use \code{set.seed} for reproducibility.
#'
#' @references
#' Elo, Laura, Filen S, Lahesmaa R, et al. Reproducibility-optimized test statistic for ranking genes in microarray studies. IEEE/ACM Trans. Comput. Biol. Bioinform. 2008; 5:423-31.\cr
#' \cr
#' Suomi T, Seyednasrollah F, Jaakkola MK, et al. ROTS: An R package for reproducibility-optimized statistical testing. PLoS Comput. Biol. 2017; 13:5.\cr
#' \cr
#' Storey JD, Bass AJ, Dabney A, et al. qvalue: Q-value estimation for false discovery rate control. 2019. \cr
#' \cr
#' Välikangas T, Suomi T, ELo LL, et al. Enhanced longitudinal differential expression detection in proteomics with robust reproducibility optimization regression. bioRxiv 2021.
#'
#'@importFrom grDevices dev.off pdf
#'@importFrom stats lm median na.omit quantile sd weighted.mean p.adjust
#'@importFrom nlme lmeControl
#'@importFrom utils tail head
#'@importFrom methods is
#'@import doRNG
#'@import rngtools
#'@import SummarizedExperiment
#'@export
#'@examples
#'#Usage of RolDE in time point aligned data without significance value estimation and 1 core
#'data("data2")
#'data("des_matrix2")
#'set.seed(1) #For reproducibility.
#'data2.res<-RolDE(data=data2, des_matrix=des_matrix2, n_cores=1, sigValSampN = 0)
RolDE<-function(data, des_matrix=NULL, aligned=TRUE, n_cores=1, model_type="auto", sigValSampN=500000, sig_adj_meth="fdr"){

        if(is(data, "SummarizedExperiment")){
                des_matrix<-as.matrix(colData(data))
                data<-as.matrix(assay(data))
        }

        if(is.null(des_matrix)){
                if(!is(data, "SummarizedExperiment")){
                        stop("The design matrix argument is NULL and the data is not a SummarizedExperiment instance. Either the data and design matrix needs
                        to be provided separately or together as a SummarizedExperiment object.")
                }
        }

        results<-tryCatch(
                {
                        RolDE_Main(data, des_matrix, aligned, min_comm_diff="auto", min_feat_obs=3, degree_RegROTS="auto", degree_PolyReg="auto", n_cores, model_type, sigValSampN, sig_adj_meth)
                },
                error=function(errorCond){
                        message("Execution halted with the following problem:")
                        message(errorCond)
                        return(NULL)
                }
        )
        if(is.null(results)){stop("\nFailure during RolDE!\n")}
        return(results)
} #end function RolDE
