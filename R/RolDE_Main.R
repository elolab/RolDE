#' Robust longitudinal Differential Expression
#'
#' Determines longitudinal differential expression between two conditions (or groups) in time point
#' aligned data or in data with non-aligned time points. A rank product from the results of three independent modules,
#' RegROTS, DiffROTS and PolyReg, is calculated to indicate the strength of differential expression of features
#' between the conditions / groups. \code{RolDE} tolerates a fair amount of missing values and is especially suitable
#' for noisy high-throughput data.
#'
#' @param data the normalized data as as a numerical matrix. Features (rows) and variables (columns) of the data must have unique identifiers.
#' @param des_matrix the design matrix for the \code{data}. Rows correspond to columns of
#' the \code{data}. Must contain four columns. First column should contain sample (column) names of the \code{data}.
#' Second column should indicate condition status (for each sample), for which longitudinal differential expression is to be examined. Third
#' column should indicate time point (time point aligned data) or time value (non-aligned time point data) for each sample. Fourth column should provide
#' the replicate (individual) information for each sample as a numeric value. Each replicate or indivdual should have a distinct number.
#' @param aligned logical; are the time points in different conditions and replicates (individuals) in the \code{data} aligned?
#' @param min_comm_diff a vector of two positive integers or string ("auto"). The minimum number of common time points for the replicates (individuals)
#' in different conditions to be compared (aligned time points) or the number of time points in the common time interval for the replicates in the different conditions
#' to be compared (non-aligned time points). The first integer refers to the minimum number of common time points for the RegROTS module
#' (aligned and non-aligned time points) and the second to DiffROTS (aligned time points). Second value needed but not used for DiffROTS when \code{aligned}
#' is set to \code{FALSE}.
#' @param min_ind_timep a positive integer. The minimum number of time points for a replicate (an individual)
#' to be included in the comparisons.
#' @param min_feat_obs a positive integer. The minimum number of non-missing obsevations a feature must have for an individual in a condition
#' to be included in the comparisons for the RegROTS module and the DiffROTS module (aligned time points).
#' @param degree_RegROTS a positive integer or string ("auto"). The degree of the polynomials used for the RegROTS module.
#' @param degree_PolyReg a positive integer or string ("auto"). The degree of the polynomials used for the PolyReg module.
#' @param B_for_ROTS a positive integer. The number of group preserving bootstraps used for ROTS method in the RegROTS and DiffROTS modules.
#' @param K_for_ROTS a positive integer or string ("auto"). The maximum size of the top list considered by ROTS.
#' @param rand_seed a positive integer seed for the random number generator
#' @param na_treat a string indicating the treatment of missing values.
#' @param n_cores a positive integer. The number of threads used for parallel computing. Only relevant if \code{doPar} is \code{TRUE}
#' @param polyregType a string indicating the type of regression to be used for the PolyReg module.
#' @param mix_degree a positive integer indicating the maximum level for which random effects should be
#' allowed. Only relevant if \code{polyregType} is "mixed".
#' @param estSigVal logical; should significance values be estimated?
#' @param sigValSampN a positive integer indicating the number of permutations for significance value
#' calculations. The overall used number will be \code{sigValSampN} * the number of rows in the \code{data}.
#' Only relevant if \code{estSigVal} is \code{TRUE}.
#' @param sig_adj_meth The multiple test hypothesis correction method for the estimated significance values. Only
#' relevant if \code{estSigVal} is \code{TRUE}.
#' @param add_refGroup logical; should a zero reference group be added for the RegROTS module? Typically does not need to be
#' altered.
#' @param iter_polyDegree logical; if a polynomial model of degree \code{degree_PolyReg} cannot be determined for a feature in the PolyReg module,
#' should RolDE try to determine lower level models for the feature?
#' @param weigh_runs logical; should the different runs in the RegROTS and DiffROTS modules be weighed according to the number of replicates/individuals
#' in the different runs when calculating internal rank products for the modules?
#'
#' @return \code{RolDE} returns a list with the following components: \cr
#' \cr
#' \code{RolDE_Results} a dataframe with the \code{RolDE} main results. Contains \code{RolDE} rank product, the estimated
#' significance values (if \code{estSigVal} is TRUE) and the multiple hypothesis adjusted estimated significance values.\cr
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
#' modules with different approaches to detecting differential expression.
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
#' are zero, no longitudinal differential expression between the two replicates
#' in different conditions exist. For a through exploration of differential
#' expression between the conditions, all possible combinations of replicates
#' in different conditions are examined.
#'
#' In the *DiffROTS* module the expression of replicates in different
#' conditions are directly compared at all time points. Again, if the expression
#' level differences at all time points are zero, no differential expression
#' between the examined replicates in different conditions exist. Similarly to
#' the RegROTS module, differential expression is examined between all possible
#' combinations of replicates in the different conditions. In non-aligned
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
#' a categorical factor within the models and by investigating the condition
#' related intercept and differences in time-related polynomial slopes at different
#' leves of condition,the average group level and the time related differences
#' in longitudinal expression between the conditions can be examined.
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
#' Parameter \code{min_comm_diff} controls how many common time points must two replicates (individuals)
#' have in different conditions to be compared. The first value controls the number of common time points
#' for the RegROTS module, while the second one controls the number of common time points for the DiffROTS module.
#' If \code{min_comm_diff} is set to "auto", \code{RolDE} will use a value of 3 for the RegROTS module and a value of
#' 1 for the DiffROTS module. Minimum values for the RegROTS and DiffROTS modules are 2 and 1, respectively.
#' In the case of data with non-aligned time points (\code{aligned} is set to \code{FALSE}), the first value of
#' \code{min_comm_diff} controls how many time values (or similar, e.g. age, temperature) must both replicates
#' in different conditions have in the common time interval to be compared. The common time interval for two replicates
#' r1 and r2 with time values t1 and t2 is defined as: \[max(min(t1,t2)),min(max(t1,t2))\]. In data with non-aligned
#' time points a value of =>1 for DiffROTS (the second value for \code{min_comm_diff}) is required but not used.
#' When \code{aligned} is \code{FALSE} an overall group comparison over all the replicates is performed
#' by the DiffROTS module.
#'
#' \code{min_ind_timep} determines the minimum number of time points/values all individuals included in the analysis
#' must have. If any individual has less time points than \code{min_ind_timep}, \code{RolDE} will throw an error. Such
#' replicates must be removed from the analysis. The default value is 3 which is also the minimum accepted value by RolDE for
#' a meaningful longitudinal analysis.
#'
#' \code{min_feat_obs} controls the number of non-missing values a feature must have for an individual in a condition to be
#' compared in the RegROTS module and the DiffROTS module (in data with aligned time points). A feature is
#' required to have at least \code{min_feat_obs} non-missing values for both replicates in the different conditions
#' to be compared. The default value used by \code{RoldE} is 3. If lowered, more missing values are allowed but
#' the analysis may become less accurate. In data with non-aligned time points, a common comparison over all the
#' replicates between the conditions is performed in the DiffROTS module and the number of allowed missing values
#' for a feature is controlled internally through other means.
#'
#' The user can control the degree of polynomials used by the RegROTS and the PolyReg modules via the
#' \code{degtree_RegROTS} and the \code{degree_PolyReg} parameters. If left to "auto", \code{RolDE} will by
#' default use as the \code{degree_RegROTS}=max(1, min(floor(median(t)/2),4)) and as the
#' \code{degree_PolyReg}=max(2, min((median(t)-1),5)), where t is a vector of the number of time points/values
#' in all the replicates.
#'
#' The parameter \code{B_for_ROTS} controls the number of group-preserving bootstraps performed
#' by ROTS when determining differential expression in the RegROTS and DiffROTS modules.
#' A larger number will lead to better estimates but will require more computational time. The default value
#' used by \code{RolDE} is 100.
#'
#' \code{K_for_ROTS} determines the largest top list size considered by ROTS when exploring
#' differential expression. The default value used by \code{RolDE} is 1/4 of the features.
#'
#' \code{na_treat} controls the treatment of missing values when calculating rank products. The default option
#' "last" will rank features with missing scores last and \code{RolDE} will deliver a final rank product
#' for all features. If \code{na_treat} is set to "keep", missing scores are retained as missing when
#' calculating internal rank products and some features might receive a missing final rank product in \code{RolDE}
#' if enough values are missing for the feature. A feature must have a valid (non missing) score from at least
#' one module to receive a non missing final rank product for \code{RolDE} if code{na_treat} is set to "keep".
#' Valid values for \code{na_treat} are"last" and "keep".
#'
#' Parallel processing can be enabled by setting the parameter \code{n_cores} as larger than the default 1 (highly recommended). With
#' parallel processing using multiple threads, the run time for \code{RolDE} can be significantly decreased. The parameter
#' \code{n_cores} controls the number of threads available for parallel processing.
#'
#' By default, \code{RolDE} uses fixed effects only regression with a common intercept and slopes for the replicates/individuals
#' in the PolyReg module and in the DiffROTS module (in data with non-aligned time points). However, the user
#' can choose to use mixed effects regression modelling when appropriate by setting the parameter \code{polyregType} as "mixed"
#' and by setting the correct level for the random effects for the replicates with the parameter \code{mix_degree}. When \code{polyregType}
#' is set to "mixed" and \code{mix_degree} is set to 0 (the default), a random intercept for each replicate will be used in the polynomial
#' regression models in the PolyReg module and the DiffROTS module (in data with non-aligned time points). Alternatively, a random effect also for
#' the linear term for each replicate can be allowed by setting the \code{mix_degree} to 1. Currently, random effects of higher polynomial degree
#' are not allowed in RolDE and the valid inputs for \code{mix_degree} are 0 and 1.
#'
#' If \code{estSigVal} is \code{TRUE}, significance values are estimated for all the features. If the interest is only in
#' ordering the features based on the strength of longitudinal differential expression between the conditions, \code{estSigVal} can
#' be set to \code{FALSE} to reduce the computational time used by \code{RolDE}. Parameter \code{sigValSampN} indicates how many permutations
#' should be performed when estimating the significance values. A larger value will lead to more accurate estimates but increases the
#' required computational time. The total number of permutataions for the significance value estimation will be approximately \code{sigValSampN}.
#' The default value used by \code{RolDE} is 500 000. The realized value of permutations might be sightly different, depending on the number of features
#' in the data. Using parallel processing greatly decreases the time needed for the significance value calculations. The estimated significance
#' values can be adjusted by any method in the \code{p.adjust} method in the \code{stats} package. Alternatively, q-values as defined by Storey et al. in the
#' Bioconductor package \code{qvalue} can be used. Valid values for \code{sig_adj_meth} are then:
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr","none", "qvalue". The default value is "fdr".
#'
#' The parameter \code{add_refGroup} controls whether a zero reference group should be added when determining differential expression in the
#' RegROTS and DiffROTS modules. The default is \code{TRUE}. If the examined data is very messy with lots of replicates / individuals and lots
#' of variation between and or/withing the replicates, setting \code{add_refGroup} to \code{FALSE} might improve the performance of \code{RolDE}.
#' This parameter will be deprecated? (It is not clear when/if this should be turned off)
#'
#' If \code{iter_polyDegree} is set to \code{TRUE} and a polynomial regression model using \code{degree_PolyReg} could not
#' be defined for a feature in the PolyReg module, \code{RolDE} iteratively decreases the degree by 1 for the given feature and tries
#' to fit a polynomial regression model of lower degree (until linear level, which is the lowest degree of model attempted to be fitted).
#' Setting \code{iter_polyDegree} to \code{TRUE} might help the PolyReg module to determine scores for some features if the number of missing
#' values is very high in the data. The default for \code{iter_polyDegree} in \code{RolDE} is \code{FALSE}.
#'
#' If the number of replicates in the ROTS runs in the RegROTS or DiffROTS modules is not equal, the unequal numbers
#' can be taken into account by weighing the runs according to the number of replicates included in each run
#' when calculating the internal rank products. This behaviour of \code{RolDE} is controlled with the parameter \code{weigh_runs}
#' which is \code{TRUE} by default.
#'
#' For more details about \code{RolDE}, see the original \code{RolDE} publication (Valikangas et al.)
#'
#' @references
#' Elo, Laura, Filen S, Lahesmaa R, et al. Reproducibility-optimized test statistic for ranking genes in microarray studies. IEEE/ACM Trans. Comput. Biol. Bioinform. 2008; 5:423-31.\cr
#' \cr
#' Suomi T, Seyednasrollah F, Jaakkola MK, et al. ROTS: An R package for reproducibility-optimized statistical testing. PLoS Comput. Biol. 2017; 13:5.\cr
#' \cr
#' Storey JD, Bass AJ, Dabney A, et al. qvalue: Q-value estimation for false discovery rate control. 2019. \cr
#' \cr
#' Valikangas T, Suomi T, Elo L.... RolDE reference (add)
#'
#'@importFrom grDevices dev.off pdf
#'@importFrom stats lm median na.omit quantile sd weighted.mean p.adjust
#'@importFrom utils head tail
#'@import knitr
#'@import printr
#'@import rmarkdown
#'@importFrom nlme lmeControl
#'@export
RolDE<-function(data, des_matrix, aligned=T, min_comm_diff="auto", min_ind_timep=3, min_feat_obs=3, degree_RegROTS="auto", degree_PolyReg="auto", B_for_ROTS=100, K_for_ROTS="auto", rand_seed=1, na_treat="last", n_cores=1, polyregType="fixed", mix_degree=0, estSigVal=T, sigValSampN=500000, sig_adj_meth="fdr", add_refGroup=T, iter_polyDegree=F, weigh_runs=T){

  if(!is.logical(aligned)){stop("Variable:aligned not logical. TRUE or FALSE needed.")}

  if(aligned){
    #The time points in the different conditions and replicates in the data are aligned.
    message("Running RolDE with time points aligned.")
    results<-tryCatch(
      {
        RolDE_Aligned(data, des_matrix, min_comm_diff, min_ind_timep, min_feat_obs, degree_RegROTS, degree_PolyReg, B_for_ROTS, K_for_ROTS, rand_seed, na_treat, n_cores, polyregType, mix_degree, estSigVal, sigValSampN, sig_adj_meth, add_refGroup, iter_polyDegree, weigh_runs)
      },
      error=function(errorCond){
        message("Execution halted with the following error:")
        message(errorCond)
        return(NULL)
      }
    )
  }else{
    #The time points in the different conditions and replicates in the data are notaligned
    message("Running RolDE with time points not aligned.")
    results<-tryCatch(
      {
        RolDE_NonAligned(data, des_matrix, min_comm_diff, min_ind_timep, min_feat_obs, degree_RegROTS, degree_PolyReg, B_for_ROTS, K_for_ROTS, rand_seed, na_treat, n_cores, polyregType, mix_degree, estSigVal, sigValSampN, sig_adj_meth, add_refGroup, iter_polyDegree, weigh_runs)
      },
      error=function(errorCond){
        message("Execution halted with the following error:")
        message(errorCond)
        return(NULL)
      }
    )
  } #end if(aligned)

  if(is.null(results)){stop(cat("\nError during RolDE!\n"))}

  return(results)
} #end function RolDE
