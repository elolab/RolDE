#' Robust longitudinal Differential Expression
#'
#' Detects longitudinal differential expression between two conditions (or groups) in time point
#' aligned data or in data with non-aligned time points. A rank product from the results of three independent modules,
#' RegROTS, DiffROTS and PolyReg, is determined to indicate the strength of differential expression of features
#' between the conditions / groups. \code{RolDE} tolerates a fair amount of missing values and is especially suitable
#' for noisy proteomics data.
#'
#' @param data the normalized data as as a numerical matrix or as a SummarizedExperiment instance. Features (rows) and variables (columns)
#' of the data must have unique identifiers. If \code{data} is a SummarizedExperiment object, the design matrix must be included in the
#' \code{colData} argument of the \code{data} object.
#' @param des_matrix the design matrix for the \code{data}. Rows correspond to columns of
#' the \code{data}. Must contain four character columns (see included example design matrices \code{des_matrix1} or \code{des_matrix3}).
#' First column should contain sample (column) names of the \code{data}. Second column should indicate condition status (for each sample),
#' for which longitudinal differential expression is to be examined. Third column should indicate time point (time point aligned data)
#' or time value (non-aligned time point data) for each sample. Fourth column should provide the replicate (individual)
#' information for each sample as a numerical value. Each replicate or indivdual should have a distinct number. If \code{data} is a
#' SummarizedExperiment object, the design matrix must be included in the \code{colData} argument of the \code{data} object.
#' @param aligned logical; are the time points in different conditions and replicates (individuals) in the \code{data} aligned? In aligned time point data,
#' the time points should be the same for each replicate (individual).
#' @param min_comm_diff a vector of two positive integers or string ("auto"). The minimum number of common time points for the replicates (individuals)
#' in different conditions to be compared (aligned time points) or the number of time points in the common time interval for the replicates (individuals)
#' in the different conditions
#' to be compared (non-aligned time points). The first integer refers to the minimum number of common time points for the RegROTS module
#' (aligned and non-aligned time points) and the second to DiffROTS (aligned time points). Second value needed but not used for DiffROTS when \code{aligned}
#' is set to \code{FALSE}.
#' @param min_feat_obs a positive integer. The minimum number of non-missing obsevations a feature must have for a replicate (individual) in a condition
#' to be included in the comparisons for the RegROTS module and the DiffROTS module (aligned time points).
#' @param degree_RegROTS a positive integer or string ("auto"). The degree of the polynomials used for the RegROTS module.
#' @param degree_PolyReg a positive integer or string ("auto"). The degree of the polynomials used for the PolyReg module.
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
#' in different conditions exist. For a through exploration of differential
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
#' Parameter \code{min_comm_diff} controls how many common time points must two replicates (individuals)
#' have in different conditions to be compared. The first value controls the number of common time points
#' for the RegROTS module, while the second one controls the number of common time points for the DiffROTS module.
#' If \code{min_comm_diff} is set to "auto", \code{RolDE} will use a value of 3 for the RegROTS module and a value of
#' 1 for the DiffROTS module. Minimum values for the RegROTS and DiffROTS modules are 2 and 1, respectively.
#' In the case of data with non-aligned time points (\code{aligned} is set to \code{FALSE}), the first value of
#' \code{min_comm_diff} controls how many time values (or similar, e.g. age, temperature) must both replicates (individuals)
#' in different conditions have in the common time interval to be compared. The common time interval for two replicates (individuals)
#' r1 and r2 with time values t1 and t2 is defined as: \[max(min(t1,t2)),min(max(t1,t2))\]. In data with non-aligned
#' time points a value of =>1 for DiffROTS (the second value for \code{min_comm_diff}) is required but not used.
#' When \code{aligned} is \code{FALSE} an overall group comparison over all the replicates (individuals) is performed
#' by the DiffROTS module.
#'
#' \code{min_feat_obs} controls the number of non-missing values a feature must have for a replicate (an individual) in a condition to be
#' compared in the RegROTS module and the DiffROTS module (in data with aligned time points). A feature is
#' required to have at least \code{min_feat_obs} non-missing values for both replicates (individuals) in the different conditions
#' to be compared. The default value used by \code{RoldE} is 3. If lowered, more missing values are allowed but
#' the analysis may become less accurate. In data with non-aligned time points, a common comparison over all the
#' replicates (individuals) between the conditions is performed in the DiffROTS module and the number of allowed missing values
#' for a feature is controlled internally through other means.
#'
#' The user can control the degree of polynomials used by the RegROTS and the PolyReg modules via the
#' \code{degtree_RegROTS} and the \code{degree_PolyReg} parameters. If left to "auto", \code{RolDE} will by
#' default use as the \code{degree_RegROTS}=max(1, min(floor(median(t)/2),4)) and as the
#' \code{degree_PolyReg}=max(2, min((median(t)-1),5)), where t is a vector of the number of time points/values
#' for all the replicates (individuals).
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
#'@import doRNG
#'@import rngtools
#'@import SummarizedExperiment
#'@importFrom methods is
#'@export
#'@examples
#'#Usage of RolDE in time point aligned data without significance value estimation and 1 core
#'data("data2")
#'data("des_matrix2")
#'set.seed(1) #For reproducibility.
#'data2.res<-RolDE_Main(data=data2, des_matrix=des_matrix2, n_cores=1, sigValSampN = 0)
RolDE_Main<-function(data, des_matrix=NULL, aligned=TRUE, min_comm_diff="auto", min_feat_obs=3, degree_RegROTS="auto", degree_PolyReg="auto", n_cores=1, model_type="auto", sigValSampN=500000, sig_adj_meth="fdr"){

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

  #Save input
  input_list<-list(data, des_matrix, aligned, min_comm_diff, min_feat_obs, degree_RegROTS, degree_PolyReg, n_cores, model_type, sigValSampN, sig_adj_meth)
  names(input_list)<-c("data", "des_matrix", "aligned", "min_comm_diff", "min_feat_obs", "degree_RegROTS", "degree_PolyReg", "n_cores", "model_type", "sigValSampN", "sig_adj_meth")

  #Validate design matrix and data
  message("Validating input")
  input_validated<-tryCatch(
    {
      validateInput(data, des_matrix, aligned, min_comm_diff, min_feat_obs, n_cores, model_type, sigValSampN, sig_adj_meth)
    },
    error=function(errorCond){
      message("Execution halted with the following error:")
      message(errorCond)
      return(NULL)
    }
  )

  #If problems during input validation
  if(is.null(input_validated)){stop(cat("\nError during input validation\n"))}

  data<-input_validated[[1]]
  des_matrix<-input_validated[[2]]
  min_comm_diff<-input_validated[[3]]
  model_type<-input_validated[[4]]
  estSigVal<-input_validated[[5]]
  message("Input validated")

  #Determine ROTS runs
  message("Preparing the Analysis")
  rots_runs<-tryCatch(
    {
      determineROTSRuns(data, des_matrix, min_comm_diff[1])
    },
    error=function(errorCond){
      message("Execution halted with the following error:")
      message(errorCond)
      return(NULL)
    }
  )

  #If problems during run determination
  if(is.null(rots_runs)){stop(cat("\nError during determining ROTS runs\n"))}
  rots_runs<-rots_runs[[1]]

  #Validate and determine the degrees.
  method_degrees<-tryCatch(
    {
      determineDegrees(des_matrix, degree_RegROTS, degree_PolyReg)
    },
    error=function(errorCond){
      message("Execution halted with the following error:")
      message(errorCond)
      return(NULL)
    }
  )

  #If problems during degree determination
  if(is.null(method_degrees)){stop(cat("\nError during degree determination\n"))}

  degree_RegROTS<-method_degrees[[1]]
  degree_PolyReg<-method_degrees[[2]]
  #Now we have everything to start running.

  #Run RegROTS
  message("Running RegROTS")
  res_regrots<-tryCatch(
    {
      RegROTS(data, des_matrix, min_feat_obs, degree_RegROTS, rots_runs, n_cores, aligned)
    },
    error=function(errorCond){
      message("Execution halted with the following error:")
      message(errorCond)
      return(NULL)
    }
  )

  #If a problem arises during RegROTS
  if(is.null(res_regrots)){stop(cat("\nError during RegROTS\n"))}

  reg_rots_pval<-res_regrots[[1]]
  regrots_weigths<-res_regrots[[3]]
  res_regrots<-res_regrots[[2]]

  #Run DiffROTS
  message("Running DiffROTS")
  if(aligned){
    #DiffROTS aligned
    res_diffrots<-tryCatch(
      {
        DiffROTSAligned(data, des_matrix, min_comm_diff[2], min_feat_obs, rots_runs, n_cores)
      },
      error=function(errorCond){
        message("Execution halted with the following error:")
        message(errorCond)
        return(NULL)
      }
    )
    #If a problem arises during DiffROTS
    if(is.null(res_diffrots)){stop(cat("\nError during DiffROTS\n"))}

    diff_rots_pval<-res_diffrots[[1]]
    diffrots_weigths<-res_diffrots[[3]]
    res_diffrots<-res_diffrots[[2]]

  }else{
    #DiffROTS non aligned
    if(model_type=="fixed"){ #Use the same parameter for DiffROTSAligned and PolyReg? Or should if be able to controlled separately? CHECK!
      res_diffrots<-tryCatch(
        {
          DiffROTSNonAligned(data, des_matrix, degree_PolyReg, n_cores)
        },
        error=function(errorCond){
          message("Execution halted with the following error:")
          message(errorCond)
          return(NULL)
        }
      ) #end trycatch DiffROTS
    }else{
      res_diffrots<-tryCatch(
        {
          DiffROTSMixNonAligned(data, des_matrix, degree_PolyReg, n_cores, model_type)
        },
        error=function(errorCond){
          message("Execution halted with the following error:")
          message(errorCond)
          return(NULL)
        }
      ) #end trycatch DiffROTS
    }

    #If a problem arises during DiffROTS, halt.
    if(is.null(res_diffrots)){stop(cat("\nError during DiffROTS\n"))}

    diffrots_weigths=NA
    diff_rots_pval<-res_diffrots[[1]]
    res_diffrots<-res_diffrots[[2]]
  }

  #Add weights to the input list
  input_list[[length(input_list)+1]]<-regrots_weigths
  input_list[[length(input_list)+1]]<-diffrots_weigths
  names(input_list)[c((length(input_list)-1):length(input_list))]<-c("regrots_weights", "diffrots_weights")

  #Run PolyReg
  message("Running PolyReg")
  if(model_type=="fixed"){
    message("Running PolyReg with fixed effects only.")
    res_polyreg<-tryCatch(
      {
        PolyReg(data, des_matrix, degree_PolyReg, n_cores)
      },
      error=function(errorCond){
        message("Execution halted with the following error:")
        message(errorCond)
        return(NULL)
      }
    )
  }else{
    message("Running PolyReg with mixed effects models.")
    res_polyreg<-tryCatch(
      {
        PolyRegMix(data, des_matrix, degree_PolyReg, n_cores, model_type)
      },
      error=function(errorCond){
        message("Execution halted with the following error:")
        message(errorCond)
        return(NULL)
      }
    ) #end trycatch PolyReg
  }

  #If a problem arises during PolyReg
  if(is.null(res_polyreg)){stop(cat("\nError during PolyReg\n"))}

  polyreg_pval<-res_polyreg[[1]]
  res_polyreg<-res_polyreg[[2]]

  #Now we should have everything to calculate the final rank product
  message("Calculating the final rank product.")
  rank_products<-tryCatch(
    {
      calculateRP(res_regrots, res_diffrots, res_polyreg)
    },
    error=function(errorCond){
      message("Execution halted with the following error:")
      message(errorCond)
      return(NULL)
    }
  )

  if(is.null(rank_products)){stop(cat("\nError during rank product calculation\n"))}

  #Estimate signigicance values
  if(estSigVal){
    message("Estimating significance values.")
    if(aligned){
      rp_with_est_sig<-tryCatch(
        {
          estSigValAlignedWeights(data, res_regrots, res_diffrots, res_polyreg, reg_rots_pval, diff_rots_pval, polyreg_pval, rank_products, n_cores, sig_adj_meth, sigValSampN, regrots_weigths, diffrots_weigths)
        },
        error=function(errorCond){
          message("Execution halted with the following error:")
          message(errorCond)
          return(NULL)
        }
      )
    }else{
      rp_with_est_sig<-tryCatch(
        {
          estSigValNonAlignedWeights(data, res_regrots, res_diffrots, res_polyreg, reg_rots_pval, diff_rots_pval, polyreg_pval, rank_products, n_cores, sig_adj_meth, sigValSampN, regrots_weigths)
        },
        error=function(errorCond){
          message("Execution halted with the following error:")
          message(errorCond)
          return(NULL)
        }
      )
    }
    fin_res<-rp_with_est_sig
    if(is.null(fin_res)){stop(cat("\nError during significance value estimation\n"))}
  } else {
    fin_res<-rank_products
  }

  ret_list<-list(fin_res, res_regrots, reg_rots_pval, res_diffrots, diff_rots_pval, res_polyreg, polyreg_pval, rots_runs, method_degrees, input_list)
  names(ret_list)<-c("RolDE_Results", "RegROTS_Results", "RegROTS_P_Values", "DiffROTS_Results", "DiffROTS_P_Values", "PolyReg_Results", "PolyReg_P_Values", "ROTS_Runs", "Method_Degrees", "Input")
  return(ret_list)
} #end function RolDE_Main
