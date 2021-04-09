#Run RolDE on data with non-aligned time points.
#Return results.
RolDE_NonAligned<-function(data, des_matrix, min_comm_diff, min_ind_timep, min_feat_obs, degree_RegROTS, degree_PolyReg, B_for_ROTS, K_for_ROTS, rand_seed, na_treat, n_cores, polyregType, mix_degree, estSigVal, sigValSampN, sig_adj_meth, add_refGroup, iter_polyDegree, weigh_runs){

  #Save input
  input_list=list(data, des_matrix, min_comm_diff, min_ind_timep, min_feat_obs, degree_RegROTS, degree_PolyReg, B_for_ROTS, K_for_ROTS, rand_seed, na_treat, n_cores, polyregType, mix_degree, estSigVal, sigValSampN, sig_adj_meth, add_refGroup, iter_polyDegree, weigh_runs)
  names(input_list)=c("data", "des_matrix", "min_comm_diff", "min_ind_timep", "min_feat_obs", "degree_RegROTS", "degree_PolyReg", "B_for_ROTS", "K_for_ROTS", "rand_seed", "na_treat", "n_cores", "polyregType", "mix_degree", "estSigVal", "sigValSampN", "sig_adj_meth", "add_refGroup", "iter_polyDegree", "weigh_runs")

  #Validate design matrix and data
  message("Validating input")
  input_validated<-tryCatch(
    {
      validateNonAlignedData(data, des_matrix, min_comm_diff, min_ind_timep, min_feat_obs, B_for_ROTS, K_for_ROTS, rand_seed, na_treat, n_cores, polyregType, mix_degree, estSigVal, sigValSampN, sig_adj_meth, add_refGroup, iter_polyDegree, weigh_runs)
    },
    error=function(errorCond){
      message("Execution halted with the following error:")
      message(errorCond)
      return(NULL)
    }
  ) #end trycatch organize aligned data

  #If problems during input validation
  if(is.null(input_validated)){stop(cat("\nError during input validation\n"))}

  data<-input_validated[[1]]
  des_matrix<-input_validated[[2]]
  test_type<-input_validated[[3]]
  min_comm_diff<-input_validated[[4]]
  message("Input validated.")

  #Determine ROTS runs
  message("Preparing the Analysis")
  rots_runs<-tryCatch(
    {
      determineROTSRunsNonAligned(data, des_matrix, test_type, min_comm_diff[1])
    },
    error=function(errorCond){
      message("Execution halted with the following error:")
      message(errorCond)
      return(NULL)
    }
  ) #end trycatch get ROTS runs

  #If problems during run determination
  if(is.null(rots_runs)){stop(cat("\nError during determining ROTS runs\n"))}
  rots_runs<-rots_runs[[1]]

  #Validate and determine the degrees.
  method_degrees<-tryCatch(
    {
      determineDegreesNonAligned(des_matrix, degree_RegROTS, degree_PolyReg)
    },
    error=function(errorCond){
      message("Execution halted with the following error:")
      message(errorCond)
      return(NULL)
    }
  ) #end trycatch define degrees

  #If problems degree determination
  if(is.null(method_degrees)){stop(cat("\nError during degree determination\n"))}

  degree_RegROTS=method_degrees[[1]]
  degree_PolyReg=method_degrees[[2]]
  #Now we have everything to start running.

  #Run RegROTS
  message("Running RegROTS")
  res_regrots<-tryCatch(
    {
      RegROTS(data, des_matrix, min_feat_obs, degree_RegROTS, rots_runs, B_for_ROTS, K_for_ROTS, rand_seed, na_treat, n_cores, add_refGroup, weigh_runs)
    },
    error=function(errorCond){
      message("Execution halted with the following error:")
      message(errorCond)
      return(NULL)
    }
  ) #end trycatch RegROTS

  #If a problem arises during RegROTS, halt.
  if(is.null(res_regrots)){stop(cat("\nError during RegROTS\n"))}

  reg_rots_pval<-res_regrots[[1]]
  if(weigh_runs){
    regrots_weigths<-res_regrots[[3]]
  }
  res_regrots<-res_regrots[[2]]

  #Run DiffROTS
  message("Running DiffROTS")
  if(polyregType=="fixed"){ #Use the same parameter for DiffROTSAligned and PolyReg? Or should if be able to controlled separately? CHECK!
    res_diffrots<-tryCatch(
      {
        DiffROTSNonAligned(data, des_matrix, degree_PolyReg, B_for_ROTS, K_for_ROTS, rand_seed, na_treat, n_cores, test_type)
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
        DiffROTSMixNonAligned(data, des_matrix, degree_PolyReg, B_for_ROTS, K_for_ROTS, rand_seed, na_treat, n_cores, test_type, mix_degree)
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

  diff_rots_pval<-res_diffrots[[1]]
  res_diffrots<-res_diffrots[[2]]

  #Run PolyReg
  message("Running PolyReg")
  if(polyregType=="fixed"){
    message("Running PolyReg with fixed effects only.")
    res_polyreg<-tryCatch(
      {
        PolyReg(data, des_matrix, degree_PolyReg, n_cores, test_type, iter_polyDegree)
      },
      error=function(errorCond){
        message("Execution halted with the following error:")
        message(errorCond)
        return(NULL)
      }
    ) #end trycatch PolyReg
  }else{
    message("Running PolyReg with mixed effects models.")
    res_polyreg<-tryCatch(
      {
        PolyRegMix(data, des_matrix, degree_PolyReg, n_cores, test_type, mix_degree, iter_polyDegree)
      },
      error=function(errorCond){
        message("Execution halted with the following error:")
        message(errorCond)
        return(NULL)
      }
    ) #end trycatch PolyReg
  } #If a too complex model (even linear mixed effects model in a paired-wise test) is used, lme will throw a warning like: In logLik.lmeStructInt(lmeSt, lmePars) : Singular precision matrix in level -1, block 1
  #Should we react somehow to this warning, or is it the user's responsibility? Given that many times in these cases, a model is fit and results are given, even if a singular fit could not be determined? CHECK!

  #If a problem arises during PolyReg, halt.
  if(is.null(res_polyreg)){stop(cat("\nError during PolyReg\n"))}

  polyreg_pval<-res_polyreg[[1]]
  res_polyreg<-res_polyreg[[2]]

  #Now we should have everything to calculate the final rank product
  message("Calculating the final rank product.")
  rank_products<-tryCatch(
    {
      calculateRP(res_regrots, res_diffrots, res_polyreg, na_treat)
    },
    error=function(errorCond){
      message("Execution halted with the following error:")
      message(errorCond)
      return(NULL)
    }
  ) #end trycatch rank products

  if(is.null(rank_products)){stop(cat("\nError during rank product calculation\n"))}

  #Estimate signigicance values
  if(estSigVal){
    message("Estimating significance values.")
    if(weigh_runs){
      rp_with_est_sig<-tryCatch(
        {
          estSigValNonAlignedWeights(data, res_regrots, res_diffrots, res_polyreg, reg_rots_pval, diff_rots_pval, polyreg_pval, rank_products, na_treat, rand_seed, n_cores, sig_adj_meth, sigValSampN, regrots_weigths)
        },
        error=function(errorCond){
          message("Execution halted with the following error:")
          message(errorCond)
          return(NULL)
        }
      ) #end trycatch estimate significance values
    }else{
      rp_with_est_sig<-tryCatch(
        {
          estSigValNonAligned(data, res_regrots, res_diffrots, res_polyreg, reg_rots_pval, diff_rots_pval, polyreg_pval, rank_products, na_treat, rand_seed, n_cores, sig_adj_meth, sigValSampN)
        },
        error=function(errorCond){
          message("Execution halted with the following error:")
          message(errorCond)
          return(NULL)
        }
      ) #end trycatch estimate significance values
    }
    fin_res<-rp_with_est_sig
    if(is.null(fin_res)){stop(cat("\nError during significance value estimation\n"))}
  } else {
    fin_res<-rank_products
  }

  ret_list<-list(fin_res, res_regrots, reg_rots_pval, res_diffrots, diff_rots_pval, res_polyreg, polyreg_pval, rots_runs, method_degrees, input_list)
  names(ret_list)<-c("RolDE_Results", "RegROTS_Results", "RegROTS_P_Values", "DiffROTS_Results", "DiffROTS_P_Values", "PolyReg_Results", "PolyReg_P_Values", "ROTS_Runs", "Method_Degrees", "Input")
  return(ret_list)
}
