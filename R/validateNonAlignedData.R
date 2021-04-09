#Organize the aligned data properly for RolDE_NonAligned and run some validity checks.

validateNonAlignedData<-function(data, des_matrix, min_comm_diff, min_ind_timep, min_feat_obs, B_for_ROTS, K_for_ROTS, rand_seed, na_treat, n_cores, polyregType, mix_degree, estSigVal, sigValSampN, sig_adj_meth, add_refGroup, iter_polyDegree, weigh_runs){

  if(K_for_ROTS!="auto"){
    if(!is.numeric(K_for_ROTS)){stop("Incorrect K_for_ROTS. Numerical value needed. See intructions.")}
  }

  if(min_comm_diff!="auto"){
    if(!is.numeric(min_comm_diff)){stop("Incorrect min_comm_diff Numerical value needed. See intructions.")}
    if(length(min_comm_diff!=2)){stop("Incorrect length for parameter min_com_diff. Must be 2 numerical values. See instructions.")}
    if(min_comm_diff[1]<2){stop("Parameter min_comm_diff for RegROTS set to <2. Longitudinal analysis not meaningful. See instructions.")}
    if(min_comm_diff[2]<1){stop("Parameter min_comm_diff for DiffROTS set to <1. Longitudinal analysis not meaningful with no common time points. See instructions.")}
  }else{
    min_comm_diff<-c(3,1)
  }

  if(!is.numeric(B_for_ROTS)){stop("Incorrect B_for_ROTS. Numerical value needed.")}
  if(B_for_ROTS<50){message("Warning! Very few permutations for ROTS. A larger value may be desired. See Instructions.")}
  if(B_for_ROTS>5000){message("Warning! Very large number of permutations for ROTS. May take a long time to finish. A value of 1000 should be sufficient. See Instructions.")}

  if(!is.numeric(rand_seed)){stop("Incorrect Rand_Seed. Numerical value needed.")}

  if(!is.numeric(min_ind_timep)){stop("Incorrect min_ind_timep Numerical value needed.")}
  if(min_ind_timep<3){stop("Too small min_ind_timep. Longitudinal analysis not meaningul. See instructions. ")}

  if(!is.numeric(min_feat_obs)){stop("Incorrect min_feat_obs. Numerical value needed.")}

  nr_timepoints=max(as.numeric(unique(des_matrix[,3])))
  if(nr_timepoints<3){stop("Too few timepoints. Longitudinal analysis not meaningful. See Instructions.")}
  if(any(min_comm_diff>nr_timepoints)){stop("Parameter min_comm_diff too large, not enough time points. See instructions.")}

  if(polyregType!="fixed" & polyregType!="mixed"){
    stop("Incorrect test type for PolyReg. Correct inputs are fixed or mixed, either required. See instructions.")
  }

  if(polyregType=="mixed"){
    if(!is.numeric(mix_degree)){stop("Incorrect mix_degree. Numerical value needed.")}
    if(mix_degree!=0 & mix_degree!=1) {
      stop("Incorrect mix_degree for PolyReg. Only 0 (intercept/baseline) or 1 (linear) level mixed effects currently allowed. See Instructions.")
    }
  }

  if(na_treat!="last" & na_treat!="keep"){
    stop("Incorrect na_treat. Correct inputs are last or keep, either required. See instructions.")
  }

  if(!is.numeric(n_cores)){stop("Incorrect n_cores. Numerical value needed.")}

  if(!is.logical(estSigVal)){stop("Variable:estSigVal not logical. TRUE or FALSE needed.")}
  if(!is.logical(add_refGroup)){stop("Variable:add_refGroup not logical. TRUE or FALSE needed.")}
  if(!is.logical(iter_polyDegree)){stop("Variable:iter_polyDegree not logical. TRUE or FALSE needed.")}
  if(!is.logical(weigh_runs)){stop("Variable:weigh_runs not logical. TRUE or FALSE needed.")}

  if(!is.numeric(sigValSampN)){stop("Incorrect sigValSampN. Numerical value needed.")}
  if(sigValSampN>10000){message("Warning! Large sigValSampN. Significance value estimation may take a long time.")}
  if(sigValSampN<50){stop("Too small sigValSampN. Please increase the number. See instructions.")}

  if(sig_adj_meth!="bonferroni" & sig_adj_meth!="holm" & sig_adj_meth!="hochberg" & sig_adj_meth!="hommel" & sig_adj_meth!="BH" & sig_adj_meth!="fdr" & sig_adj_meth!="BY" & sig_adj_meth!="qvalue"){
    stop("Incorrect sig_adj_meth. Correct inputs are: bonferroni, holm, hochberg, hommel, BH, fdr, BY or qvalue. See instructions.")
  }

  #Validation checks for the design matrix
  #Right numner of columns
  if(ncol(des_matrix)!=4){
    stop("Wrong number of columns in the design matrix. Four columns needed, see instructions.")
  }

  #All samples accounted for
  if(!all(colnames(data)==as.character(des_matrix[,1]))){
    stop("Sample names not found. Sample names needed in the first column of the design matrix.")
  }

  #Number of unique conditions is correct
  if(length(unique(as.character(des_matrix[,2])))!=2){
    stop("Number of conditions wrong in the design matrix. Two conditions needed and allowed.")
  }

  unique_conditions=unique(as.character(des_matrix[,2]))

  if(!all(as.character(des_matrix[,2])%in%unique_conditions)){
    stop("Wrong/missing Condition for some samples in the design matrix. All samples must have a valid condition. See instructions.")
  }

  unique_timepoints<-unique(as.numeric(as.character(des_matrix[,3])))

  if(!all(as.numeric(as.character(des_matrix[,3]))%in%unique_timepoints)){
    stop("Wrong/missing Timepoint information for some samples in the design matrix. All samples must have a valid timepoint. See instructions.")
  }

  unique_individuals<-unique(as.numeric(as.character(des_matrix[,4])))

  if(!all(as.numeric(as.character(des_matrix[,4]))%in%unique_individuals)){
    stop("Wrong/missing individual ID for some samples in the design matrix. All samples must have a valid individual ID. See instructions.")
  }

  #check that each individual has enough timepoints
  individual_timepoints<-numeric(length(unique_individuals))
  for(ind in 1:length(unique_individuals)){
    ind_locs<-which(as.numeric(as.character(des_matrix[,4]))==unique_individuals[ind])
    individual_timepoints[ind]<-length(unique(as.numeric(as.character(des_matrix[ind_locs,3]))))
  }

  if(any(individual_timepoints<min_ind_timep)){
    stop("Some individual have too few timepoints in the design matrix. Please remove such individuals or change Variable:min_ind_timep. See instructions.")
  }

  #Determine that each condition has enough individuals. Enough timepoints per individual already checked.
  individuals_in_condition<-numeric(2)

  for(cond in 1:2){
    cond_locs<-which(as.character(des_matrix[,2])==unique_conditions[cond])
    individuals_in_condition[cond]<-length(unique(as.numeric(as.character(des_matrix[cond_locs,4]))))
  } #end for loop

  #Design matrix is good, proceed
  des_matrix=data.frame(des_matrix)
  #colnames(des_matrix)=c("Sample Names", "Condition", "Timepoint", "Individual", "Block")
  colnames(des_matrix)=c("Sample Names", "Condition", "Timepoint", "Individual")

  des_matrix[,1]<-as.character(des_matrix[,1])
  des_matrix[,2]<-as.factor(as.character(des_matrix[,2]))
  des_matrix[,3]<-as.numeric(as.character(des_matrix[,3]))
  des_matrix[,4]<-as.numeric(as.character(des_matrix[,4]))

  #Validate data
  if(!is.matrix(data)){
    stop("Data is not matrix. A numerical matrix required, see instuctions.")
  }

  if(!is.numeric(data)){
    stop("Data is not a numerical matrix. A numerical matrix required, see instuctions.")
  }

  test_type<-"unpaired"
  return_list<-list(data, des_matrix, test_type, min_comm_diff)
  names(return_list)<-c("Data", "Design_Matrix", "Test_Type", "Min_Com_Diff")
  return(return_list)
}
