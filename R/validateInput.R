#Organize the data properly and run some validity checks.
validateInput<-function(data, des_matrix, aligned, min_comm_diff, min_feat_obs, n_cores, model_type, sigValSampN, sig_adj_meth){

  if(!is.logical(aligned)){stop("Variable:aligned not logical. TRUE or FALSE needed.")}

  if(length(min_comm_diff)==2){
    if(!is.numeric(min_comm_diff)){stop("Incorrect min_comm_diff. Numerical values (two) needed. See documentation.")}
    if(min_comm_diff[1]<2){stop("Parameter min_comm_diff for RegROTS set to <2. Longitudinal analysis not meaningful. See documentation.")}
    if(min_comm_diff[2]<1){stop("Parameter min_comm_diff for DiffROTS set to <1. Longitudinal analysis not meaningful with no common time points. See documentation.")}
  }else{
    if(min_comm_diff!="auto"){
      stop("Incorrect min_comm_diff. Must be 2 numerical values or auto. See documentation.")
    }else{
      min_comm_diff<-c(3,1)
    }
  }

  if(!is.numeric(min_feat_obs)){stop("Incorrect min_feat_obs. Numerical value needed.")}

  if(model_type!="auto" & model_type!="fixed" & model_type!="mixed0" & model_type!="mixed1"){
    stop("Incorrect value for model_type. Allowed values for the parameter are: auto, fixed, mixed0 or mixed1.")
  }
  if(model_type=="auto"){
    if(aligned){
      model_type<-"fixed"
    }else{
      model_type<-"mixed0"
    }
  }

  if(!is.numeric(n_cores)){stop("Incorrect n_cores. Numerical value needed.")}

  if(!is.numeric(sigValSampN)){stop("Incorrect sigValSampN. Numerical value needed.")}
  if(sigValSampN>10000000){message("Caution! Large sigValSampN. Significance value estimation may take a long time.")}

  if(sigValSampN==0){
    estSigVal<-FALSE
  }else{
    if(sigValSampN<100000){stop("Too small sigValSampN. Please increase the number.")}
    estSigVal<-TRUE
  }

  if(sig_adj_meth!="bonferroni" & sig_adj_meth!="holm" & sig_adj_meth!="hochberg" & sig_adj_meth!="hommel" & sig_adj_meth!="BH" & sig_adj_meth!="fdr" & sig_adj_meth!="BY" & sig_adj_meth!="qvalue"){
    stop("Incorrect sig_adj_meth. Correct inputs are: bonferroni, holm, hochberg, hommel, BH, fdr, BY or qvalue. See documentation")
  }

  #Validation checks for the design matrix
  #Right numner of columns
  if(ncol(des_matrix)!=4){
    stop("Wrong number of columns in the design matrix. Four columns needed, see documentation.")
  }

  #All samples accounted for
  if(!all(colnames(data)==as.character(des_matrix[,1]))){
    stop("Sample names not found. Sample names needed in the first column of the design matrix.")
  }

  #Number of unique conditions is correct
  if(length(unique(as.character(des_matrix[,2])))!=2){
    stop("Number of conditions wrong in the design matrix. Two conditions needed and allowed.")
  }

  unique_conditions<-unique(as.character(des_matrix[,2]))

  if(!all(as.character(des_matrix[,2])%in%unique_conditions)){
    stop("Wrong/missing Condition for some samples in the design matrix. All samples must have valid condition information.")
  }

  unique_timepoints<-unique(as.numeric(as.character(des_matrix[,3])))

  if(!all(as.numeric(as.character(des_matrix[,3]))%in%unique_timepoints)){
    stop("Wrong/missing Timepoint information for some samples in the design matrix. All samples must have valid timepoint information.")
  }

  unique_individuals<-unique(as.numeric(as.character(des_matrix[,4])))

  if(!all(as.numeric(as.character(des_matrix[,4]))%in%unique_individuals)){
    stop("Wrong/missing individual ID for some samples in the design matrix. All samples must have valid individual IDs.")
  }

  #check that each individual has enough timepoints
  individual_timepoints<-unlist(lapply(unique_individuals, function(x) {
    ind_locs<-which(as.numeric(as.character(des_matrix[,4]))==x)
    length(unique(as.numeric(as.character(des_matrix[ind_locs,3]))))
  }))

  if(any(individual_timepoints<3)){
    stop("Some individual have less than three time points in the design matrix. Please remove such individuals, longitudinal analysis not meaningful.")
  }

  #Determine that each condition has enough individuals. Enough timepoints per individual already checked.
  individuals_in_condition<-unlist(lapply(unique_conditions, function(x) {
    cond_locs<-which(as.character(des_matrix[,2])==x)
    length(unique(as.numeric(as.character(des_matrix[cond_locs,4]))))
  }))

  if(any(individuals_in_condition<2)){
    stop("Either condition does not have at least 2 replicates (individuals). Differential expression analysis not meaningful.")
  }

  #Validate that if timepoints should be aligned, they truly are
  if(aligned){
    if(!all(sort(unique(des_matrix[which(des_matrix[,2]%in%unique_conditions[1]),3]))==
            sort(unique(des_matrix[which(des_matrix[,2]%in%unique_conditions[2]),3])))){
      stop("Timepoints indicated as aligned but different timepoints in conditions")
    }
  }

  #Design matrix is good, proceed
  des_matrix<-data.frame(des_matrix, stringsAsFactors = FALSE)
  colnames(des_matrix)<-c("Sample Names", "Condition", "Timepoint", "Individual")
  des_matrix[,1]<-as.character(des_matrix[,1])
  des_matrix[,2]<-as.factor(as.character(des_matrix[,2]))
  des_matrix[,3]<-as.numeric(as.character(des_matrix[,3]))
  des_matrix[,4]<-as.numeric(as.character(des_matrix[,4]))

  #Validate data
  if(!is.matrix(data)){
    stop("Data is not matrix. A numerical matrix required, see documentation.")
  }

  if(!is.numeric(data)){
    stop("Data is not a numerical matrix. A numerical matrix required, see documentation.")
  }
  return_list<-list(data, des_matrix, min_comm_diff, model_type, estSigVal)
  names(return_list)<-c("Data", "Design_Matrix", "Min_Comm_Diff", "Model_Type", "EstSigVal")
  return(return_list)
}
