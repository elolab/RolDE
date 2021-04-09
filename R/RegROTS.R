RegROTS<-function(data, des_matrix, min_feat_obs, degree_RegROTS, rots_runs, B_for_ROTS, K_for_ROTS, rand_seed, na_treat, n_cores, add_refGroup, weigh_runs){

  degree<-degree_RegROTS
  unique_conditions<-unique(as.character(des_matrix[,2]))

  #Determine case and control. Control is always the first condition encountered in the design matrix.
  #Doesn't really matter, as there is no fold change.
  control<-unique_conditions[1]
  case<-unique_conditions[2]

  #The main function to compare coefficients
  getCoef<-function(r1, r2, time1, time2, degree, thres_feat_both_cond, min_feat_obs) { #a function to fit two models for strains separately and get the differences of coefficients

    if( sum(!is.na(r1))<min_feat_obs | sum(!is.na(r2))<min_feat_obs | sum(is.na(c(r1,r2)))>thres_feat_both_cond ){
      diffs<-rep(NA, degree+1)
    } else {

      if(any(is.na(r1))){
        time1<-time1[-which(is.na(r1))]
        r1<-r1[-which(is.na(r1))]
      }

      if(any(is.na(r2))){
        time2<-time2[-which(is.na(r2))]
        r2<-r2[-which(is.na(r2))]
      }

      if((length(time1)-1)>=degree){d1<-degree}else{d1<-length(time1)-1}
      if((length(time2)-1)>=degree){d2<-degree}else{d2<-length(time2)-1}

      mod1<-tryCatch({
        lm(as.numeric(r1)~poly(time1,(d1)))
      },error = function(e) {
        NULL
      })

      mod2<-tryCatch({
        lm(as.numeric(r2)~poly(time2,(d2)))
      },error = function(e) {
        NULL
      })

      if(is.null(mod1) | is.null(mod2)){
        diffs<-rep(NA, degree+1)
      } else {
        coeffs1<-mod1$coefficients
        coeffs2<-mod2$coefficients
        diffs<-numeric(degree+1)
        diffs[1:(min(d1,d2, na.rm = T)+1)]<-coeffs1[1:(min(d1,d2, na.rm = T)+1)]-coeffs2[1:(min(d1,d2, na.rm = T)+1)]
        diffs[diffs==0]<-NA
      }
    }
    return(diffs)
  }
  run<-NULL
  rots_res_frame=matrix(nrow=nrow(data), ncol = length(rots_runs))

  #Parallel run - always, in case of sequential, the backend will be registered only with 1 thread (n_cores=1)
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  rots_res_frame <- foreach::foreach(run=1:length(rots_runs), .export = c("fillGaps_New", "fillGapsAll_New"), .packages = c("ROTS", "doParallel", "matrixStats", "foreach"), .combine=cbind) %dopar% {
    #Determine groups for ROTS
    run_comparisons<-rots_runs[[run]]
    replicate_comparisons<-ncol(run_comparisons)
    groups_for_rots<-sort(rep(seq(1:(degree+1)),replicate_comparisons))

    #Get the sample locations (columns) for the individuals in the comparisons for the current run.
    control_locs<-list()
    case_locs<-list()

    #Store the timepoints for each individual in each comparison here
    times_control<-list()
    times_case<-list()

    for(comp in 1:ncol(run_comparisons)){ #Case stays always the same, control changes. Shouldn't make a large difference?
      control_locs[[comp]]<-which(des_matrix$Individual%in%run_comparisons[1,comp])
      case_locs[[comp]]<-which(des_matrix$Individual%in%run_comparisons[2,comp])

      times_control[[comp]]=des_matrix$Timepoint[control_locs[[comp]]]
      times_case[[comp]]=des_matrix$Timepoint[case_locs[[comp]]]
    }

    #A data frame to save the coefficient differences in
    coef_frame<-data.frame(matrix(nrow = nrow(data), ncol = replicate_comparisons*(degree+1)))
    rownames(coef_frame)<-rownames(data)
    colnames(coef_frame)<-groups_for_rots

    for(r in 1:nrow(data)) {

      res_row<-numeric(ncol(coef_frame))

      for(comp in 1:ncol(run_comparisons)){ #is there a better way than another loop?

        cont_loc<-control_locs[[comp]]
        case_loc<-case_locs[[comp]]

        row_control<-data[r,cont_loc]
        row_case<-data[r,case_loc]

        time_control<-times_control[[comp]]
        time_case<-times_case[[comp]]

        thres_feat_both_cond<-(length(time_control)+length(time_case))-4 #Interna threshold. Should be able to be controlled or should be removed?

        res_vals<-getCoef(r1 = row_case, r2 = row_control, time1 = time_case, time2 = time_control, degree = degree, thres_feat_both_cond = thres_feat_both_cond, min_feat_obs = min_feat_obs)
        res_row[seq(from=comp, to=length(res_row), by=replicate_comparisons)]=res_vals

      } #end comparison loop
      coef_frame[r,]<-res_row
    } #end getting coefficients for rows

    #Scale the coefficients differences within each group to acquire similar distributions for correct testing in ROTS
    coef_scaled<-data.frame(matrix(nrow = nrow(coef_frame), ncol = ncol(coef_frame)))
    rownames(coef_scaled)<-rownames(coef_frame)
    colnames(coef_scaled)<-colnames(coef_frame)

    for(group in 1:max(groups_for_rots)){
      data_temp<-coef_frame[, which(groups_for_rots==group)]
      data_temp<-as.numeric(unlist(data_temp))
      data_temp<-data_temp/sd(data_temp, na.rm = T)
      coef_scaled[,which(groups_for_rots==group)]<-data_temp
    }

    #Remove features with noth enough values in the first two groups, intercept (expression) level
    #differences and differences in longitudinal linear patterns of expression. Demand at least those.
    #Require at least two values for group, 2 required for ROTS.
    rem_feat<-numeric(0)
    for(group in 1:2){ #so we need to remove those that have too many nas in group 1 or 2
      rem_feat<-c(rem_feat,which(apply(coef_scaled[,which(groups_for_rots==group)], 1, function(x) (sum(!is.na(x))<2))))
    }

    rem_feat<-unique(rem_feat)

    #Save and remove the possible removed features
    if(length(rem_feat)>0){
      rem_names<-rownames(coef_scaled)[rem_feat] #save these and input back into the results
      coef_scaled<-coef_scaled[-rem_feat,]
    }

    #Make sure that each group has enough values for ROTS. While each protein is required
    #to have enough numerical difference values in the first two groups in the prior step,
    #the rest of the group might have insufficient number of values. Furthermore,
    #each feature is required to have the same number of groups. So in each group
    #for each feature, if not at least two difference values per group exist,
    #some random values are imputed from the total distribution of differences over
    #all features and groups. Since most of the features are not typically
    #differentially expressed, most of these values are near zero. Now uses all values,
    #could also be values between quantiles 0.25-0.75, so more pronounced zero values. Change if needed?
    all_num_values<-as.numeric(unlist(coef_scaled))
    all_num_values<-na.omit(all_num_values)

    #get enough quantile values used for imputation
    all_quant_vals<-quantile(all_num_values, seq(from=0.001, to=1, by=0.001))

    #or more tightly around the zero. Will bias the results?
    #all_quant_vals<-quantile(all_num_values, seq(from=0.2, to=0.8, by=0.0005))

    #make the imputation repeatable
    for(r in 1:nrow(coef_scaled)){
      row<-coef_scaled[r,]
      seed_for_imp<-rand_seed+r
      coef_scaled[r,]<-fillGaps_New(all_quant_vals = all_quant_vals, row = row, groups_for_rots = groups_for_rots, rand_seed = seed_for_imp)
    }

    #non-repeatble, but slightly faster, without seed.
    #coef_scaled<-t(apply(coef_scaled, 1, function (x) fillGaps_New(all_quant_vals = all_quant_vals, row = x, groups_for_rots = groups_for_rots)))

    #Add a zeroish reference group to compare to.
    if(add_refGroup){
      ref_group<-data.frame(matrix(nrow = nrow(coef_scaled), ncol = length(which(groups_for_rots==1))))

      #Fill compeletely in with random values from the distribution of all differences
      for(r in 1:nrow(ref_group)){
        row<-ref_group[r,]
        seed_for_imp<-rand_seed+r
        ref_group[r,]<-fillGapsAll_New(all_quant_vals = all_quant_vals, row = row, rand_seed = seed_for_imp)
      }

      #Non-repeatble but slightly faster
      #ref_group=t(apply(ref_group, 1, function (x) fillGapsAll_New(all_quant_vals = all_quant_vals, row = x)))

      #Make a new dataset and update groups for ROTS
      coef_new<-cbind(coef_scaled, ref_group)
      groups_for_rots<-c(groups_for_rots, rep((max(groups_for_rots)+1), length(which(groups_for_rots==1))))
      colnames(coef_new)<-groups_for_rots
    }else{
      coef_new<-coef_scaled
      colnames(coef_new)<-groups_for_rots
    }

    #Perform multigroup ROTS

    #determine K
    if(K_for_ROTS=="auto"){
      K_for_ROTS<-nrow(coef_new)/4
    }
    suppressMessages(expr = {rots_out<-ROTS::ROTS(data = coef_new, groups = groups_for_rots, B = B_for_ROTS, K = K_for_ROTS, paired = F, seed = rand_seed, progress = F)})
    rots_frame<-data.frame(d=rots_out$d, p=rots_out$pvalue)

    res_mat<-matrix(nrow = nrow(data), ncol = 1)
    rownames(res_mat)<-rownames(data)

    res_mat[match(rownames(rots_frame), rownames(res_mat)),1]<-rots_frame$p
    res_mat
  } #end one ROTS run
  parallel::stopCluster(cl)
  rownames(rots_res_frame)<-rownames(data)

  #Calculate RegROTS rank product from the different ROTS runs
  #Allow control for ties.method? Currently no.
  if(na_treat=="last"){
    ranks<-apply(rots_res_frame, 2, function(x) rank(x, na.last = T, ties.method = "average"))
  }else if(na_treat=="keep"){
    ranks<-apply(rots_res_frame, 2, function(x) rank(x, na.last = "keep", ties.method = "average"))
  }

  if(weigh_runs){
    nr_comps<-as.numeric(unlist(lapply(rots_runs, ncol)))
    max_comps<-max(nr_comps)
    weights_runs<-nr_comps/max_comps
    rank_prods<-apply(ranks, 1, function(x) exp(weighted.mean(log(x), weights_runs, na.rm=T)))
  }else{
    rank_prods<-apply(ranks, 1, function(x) exp(mean(log(x), na.rm=T)))
  }

  names(rank_prods)<-rownames(rots_res_frame)

  fin_res<-cbind(id=names(rank_prods), rp=rank_prods)
  if(any(is.nan(fin_res[,2]))){fin_res[which(is.nan(fin_res[,2])),2]=NA}
  if(any(fin_res[,2]=="NaN", na.rm = T)){fin_res[which(fin_res[,2]=="NaN"),2]=NA}
  fin_res<-data.frame(fin_res, stringsAsFactors = F)
  fin_res[,1]<-as.character(fin_res[,1])
  fin_res[,2]<-as.numeric(as.character(fin_res[,2]))

  if(all(is.na(fin_res[,2]))){stop("Unkown error during RegROTS.")}

  #make the return list
  if(weigh_runs){
    ret_list<-list(rots_res_frame, fin_res, weights_runs)
    names(ret_list)<-c("p-values", "rank products", "weights for runs")
  }else{
    ret_list<-list(rots_res_frame, fin_res)
    names(ret_list)<-c("p-values", "rank products")
  }

  return(ret_list)
}
