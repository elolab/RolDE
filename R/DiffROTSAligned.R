DiffROTSAligned<-function(data, des_matrix, min_comm_diff, min_feat_obs, rots_runs, n_cores){

  unique_conditions<-unique(as.character(des_matrix[,2]))

  #Determine the minimum numnber of common timepoints required. In aligned data currently only 1.
  if(min_comm_diff=="auto"){min_comm_diff<-1}

  #Determine case and control. Control is always the first condition encountered in the design matrix. Need to be changed?
  control<-unique_conditions[1]
  case<-unique_conditions[2]

  time_points<-unique(as.numeric(as.character(des_matrix$Timepoint)))
  nr_timepoints<-max(unique(as.numeric(as.character(des_matrix$Timepoint))))

  #The main function to compare expression
  getDiff<-function(r1, r2, time1, time2, time_points){

    times<-numeric(length(time_points))
    names(times)<-as.character(time_points)
    names(r1)<-as.character(time1)
    names(r2)<-as.character(time2)

    if( sum(!is.na(r1))<min_feat_obs | sum(!is.na(r2))<min_feat_obs){
      diffs<-rep(NA, length(times))
    } else {
      int_times<-intersect(time1, time2)
      if(length(int_times)<min_comm_diff){ #do not compare
        diffs<-rep(NA, length(times))
      }else{
        times[as.character(int_times)]<-r1[as.character(int_times)]-r2[as.character(int_times)]
        times[times==0]<-NA
        diffs<-times
      }
    }
    return(diffs)
  } #End compare expression

  run<-NULL

  rots_res_frame<-matrix(nrow=nrow(data), ncol = length(rots_runs))
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  rots_res_frame <- foreach::foreach(run=seq_len(length(rots_runs)), .export = c("fillGaps_New", "fillGapsAll_New"), .packages = c("ROTS", "doParallel", "matrixStats", "foreach"), .combine=cbind) %dorng% {

    #Determine groups for ROTS
    run_comparisons<-rots_runs[[run]]
    replicate_comparisons<-ncol(run_comparisons)
    groups_for_rots<-sort(rep(seq_len(nr_timepoints),replicate_comparisons))

    #Get the sample locations (columns) for the individuals in the comparisons for the current run.
    control_locs<-as.list(as.data.frame(apply(run_comparisons, 2, function(x){which(des_matrix$Individual%in%x[1])})))
    case_locs<-as.list(as.data.frame(apply(run_comparisons, 2, function(x){which(des_matrix$Individual%in%x[2])})))

    #Store the timepoints for each individual in each comparison here
    times_control<-lapply(control_locs, function(x) {des_matrix$Timepoint[x]})
    times_case<-lapply(case_locs, function(x) {des_matrix$Timepoint[x]})

    #A dataframe to store the expression level differences in.
    diff_frame<-data.frame(matrix(nrow = nrow(data), ncol = replicate_comparisons*nr_timepoints))
    rownames(diff_frame)<-rownames(data)
    colnames(diff_frame)<-groups_for_rots

    for(r in seq_len(nrow(data))) {

      res_row<-numeric(ncol(diff_frame))

      for(comp in seq_len(ncol(run_comparisons))){

        cont_loc<-control_locs[[comp]]
        case_loc<-case_locs[[comp]]

        row_control<-data[r,cont_loc]
        row_case<-data[r,case_loc]

        time_control<-times_control[[comp]]
        time_case<-times_case[[comp]]

        res_vals<-getDiff(r1 = row_case, r2 = row_control, time1 = time_case, time2 = time_control, time_points = time_points)
        res_row[seq(from=comp, to=length(res_row), by=replicate_comparisons)]<-res_vals

      } #end comparison loop
      diff_frame[r,]<-res_row
    }

    #Check how many non-missing difference values each feature has in each group.
    na_values<-t(apply(diff_frame, 1, function(x){
      unlist(lapply(unique(groups_for_rots), function(z){
        vals<-as.numeric(x[which(groups_for_rots==z)])
        if(length(which(!is.na(vals)))>=2){1}else{0}
      }))
    }))

    rem_feat<-which(rowSums(na_values)<1) #Check if a feature has not a single valid group for ROTS. Then remove it.

    #Remove such features
    if(length(rem_feat)>0){
      rem_names<-rownames(diff_frame)[rem_feat] #save these and input back into the results
      diff_frame<-diff_frame[-rem_feat,]
    }


    #Each feature is required to have the same number of groups. So in each group
    #for each feature, if not at least two difference values per group exist,
    #some random values are imputed from the total distribution of differences over
    #all features and groups. Since most of the features are not typically
    #differentially expressed, most of these values are near zero.
    all_num_values<-as.numeric(unlist(diff_frame))
    all_num_values<-na.omit(all_num_values)

    #get enough quantile values used for imputation
    all_quant_vals<-quantile(all_num_values, seq(from=0.001, to=1, by=0.001))

    diff_frame<-data.frame(t(apply(diff_frame, 1, function(x) {fillGaps_New(all_quant_vals = all_quant_vals, row = x, groups_for_rots = groups_for_rots)})), check.names = FALSE)

    #Add a zeroish reference group to compare to.
    ref_group<-data.frame(matrix(nrow = nrow(diff_frame), ncol = length(which(groups_for_rots==1))))

    #Fill compeletely in with random values from the distribution of all differences
    ref_group<-data.frame(t(apply(ref_group, 1, function(x) {fillGapsAll_New(all_quant_vals = all_quant_vals, row = x)})), check.names = FALSE)

    #Make a new dataset and update groups for ROTS
    diff_new<-cbind(diff_frame, ref_group)
    groups_for_rots<-c(groups_for_rots, rep((max(groups_for_rots)+1), length(which(groups_for_rots==1))))
    colnames(diff_new)<-groups_for_rots

    #Perform multigroup ROTS
    suppressMessages(expr = {rots_out<-ROTS::ROTS(data = diff_new, groups = groups_for_rots, B = 100, K = nrow(diff_new)/4, paired = FALSE, progress = FALSE)})
    rots_frame<-data.frame(d=rots_out$d, p=rots_out$pvalue)

    res_mat<-matrix(nrow = nrow(data), ncol = 1)
    rownames(res_mat)<-rownames(data)

    res_mat[match(rownames(rots_frame), rownames(res_mat)),1]<-rots_frame$p
    res_mat
  } #end one ROTS run
  parallel::stopCluster(cl)
  rownames(rots_res_frame)<-rownames(data)

  ranks<-apply(rots_res_frame, 2, function(x) rank(x, na.last = TRUE, ties.method = "average"))

  nr_comps<-as.numeric(unlist(lapply(rots_runs, ncol)))
  max_comps<-max(nr_comps)
  weights_runs<-nr_comps/max_comps
  rank_prods<-apply(ranks, 1, function(x) exp(weighted.mean(log(x), weights_runs, na.rm=TRUE)))

  names(rank_prods)<-rownames(rots_res_frame)

  fin_res<-cbind(id=names(rank_prods), rp=rank_prods)
  if(any(is.nan(fin_res[,2]))){fin_res[which(is.nan(fin_res[,2])),2]<-NA}
  if(any(fin_res[,2]=="NaN", na.rm = TRUE)){fin_res[which(fin_res[,2]=="NaN"),2]<-NA}
  fin_res<-data.frame(fin_res, stringsAsFactors = FALSE)
  fin_res[,1]<-as.character(fin_res[,1])
  fin_res[,2]<-as.numeric(as.character(fin_res[,2]))

  if(all(is.na(fin_res[,2]))){stop("Unkown error during DiffROTS.")}

  #make the return list
  ret_list<-list(rots_res_frame, fin_res, weights_runs)
  names(ret_list)<-c("p-values", "rank products", "weights for runs")

  return(ret_list)
}
