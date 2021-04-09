determineROTSRunsAligned<-function(data, des_matrix, test_type, min_comm_diff){

  unique_conditions<-unique(as.character(des_matrix[,2]))
  #rots_runs<-list()

  if(test_type=="unpaired"){
    inds_condition1<-unique(des_matrix$Individual[which(des_matrix$Condition==unique_conditions[1])])
    inds_condition2<-unique(des_matrix$Individual[which(des_matrix$Condition==unique_conditions[2])])

    #determine which individual in different conditions have at least min_comm_diff timepoints in their common range.
    combs=matrix(nrow = 0, ncol = 2)
    colnames(combs)=unique_conditions

    for(indc1 in 1:length(inds_condition1)){
      temp_rows1=which(des_matrix$Individual%in%inds_condition1[indc1])
      temp_time1=as.numeric(as.character(des_matrix$Timepoint[temp_rows1]))
      #now look for individuals with at least 3 timepoints during the same age time in the case group
      for(indc2 in 1:length(inds_condition2)){
        temp_rows2=which(des_matrix$Individual%in%inds_condition2[indc2])
        temp_time2=as.numeric(as.character(des_matrix$Timepoint[temp_rows2]))
        temp_common_time=getLengthCommonTimepoints(temp_time1 = temp_time1, temp_time2 = temp_time2)
        if(temp_common_time[1]>=min_comm_diff & temp_common_time[2]>=min_comm_diff){
          combs=rbind(combs, c(inds_condition1[indc1],inds_condition2[indc2]))
        }
      }
    }

    #determine the number of ROTS runs, the same as the numnber of individuals as in the condition with
    #more individuals
    nr_ROTS_runs<-max(length(inds_condition1), length(inds_condition2))

    if(length(inds_condition1)==length(inds_condition2)){
      min_cond<-1
    }else if(length(inds_condition1)<length(inds_condition2)){
      min_cond<-1
      diff_val<-length(inds_condition2)-length(inds_condition1)
      inds_condition1<-c(inds_condition1,rep(NA, diff_val))
      #inds_condition1<-c(inds_condition1,NA)
    }else{
      min_cond<-2
      diff_val<-length(inds_condition1)-length(inds_condition2)
      inds_condition2<-c(inds_condition2,rep(NA, diff_val))
      #inds_condition2<-c(inds_condition2,NA)
    }

    all_pos_combs<-matrix(nrow=2, ncol = 0)
    temp_mat<-rbind(inds_condition1, inds_condition2)
    for(run in 1:nr_ROTS_runs){
      row1<-temp_mat[1,]
      row2<-temp_mat[2,]
      locs_na<-which(is.na(temp_mat), arr.ind <- T)
      locs_na<-as.numeric(locs_na[,2])
      if(length(locs_na)>0){temp_mat<-temp_mat[,-c(locs_na)]}
      if(!is.matrix(temp_mat)){temp_mat=as.matrix(temp_mat)} #if we have only one comparison
      all_pos_combs<-cbind(all_pos_combs, temp_mat)

      if(min_cond==1){
        row1=shifter(x <- row1, n<-1)
      }else{
        row2<-shifter(x <- row2, n<-1)
      }
      temp_mat<-rbind(row1, row2)
    }

    #validate which of all the possible combinations are possible to be realized based on enough timepoints in common range
    all_combs_one<-character(ncol(all_pos_combs))
    for(co in 1:ncol(all_pos_combs)){all_combs_one[co]<-paste(all_pos_combs[,co], collapse = ",")}

    pos_combs_one<-character(nrow(combs))
    for(r in 1:nrow(combs)){pos_combs_one[r]<-paste(combs[r,], collapse = ",")}

    #realized combinations from all possible combinations
    real_combs<-character()
    for(com in 1:length(all_combs_one)){if(all_combs_one[com]%in%pos_combs_one){real_combs<-c(real_combs,all_combs_one[com])}}

    #divide into rots runs, make runs as equally sized as possible.
    rots_runs<-getUniqueCombRunsNew(real_combs = real_combs, unique_conditions = unique_conditions)

    col_check<-as.numeric(unlist(lapply(rots_runs, ncol)))
    if(any(col_check<2)){
      error_msg<-paste("Inbalanced block design for RegROTS. Unable to Divide the comparisons into ROTS runs. Please remove or average samples. See instructions.", sep="")
      stop(error_msg)
    }
  }
  #Have removed the old paired and block design compile runs. If needed, see older versions.
  check1<-any(unlist(lapply(rots_runs, function(x) any(duplicated(x[1,])))))
  check2<-any(unlist(lapply(rots_runs, function(x) any(duplicated(x[2,])))))

  if(check1 | check2){stop("Not unique combinations in ROTS runs.")}

  ret_list<-list(rots_runs)
  names(ret_list)<-c("ROTS runs")
  return(ret_list)
}
