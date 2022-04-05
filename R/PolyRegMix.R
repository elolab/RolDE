PolyRegMix<-function(data, des_matrix, degree_PolyReg, n_cores, model_type){

  unique_conditions<-unique(as.character(des_matrix[,2]))

  #Determine case and control. Control is always the first condition encountered in the design matrix. Need to be changed?
  control<-unique_conditions[1]
  case<-unique_conditions[2]

  nr_timepoints<-max(unique(as.numeric(as.character(des_matrix$Timepoint))))

  res_frame<-matrix(nrow=nrow(data), ncol=(3+degree_PolyReg))
  r<-NULL

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  res_frame <- foreach::foreach(r=seq_len(nrow(data)), .packages = c("nlme", "foreach"), .combine=rbind) %dorng% {

    temp_data<-matrix(nrow = nrow(des_matrix), ncol = 4)
    temp_data<-data.frame(temp_data, stringsAsFactors = FALSE)

    colnames(temp_data)<-c("Intensity", "Time", "Condition", "Individual")
    temp_data[,1]<-as.numeric(as.character(data[r,]))
    temp_data[,2]<-as.numeric(as.character(des_matrix$Timepoint))
    temp_data[,3]<-as.factor(des_matrix$Condition)
    temp_data[,4]<-as.factor(as.numeric(as.character(des_matrix$Individual)))
    rownames(temp_data)<-colnames(data)

    if(model_type=="mixed0"){ #Baseline or intercept level mixed effects
      mod1 <- tryCatch({
        nlme::lme(Intensity ~ poly(Time, degree = get("degree_PolyReg"))*Condition, random = ~ 1 | Individual,
                  data = temp_data, na.action = na.omit, method = "REML",control = lmeControl(opt = 'optim')) #optim chosen as the optimized, as based on experience, it seems to converge more often than the default.
      }, error = function(e) {
        NULL
      })

      if(!is.null(mod1) && (length(grep("Condition", rownames(summary(mod1)$tTable)))>0)){
        sum<-summary(mod1)
        sum<-sum$tTable
        locs<-grep("Condition", rownames(sum))
        if(length(locs)==1){
          sum<-t(data.frame(sum[locs,], stringsAsFactors = FALSE))
        } else {
          sum<-sum[locs,]
        }
        sum<-data.frame(sum, stringsAsFactors = FALSE)
        all_p_vals<-as.numeric(sum$p.value)
        min_p<-min(as.numeric(sum$p.value), na.rm = TRUE)
        if(is.infinite(min_p)){min_p<-NA}
        temp_p_vals<-rep(NA, (degree_PolyReg+1))
        temp_p_vals[seq_len(length(all_p_vals))]<-all_p_vals
        all_p_vals<-temp_p_vals
        row_val<-c(rownames(data)[r], min_p, all_p_vals)
      } else {
        row_val<-c(rownames(data)[r], NA, rep(NA, (degree_PolyReg+1)))
      }
    } else {
      mod1 <- tryCatch({
        nlme::lme(Intensity ~ poly(Time, degree = get("degree_PolyReg"))*Condition, random = ~ Time | Individual,
                  data = temp_data, na.action = na.omit, method = "REML",control = lmeControl(opt = 'optim', msMaxIter=500)) #For more complex models, increased the maximum number of iterations for the model to converge. Default is 50.
      }, error = function(e) {
        NULL
      })

      if(!is.null(mod1) && (length(grep("Condition", rownames(summary(mod1)$tTable)))>0)){
        sum<-summary(mod1)
        sum<-sum$tTable
        locs<-grep("Condition", rownames(sum))
        if(length(locs)==1){
          sum<-t(data.frame(sum[locs,], stringsAsFactors = FALSE))
        } else {
          sum<-sum[locs,]
        }
        sum<-data.frame(sum, stringsAsFactors = FALSE)
        all_p_vals<-as.numeric(sum$p.value)
        min_p<-min(as.numeric(sum$p.value), na.rm = TRUE)
        if(is.infinite(min_p)){min_p<-NA}
        temp_p_vals<-rep(NA, (degree_PolyReg+1))
        temp_p_vals[seq_len(length(all_p_vals))]<-all_p_vals
        all_p_vals<-temp_p_vals
        row_val<-c(rownames(data)[r], min_p, all_p_vals)
      } else {
        row_val<-c(rownames(data)[r], NA, rep(NA, (degree_PolyReg+1)))
      }
    }
    row_val
  }
  parallel::stopCluster(cl)

  all_cond_pvals<-res_frame[,c(3:ncol(res_frame))]
  if(any(is.nan(all_cond_pvals))){all_cond_pvals[which(is.nan(all_cond_pvals))]<-NA}
  if(any(all_cond_pvals=="NaN", na.rm = TRUE)){all_cond_pvals[which(all_cond_pvals=="NaN")]<-NA}

  rownames(all_cond_pvals)<-rownames(data)
  all_cond_pvals<-data.frame(all_cond_pvals, stringsAsFactors = FALSE)
  #for(c in seq_len(ncol(all_cond_pvals))){all_cond_pvals[,c]<-as.numeric(as.character(all_cond_pvals[,c]))}
  all_cond_pvals<-data.frame(apply(all_cond_pvals, 2, function(x) {as.numeric(as.character(x))}), check.names = FALSE)
  colnames(all_cond_pvals)<-c("intercept", paste("degree", seq_len(degree_PolyReg)))

  res_frame<-res_frame[,c(seq_len(2))]
  res_frame<-data.frame(res_frame, stringsAsFactors = FALSE)
  rownames(res_frame)<-rownames(data)
  colnames(res_frame)<-c("id","rep p-value")
  res_frame[,1]<-as.character(res_frame[,1])
  res_frame[,2]<-as.numeric(as.character(res_frame[,2]))

  if(all(is.na(res_frame[,2]))){stop("Unkown error during Polyreg.")}

  ret_list<-list(all_cond_pvals, res_frame)
  names(ret_list)<-c("all cond p-values", "rep p-values")

  return(ret_list)
}
