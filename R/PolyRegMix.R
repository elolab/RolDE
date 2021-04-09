PolyRegMix<-function(data, des_matrix, degree_PolyReg, n_cores, test_type, mix_degree, iter_polyDegree){

  message("Running PolyReg with mixed effects. If the results seem poor, try lower-level mixed effects or fixed only effect models. See vignette.")
  message("Especially missing values might pose a problem for the mixed-effects models.")

  unique_conditions<-unique(as.character(des_matrix[,2]))

  #Determine case and control. Control is always the first condition encountered in the design matrix. Need to be changed?
  #Doesn't really matter, as there is no fold change.
  control<-unique_conditions[1]
  case<-unique_conditions[2]

  nr_timepoints<-max(unique(as.numeric(as.character(des_matrix$Timepoint))))

  res_frame<-matrix(nrow=nrow(data), ncol=(3+degree_PolyReg))
  r<-NULL

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  res_frame <- foreach::foreach(r=1:nrow(data), .export = c("degree_PolyReg"), .packages = c("nlme", "foreach"), .combine=rbind) %dopar% {

    temp_data<-matrix(nrow = nrow(des_matrix), ncol = 5)
    temp_data<-data.frame(temp_data, stringsAsFactors = F)

    colnames(temp_data)<-c("Intensity", "Time", "Condition", "Block", "Individual")
    temp_data[,1]<-as.numeric(as.character(data[r,]))
    temp_data[,2]<-as.numeric(as.character(des_matrix$Timepoint))
    temp_data[,3]<-as.factor(des_matrix$Condition)
    temp_data[,4]<-as.factor(as.numeric(as.character(des_matrix$Block)))
    temp_data[,5]<-as.factor(as.numeric(as.character(des_matrix$Individual)))
    rownames(temp_data)<-colnames(data)


    if(mix_degree==0){ #Baseline or intercept level mixed effects
      # message("Running PolyReg with baseline mixed effects. If the results seem poor, try fixed only effects. See instructions.") #Remove or not?
      if(test_type=="unpaired"){
        mod1 = tryCatch({
          nlme::lme(Intensity ~ poly(Time, degree = get("degree_PolyReg"))*Condition, random = ~ 1 | Individual,
                    data = temp_data, na.action = na.omit, method = "REML",control = lmeControl(opt = 'optim')) #optim chosen as the optimized, as based on experience, it seems to converge more often than the default.
        }, error = function(e) {
          NULL
        })

        if(is.null(mod1) && iter_polyDegree){
          degree_try<-degree_PolyReg-1
          while(degree_try>0 && is.null(mod1)){
            mod1<-tryCatch({
              nlme::lme(Intensity ~ poly(Time, degree = get("degree_try"))*Condition, random = ~ 1 | Individual,
                        data = temp_data, na.action = na.omit, method = "REML",control = lmeControl(opt = 'optim'))
            },error = function(e) {
              NULL
            })
            degree_try<-degree_try-1
          }
        }

        if(!is.null(mod1) && (length(grep("Condition", rownames(summary(mod1)$tTable)))>0)){
          sum<-summary(mod1)
          sum<-sum$tTable
          locs<-grep("Condition", rownames(sum))
          if(length(locs)==1){
            sum<-t(data.frame(sum[locs,], stringsAsFactors = F))
          } else {
            sum<-sum[locs,]
          }
          sum<-data.frame(sum, stringsAsFactors = F)
          all_p_vals<-as.numeric(sum$p.value)
          min_p<-min(as.numeric(sum$p.value), na.rm = T)
          if(is.infinite(min_p)){min_p<-NA}
          temp_p_vals<-rep(NA, (degree_PolyReg+1))
          temp_p_vals[1:length(all_p_vals)]<-all_p_vals
          all_p_vals<-temp_p_vals
          row_val<-c(rownames(data)[r], min_p, all_p_vals)
        } else {
          row_val<-c(rownames(data)[r], NA, rep(NA, (degree_PolyReg+1)))
        }
      }
    } else {
      #message("Running PolyReg with linear mixed effects. If the results seem poor, try lower-level mixed effects of fixed only effects. See instructions.")
      if(test_type=="unpaired"){
        mod1 = tryCatch({
          nlme::lme(Intensity ~ poly(Time, degree = get("degree_PolyReg"))*Condition, random = ~ Time | Individual,
                    data = temp_data, na.action = na.omit, method = "REML",control = lmeControl(opt = 'optim', msMaxIter=500)) #For more complex models, increased the maximum number of iterations for the model to converge. Default is 50.
        }, error = function(e) {
          NULL
        })

        if(is.null(mod1) && iter_polyDegree){
          degree_try<-degree_PolyReg-1
          while(degree_try>0 && is.null(mod1)){
            mod1<-tryCatch({
              nlme::lme(Intensity ~ poly(Time, degree = get("degree_try"))*Condition, random = ~ Time | Individual,
                        data = temp_data, na.action = na.omit, method = "REML",control = lmeControl(opt = 'optim', msMaxIter=500))
            },error = function(e) {
              NULL
            })
            degree_try<-degree_try-1
          }
        }

        if(!is.null(mod1) && (length(grep("Condition", rownames(summary(mod1)$tTable)))>0)){
          sum<-summary(mod1)
          sum<-sum$tTable
          locs<-grep("Condition", rownames(sum))
          if(length(locs)==1){
            sum<-t(data.frame(sum[locs,], stringsAsFactors = F))
          } else {
            sum<-sum[locs,]
          }
          sum<-data.frame(sum, stringsAsFactors = F)
          all_p_vals<-as.numeric(sum$p.value)
          min_p<-min(as.numeric(sum$p.value), na.rm = T)
          if(is.infinite(min_p)){min_p<-NA}
          temp_p_vals<-rep(NA, (degree_PolyReg+1))
          temp_p_vals[1:length(all_p_vals)]<-all_p_vals
          all_p_vals<-temp_p_vals
          row_val<-c(rownames(data)[r], min_p, all_p_vals)
        } else {
          row_val<-c(rownames(data)[r], NA, rep(NA, (degree_PolyReg+1)))
        }
      }
    }
    row_val
  }
  parallel::stopCluster(cl)

  all_cond_pvals<-res_frame[,(3:ncol(res_frame))]
  if(any(is.nan(all_cond_pvals))){all_cond_pvals[which(is.nan(all_cond_pvals))]=NA}
  if(any(all_cond_pvals=="NaN", na.rm = T)){all_cond_pvals[which(all_cond_pvals=="NaN")]=NA}

  rownames(all_cond_pvals)<-rownames(data)
  all_cond_pvals=data.frame(all_cond_pvals, stringsAsFactors = F)
  for(c in 1:ncol(all_cond_pvals)){all_cond_pvals[,c]<-as.numeric(as.character(all_cond_pvals[,c]))}
  colnames(all_cond_pvals)=c("intercept", paste("degree", seq(1:degree_PolyReg)))

  res_frame=res_frame[,c(1:2)]
  res_frame<-data.frame(res_frame, stringsAsFactors = F)
  rownames(res_frame)<-rownames(data)
  colnames(res_frame)<-c("id","rep p-value")
  res_frame[,1]<-as.character(res_frame[,1])
  res_frame[,2]<-as.numeric(as.character(res_frame[,2]))

  if(all(is.na(res_frame[,2]))){stop("Unkown error during Polyreg.")}

  ret_list=list(all_cond_pvals, res_frame)
  names(ret_list)<-c("all cond p-values", "rep p-values")

  return(ret_list)
}
