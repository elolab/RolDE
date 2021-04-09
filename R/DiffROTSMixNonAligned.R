DiffROTSMixNonAligned<-function(data, des_matrix, degree_PolyReg, B_for_ROTS, K_for_ROTS, rand_seed, na_treat, doPar, n_cores, test_type, mix_degree){

  message("Running DiffROTS with mixed effects. If the results seem poor, try lower-level mixed effects or models with fixed only effects. See vignette.")
  message("Especially missing values might pose a problem for the mixed-effects models.")

  unique_conditions<-unique(as.character(des_matrix[,2]))

  #Determine case and control. Control is always the first condition encountered in the design matrix. Need to be changed?
  #Doesn't really matter, as there is no fold change.
  control<-unique_conditions[1]
  case<-unique_conditions[2]

  rots_frame=matrix(nrow = nrow(data), ncol = degree_PolyReg+1)
  d<-NULL
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  rots_frame <- foreach::foreach(d=0:degree_PolyReg, .packages = c("ROTS", "nlme", "foreach"), .export = ("degree_PolyReg"), .combine=cbind) %dopar% {

    resid_frame<-matrix(nrow = nrow(data), ncol = ncol(data))
    rownames(resid_frame)<-rownames(data)
    colnames(resid_frame)<-colnames(data)

    for(r in 1:nrow(data)){

      temp_data<-matrix(nrow = nrow(des_matrix), ncol = 5)
      temp_data<-data.frame(temp_data, stringsAsFactors = F)

      colnames(temp_data)<-c("Intensity", "Time", "Condition", "Block", "Individual")
      temp_data[,1]<-as.numeric(as.character(data[r,]))
      temp_data[,2]<-as.numeric(as.character(des_matrix$Timepoint))
      temp_data[,3]<-as.factor(des_matrix$Condition)
      temp_data[,4]<-as.factor(as.numeric(as.character(des_matrix$Block)))
      temp_data[,5]<-as.factor(as.numeric(as.character(des_matrix$Individual)))
      rownames(temp_data)<-colnames(data)

      if(mix_degree==0){
        if(test_type=="unpaired"){
          if(d==0){ #intercept only model
            mod1<-tryCatch({
              nlme::lme(Intensity~1, random = ~ 1 | Individual, data=temp_data, na.action = na.omit, method = "REML",control = lmeControl(opt = 'optim'))
            },error = function(e) {
              NULL
            })
          }else if(d==1){
            mod1<-tryCatch({
              nlme::lme(Intensity~Time, random = ~ 1 | Individual, data=temp_data, na.action = na.omit, method = "REML",control = lmeControl(opt = 'optim'))
            },error = function(e) {
              NULL
            })
          }else {
            mod1<-tryCatch({
              nlme::lme(Intensity~poly(Time, degree = get("degree_PolyReg")), random = ~ 1 | Individual, data=temp_data, na.action = na.omit, method = "REML",control = lmeControl(opt = 'optim'))
            },error = function(e) {
              NULL
            })
          }
        }
      }else{ #mix_degree = 1
        if(test_type=="unpaired"){
          if(d==0){ #intercept only model
            mod1<-tryCatch({
              nlme::lme(Intensity~1, random = ~ Time | Individual, data=temp_data, na.action = na.omit, method = "REML",control = lmeControl(opt = 'optim', msMaxIter=500))
            },error = function(e) {
              NULL
            })
          }else if(d==1){
            mod1<-tryCatch({
              nlme::lme(Intensity~Time, random = ~ Time | Individual, data=temp_data, na.action = na.omit, method = "REML",control = lmeControl(opt = 'optim', msMaxIter=500))
            },error = function(e) {
              NULL
            })
          }else {
            mod1<-tryCatch({
              nlme::lme(Intensity~poly(Time, degree = get("degree_PolyReg")), random = ~ Time | Individual, data=temp_data, na.action = na.omit, method = "REML",control = lmeControl(opt = 'optim', msMaxIter=500))
            },error = function(e) {
              NULL
            })
          }
        }
      }
      if(!is.null(mod1)){
        resid=mod1$residuals[,1] #use fixed effects residuals, mixed effects residuals don't seem to work good. And random effects used already in the estimation of the fixed effects as well.
        resid_frame[r,match(names(resid), colnames(resid_frame))]=resid
      }
    } #end for loop get resid frame
    ###DO ROTS
    groups_for_rots<-rep(0, ncol(resid_frame))
    groups_for_rots[which(des_matrix$Condition==case)]<-1

    ###remove features with too few observations
    nas1=apply(resid_frame[,which(groups_for_rots==0)], 1, function(x) length(which(!is.na(x))))
    nas2=apply(resid_frame[,which(groups_for_rots==1)], 1, function(x) length(which(!is.na(x))))
    rems1=which(nas1<2) #Should this be larger than 2? Seems pretty tolerant.
    rems2=which(nas2<2)
    rems=c(rems1, rems2)
    rems=unique(rems)
    if(length(rems)>0){
      resid_frame=resid_frame[-rems,]
    }

    if(K_for_ROTS=="auto"){
      used_K_for_ROTS<-nrow(resid_frame)/4
    } else{
      used_K_for_ROTS<-K_for_ROTS
    }

    suppressMessages(expr <- {rots_out<-ROTS::ROTS(data = resid_frame, groups = groups_for_rots, B = B_for_ROTS, K = used_K_for_ROTS, paired = F, seed = rand_seed, progress = F)})
    rots_res<-data.frame(d=rots_out$d, p=rots_out$pvalue)

    res_mat<-matrix(nrow = nrow(data), ncol = 1)
    rownames(res_mat)<-rownames(data)

    res_mat[match(rownames(rots_res), rownames(res_mat)),1]<-as.numeric(rots_res$p)
    res_mat
  }
  parallel::stopCluster(cl)
  rownames(rots_frame)=rownames(data)
  colnames(rots_frame)=c(0, seq(1:degree_PolyReg))

  #Combine the result from different degrees
  combin_fun<-"min" #minimum the default. Should alternatives be allowed? Rank product? Not allowed right now. CHECK!
  c_func<-get(combin_fun) #minimum now
  fin_p<-apply(rots_frame, 1, function(x) c_func(x, na.rm = T))
  fin_p[is.infinite(fin_p)]<-NA
  fin_res<-cbind(id=names(fin_p), p=as.numeric(fin_p))
  fin_res<-data.frame(fin_res, stringsAsFactors = F)
  rownames(fin_res)<-rownames(data)
  colnames(fin_res)<-c("id","rep p-value")
  fin_res[,1]<-as.character(fin_res[,1])
  fin_res[,2]<-as.numeric(as.character(fin_res[,2]))

  if(all(is.na(fin_res[,2]))){stop("Unkown error during DiffROTS.")} #Problem is that we have a lot of zeros with ROTS this way.
  ret_list=list(rots_frame, fin_res)
  names(ret_list)<-c("all cond p-values", "rep p-values")
  return(ret_list)
}
