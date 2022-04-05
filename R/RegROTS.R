RegROTS<-function(data, des_matrix, min_feat_obs, degree_RegROTS, rots_runs, n_cores, aligned){

  degree<-degree_RegROTS
  unique_conditions<-unique(as.character(des_matrix[,2]))

  #Determine case and control. Control is always the first condition encountered in the design matrix.
  #Doesn't really matter.
  control<-unique_conditions[1]
  case<-unique_conditions[2]

  #Function to get orthogonal polynomials in case of nonaligned timepoints
  my_poly <- function (x, degree = 1) {
    ## check feasibility
    if (length(unique(x)) < degree)
      stop("insufficient unique data points for specified degree!")
    ## centring covariates (so that `x` is orthogonal to intercept)
    centre <- mean(x)
    x <- x - centre
    beta <- alpha <- norm2 <- numeric(degree)
    ## repeat first polynomial `x` on all columns to initialize design matrix X
    X <- matrix(x, nrow = length(x), ncol = degree, dimnames = list(NULL, seq_len(degree)))
    ## compute alpha[1] and beta[1]
    norm2[1] <- new_norm <- drop(crossprod(x))
    alpha[1] <- sum(x ^ 3) / new_norm
    beta[1] <- new_norm / length(x)
    if (degree > 1L) {
      old_norm <- new_norm
      ## second polynomial
      X[, 2] <- Xi <- (x - alpha[1]) * X[, 1] - beta[1]
      norm2[2] <- new_norm <- drop(crossprod(Xi))
      alpha[2] <- drop(crossprod(Xi * Xi, x)) / new_norm
      beta[2] <- new_norm / old_norm
      old_norm <- new_norm
      ## further polynomials obtained from recursion
      i <- 3
      while (i <= degree) {
        X[, i] <- Xi <- (x - alpha[i - 1]) * X[, i - 1] - beta[i - 1] * X[, i - 2]
        norm2[i] <- new_norm <- drop(crossprod(Xi))
        alpha[i] <- drop(crossprod(Xi * Xi, x)) / new_norm
        beta[i] <- new_norm / old_norm
        old_norm <- new_norm
        i <- i + 1
      }
    }
    ## add attributes and return
    attr(X, "coefs") <- list(centre = centre, alpha = alpha[-degree], beta = beta[-degree])
    X
  }

  #The main function to compare coefficients
  getCoef<-function(r1, r2, time1, time2, degree, thres_feat_both_cond, min_feat_obs, aligned) {

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

      if(aligned){
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

      }else{
        poly.time1<-my_poly(x = time1, degree = d1)
        poly.time2<-my_poly(x = time2, degree = d2)

        mod1<-tryCatch({
          lm(as.numeric(r1)~poly.time1)
        },error = function(e) {
          NULL
        })

        mod2<-tryCatch({
          lm(as.numeric(r2)~poly.time2)
        },error = function(e) {
          NULL
        })

      }

      if(is.null(mod1) | is.null(mod2)){
        diffs<-rep(NA, degree+1)
      } else {
        coeffs1<-mod1$coefficients
        coeffs2<-mod2$coefficients
        diffs<-numeric(degree+1)
        diffs[seq_len(min(d1,d2, na.rm = TRUE)+1)]<-coeffs1[seq_len(min(d1,d2, na.rm = TRUE)+1)]-coeffs2[seq_len(min(d1,d2, na.rm = TRUE)+1)]
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
  rots_res_frame <- foreach::foreach(run=seq_len(length(rots_runs)), .export = c("fillGaps_New", "fillGapsAll_New"), .packages = c("ROTS", "doParallel", "matrixStats", "foreach"), .combine=cbind) %dorng% {
    #Determine groups for ROTS
    run_comparisons<-rots_runs[[run]]
    replicate_comparisons<-ncol(run_comparisons)
    groups_for_rots<-sort(rep(seq_len(degree+1),replicate_comparisons))

    #Get the sample locations (columns) for the individuals in the comparisons for the current run.
    control_locs<-as.list(as.data.frame(apply(run_comparisons, 2, function(x){which(des_matrix$Individual%in%x[1])})))
    case_locs<-as.list(as.data.frame(apply(run_comparisons, 2, function(x){which(des_matrix$Individual%in%x[2])})))

    #Store the timepoints for each individual in each comparison here
    times_control<-lapply(control_locs, function(x) {des_matrix$Timepoint[x]})
    times_case<-lapply(case_locs, function(x) {des_matrix$Timepoint[x]})

    #A data frame to save the coefficient differences in.
    coef_frame<-data.frame(matrix(nrow = nrow(data), ncol = replicate_comparisons*(degree+1)))
    rownames(coef_frame)<-rownames(data)
    colnames(coef_frame)<-groups_for_rots

    for(r in seq_len(nrow(data))) {

      res_row<-numeric(ncol(coef_frame))

      for(comp in seq_len(ncol(run_comparisons))){

        cont_loc<-control_locs[[comp]]
        case_loc<-case_locs[[comp]]

        row_control<-data[r,cont_loc]
        row_case<-data[r,case_loc]

        time_control<-times_control[[comp]]
        time_case<-times_case[[comp]]

        thres_feat_both_cond<-(length(time_control)+length(time_case))-4 #Internal threshold. Guarantees the filtering away extreme scenarios. Doesn't really play a role? To be removed in the future?

        res_vals<-getCoef(r1 = row_case, r2 = row_control, time1 = time_case, time2 = time_control, degree = degree, thres_feat_both_cond = thres_feat_both_cond, min_feat_obs = min_feat_obs, aligned = aligned)
        res_row[seq(from=comp, to=length(res_row), by=replicate_comparisons)]=res_vals

      } #end comparison loop
      coef_frame[r,]<-res_row
    } #end getting coefficients for rows

    #Scale the coefficients differences within each group to acquire similar distributions for correct testing in ROTS
    coef_temp<-lapply(unique(groups_for_rots), function(x){
      data_temp1<-coef_frame[, which(groups_for_rots==x)]
      data_temp2<-as.numeric(unlist(data_temp1))
      data_temp2<-data_temp2/sd(data_temp2, na.rm = TRUE)
      data_temp1[,c(seq_len(ncol(data_temp1)))]<-data_temp2
      data_temp1
    })

    coef_scaled=do.call(cbind, coef_temp)
    rownames(coef_scaled)<-rownames(coef_frame)
    colnames(coef_scaled)<-colnames(coef_frame)

    #Remove features with noth enough values in the first two groups, intercept (expression) level
    #differences and differences in longitudinal linear patterns of expression.
    #Require at least two values for group, 2 required for ROTS.
    rem_feat<-unique(unlist(lapply(seq_len(2), function(z){
      which(apply(coef_scaled[,which(groups_for_rots==z)], 1, function(x) (sum(!is.na(x))<2)))
    })))

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
    #could also be values between quantiles 0.25-0.75, so more pronounced zero values. Change if needed.
    all_num_values<-as.numeric(unlist(coef_scaled))
    all_num_values<-na.omit(all_num_values)

    #get enough quantile values used for imputation
    all_quant_vals<-quantile(all_num_values, seq(from=0.001, to=1, by=0.001))

    coef_scaled<-data.frame(t(apply(coef_scaled, 1, function(x) {fillGaps_New(all_quant_vals = all_quant_vals, row = x, groups_for_rots = groups_for_rots)})), check.names = FALSE)

    #Add a zeroish reference group to compare to.
    ref_group<-data.frame(matrix(nrow = nrow(coef_scaled), ncol = length(which(groups_for_rots==1))))

    #Fill compeletely in with random values from the distribution of all differences
    ref_group<-data.frame(t(apply(ref_group, 1, function(x) {fillGapsAll_New(all_quant_vals = all_quant_vals, row = x)})), check.names = FALSE)

    #Make a new dataset and update groups for ROTS
    coef_new<-cbind(coef_scaled, ref_group)
    groups_for_rots<-c(groups_for_rots, rep((max(groups_for_rots)+1), length(which(groups_for_rots==1))))
    colnames(coef_new)<-groups_for_rots

    #Perform multigroup ROTS
    suppressMessages(expr = {rots_out<-ROTS::ROTS(data = coef_new, groups = groups_for_rots, B = 100, K = nrow(coef_new)/4, paired = FALSE, progress = FALSE)})
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

  ranks<-apply(rots_res_frame, 2, function(x) rank(x, na.last = TRUE, ties.method = "average"))

  nr_comps<-as.numeric(unlist(lapply(rots_runs, ncol)))
  max_comps<-max(nr_comps)
  weights_runs<-nr_comps/max_comps
  rank_prods<-apply(ranks, 1, function(x) exp(weighted.mean(log(x), weights_runs, na.rm=TRUE)))

  names(rank_prods)<-rownames(rots_res_frame)

  fin_res<-cbind(id=names(rank_prods), rp=rank_prods)
  if(any(is.nan(fin_res[,2]))){fin_res[which(is.nan(fin_res[,2])),2]=NA}
  if(any(fin_res[,2]=="NaN", na.rm = TRUE)){fin_res[which(fin_res[,2]=="NaN"),2]=NA}
  fin_res<-data.frame(fin_res, stringsAsFactors = FALSE)
  fin_res[,1]<-as.character(fin_res[,1])
  fin_res[,2]<-as.numeric(as.character(fin_res[,2]))

  if(all(is.na(fin_res[,2]))){stop("Unkown error during RegROTS.")}

  #make the return list
  ret_list<-list(rots_res_frame, fin_res, weights_runs)
  names(ret_list)<-c("p-values", "rank products", "weights for runs")
  return(ret_list)
}
