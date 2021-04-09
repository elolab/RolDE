#Functions needed in the process by RolDE

#to help with foreach
`%dopar%` <- foreach::`%dopar%`
`%do%` <- foreach::`%do%`

#shifter helper function
shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

#A Function to divide the individual into ROTS runs to preserve the proper degrees of
#freedom for statistical testing, i.e. so that each individual is used only once
#in each run. Make the sizes of the runs as balanced as possible.
getUniqueCombRunsNew<-function(real_combs, unique_conditions){ #Still corrected 23.9.2020
  test_mat<-matrix(nrow = length(real_combs), ncol=2)
  for(i in 1:nrow(test_mat)){
    temp<-as.numeric(unlist(strsplit(real_combs[i], ",")[[1]]))
    test_mat[i,]<-temp
  }

  real_ind_cond1<-unique(test_mat[,1])
  real_ind_cond2<-unique(test_mat[,2])

  lengths1<-numeric(length(real_ind_cond1))
  for(i in 1:length(lengths1)){
    lengths1[i]<-length(which(test_mat[,1]==real_ind_cond1[i]))
  }

  names(lengths1)<-real_ind_cond1

  lengths2<-numeric(length(real_ind_cond2))
  for(i in 1:length(lengths2)){
    lengths2[i]<-length(which(test_mat[,2]==real_ind_cond2[i]))
  }

  names(lengths2)<-real_ind_cond2

  min_comps<-max(lengths1, lengths2)

  temp_list<-list()
  for(i in 1:min_comps){
    temp_list[[i]]<-matrix(nrow=2, ncol=0)
  }

  used_inds=list()
  for(i in 1:min_comps){
    temp_vals=matrix(nrow=2, ncol=0)
    used_inds[[i]]=temp_vals
  }

  all_lengths<-c(lengths1, lengths2)
  all_lengths<-sort(all_lengths, decreasing = T)

  start_ind<-1

  for(ind in 1:length(all_lengths)){
    if(names(all_lengths[ind])%in%real_ind_cond1){
      ind_cond=1
    }else{ind_cond=2}

    locs_ind<-which(test_mat[,ind_cond]%in%names(all_lengths[ind]))
    ind_mat<-test_mat[locs_ind,]
    if(!is.matrix(ind_mat)){ind_mat<-t(as.matrix(ind_mat))}
    found<-logical(nrow(ind_mat))

    for(j in 1:nrow(ind_mat)){
      for(i in 1:length(temp_list)){
        temp_mat<-temp_list[[i]]
        if(ncol(temp_mat)==0){next}
        for(co in 1:ncol(temp_mat)){
          if(all(ind_mat[j,]==temp_mat[,co]))
            found[j]=T
        }
      }
    }

    ind_mat<-ind_mat[!found,]
    if(!is.matrix(ind_mat)){ind_mat<-t(as.matrix(ind_mat))}
    if(nrow(ind_mat)==0){next}
    continue=TRUE
    start_time=Sys.time()
    stop_ind<-0
    while(continue){
      for(i in start_ind:length(temp_list)){
        temp_mat<-temp_list[[i]]
        if(!is.matrix(ind_mat)){ind_mat<-t(as.matrix(ind_mat))}
        if(nrow(ind_mat)==0){next}
        if(ncol(temp_mat)==0){
          used_inds[[i]]=cbind(used_inds[[i]], ind_mat[1,])
          temp_mat<-cbind(temp_mat, ind_mat[1,])
          ind_mat<-ind_mat[-1,]
          if(!is.matrix(ind_mat)){ind_mat<-t(as.matrix(ind_mat))}
          if(nrow(ind_mat)==0){
            stop_ind<-i
          }
        }else{
          if(!ind_mat[1,1]%in%used_inds[[i]][1,] & !ind_mat[1,2]%in%used_inds[[i]][2,] ){
            used_inds[[i]]=cbind(used_inds[[i]], ind_mat[1,])
            temp_mat<-cbind(temp_mat, ind_mat[1,])
            ind_mat<-ind_mat[-1,]
            if(!is.matrix(ind_mat)){ind_mat<-t(as.matrix(ind_mat))}
            if(nrow(ind_mat)==0){
              stop_ind<-i
            }
          }

        }

        temp_list[[i]]<-temp_mat
      }

      if(stop_ind!=0){
        if(stop_ind==length(temp_list)){
          start_ind<-1
        }else{
          start_ind<-stop_ind+1
        }
      }else{
        start_ind<-1
      }

      if(nrow(ind_mat)==0){continue=FALSE}
      cur_time=Sys.time()
      if(as.numeric(cur_time-start_time)>300){stop("Error during ROTS run assignment. Could not assign runs.")}
    }
  }

  for(i in 1:length(temp_list)){
    rownames(temp_list[[i]])<-unique_conditions
  }

  return(temp_list)
}

#Function to fill in missing values, only fill if a row has less than 2 missing values per group.
#Do not impute more than two values per group for a row.
fillGaps_New<-function(all_quant_vals, row, groups_for_rots, rand_seed){
  for(group in 1:max(groups_for_rots)){
    locs_na<-which(is.na(row[which(groups_for_rots==group)]))
    locs_nonna<-which(!is.na(row[which(groups_for_rots==group)]))
    if(length(locs_nonna)<2){
      imp.num<-2-length(locs_nonna)
      #set.seed(rand_seed)
      vals<-sample(1:1000, length(imp.num), replace = T)
      vals<-all_quant_vals[vals]
      row[which(groups_for_rots==group)][locs_na][1:imp.num]<-vals
    }
  }
  return(row)
}

#Fill in all gaps in in a row.
fillGapsAll_New<-function(all_quant_vals, row, rand_seed){
  locs_na<-which(is.na(row))
  set.seed(rand_seed)
  vals<-sample(1:1000, length(locs_na), replace = T)
  vals<-all_quant_vals[vals]
  row[locs_na]=vals
  return(row)
}


#Get a rank for a simulated p value based on the experimental data.
#p is the pvalues from the experiment, e.g, from a single rots run, a vector.
#the ranks within the vector need to be attached as names
#s is a single simulated p-value
getRank<-function(p,s){
  if(is.na(s)){
    sr<-NA
  }else{
    loc<-which(abs(p - s) == min(abs(p - s), na.rm = T))
    sr<-as.numeric(names(p)[loc])
    if(length(sr)>1){
      sr=sample(c(sr),1)
    }
  }
  return(sr)
}

#Get simulated internal rank products for RegROTS, DiffROTS
getSimRPs<-function(pvals, sim_p, rand_seed, na_treat){
  pvals=rbind(seq(1:ncol(pvals)), pvals)
  est_ranks<-apply(pvals, 2, function(x){
    us.seed=rand_seed+as.numeric(x[1])
    x=x[-1]
    set.seed(us.seed)
    locs_non_na<-which(!is.na(x))
    sim_sample<-sample(x = sim_p, size = length(locs_non_na), replace = T)
    rank_org<-rank(x, na.last = T, ties.method = "random")
    names(x)<-rank_org
    sim_sample<-sort(sim_sample)
    sim_sample<-sim_sample[rank_org]
    sr<-as.numeric(unlist(lapply(sim_sample, function(y) getRank(p = x, s = y))))
    if(na_treat=="last"){ #add ranks for NA values, else keep them as such
      locs_na<-which(is.na(x))
      if(length(locs_na)>0){
        max_sr<-max(sr[-locs_na], na.rm = T)
        sr[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
      }
    }
    return(sr)
  })
  #sim_rp=apply(est_ranks, 1, function(x) exp(mean(log(x), na.rm=T)))
  sim_rp=exp(rowMeans(log(est_ranks), na.rm = T))
  if(any(is.nan(sim_rp))){sim_rp[which(is.nan(sim_rp))]=NA}
  return(sim_rp)
}
#Get simulated internal rank products for RegROTS, DiffROTS with weights
getSimRPsWeights<-function(pvals, sim_p, rand_seed, na_treat, weights){
  pvals=rbind(seq(1:ncol(pvals)), pvals)
  est_ranks<-apply(pvals, 2, function(x){
    us.seed=rand_seed+as.numeric(x[1])
    x=x[-1]
    set.seed(us.seed)
    locs_non_na<-which(!is.na(x))
    sim_sample<-sample(x = sim_p, size = length(locs_non_na), replace = T)
    rank_org<-rank(x, na.last = T, ties.method = "random")
    names(x)<-rank_org
    sim_sample<-sort(sim_sample)
    sim_sample<-sim_sample[rank_org]
    sr<-as.numeric(unlist(lapply(sim_sample, function(y) getRank(p = x, s = y))))
    if(na_treat=="last"){ #add ranks for NA values, else keep them as such
      locs_na<-which(is.na(x))
      if(length(locs_na)>0){
        max_sr<-max(sr[-locs_na], na.rm = T)
        sr[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
      }
    }
    return(sr)
  })
  #sim_rp=apply(est_ranks, 1, function(x) exp(mean(log(x), na.rm=T)))
  sim_rp=exp(matrixStats::rowWeightedMeans(x=log(est_ranks),w = weights,  na.rm = T))
  if(any(is.nan(sim_rp))){sim_rp[which(is.nan(sim_rp))]=NA}
  return(sim_rp)
}

#Get Simulated P-values for non aligned timepoint DiffROTS
getSimPDif<-function(pvals, sim_p, rand_seed, na_treat){
  pvals=rbind(seq(1:ncol(pvals)), pvals)
  est_pvals<-apply(pvals, 2, function(x){
    us.seed=rand_seed+as.numeric(x[1])
    x=x[-1]
    set.seed(us.seed)
    locs_non_na<-which(!is.na(x))
    sim_sample<-sample(x = sim_p, size = length(locs_non_na), replace = T)
    rank_org<-rank(x, na.last = T, ties.method = "random")
    names(x)<-rank_org
    sim_sample<-sort(sim_sample)
    sim_sample<-sim_sample[rank_org]
    })
  sim_pval=apply(est_pvals,1,min,na.rm=T) #now built for minimum, if alternatives will be provided, this needs to be updated as well. CHECK!
  if(any(is.infinite(sim_pval))){sim_pval[which(is.infinite(sim_pval))]=NA}
  return(sim_pval)
}

#Get Simulated P-values for PolyReg
getSPPoly<-function(all_p_vals, sim_p, ord_poly, all_cond_pvals, rand_seed){
  set.seed(rand_seed)
  est_p<-sample(x = sim_p, size = length(na.omit(all_p_vals)), replace = T)
  est_p<-sort(est_p)
  est_p<-est_p[ord_poly] #If there are NAs, they get put back in at the right indexes
  #if(any(is.na(all_p_vals))){est_p[which(is.na(all_p_vals))]=NA}
  all_p_sim<-matrix(est_p,nrow = nrow(all_cond_pvals),ncol = ncol(all_cond_pvals))
  p_sim<-apply(all_p_sim, 1, min, na.rm=T)
  if(any(is.infinite(p_sim))){p_sim[which(is.infinite(p_sim))]=NA}
  return(p_sim)
}

#Get simulated overall Rank products.
getSimRankProds<-function(reg_rots_pval, diff_rots_pval, polyreg_pval, all_poly_pvals, ord_poly, sim_p, rp_reg, rp_diff, p_poly, rand_seed, na_treat, sigValSampN){

  #always repeatable seeds, but different for each individual sample.
  seed1<-rand_seed
  seed2<-rand_seed+(sigValSampN+1)
  seed3<-rand_seed+((sigValSampN*2)+1)

  sim_rp_reg<-getSimRPs(pvals = reg_rots_pval, sim_p = sim_p, rand_seed = seed1, na_treat = na_treat)
  sim_rp_diff<-getSimRPs(pvals = diff_rots_pval, sim_p = sim_p, rand_seed = seed2, na_treat = na_treat)

  sim_p_poly<-getSPPoly(all_p_vals = all_poly_pvals, sim_p = sim_p, ord_poly = ord_poly, all_cond_pvals = polyreg_pval, seed3)

  sim_r1<-as.numeric(unlist(lapply(sim_rp_reg, function(x) getRank(p = rp_reg, s = x))))
  sim_r2<-as.numeric(unlist(lapply(sim_rp_diff, function(x) getRank(p = rp_diff, s = x))))
  sim_r3<-as.numeric(unlist(lapply(sim_p_poly, function(x) getRank(p = p_poly, s = x))))

  if(na_treat=="last"){ #add ranks for NA values, else keep them as such
    locs_na<-which(is.na(sim_r1))
    if(length(locs_na)>0){
      max_sr<-max(sim_r1[-locs_na], na.rm = T)
      sim_r1[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
    }

    locs_na<-which(is.na(sim_r2))
    if(length(locs_na)>0){
      max_sr<-max(sim_r2[-locs_na], na.rm = T)
      sim_r2[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
    }

    locs_na<-which(is.na(sim_r3))
    if(length(locs_na)>0){
      max_sr<-max(sim_r3[-locs_na], na.rm = T)
      sim_r3[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
    }
  }

  sim_ranks=cbind(sim_r1, sim_r2, sim_r3)
 # sim_rank_prods<-apply(sim_ranks, 1, function(x) exp(mean(log(x), na.rm=T)))
  sim_rank_prods<-exp(rowMeans(log(sim_ranks), na.rm = T))
}

#Get simulated overall Rank products when weighing the runs.
getSimRankProdsWeights<-function(reg_rots_pval, diff_rots_pval, polyreg_pval, all_poly_pvals, ord_poly, sim_p, rp_reg, rp_diff, p_poly, rand_seed, na_treat, regrots_weigths, diffrots_weigths, sigValSampN){

  #always repeatable seeds, but different for each individual sample.
  seed1<-rand_seed
  seed2<-rand_seed+(sigValSampN+1)
  seed3<-rand_seed+((sigValSampN*2)+1)

  sim_rp_reg<-getSimRPsWeights(pvals = reg_rots_pval, sim_p = sim_p, rand_seed = seed1, na_treat = na_treat, weights = regrots_weigths)
  sim_rp_diff<-getSimRPsWeights(pvals = diff_rots_pval, sim_p = sim_p, rand_seed = seed2, na_treat = na_treat, weights = diffrots_weigths)

  sim_p_poly<-getSPPoly(all_p_vals = all_poly_pvals, sim_p = sim_p, ord_poly = ord_poly, all_cond_pvals = polyreg_pval, seed3)

  sim_r1<-as.numeric(unlist(lapply(sim_rp_reg, function(x) getRank(p = rp_reg, s = x))))
  sim_r2<-as.numeric(unlist(lapply(sim_rp_diff, function(x) getRank(p = rp_diff, s = x))))
  sim_r3<-as.numeric(unlist(lapply(sim_p_poly, function(x) getRank(p = p_poly, s = x))))

  if(na_treat=="last"){ #add ranks for NA values, else keep them as such
    locs_na<-which(is.na(sim_r1))
    if(length(locs_na)>0){
      max_sr<-max(sim_r1[-locs_na], na.rm = T)
      sim_r1[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
    }

    locs_na<-which(is.na(sim_r2))
    if(length(locs_na)>0){
      max_sr<-max(sim_r2[-locs_na], na.rm = T)
      sim_r2[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
    }

    locs_na<-which(is.na(sim_r3))
    if(length(locs_na)>0){
      max_sr<-max(sim_r3[-locs_na], na.rm = T)
      sim_r3[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
    }
  }

  sim_ranks=cbind(sim_r1, sim_r2, sim_r3)
  # sim_rank_prods<-apply(sim_ranks, 1, function(x) exp(mean(log(x), na.rm=T)))
  sim_rank_prods<-exp(rowMeans(log(sim_ranks), na.rm = T))
}

#Get simulated rank products in non-aligned timepoint data.
getSimRankProdsNonAligned<-function(reg_rots_pval, diff_rots_pval, polyreg_pval, all_poly_pvals, ord_poly, sim_p, rp_reg, p_diff, p_poly, rand_seed, na_treat, sigValSampN){

  #always repeatable seeds, but different for each individual sample.
  seed1<-rand_seed
  seed2<-rand_seed+(sigValSampN+1)
  seed3<-rand_seed+((sigValSampN*2)+1)

  sim_rp_reg<-getSimRPs(pvals = reg_rots_pval, sim_p = sim_p, rand_seed = seed1, na_treat = na_treat)
  sim_p_diff<-getSimPDif(pvals = diff_rots_pval, sim_p = sim_p, rand_seed = seed2, na_treat = na_treat)

  sim_p_poly<-getSPPoly(all_p_vals = all_poly_pvals, sim_p = sim_p, ord_poly = ord_poly, all_cond_pvals = polyreg_pval, seed3)

  sim_r1<-as.numeric(unlist(lapply(sim_rp_reg, function(x) getRank(p = rp_reg, s = x))))
  sim_r2<-as.numeric(unlist(lapply(sim_p_diff, function(x) getRank(p = p_diff, s = x))))
  sim_r3<-as.numeric(unlist(lapply(sim_p_poly, function(x) getRank(p = p_poly, s = x))))

  if(na_treat=="last"){ #add ranks for NA values, else keep them as such
    locs_na<-which(is.na(sim_r1))
    if(length(locs_na)>0){
      max_sr<-max(sim_r1[-locs_na], na.rm = T)
      sim_r1[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
    }

    locs_na<-which(is.na(sim_r2))
    if(length(locs_na)>0){
      max_sr<-max(sim_r2[-locs_na], na.rm = T)
      sim_r2[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
    }

    locs_na<-which(is.na(sim_r3))
    if(length(locs_na)>0){
      max_sr<-max(sim_r3[-locs_na], na.rm = T)
      sim_r3[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
    }
  }

  sim_ranks=cbind(sim_r1, sim_r2, sim_r3)
  # sim_rank_prods<-apply(sim_ranks, 1, function(x) exp(mean(log(x), na.rm=T)))
  sim_rank_prods<-exp(rowMeans(log(sim_ranks), na.rm = T))
}

#Get simulated rank products in non-aligned timepoint data with weights.
getSimRankProdsNonAlignedWeights<-function(reg_rots_pval, diff_rots_pval, polyreg_pval, all_poly_pvals, ord_poly, sim_p, rp_reg, p_diff, p_poly, rand_seed, na_treat, regrots_weigths, sigValSampN){

  #always repeatable seeds, but different for each individual sample.
  seed1<-rand_seed
  seed2<-rand_seed+(sigValSampN+1)
  seed3<-rand_seed+((sigValSampN*2)+1)

  sim_rp_reg<-getSimRPsWeights(pvals = reg_rots_pval, sim_p = sim_p, rand_seed = seed1, na_treat = na_treat, weights = regrots_weigths)
  sim_p_diff<-getSimPDif(pvals = diff_rots_pval, sim_p = sim_p, rand_seed = seed2, na_treat = na_treat)

  sim_p_poly<-getSPPoly(all_p_vals = all_poly_pvals, sim_p = sim_p, ord_poly = ord_poly, all_cond_pvals = polyreg_pval, seed3)

  sim_r1<-as.numeric(unlist(lapply(sim_rp_reg, function(x) getRank(p = rp_reg, s = x))))
  sim_r2<-as.numeric(unlist(lapply(sim_p_diff, function(x) getRank(p = p_diff, s = x))))
  sim_r3<-as.numeric(unlist(lapply(sim_p_poly, function(x) getRank(p = p_poly, s = x))))

  if(na_treat=="last"){ #add ranks for NA values, else keep them as such
    locs_na<-which(is.na(sim_r1))
    if(length(locs_na)>0){
      max_sr<-max(sim_r1[-locs_na], na.rm = T)
      sim_r1[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
    }

    locs_na<-which(is.na(sim_r2))
    if(length(locs_na)>0){
      max_sr<-max(sim_r2[-locs_na], na.rm = T)
      sim_r2[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
    }

    locs_na<-which(is.na(sim_r3))
    if(length(locs_na)>0){
      max_sr<-max(sim_r3[-locs_na], na.rm = T)
      sim_r3[locs_na]=seq(from=(max_sr+1), length.out = length(locs_na))
    }
  }

  sim_ranks=cbind(sim_r1, sim_r2, sim_r3)
  # sim_rank_prods<-apply(sim_ranks, 1, function(x) exp(mean(log(x), na.rm=T)))
  sim_rank_prods<-exp(rowMeans(log(sim_ranks), na.rm = T))
}

#Get length of common timeinterval for two time lines.
getLengthCommonTimepoints<-function(temp_time1, temp_time2){
  if(min(temp_time2)>(min(temp_time1))){global_min<-min(temp_time2)}else{global_min<-min(temp_time1)}
  if(max(temp_time2)<(max(temp_time1))){global_max<-max(temp_time2)}else{global_max<-max(temp_time1)}
  l_timeps<-c(length(which(temp_time1>=global_min & temp_time1<=global_max)), length(which(temp_time2>=global_min & temp_time2<=global_max)))
  return(l_timeps)
}

#Get timepoints in common time interval for two time lines.
getCommonTimepoints<-function(temp_time1, temp_time2){
  if(min(temp_time2)>(min(temp_time1))){global_min<-min(temp_time2)}else{global_min<-min(temp_time1)}
  if(max(temp_time2)<(max(temp_time1))){global_max<-max(temp_time2)}else{global_max<-max(temp_time1)}
  new_time1<-temp_time1[which(temp_time1>=global_min & temp_time1<=global_max)]
  new_time2<-temp_time2[which(temp_time2>=global_min & temp_time2<=global_max)]
  return(list(new_time1, new_time2))
}

#Get common min and max for two time lines.
getCommonRange<-function(temp_time1, temp_time2){
  if(min(temp_time2)>(min(temp_time1))){global_min<-min(temp_time2)}else{global_min<-min(temp_time1)}
  if(max(temp_time2)<(max(temp_time1))){global_max<-max(temp_time2)}else{global_max<-max(temp_time1)}
  return(list(global_min, global_max))
}

