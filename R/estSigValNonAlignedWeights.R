estSigValNonAlignedWeights<-function(data, res_regrots, res_diffrots, res_polyreg, reg_rots_pval, diff_rots_pval, polyreg_pval, rank_products, na_treat, rand_seed, n_cores, sig_adj_meth, sigValSampN, regrots_weigths){

  #Organize
  res_regrots<-res_regrots[rownames(rank_products),]
  res_diffrots<-res_diffrots[rownames(rank_products),]
  res_polyreg<-res_polyreg[rownames(rank_products),]
  reg_rots_pval<-reg_rots_pval[rownames(rank_products),]
  diff_rots_pval<-diff_rots_pval[rownames(rank_products),]
  polyreg_pval<-polyreg_pval[rownames(rank_products),]

  if(!is.matrix(reg_rots_pval) & is.numeric(reg_rots_pval)){
    reg_rots_pval=as.matrix(reg_rots_pval)
  }

  if(!is.matrix(diff_rots_pval) & is.numeric(diff_rots_pval)){
    diff_rots_pval=as.matrix(diff_rots_pval)
  }

  if(!is.matrix(polyreg_pval) & is.numeric(polyreg_pval)){
    polyreg_pval=as.matrix(polyreg_pval)
  }

  rp_reg<-as.numeric(res_regrots[,2])
  p_diff<-as.numeric(res_diffrots[,2])
  p_poly<-as.numeric(res_polyreg[,2])

  rank_reg<-rank(rp_reg,ties.method = "random", na.last = T) #CHECK ties.method
  rank_diff<-rank(p_diff,ties.method = "random", na.last = T)
  rank_poly<-rank(p_poly,ties.method = "random", na.last = T)

  names(rp_reg)<-rank_reg
  names(p_diff)<-rank_diff
  names(p_poly)<-rank_poly

  #create a uniform distribution of enough p-values
  sim_p<-seq(from=0, to=1, length.out = nrow(rank_products)*1000) #Should be enough?
  all_poly_pvals=as.numeric(unlist(polyreg_pval)) #uniform if random data
  ord_poly=rank(all_poly_pvals, na.last = T, ties.method = "random")

  sigValSampN<-ceiling(sigValSampN/nrow(rank_products))
  sim_rps<-matrix(nrow = nrow(data), ncol = sigValSampN)
  co<-NULL

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  sim_rps <- foreach::foreach(co=1:ncol(sim_rps), .export=c("getSimRankProdsNonAligned", "getSPPoly", "getSimPDif", "getSimRPs", "getRank"), .packages = c("foreach"), .combine=cbind) %dopar% {
    used_seed<-rand_seed+co
    col_Val<-getSimRankProdsNonAlignedWeights(reg_rots_pval, diff_rots_pval, polyreg_pval, all_poly_pvals, ord_poly, sim_p, rp_reg, p_diff, p_poly, rand_seed, na_treat, regrots_weigths, sigValSampN)
    col_Val
  }
  parallel::stopCluster(cl)

  #estimate significance values. If the data is random, will result in a uniform p-value
  #distribution.
  all_simrps=as.numeric(unlist(sim_rps))
  rank_prods<-as.numeric(as.character(rank_products[,2]))
  estSigVals=numeric(length(rank_prods))
  for(i in 1:length(estSigVals)){
    if(is.na(rank_prods[i])){
      estSigVals[i]=NA
    }else{
      estSigVals[i]=length(which(all_simrps<=rank_prods[i]))/length(na.omit(all_simrps))
    }
  }

  if(sig_adj_meth=="qvalue"){
    estAdjSigVal<-qvalue::qvalue(estSigVals)
    estAdjSigVal<-as.numeric(estAdjSigVal$qvalues)
  }else{
    estAdjSigVal<-p.adjust(p = estSigVals, method = sig_adj_meth)
  }

  fin_res<-cbind(rank_products, estSigVals, estAdjSigVal)
  fin_res<-data.frame(fin_res, stringsAsFactors = F)
  colnames(fin_res)=c("Feature ID", "RolDE Rank Product", "Estimated Significance Value", "Adjusted Estimated Significance Value")

  fin_res[,1]=as.character(fin_res[,1])
  fin_res[,2]=as.numeric(as.character(fin_res[,2]))
  fin_res[,3]=as.numeric(as.character(fin_res[,3]))
  fin_res[,4]=as.numeric(as.character(fin_res[,4]))

  if(all(is.na(fin_res[,3]))){stop("Unkown error during estimating significance values")}

  return(fin_res)
}
