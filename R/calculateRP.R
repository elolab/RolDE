calculateRP<-function(res_regrots, res_diffrots, res_polyreg){

  res_diffrots<-res_diffrots[rownames(res_regrots),] #Make sure the order is the same for all results.
  res_polyreg<-res_polyreg[rownames(res_regrots),]

  r1<-rank(as.numeric(res_regrots[,2]),ties.method = "average",  na.last = TRUE)
  r2<-rank(as.numeric(res_diffrots[,2]),ties.method = "average", na.last = TRUE)
  r3<-rank(as.numeric(res_polyreg[,2]),ties.method = "average",  na.last = TRUE)

  ranks<-cbind(r1, r2, r3)
  rownames(ranks)<-rownames(res_regrots)

  rank_prods=exp(rowMeans(log(ranks), na.rm = TRUE))

  if(any(is.nan(rank_prods))){rank_prods[is.nan(rank_prods)]<-NA}
  if(any(rank_prods=="NaN", na.rm = TRUE)){rank_prods[which(rank_prods=="NaN")]=NA}

  fin_res<-cbind(id=names(rank_prods), rp=rank_prods)
  fin_res<-data.frame(fin_res, stringsAsFactors = FALSE)
  fin_res[,1]<-as.character(fin_res[,1])
  fin_res[,2]<-as.numeric(as.character(fin_res[,2]))

  if(all(is.na(fin_res[,2]))){stop("Unkown error during Rank Product Calculation.")}

  return(fin_res)
}
