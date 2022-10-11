#Functions needed in the process by RolDE

#to help with foreach
`%dopar%` <- foreach::`%dopar%`
`%do%` <- foreach::`%do%`
`%dorng%` <- doRNG::`%dorng%`

#shifter helper function
shifter <- function(x, n = 1) {
        if (n == 0) x else c(tail(x, -n), head(x, n))
}

#A Function to divide the individual into ROTS runs to preserve the proper degrees of
#freedom for statistical testing, i.e. so that each individual is used only once
#in each run. Make the sizes of the runs as balanced as possible.

#Renewed 4.10.2022
getUniqueCombRunsNew<-function(real_combs, unique_conditions){
        test_mat<-lapply(real_combs, function(x){
                temp<-as.numeric(unlist(strsplit(x, ",")[[1]]))
        })
        test_mat<-matrix(unlist(test_mat), ncol=2, byrow = TRUE)
        
        real_ind_cond1<-unique(test_mat[,1])
        real_ind_cond2<-unique(test_mat[,2])
        
        lengths1<-unlist(lapply(real_ind_cond1, function(x){length(which(test_mat[,1]==x))}))
        names(lengths1)<-real_ind_cond1
        
        lengths2<-unlist(lapply(real_ind_cond2, function(x){length(which(test_mat[,2]==x))}))
        names(lengths2)<-real_ind_cond2
        
        min_comps<-max(lengths1, lengths2)
        
        temp_list<-lapply(seq_len(min_comps), function(x){matrix(nrow=2, ncol=0)})
        
        placed<-logical(nrow(test_mat))
        start_time<-Sys.time() #make sure doesn't stay in eternal loop
        while(!all(placed)){ # make a first run list, this should always converge?
                for(i in seq_len(nrow(test_mat))){
                        if(placed[i]){next}
                        temp_comp<-test_mat[i,]
                        for(j in seq_len(length(temp_list))){
                                temp_mat<-temp_list[[j]]
                                if(placed[i]){next}
                                if(ncol(temp_mat)==0){
                                        temp_mat<-cbind(temp_mat, temp_comp)
                                        temp_list[[j]]<-temp_mat
                                        placed[i]<-TRUE
                                }else{
                                        if(!temp_comp[1]%in%temp_mat[1,] & !temp_comp[2]%in%temp_mat[2,]){
                                                temp_mat<-cbind(temp_mat, temp_comp)
                                                temp_list[[j]]<-temp_mat
                                                placed[i]<-TRUE
                                        }
                                }
                        }
                }
                cur_time<-Sys.time()
                if(as.numeric(cur_time)-as.numeric(start_time)>300){stop("Failure during ROTS run assignment. Could not assign runs.")}
        }
        
        #Now try to balance runs if possible.
        converged<-FALSE
        l1<-lengths(temp_list)
        max_l<-max(l1)
        min_l<-min(l1)
        
        org_diff<-max_l-min_l
        org_diff_diffs<-numeric()
        
        start_time<-Sys.time() 
        while(!converged){
                changed<-FALSE
                l1<-lengths(temp_list)
                max_l_loc<-which.max(l1)[1]
        
                temp_mat<-temp_list[[max_l_loc]]
                ords<-rank(l1, ties.method = "first")
                names(ords)<-seq_len(length(temp_list))
                ords<-sort(ords)
                
                for(i in seq_len(length(ords))){
                        temp_mat_compare<-temp_list[[as.numeric(names(ords)[i])]] #start filling from the smallest end
                        pos_moves<-apply(temp_mat, 2, function(x) {!x[1]%in%temp_mat_compare[1,]&!x[2]%in%temp_mat_compare[2,]})
                        if(any(pos_moves)){
                                sel_move<-which(pos_moves)[1] #always just move one
                                move_comp<-temp_mat[,sel_move]
                                temp_mat<-temp_mat[,-sel_move]
                                temp_mat_compare<-cbind(temp_mat_compare, move_comp)
                                temp_list[[max_l_loc]]<-temp_mat
                                temp_list[[as.numeric(names(ords)[i])]]<-temp_mat_compare
                                changed<-TRUE
                                break #break from the for loop, only take one element in each round
                        }
                }
                diff2<-max(lengths(temp_list))-min(lengths(temp_list))
                org_diff_diffs<-c(org_diff_diffs, org_diff-diff2)
                
                if(!changed){converged<-TRUE}
                
                if(length(org_diff_diffs)>50){ #look at least 50 rounds
                        if(all(org_diff_diffs[c((length(org_diff_diffs)-49):(length(org_diff_diffs)))]==org_diff_diffs[length(org_diff_diffs)])){
                                converged<-TRUE
                        }
                }
                cur_time<-Sys.time()
                if(as.numeric(cur_time)-as.numeric(start_time)>300){stop("Failure during ROTS run assignment. Could not assign runs.")}
        }
        
        temp_list<-lapply(temp_list, function(x) {
                rownames(x)<-unique_conditions
                x
        })
        
        return(temp_list)
}

#Function to fill in missing values, only fill if a row has less than 2 missing values per group.
#Do not impute more than two values per group for a row.
fillGaps_New<-function(all_num_values, row, groups_for_rots){
        row2<-unlist(lapply(seq_len(max(groups_for_rots)), function(x){
                locs_na<-which(is.na(row[which(groups_for_rots==x)]))
                locs_nonna<-which(!is.na(row[which(groups_for_rots==x)]))

                row_vals<-row[which(groups_for_rots==x)]
                if(length(locs_nonna)<2){
                        imp.num<-2-length(locs_nonna)
                        vals<-sample(all_num_values, length(imp.num), replace = TRUE)
                        row_vals[locs_na][seq_len(imp.num)]<-vals
                }
                row_vals
        }))
        return(row2)
}

#Fill in all gaps in in a row.
fillGapsAll_New<-function(all_num_values, row){
        locs_na<-which(is.na(row))
        vals<-sample(all_num_values, length(locs_na), replace = TRUE)
        row[locs_na]<-vals
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
                loc<-which(abs(p - s) == min(abs(p - s), na.rm = TRUE))
                sr<-as.numeric(names(p)[loc])
                if(length(sr)>1){
                        sr<-sample(c(sr),1)
                }
        }
        return(sr)
}

#Get simulated internal rank products for RegROTS, DiffROTS with weights
getSimRPsWeights<-function(pvals, sim_p, weights){
        est_ranks<-apply(pvals, 2, function(x){
                locs_non_na<-which(!is.na(x))
                sim_sample<-sample(x = sim_p, size = length(locs_non_na), replace = TRUE)
                rank_org<-rank(x, na.last = TRUE, ties.method = "random")
                names(x)<-rank_org
                sim_sample<-sort(sim_sample)
                sim_sample<-sim_sample[rank_org]
                sr<-as.numeric(unlist(lapply(sim_sample, function(y) getRank(p = x, s = y))))
                locs_na<-which(is.na(x))
                if(length(locs_na)>0){
                        max_sr<-max(sr[-locs_na], na.rm = TRUE)
                        sr[locs_na]<-seq(from=(max_sr+1), length.out = length(locs_na))
                }
                return(sr)
        })
        #sim_rp=apply(est_ranks, 1, function(x) exp(mean(log(x), na.rm=T)))
        sim_rp<-exp(matrixStats::rowWeightedMeans(x=log(est_ranks),w = weights,  na.rm = TRUE))
        if(any(is.nan(sim_rp))){sim_rp[which(is.nan(sim_rp))]<-NA}
        return(sim_rp)
}

#Get Simulated P-values for non aligned timepoint DiffROTS
getSimPDif<-function(pvals, sim_p){
        est_pvals<-apply(pvals, 2, function(x){
                locs_non_na<-which(!is.na(x))
                sim_sample<-sample(x = sim_p, size = length(locs_non_na), replace = TRUE)
                rank_org<-rank(x, na.last = TRUE, ties.method = "random")
                names(x)<-rank_org
                sim_sample<-sort(sim_sample)
                sim_sample<-sim_sample[rank_org]
        })
        sim_pval<-apply(est_pvals,1,min,na.rm=TRUE) #now built for minimum.
        if(any(is.infinite(sim_pval))){sim_pval[which(is.infinite(sim_pval))]<-NA}
        return(sim_pval)
}

#Get Simulated P-values for PolyReg
getSPPoly<-function(all_p_vals, sim_p, ord_poly, all_cond_pvals){
        est_p<-sample(x = sim_p, size = length(na.omit(all_p_vals)), replace = TRUE)
        est_p<-sort(est_p)
        est_p<-est_p[ord_poly] #If there are NAs, they get put back in at the right indexes
        #if(any(is.na(all_p_vals))){est_p[which(is.na(all_p_vals))]=NA}
        all_p_sim<-matrix(est_p,nrow = nrow(all_cond_pvals),ncol = ncol(all_cond_pvals))
        p_sim<-apply(all_p_sim, 1, min, na.rm=TRUE)
        if(any(is.infinite(p_sim))){p_sim[which(is.infinite(p_sim))]<-NA}
        return(p_sim)
}


#Get simulated overall Rank products when weighing the runs.
getSimRankProdsWeights<-function(reg_rots_pval, diff_rots_pval, polyreg_pval, all_poly_pvals, ord_poly, sim_p, rp_reg, rp_diff, p_poly, regrots_weigths, diffrots_weigths, sigValSampN){

        sim_rp_reg<-getSimRPsWeights(pvals = reg_rots_pval, sim_p = sim_p, weights = regrots_weigths)
        sim_rp_diff<-getSimRPsWeights(pvals = diff_rots_pval, sim_p = sim_p, weights = diffrots_weigths)

        sim_p_poly<-getSPPoly(all_p_vals = all_poly_pvals, sim_p = sim_p, ord_poly = ord_poly, all_cond_pvals = polyreg_pval)

        sim_r1<-as.numeric(unlist(lapply(sim_rp_reg, function(x) getRank(p = rp_reg, s = x))))
        sim_r2<-as.numeric(unlist(lapply(sim_rp_diff, function(x) getRank(p = rp_diff, s = x))))
        sim_r3<-as.numeric(unlist(lapply(sim_p_poly, function(x) getRank(p = p_poly, s = x))))

        locs_na<-which(is.na(sim_r1))
        if(length(locs_na)>0){
                max_sr<-max(sim_r1[-locs_na], na.rm = TRUE)
                sim_r1[locs_na]<-seq(from=(max_sr+1), length.out = length(locs_na))
        }

        locs_na<-which(is.na(sim_r2))
        if(length(locs_na)>0){
                max_sr<-max(sim_r2[-locs_na], na.rm = TRUE)
                sim_r2[locs_na]<-seq(from=(max_sr+1), length.out = length(locs_na))
        }

        locs_na<-which(is.na(sim_r3))
        if(length(locs_na)>0){
                max_sr<-max(sim_r3[-locs_na], na.rm = TRUE)
                sim_r3[locs_na]<-seq(from=(max_sr+1), length.out = length(locs_na))
        }

        sim_ranks<-cbind(sim_r1, sim_r2, sim_r3)
        # sim_rank_prods<-apply(sim_ranks, 1, function(x) exp(mean(log(x), na.rm=T)))
        sim_rank_prods<-exp(rowMeans(log(sim_ranks), na.rm = TRUE))
}

#Get simulated rank products in non-aligned timepoint data with weights.
getSimRankProdsNonAlignedWeights<-function(reg_rots_pval, diff_rots_pval, polyreg_pval, all_poly_pvals, ord_poly, sim_p, rp_reg, p_diff, p_poly, regrots_weigths, sigValSampN){

        sim_rp_reg<-getSimRPsWeights(pvals = reg_rots_pval, sim_p = sim_p, weights = regrots_weigths)
        sim_p_diff<-getSimPDif(pvals = diff_rots_pval, sim_p = sim_p)

        sim_p_poly<-getSPPoly(all_p_vals = all_poly_pvals, sim_p = sim_p, ord_poly = ord_poly, all_cond_pvals = polyreg_pval)

        sim_r1<-as.numeric(unlist(lapply(sim_rp_reg, function(x) getRank(p = rp_reg, s = x))))
        sim_r2<-as.numeric(unlist(lapply(sim_p_diff, function(x) getRank(p = p_diff, s = x))))
        sim_r3<-as.numeric(unlist(lapply(sim_p_poly, function(x) getRank(p = p_poly, s = x))))

        locs_na<-which(is.na(sim_r1))
        if(length(locs_na)>0){
                max_sr<-max(sim_r1[-locs_na], na.rm = TRUE)
                sim_r1[locs_na]<-seq(from=(max_sr+1), length.out = length(locs_na))
        }

        locs_na<-which(is.na(sim_r2))
        if(length(locs_na)>0){
                max_sr<-max(sim_r2[-locs_na], na.rm = TRUE)
                sim_r2[locs_na]<-seq(from=(max_sr+1), length.out = length(locs_na))
        }

        locs_na<-which(is.na(sim_r3))
        if(length(locs_na)>0){
                max_sr<-max(sim_r3[-locs_na], na.rm = TRUE)
                sim_r3[locs_na]<-seq(from=(max_sr+1), length.out = length(locs_na))
        }

        sim_ranks<-cbind(sim_r1, sim_r2, sim_r3)
        # sim_rank_prods<-apply(sim_ranks, 1, function(x) exp(mean(log(x), na.rm=T)))
        sim_rank_prods<-exp(rowMeans(log(sim_ranks), na.rm = TRUE))
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

