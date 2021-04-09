determineDegreesNonAligned<-function(des_matrix, degree_RegROTS, degree_PolyReg){

  #Require automatic determination or numerical input
  if(degree_RegROTS!="auto"){
    if(!is.numeric(degree_RegROTS)){stop("Incorrect degree for RegROTS. Numerical value needed.")}
  }

  if(degree_PolyReg!="auto"){
    if(!is.numeric(degree_PolyReg)){stop("Incorrect degree for PolyReg. Numerical value needed.")}
  }

  unique_individuals<-unique(as.numeric(as.character(des_matrix[,4])))

  individual_timepoints<-numeric(length(unique_individuals))
  for(ind in 1:length(unique_individuals)){
    ind_locs<-which(as.numeric(as.character(des_matrix[,4]))==unique_individuals[ind])
    individual_timepoints[ind]<-length(unique(as.numeric(as.character(des_matrix[ind_locs,3]))))
  }

  if(degree_PolyReg=="auto"){
    degree_PolyReg=min((median(as.numeric(individual_timepoints))-1), 5) #Should ever be over 5?
    if(degree_PolyReg<2){degree_PolyReg=2}
  }

  if(degree_RegROTS=="auto"){
    #degree_RegROTS=min(floor(min(as.numeric(individual_timepoints))/2), 4) #is good? Check!
    degree_RegROTS=min(floor(median(as.numeric(individual_timepoints))/2), 4) #is good? Check!
    if(degree_RegROTS<1){degree_RegROTS=1}
  }

  if(degree_RegROTS>=degree_PolyReg){
    message("Warning! The degree for RegROTS as large or larger than the degree for PolyReg. Smaller degree recommended for RegROTS than for PolyReg.") #CHECK! Should we stop operation?
  }
  return_list=list(degree_RegROTS, degree_PolyReg)
  names(return_list)=c("Degree_RegROTS", "Degree_PolyReg")
  return(return_list)
}


