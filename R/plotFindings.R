#' Plot RolDE results
#'
#' Plot the findings from longitudinal differential expression analysis with RolDE.
#'
#' @param file_name a string indicating the file name in which the results should be plotted. Should have a ".pdf" extension. Default is NULL, no file is created.
#' @param RolDE_res the RolDE result object.
#' @param top_n an integer or a vector of integers indicating what top differentially expressed features should be plotted. If \code{top_n} is a single number, the \code{top_n} most
#' differentially expressed feature will be plotted (e.g \code{top_n}=1 will plot the most differentially expressed feature). If \code{top_n} is a vector of numbers,
#' the differentially expressed features corresponding to top detections within the given range will be plotted (e.g. \code{top_n}=seq(1:50) will plot the top 50 differentially expressed features).
#' If more than one feature will be plotted, it is advisable to define a suitable file name in \code{file_name}.
#' @param col1 a string indicating which color should be used for replicates / individualsin in condition 1. The default is blue.
#' @param col2 a string indicating which color should be used for replicates / individualsin in condition 2. The default is red.
#' @return \code{plotFindings} Plots the results from the RolDE object.

#' @details The function plots the longitudinal expression of the top RolDE findings. The function can plot either the expression of a single finding
#' or multiple top findings as indicated by the \code{top_n}. The findings can be plotted into a pdf file as indicated by the \code{file_name}.
#' The given \code{file_name} should have a ".pdf" extension.
#' @importFrom graphics legend lines mtext plot points
#' @export
#' @examples
#' data("res3")
#'#Plotting the most DE finding. DE results are in the res3 object.
#'plotFindings(file_name = NULL, RolDE_res = res3, top_n = 1)

plotFindings<-function(file_name=NULL, RolDE_res, top_n, col1="blue", col2="red"){

  #Data is in the Rolde object
  results=RolDE_res$RolDE_Results
  results=results[order(as.numeric(results$`RolDE Rank Product`)),]
  inputs=RolDE_res$Input
  data_plot=inputs$data
  des_plot=inputs$des_matrix

  des_plot=data.frame(des_plot, stringsAsFactors = FALSE)
  colnames(des_plot)=c("Sample", "Condition", "Time", "Individual")
  des_plot[,1]<-as.character(des_plot[,1])
  des_plot[,2]<-as.factor(as.character(des_plot[,2]))
  des_plot[,3]<-as.numeric(as.character(des_plot[,3]))
  des_plot[,4]<-as.numeric(as.character(des_plot[,4]))

  #Now, we want to plot the longitudinal expression of indivduals in condition 1
  #and the longitudinal expression of individual in condition 2
  uniq_conds=unique(des_plot$Condition)
  locs_cond1=which(des_plot$Condition==uniq_conds[1])
  locs_cond2=which(des_plot$Condition==uniq_conds[2])

  des_plot_cond1=des_plot[locs_cond1,]
  des_plot_cond2=des_plot[locs_cond2,]

  features_to_plot=results$`Feature ID`[top_n]

  range_time=range(as.numeric(as.character(des_plot$Time)), na.rm = TRUE)
  range_time=c((range_time[1]-range_time[1]*0.1), (range_time[2]+range_time[1]*0.1))
  #Need to be adjusted?

  #plot each top feature
  if(!is.null(file_name)){
    pdf(file=file_name, width = 11.7, height = 8.3,onefile = TRUE)
  }

  for(j in seq_len(length(features_to_plot))){
    range_exprs=range(data_plot[features_to_plot[j],c(des_plot$Sample)], na.rm = TRUE)
    range_exprs=c((range_exprs[1]-range_exprs[1]*0.1), (range_exprs[2]+range_exprs[1]*0.1))
    plot_name=features_to_plot[j]
    de_info=paste("RolDE rank product: ", round(results$`RolDE Rank Product`[j],3), ", estimated significance value (ESV): ", round(results$`Estimated Significance Value`[j],3), ", multiple hypothesis adjusted ESV: ", round(results$`Adjusted Estimated Significance Value`[j],3), sep = "")
    #Plot by conditions
    #Condition 1
    uniq_inds_cond1=unique(des_plot_cond1$Replicate)
    for(i in seq_len(length(uniq_inds_cond1))){
      inds_samp=which(des_plot_cond1$Replicate==uniq_inds_cond1[i])
      temp_des=des_plot_cond1[inds_samp,]
      #order by time
      temp_des=temp_des[order(as.numeric(temp_des$Time)),]

      x=as.numeric(data_plot[features_to_plot[j], c(temp_des$Sample)])
      pchs=rep(19, length(x))
      col.pchs=rep(col1, length(x))
      if(any(is.na(x))){
        pchs[which(is.na(x))]=1
        col.pchs[which(is.na(x))]="black"
        x[which(is.na(x))]=mean(x, na.rm=TRUE)
      }
      y=as.numeric(as.character(temp_des$Time))

      #plot
      if(i==1){
        plot(y, x, xlim = range_time, ylim = range_exprs, pch=pchs, col=col.pchs, ylab = "Expression", xlab="Time / Longitudinal variable" , main = plot_name)
        lines(y,x, col=col1)
        mtext(3, text = de_info)
        legend("bottomright", legend = c(uniq_conds[1], uniq_conds[2], "Missing value - replicate / individua mean imputed for plot"), col=c(col1, col2, 1), pch = c(19,19,1), bty = "n")
      }else{
        points(y,x, pch=pchs, col=col.pchs)
        lines(y,x, col=col1)
      }
    }

    #Condition2
    uniq_inds_cond2=unique(des_plot_cond2$Replicate)
    for(i in seq_len(length(uniq_inds_cond2))){
      inds_samp=which(des_plot_cond2$Replicate==uniq_inds_cond2[i])
      temp_des=des_plot_cond2[inds_samp,]
      #order by time
      temp_des=temp_des[order(as.numeric(temp_des$Time)),]

      x=as.numeric(data_plot[features_to_plot[j], c(temp_des$Sample)])
      pchs=rep(19, length(x))
      col.pchs=rep(col2, length(x))
      if(any(is.na(x))){
        pchs[which(is.na(x))]=1
        col.pchs[which(is.na(x))]="black"
        x[which(is.na(x))]=mean(x, na.rm=TRUE)
      }
      y=as.numeric(as.character(temp_des$Time))
      points(y,x, pch=pchs, col=col.pchs)
      lines(y,x, col=col2)
    }
  }

  if(!is.null(file_name)){
    dev.off()
  }
}

