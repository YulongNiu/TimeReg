
##' Plot time-course differential changes data of regulators and regulon.
##'
##' Use ggplot2 system to plot a cluster of regulators and regulon with 
##' differential changes data. The cluster can include multiple regulators (repressors or activators). 
##' 
##' @title Plot time-couse regulator and regulons
##' @param mergedRegData a data-frame, the first column is the regulator/regulon type, the second column is gene name, the third column is time points, and the fourth column is the differential changes.
##' @param sigCutoff a numeric vector with length 2 to define the low and high threshold of significance, for example c(-1, 1)
##' @param maxDiff the maximum difference
##' @param minDiff the minimum difference
##' @param hasRegulator whether to plot regulors, default value is TRUE. If FALSE is set, it can be also used to plot transcription units. 
##' @return a ggplot object 
##' @examples
##' data(smuMergedReg)
##' data(smuHeatFC)
##' maxFC <- max(smuHeatFC)
##' minFC <- min(smuHeatFC)
##' 
##' # plot Hrc regulons and regulators
##' HrcAReg <- smuMergedReg$SMU_80
##' RegTimePlot(HrcAReg, sigCutoff = c(log2(1.5), -log2(1.5)), maxDiff = maxFC, minDiff = minFC)
##'
##' # only plot CtsR regulons
##' CtsRReg <- smuMergedReg$SMU_2030
##' CtsRReg <- CtsRReg[CtsRReg$Type == 'regulon', ]
##' RegTimePlot(CtsRReg, sigCutoff = c(log2(1.5), -log2(1.5)), maxDiff = maxFC, minDiff = minFC, hasRegulator = FALSE)
##' 
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom ggplot2 ggplot aes scale_y_continuous geom_line facet_wrap geom_hline geom_rect scale_fill_manual
##' @export
##'
RegTimePlot <- function(mergedRegData, sigCutoff, maxDiff, minDiff, hasRegulator = TRUE){


  # coerce column names
  colnames(mergedRegData) <- c('Type', 'GeneName', 'TimePoint', 'Diff')
  
  # up- and down-threshold
  diffRect <- data.frame(SigType = c('over-expressed', 'under-expressed'), Start = sigCutoff, End = c(Inf, -Inf))

  # plot basic line image
  plotReg <- ggplot(mergedRegData, aes(x = TimePoint, y = Diff, group = GeneName, colour = GeneName)) +
    scale_y_continuous(limits = c(minDiff, maxDiff)) +
      geom_line() +
        geom_hline(yintercept = sigCutoff, colour = c('darkred', 'darkgreen'), alpha = 0.4, linetype = 'longdash') +
          geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = Start, ymax = End, fill = SigType), data = diffRect, inherit.aes = FALSE, alpha = 0.1) +
            scale_fill_manual(values = c('red', 'green'))

  if (hasRegulator) {
    plotReg <- plotReg + facet_wrap(~Type)
  } else {}

  return(plotReg)
  
}


##' Streptococcus mutans heat stress time-course data of regulators and regulons
##'
##' A list containing the S.mutans time-course fold change data of regulators and regulons.
##' The S.mutans were transferred from 37C to 42C for 5 minutes, 10 minutes, 15 minutes, 30 minutes, 45 minutes, and 60 minutes. Fold changes were calculated by comparing the heat stress probe indensity with 37C.
##' Regulator and regulons data were from RegPrecise database
##' @docType data
##' @name smuMergedReg
##' @format A data frame
##' @source GEO Database GSE59302
##' @source \url{http://regprecise.lbl.gov/RegPrecise/genome.jsp?genome_id=37}
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' 
NULL
