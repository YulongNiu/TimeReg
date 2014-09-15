##' Transfer data-frame regulon data to list.
##'
##' The regulon data file can be stored in a .txt or .csv file, which are read in and stored as data-frame. The row names are used as regulons' name, and the same regulons are arranged in into one column.
##' @title Data-frame to regulons list
##' @param dfFile a data-frame regulons
##' @return A list of regulons
##' @examples
##' data(smuRegulons)
##' smuRegulons <- df2regulon(smuRegulons)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
df2regulon <- function(dfFile) {
  
  clusterList <- lapply(dfFile, function(x) {
    clusterGene <- as.character(x)
    # remove length == 0
    clusterGene <- clusterGene[!(clusterGene == '')]
    clusterMat <- cbind(rep('regulon', length(clusterGene)), clusterGene)
    colnames(clusterMat) <- NULL
    
    return(clusterMat)
  })

  return(clusterList)
} 


##' Transfer data-frame regulor data to list.
##'
##' The regulon data file can be stored in a .txt or .csv file, which are read in and stored as data-frame. The row names are used as regulors' name, and the same regulors are arranged in into one column. There are three types of regulors, "activator", "repressor", and "regulatorFNS". "regulatorFNS" represents regulators whose function is not sure, maybe activator or repressor. The regulator data is separated by "|", for example "repressor|SMU_1995c".
##' @title Data-frame to regulons list
##' @param dfFile a data-frame regulors
##' @return A list of regulors
##' @examples
##' data(smuRegulators)
##' smuRegulators <- df2regulator(smuRegulators)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
df2regulator <- function(dfFile){
  
  clusterList <- lapply(dfFile, function(x) {
    clusterGene <- as.character(x)
    # remove length == 0
    clusterGene <- clusterGene[!(clusterGene == '')]

    clusterMat <- do.call(rbind, strsplit(clusterGene, split = '|', fixed = TRUE))
    return(clusterMat)
  })

  return(clusterList)

}



##' Transfer list to proper list data
##'
##' Transcription unit data are commonly stored in a list. This function convert the TU data to proper list format.
##' @title List to TU list
##' @param lsFile TU list data. In each TU, a character vector contains all the genes in this TU.
##' @return A list of TU
##' @examples
##' data(smuTU)
##' smuTU <- ls2TU(smuTU)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
ls2TU <- function(lsFile) {

  TUList <- lapply(lsFile, function(x) {
    len <- length(x)
    x <- cbind(rep('strucgene', len), x)
    colnames(x) <- NULL
    return(x)
  })

  return(TUList)
}




##' Melt regulator data, regulon data, and the degree of differential changes into a list.
##'
##' For a time-course experiment, we want to see the gene expression patterns of the regulators/regulons. "MeltReg" function provides a way to melt regulator, regulon, and degree of changes. The corresponding regulator and regulon data are recognized by the name.
##' @title Melt regulator, regulon, and degree of changes data
##' @param regulatorData regulator data.
##' @param regulonData regulon data.
##' @param diffData degree of differential changes data, which can be fold change data, q-PCR data. The row names of "diffData" should be in the same fomat as the second column in "regulatorData" and "regulonData".
##' @return a list of merged data.
##' @examples
##' data(smuRegulators)
##' smuRegulators <- df2regulator(smuRegulators)
##' data(smuRegulons)
##' smuRegulons <- df2regulon(smuRegulons)
##' data(smuHeatFCb)
##' smuMergedReg <- MeltReg(smuRegulators, smuRegulons, smuHeatFC)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom reshape2 melt
##' @export
##' 
MeltReg <- function(regulatorData, regulonData, diffData) {

  # common cluster names
  regulatorNames <- names(regulatorData)
  regulonNames <- names(regulonData)
  comClustNames <- intersect(regulatorNames, regulonNames)
  comRegulator <- regulatorData[match(comClustNames, regulatorNames)]
  comRegulon <- regulonData[match(comClustNames, regulonNames)]

  # merge common regulators, regulons, and diffData
  mergedReg <- list()
  for(i in 1:length(comClustNames)) {
    mergedReg[[i]] <- rbind(comRegulator[[i]], comRegulon[[i]])
    geneNames <- mergedReg[[i]][, 2]
    diffReg <- diffData[match(geneNames, rownames(diffData)), ]
    mergedReg[[i]] <- cbind(mergedReg[[i]], diffReg)
    colnames(mergedReg[[i]])[1:2] <- c('Type', 'GeneName')
    
    # melt data for ggplot2
    mergedReg[[i]] <- melt(mergedReg[[i]], id.vars = c('Type', 'GeneName'),
                           variable.name = 'TimePoint',
                           value.name = 'Diff')
  }

  names(mergedReg) <- comClustNames

  
  return(mergedReg)
}



##' Melt transcriptional units data and the degree of differential changes into a list.
##'
##' In prokarytic genomes, the gene expression patterns are closely related to transcrional units (TU). "mergeTU" provides a way to merge TU and degree of changes.
##' @title Melt TU and degree of changes data
##' @param TUData formatted TU list data. The structural genes in TU should be a subset of row names of "diffData". If not, NA will be returned.
##' @inheritParams MeltReg
##' @return A list merged data
##' @examples
##' data(smuTU)
##' smuTU <- ls2TU(smuTU)
##' mergedTU <- MeltTU(smuTU, smuHeatFC)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom reshape2 melt
##' @export
##' 
MeltTU <- function(TUData, diffData) {

  mergedTU <- lapply(TUData, function(x) {
    
    diffReg <- diffData[match(x[, 2], rownames(diffData)), ]
    mergedTUEach <- cbind(x, diffReg)
    colnames(mergedTUEach)[1:2] <- c('Type', 'GeneName')

    # melt data for ggplot2
    meltTUEach <- melt(mergedTUEach, id.vars = c('Type', 'GeneName'),
                           variable.name = 'TimePoint',
                           value.name = 'Diff')

    return(meltTUEach)
    
  })

  return(mergedTU)
}







##' Streptococcus mutans regulons.
##'
##' A data-frame containing the S.mutans regulons, and the dataset is read from
##' an csv file by read.csv(csvFilePath, header = TRUE). 
##' Each column is the same regulons which maybe regulated by several regulators. 
##' The first row is the regulons' name.
##' @docType data
##' @name smuRegulons
##' @format A data frame
##' @source \url{http://regprecise.lbl.gov/RegPrecise/genome.jsp?genome_id=37}
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' 
NULL


##' Streptococcus mutans regulators.
##'
##' A data-frame containing the S.mutans regulators, and the dataset is read from
##' an csv file by read.csv(csvFilePath, header = TRUE). 
##' Each column is the cluster of regulators which may regulate a cluster of regulons. 
##' The first row is the name of regulator clusters.
##' @docType data
##' @name smuRegulators
##' @format A data frame
##' @source \url{http://regprecise.lbl.gov/RegPrecise/genome.jsp?genome_id=37}
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' 
NULL



##' Streptococcus mutans heat stress time-course microarray data
##'
##' A data-frame containing the S.mutans time-course microarray data. The S.mutans
##' were transferred from 37C to 42C for 5 minutes, 10 minutes, 15 minutes, 30 minutes, 45 minutes, and 60 minutes. Fold changes were calculated by comparing the heat stress probe indensity with 37C. 
##' @docType data
##' @name smuHeatFC
##' @format A data frame
##' @source GEO Database GSE59302
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' 
NULL



##' Streptococcus mutans transcriptional units data
##'
##' A list containing the S.mutans transcriptional units data from "Prokaryotic Operon DataBase". 
##' @docType data
##' @name smuTU
##' @format A list
##' @source \url{http://operons.ibt.unam.mx/OperonPredictor/}
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' 
NULL


