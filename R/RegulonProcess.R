##' Transfer data-frame regulon data to list 
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
    return(clusterGene)
  })

  return(clusterList)
} 


##' Transfer data-frame regulor data to list.
##'
##' The regulon data file can be stored in a .txt or .csv file, which are read in and stored as data-frame. The row names are used as regulors' name, and the same regulors are arranged in into one column. There are three types of regulors, "activator", "repressor", and "notsure". The regulator data is sperated by "|", for example "repressor|SMU_1995c"
##' @title Data-frame to regulons list
##' @param dfFile a data-frame regulors
##' @return A list of regulors
##' @examples
##' data(smuRegulors)
##' smuRegulors <- df2regulon(smuRegulors)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
df2regulor <- function(dfFile){
  
  clusterList <- lapply(dfFile, function(x) {
    clusterGene <- as.character(x)
    # remove length == 0
    clusterGene <- clusterGene[!(clusterGene == '')]

    clusterMat <- do.call(rbind, strsplit(clusterGene, split = '|', fixed = TRUE))
    return(clusterMat)
  })

  return(clusterList)

}


