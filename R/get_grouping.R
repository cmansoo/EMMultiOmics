#' Get Gene Groupings from Functional classification text file
#' @export
get_grouping <- function(full_gene_list, file_dir){
  K0 = length(full_gene_list)
  functionalClassification <- read.delim(file_dir, header=FALSE, comment.char="#") 
  breakPoint = which(startsWith(as.character(functionalClassification[,1]),'Gene Group'))
  gene_names = unique(functionalClassification[-c(breakPoint,breakPoint+1),1])
  grouping = matrix(0,length(gene_names),length(breakPoint))
  rownames(grouping) = gene_names
  for (i in 1:(length(breakPoint)-1)){
    genes = functionalClassification[(breakPoint[i]+2):(breakPoint[i+1]-1),1]
    grouping[as.character(genes),i]=1
  }
  genes = functionalClassification[(tail(breakPoint,1)+2):nrow(functionalClassification),1]
  grouping[as.character(genes),i+1] = 1
  R = ncol(grouping)+1
  gene_group = rep(R,K0)
  for (k in 1:K0){
    genename = full_gene_list[k]
    if (genename %in% rownames(grouping)){
      gene_group[k] = which(grouping[genename,]>0)[1]
    }
  }
  names(gene_group) = full_gene_list
  return(gene_group)
}