# library(cBioPortalData)
# > library(cBioPortalData)
# Loading required package: AnVIL
# Loading required package: dplyr
# 
# Attaching package: ‘dplyr’
# 
# The following objects are masked from ‘package:stats’:
#   
#   filter, lag
# 
# The following objects are masked from ‘package:base’:
#   
#   intersect, setdiff, setequal, union
# 
# Loading required package: MultiAssayExperiment
# Loading required package: SummarizedExperiment
# Loading required package: MatrixGenerics
# Loading required package: matrixStats
# 
# Attaching package: ‘matrixStats’
# 
# The following object is masked from ‘package:dplyr’:
#   
#   count
# 
# 
# Attaching package: ‘MatrixGenerics’
# 
# The following objects are masked from ‘package:matrixStats’:
#   
#   colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse, colCounts, colCummaxs,
# colCummins, colCumprods, colCumsums, colDiffs, colIQRDiffs, colIQRs, colLogSumExps,
# colMadDiffs, colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
# colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds, colSums2,
# colTabulates, colVarDiffs, colVars, colWeightedMads, colWeightedMeans,
# colWeightedMedians, colWeightedSds, colWeightedVars, rowAlls, rowAnyNAs, rowAnys,
# rowAvgsPerColSet, rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
# rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps, rowMadDiffs, rowMads,
# rowMaxs, rowMeans2, rowMedians, rowMins, rowOrderStats, rowProds, rowQuantiles,
# rowRanges, rowRanks, rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs,
# rowVars, rowWeightedMads, rowWeightedMeans, rowWeightedMedians, rowWeightedSds,
# rowWeightedVars
# 
# Loading required package: GenomicRanges
# Loading required package: stats4
# Loading required package: BiocGenerics
# 
# Attaching package: ‘BiocGenerics’
# 
# The following objects are masked from ‘package:dplyr’:
#   
#   combine, intersect, setdiff, union
# 
# The following objects are masked from ‘package:stats’:
#   
#   IQR, mad, sd, var, xtabs
# 
# The following objects are masked from ‘package:base’:
#   
#   anyDuplicated, aperm, append, as.data.frame, basename, cbind, colnames, dirname,
# do.call, duplicated, eval, evalq, Filter, Find, get, grep, grepl, intersect,
# is.unsorted, lapply, Map, mapply, match, mget, order, paste, pmax, pmax.int, pmin,
# pmin.int, Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort, table,
# tapply, union, unique, unsplit, which.max, which.min
# 
# Loading required package: S4Vectors
# 
# Attaching package: ‘S4Vectors’
# 
# The following objects are masked from ‘package:dplyr’:
#   
#   first, rename
# 
# The following object is masked from ‘package:utils’:
#   
#   findMatches
# 
# The following objects are masked from ‘package:base’:
#   
#   expand.grid, I, unname
# 
# Loading required package: IRanges
# 
# Attaching package: ‘IRanges’
# 
# The following objects are masked from ‘package:dplyr’:
#   
#   collapse, desc, slice
# 
# The following object is masked from ‘package:grDevices’:
#   
#   windows
# 
# Loading required package: GenomeInfoDb
# Loading required package: Biobase
# Welcome to Bioconductor
# 
# Vignettes contain introductory material; view with 'browseVignettes()'. To cite
# Bioconductor, see 'citation("Biobase")', and for packages 'citation("pkgname")'.
# 
# 
# Attaching package: ‘Biobase’
# 
# The following object is masked from ‘package:MatrixGenerics’:
#   
#   rowMedians
# 
# The following objects are masked from ‘package:matrixStats’:
#   
#   anyMissing, rowMedians
# 
# Warning messages:
#   1: package ‘AnVIL’ was built under R version 4.3.3 
# 2: package ‘dplyr’ was built under R version 4.3.3 
# 3: package ‘matrixStats’ was built under R version 4.3.3 
# 4: package ‘GenomeInfoDb’ was built under R version 4.3.3 
# library(survival)


# get data
mydat <- cBioPortalData::cBioDataPack("gbm_tcga_pub")

# expression_med <- mydat[[4L]] # what is this fourth element of the list?
expression_med <- mydat[["methylation_hm27"]]
#  [4] methylation_hm27: SummarizedExperiment with 6 rows and 58 columns <- this is it?

# get expression vector
# SummarizedExperiment::assay(expression_med) |> as.data.frame() |> View()
exp_med <- SummarizedExperiment::assay(expression_med)
exp_med <- exp_med[complete.cases(exp_med),] # it went from matrix to flattened vector 9020 *58, are they just squentially appended then? i.e. c(col1, col2, ...)

exp_med_df <- data.frame(SummarizedExperiment::assay(expression_med)) # for demo
exp_med_df <- exp_med_df[stats::complete.cases(exp_med_df), ] |> unlist() # for demo

#get miRNA profile
# expression_miRNA = mydat[[6L]]
expression_miRNA <- mydat[["mirna_zscores"]]
# expression_miRNA@colData
# expression_miRNA@assays
# ?SummarizedExperiment::assay
exp_miRNA <- SummarizedExperiment::assay(expression_miRNA)
exp_miRNA <- exp_miRNA[stats::complete.cases(exp_miRNA),]
# exp_miRNA

#get methylation profile
mydat@ExperimentList
mydat@metadata
expression_methy <- mydat[["mrna_agilent_microarray_zscores_ref_diploid_samples"]]
exp_methy <- SummarizedExperiment::assay(expression_methy)
exp_methy <- 
  exp_methy[stats::complete.cases(exp_methy), 
            base::setdiff(colnames(exp_methy), c("TCGA-02-0001-01", "TCGA-02-0003-01"))]

# expression_med@NAMES

# init cbioportal API
# cBioPortalData::cBioPortal()
# equivalent
cbio <- cBioPortalData::cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/v2/api-docs",
  token = character()
)
# cbio@api
# cbio@api_header

# get molecular profile
mols <- cBioPortalData::molecularProfiles(
  api=cbio,
  studyId="gbm_tcga_pub"
)
# mols[["molecularProfileId"]]


# get clinical data
clinicalDat <- cBioPortalData::clinicalData(cbio, "gbm_tcga_pub")
clinicalDat2013 <- cBioPortalData::clinicalData(cbio,"gbm_tcga_pub2013")
# copy age
clinicalDat <- merge(
  x=clinicalDat,
  y=clinicalDat2013[, c("patientId", "AGE")],
  by="patientId",
  all.x=TRUE
) 

clinicalDF <- clinicalDat[, c("patientId", "AGE","OS_MONTHS","OS_STATUS","PRIOR_GLIOMA",
                              "SEX","KARNOFSKY_PERFORMANCE_SCORE","PRETREATMENT_HISTORY")]

# extract desired patient records
clinicalIDs <- clinicalDat$patientId[complete.cases(clinicalDF)]
myIDs <- intersect(paste0(clinicalIDs,'-01'), colnames(exp_med))
finalIDs <- intersect(myIDs, colnames(exp_methy))
myexp_med <- exp_med[, finalIDs]
# myexp_miRNA <- exp_miRNA[, finalIDs] 
myexp_methy <- exp_methy[, finalIDs]
myclinical <- clinicalDF[paste0(clinicalDF$patientId,'-01') %in% finalIDs, ]

# format columns
myclinical <- transform(
  myclinical,
  AGE = as.numeric(AGE),
  OS_MONTHS = as.numeric(OS_MONTHS),
  OS_STATUS = ifelse(OS_STATUS == "1:DECEASED", 1, 0),
  PRIOR_GLIOMA = ifelse(PRIOR_GLIOMA=="NO", 0, 1),
  SEX = ifelse(SEX == "Male", 0 , 1),
  KARNOFSKY_PERFORMANCE_SCORE = as.numeric(KARNOFSKY_PERFORMANCE_SCORE),
  PRETREATMENT_HISTORY = ifelse(PRETREATMENT_HISTORY == "NO", 0, 1)
)
Y <- myclinical$OS_MONTHS
Delta <- myclinical$OS_STATUS

# initialize vector for survival
p_surv <- matrix(nrow=nrow(myexp_med), ncol=1)
p_surv <- sapply(
  1:nrow(myexp_med),
  function(i){
    survregWeibull <- survival::survreg(survival::Surv(OS_MONTHS, OS_STATUS) ~
                                          AGE + myexp_med[i, ], 
                                        data=myclinical,
                                        dist="weibull")
    s <- summary(survregWeibull)
    # print(row.names(s$table))
    return(s$table["myexp_med[i, ]", "p"])
  }
)

G <- t(myexp_med[order(p_surv)[1:1000], ])
C <- as.matrix(myclinical[, c("AGE", "PRIOR_GLIOMA", "SEX", "PRETREATMENT_HISTORY")])



###### PICKUP FROM LINE 312 #use AFT model to select 1000 genes having strongest association with survival
# M.raw giving nothing? -> some how Hao got it to work and we have GBM data file -> use that for now
#use AFT model to select 1000 genes having strongest association with survival
C = as.matrix(myclinical[,c('AGE','PRIOR_GLIOMA','SEX','PRETREATMENT_HISTORY')])
M.raw=t(myexp_methy[colnames(G)[colnames(G)%in%rownames(myexp_methy)],])
N=nrow(M.raw)
if (d==3){
  M = matrix(0,nrow=N,ncol = 3*ncol(M.raw))
  for (j in 1:ncol(M.raw)) {
    M[,(3*(j-1)+1):(3*j)]=ns(M.raw[,j],df=3)
  }
}else if (d==2){
  M = matrix(0,nrow=N,ncol = 2*ncol(M.raw))
  for (j in 1:ncol(M.raw)) {
    M[,(2*(j-1)+1)]=M.raw[,j]
    M[,(2*j)]=M.raw[,j]^2
  }
}else{
  M=M.raw}



# orig fx
get_GBM_data = function(d=1){
  if (!d%in%c(1,2,3)){
    stop('degree of mechanistic model should be 1, 2 or 3')
  }
  mydat = cBioDataPack("gbm_tcga_pub")
  #preprocess data from cbioportal
  #caveat: the data are updated irregularly, modification might be needed
  #get expression
  expression_med = mydat[[4L]]
  exp_med = assay(expression_med)
  exp_med = exp_med[complete.cases(exp_med),]
  #get miRNA profile
  expression_miRNA = mydat[[6L]]
  exp_miRNA = assay(expression_miRNA)
  exp_miRNA = exp_miRNA[complete.cases(exp_miRNA),]
  #get methylation profile
  expression_methy = mydat[[9L]]
  exp_methy = assay(expression_methy)
  exp_methy = exp_methy[complete.cases(exp_methy),-c(1,2)]
  #expression_med@NAMES
  (cbio <- cBioPortal())
  mols <- molecularProfiles(cbio, "gbm_tcga_pub")
  #mols[["molecularProfileId"]]
  #get clinical data
  clinicalDat = clinicalData(cbio,"gbm_tcga_pub")
  clinicalDat2013 = clinicalData(cbio,"gbm_tcga_pub2013")
  clinicalDat$AGE = clinicalDat2013[clinicalDat2013$patientId%in%clinicalDat$patientId,]$AGE
  clinicalDF = data.frame(clinicalDat[,c("AGE","OS_MONTHS","OS_STATUS","PRIOR_GLIOMA",
                                         "SEX","KARNOFSKY_PERFORMANCE_SCORE","PRETREATMENT_HISTORY")])
  clinicalDF$AGE = as.vector(clinicalDat$AGE)
  rownames(clinicalDF) = clinicalDat$patientId
  clinicalIDs = clinicalDat$patientId[complete.cases(clinicalDF)]
  myIDs = intersect(paste0(clinicalIDs,'-01'),colnames(exp_med))
  finalIDs = intersect(myIDs, colnames(exp_methy))
  myexp_med = exp_med[,finalIDs]
  #myexp_miRNA = exp_miRNA[,finalIDs]
  myexp_methy = exp_methy[,finalIDs]
  myclinical = clinicalDF[paste0(rownames(clinicalDF),'-01')%in%finalIDs,]
  myclinical$AGE = as.numeric(myclinical$AGE)
  myclinical$OS_MONTHS = as.numeric(myclinical$OS_MONTHS)
  myclinical$OS_STATUS = ifelse(myclinical$OS_STATUS=="1:DECEASED",1,0)
  myclinical$PRIOR_GLIOMA = ifelse(myclinical$PRIOR_GLIOMA=="NO",0,1)
  myclinical$SEX = ifelse(myclinical$SEX=="Male",0,1)
  myclinical$KARNOFSKY_PERFORMANCE_SCORE = as.numeric(myclinical$KARNOFSKY_PERFORMANCE_SCORE)
  myclinical$PRETREATMENT_HISTORY = ifelse(myclinical$PRETREATMENT_HISTORY=="NO",0,1)
  Y=myclinical$OS_MONTHS
  Delta = myclinical$OS_STATUS
  #use AFT model to select 1000 genes having strongest association with survival
  p.surv = matrix(nrow=nrow(myexp_med),ncol=1)
  for (i in 1:nrow(myexp_med)){
    survregWeibull <- survreg(Surv(OS_MONTHS, OS_STATUS) ~ AGE+myexp_med[i,],
                              myclinical, dist = "weibull")
    s = summary(survregWeibull)
    p.surv[i] = s$table[3,4]
  }
  G=t(myexp_med[order(p.surv)[1:1000],])
  C = as.matrix(myclinical[,c('AGE','PRIOR_GLIOMA','SEX','PRETREATMENT_HISTORY')])
  M.raw=t(myexp_methy[colnames(G)[colnames(G)%in%rownames(myexp_methy)],])
  N=nrow(M.raw)
  if (d==3){
    M = matrix(0,nrow=N,ncol = 3*ncol(M.raw))
    for (j in 1:ncol(M.raw)) {
      M[,(3*(j-1)+1):(3*j)]=ns(M.raw[,j],df=3)
    }
  }else if (d==2){
    M = matrix(0,nrow=N,ncol = 2*ncol(M.raw))
    for (j in 1:ncol(M.raw)) {
      M[,(2*(j-1)+1)]=M.raw[,j]
      M[,(2*j)]=M.raw[,j]^2
    }
  }else{
    M=M.raw}
  return(list(G=G,C=C,M=M,Y=Y,Delta=Delta))
}

get_grouping = function(full_gene_list,file_dir="functional_classification.txt"){
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

GBM_data = get_GBM_data(d=2)
gene_group = get_grouping(colnames(GBM_data$G), "functional_classification.txt")
save(GBM_data, gene_group, file='GBM_data2.RData')
