#' g
#' 
#' 
#' 
#' 
#' 
get_grouping <- function(){
  
}



URLencode()
u1 <- "https://davidbioinformatics.nih.gov/api.jsp?type=ENTREZ_GENE_ID&ids=2919,6347,6348,6364&tool=gene2gene"
u1 |> URLencode()
u1_encoded <- URLencode(u1)


resp <- httr::GET(u1_encoded)
resp$status_code
resp$content |> rawToChar()

u2 <- "https://davidbioinformatics.nih.gov/annotationReport.jsp" |> URLencode()
resp <- httr::GET(u2)


u3 <- "https://davidbioinformatics.nih.gov/data/download/tr_F665BB82C3C71726015768461.txt" |> URLencode()
resp <- httr::GET(u3)
resp$content |> rawToChar()



# e.g. URL
# https://davidbioinformatics.nih.gov/api.jsp?type=xxxxx&ids=XXXXX,XXXXX,XXXXXX,&tool=xxxx&annot=xxxxx,xxxxxx,xxxxx,
# example ID's (we want GENE SYMBOL)
# read from geneNames.txt
# use is supposed to provide this geneNames
list.files("ignore_dev")
gene_names <- read.table("ignore_dev/geneNames.txt") |> unlist(use.names=FALSE)
gene_names
# collapse into string
# the API has limitations:

# DAVID APIs (alpha) allow other bioinformatics web sites to directly link to DAVID tools and functions ONLY for light-duty jobs (i.e. a gene list with no more than 400 genes).
# DAVID APIs (alpha) are not for high-throughput or large gene list jobs, such as: a job for a gene list with more than 500 genes; try to loop DAVID data through scripts for hundreds/thousands of gene lists. For big jobs like above, you need to do them in a regular way through manual submission form on DAVID web site, or to download DAVID Knowledgebase to setup in-house analysis engines, or to contact DAVID team for alternative automatic solutions to meet you specific situations.
# URL has a charactor size limitation (<= 2048 characters in total), i.e., the very large gene list may not be able to completely passed by URL.
# No more than 200 hits in a day from one computer.
# 10 seconds interval between hits.
# DAVID Team reserves right to suspend any improper uses of DAVID APIs without notice.


# we need to break them into vectors, every 400 genes
gene_names |> matrix(nrow=400) |> View()
# if we don't have enough records to break it down, fill the last column with NAs?
extra_len <- length(gene_names) %% 400
gene_names_mat <- append(gene_names, rep(NA, extra_len)) |> matrix(nrow=400)
gene_names_mat

# repeat for every column in the matrix?

# try with the whole list and see what errors we get?
# gene_names_new <- gene_names |> paste0(collapse=",")
# try first 400
gene_names_new <- gene_names_mat[,1] |> paste0(collapse=",")
# make even smaller, first 10
gene_names_new <- gene_names_mat[1:10, 1] |> paste0(collapse=",")
gene_names_new

# make URL
base_url <- "https://davidbioinformatics.nih.gov/api.jsp?"
my_type <- "type=OFFICIAL_GENE_SYMBOL"
my_ids <- paste0("&ids=", gene_names_new)
my_tool <- "&tool=gene2gene"
request_url <- paste0(base_url, my_type, my_ids, my_tool) |> URLencode()
request_url

# make request?
resp <- httr::GET(request_url)
resp$content |> rawToChar()

### API doesnt seem to work well because even small number of search gives a lot of irrelevant genes
### and there is no way to specify species as homo sapiens
# 1. preferable is to use the manual download
# 2. try R webs service if possible?
# 3. try python web service
