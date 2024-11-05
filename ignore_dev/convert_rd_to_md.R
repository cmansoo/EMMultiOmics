library(Rd2md)
# setwd("C:/Users/mc2614/Desktop/mcho/EMMultiOmics/")
# rdf <- read_rdfile(path = "man/multiOmics.Rd")
# as_markdown(rdf)
# Rd2markdown("man/multiOmics.Rd", "test.md")
rd_files <- list.files("man")

md_files <- gsub("Rd$", "md", rd_files)


mapply(function(r, m) {
  r <- paste0("man/", r)
  m <- paste0("ignore_dev/md/", m)
  message("r: ", r)
  message("m: ", m)
  Rd2markdown(r, m)
},
rd_files, md_files)
