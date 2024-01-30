
library("tidyr")
library("dplyr")
library("stringr")
library("purrr")

# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if an argument was passed
if (length(args) > 0) {
  # Set the working directory to the argument passed
  setwd(args[1])
}

#-----------------------------------------------------------------
#        Defining the gffRead function for loading gff data       
#-----------------------------------------------------------------
gffRead <- function(gffFile){
  gff <- read.delim(gffFile, header=F, comment.char="#")
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  gff <- gff[complete.cases(gff), ]
  
  gff <- gff %>%
    filter(feature == "CDS") %>%
    separate(attributes, c("prokkaID", "attributes"), sep = ";", extra = "merge") %>%
    separate(prokkaID, c("ID1", "ProkkaNO"),
             sep = "_(?!.*_)")
  
  return(gff)}

#-----------------------------------------------------------------
#          Using the gffRead function to import gff files         
#-----------------------------------------------------------------
gff <- tibble(
  file = list.files("./data/prokka/", recursive = TRUE, pattern = ".gff"),
  location = list.files("./data/prokka/", full.names = TRUE, recursive = TRUE, pattern = ".gff"),
  ID = str_remove(file, ".gff"),
  df = map(location, gffRead)
  ) %>% 
  unnest(df)

#-----------------------------------------------------------------
#  Removing columns redundant to psiblast data in gff dataframe   
#-----------------------------------------------------------------
gff <-  gff %>% mutate(ID = substr(ID1, 4, nchar(ID1))) %>% 
  select(!c("score", "feature", "file", "location", 
                          "frame", "attributes",  "source", "ID1"))

saveRDS(gff, file="./data/gff.rds")


