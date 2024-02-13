
library("data.table")
library("dplyr")
library("tidyr")
library("stringr")
library("ggplot2")
library("gggenes")
library("ggtext")
library("glue")
library("readxl")
library("ggnewscale")
library("here")
`%ni%` <- Negate(`%in%`)

# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if an argument was passed
if (length(args) > 0) {
  # Set the working directory to the argument passed
  setwd(args[1])
  if (args[2] == "TRUE") {
    IPS = TRUE
  }
  if (args[2] != "TRUE") {
    IPS = FALSE
  }

}

##############################################
# Plotting of operons from results
# Remember to run ips before this
##############################################
source("./scripts/plot_operon.R")

for (file in list.dirs("./data/output_proximity_filtration/fasta_output", full.names=FALSE)) {
    if (file %ni% c("", "cellulose_All")) {
    plot_operon(file, article_plot_domain = IPS)
    }

    if (file %in% c("cellulose_All")) {
        plot_operon("cellulose", same_database = file, article_plot_domain = IPS)
    }
}

