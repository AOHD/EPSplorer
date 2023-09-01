.libPaths( c( "/user_data/ahd/r_v.4.0.0", .libPaths() ) )

# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if an argument was passed
if (length(args) > 0) {
  # Set the working directory to the argument passed
  setwd(args[1])
}

`%ni%` <- Negate(`%in%`)

##############################################
# Plotting of operons from results
# Remember to run ips before this
##############################################
source("./scripts/plot_operon.R")

for (file in list.dirs("./data/output_proximity_filtration/fasta_output", full.names=FALSE)) {
    if (file %ni% c("", "cellulose", "celluloseI", "celluloseII", "celluloseIII", "cellulose_Ac", "cellulose_NA")) {
    plot_operon(file)
    }

    if (file %in% c("cellulose", "celluloseI", "celluloseII", "celluloseIII", "cellulose_Ac", "cellulose_NA")) {
        plot_operon("cellulose", same_database = file)
    }
}

