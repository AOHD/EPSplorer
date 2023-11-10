
# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if an argument was passed
if (length(args) > 0) {
  # Set the working directory to the argument passed
  setwd(args[1])
}

library("data.table")
library("readxl")
library("glue")
library("gggenes")
#library("ggtext")
library("broom")
library("seqinr")
library("dplyr")
library("tidyr")
#library("tidyverse")
`%ni%` <- Negate(`%in%`)

##----------------------------------------------------------------
##  Loading MAG statistic, Prokka annotations and query metadata  
##----------------------------------------------------------------
# Statistics on the HQ-MAGs from singleton et.al 2021
magstats. <- fread("./magstats.tsv", sep = "\t") %>% rename(ID = bin)

# File with metadata on the query genes, e.g. function
query_metadata. <- excel_sheets("./Query_figur.xlsx") %>%
  sapply(function(X) read_xlsx("./Query_figur.xlsx", sheet = X, skip = 1), USE.NAMES = T) %>% 
  lapply(as.data.frame) %>% 
  `[`(!(names(.) %in% c("Abbreveations", "HA_S_pyogenes"))) %>%
  rbindlist(fill = TRUE) %>%
  mutate(Polysaccheride = ifelse(Polysaccheride == "S88","s88",Polysaccheride),
         Psiblast = ifelse(Psiblast == "S88","s88",Psiblast),
         Function = ifelse(Function == "NA", "UNKNOWN", Function)) 

# Additional information of the prokka annotations
gff. <- readRDS("./data/gff.rds") %>%
  filter(!(ProkkaNO == "units")) %>% filter(!grepl("\\.|-", ProkkaNO)) %>%
  group_by(ID) %>% distinct(ProkkaNO, .keep_all = T) %>%
  mutate(Target_label = paste(ID, ProkkaNO, sep = "_"),
         ProkkaNO = as.numeric(ProkkaNO))


##############################################
# The main proximity filtration
##############################################
source("./scripts/proximity_filtration.R")
# Default: min_genes=2 and perc_id=20
proximity_filtration("alginate", 
                     min_genes = 6)
proximity_filtration("psl", 
                     min_genes = 7)
proximity_filtration("pel_merged",
                     min_genes = 6)
proximity_filtration("cellulose",
                     same_database = "celluloseI",
                     min_genes = 2,
                     essential_genes = "bcsD")
proximity_filtration("cellulose",
                     same_database = "celluloseII", 
                     min_genes = 2,
                     essential_genes = c("bcsE", "bcsG"))
proximity_filtration("cellulose", 
                     same_database = "celluloseIII",
                     min_genes = 2,
                     essential_genes = c("bcsK"))
proximity_filtration("cellulose", 
                     same_database = "cellulose_Ac",
                     min_genes = 2,
                     essential_genes = c("algF", "algI", "algJ", "algX"))
proximity_filtration("cellulose", 
                     same_database = "cellulose_NA",
                     min_genes = 2,
                     exclude_gene = c("algF", "algI", "algJ", "algX", "bcsK", "bcsE", "bcsG", "bcsD"))
proximity_filtration("succinoglycan",
                     min_genes = 10)
proximity_filtration("xanthan",
                     min_genes = 6)
proximity_filtration("curdlan",
                     min_genes = 2)
proximity_filtration("pnag_pga",
                     min_genes = 3)
proximity_filtration("pnag_ica",
                     min_genes = 3)
proximity_filtration("B_subtilis_EPS", 
                     min_genes = 3)
proximity_filtration("diutan",
                     min_genes = 10)
proximity_filtration("s88", 
                     min_genes = 9)
proximity_filtration("HA_Pasteurella",
                     min_genes = 5,
                     perc_id = 20)
proximity_filtration("HA_streptococcus", 
                     min_genes = 3, 
                     perc_id = 20)
proximity_filtration("acetan", 
                     min_genes = 8)
proximity_filtration("amylovoran", 
                     min_genes = 6)
proximity_filtration("ColA", 
                     min_genes = 9)
proximity_filtration("gellan1", 
                     min_genes = 2)
proximity_filtration("gellan2", 
                     min_genes = 9)
proximity_filtration("salecan", 
                     min_genes = 6)
proximity_filtration("rhizobium_eps", 
                     min_genes = 9)
proximity_filtration("stewartan", 
                     min_genes = 6)
proximity_filtration("vps", 
                     min_genes = 9)
proximity_filtration("burkholderia_eps", 
                     min_genes = 6)
proximity_filtration("levan", 
                     min_genes = 2)
proximity_filtration("synechan", 
                     min_genes = 10)
proximity_filtration("methanolan", 
                     min_genes = 11)
proximity_filtration("galactoglucan", 
                     min_genes = 9)

proximity_filtration("B_fragilis_PS_A", 
                     min_genes = 7)
proximity_filtration("B_fragilis_PS_B", 
                     min_genes = 10)
proximity_filtration("GG", 
                     min_genes = 8)
proximity_filtration("B_pseudomallei_EPS", 
                     min_genes = 10)
proximity_filtration("phosphonoglycan", 
                     min_genes = 10)
proximity_filtration("E_faecalis_PS", 
                     min_genes = 8)
proximity_filtration("glucorhamnan", 
                     min_genes = 10)
proximity_filtration("L_plantarum_HePS", 
                     min_genes = 10)
proximity_filtration("L_johnsonii_ATCC_33200_EPS_A", 
                     min_genes = 7)
proximity_filtration("EPS273", 
                     min_genes = 8)
proximity_filtration("L_johnsonii_ATCC_11506_EPS_B", 
                     min_genes = 10)
proximity_filtration("emulsan", 
                     min_genes = 10)
proximity_filtration("L_johnsonii_ATCC_2767_EPS_C", 
                     min_genes = 10)
proximity_filtration("L_lactis_EPS", 
                     min_genes = 8)
proximity_filtration("cepacian", 
                     min_genes = 10)
