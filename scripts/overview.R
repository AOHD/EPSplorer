
library(tidyr)
library(dplyr) 
library(stringr)
library(ggupset)
library(ggplot2)
library(ggbeeswarm)
library(openxlsx)
`%ni%` <- Negate(`%in%`)

# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if an argument was passed
if (length(args) > 0) {
  # Set the working directory to the argument passed
  setwd(args[1])
}

database_all <- data.frame(
  Target_label = as.character(),
  Query_label = as.character(),
  Psiblast = as.character(),
  operon = as.integer(),
  midas4_tax = as.character()
)


for (file in list.files("./data/output_proximity_filtration/psi_operon_full")) {
 temp <- read.csv2(paste0("./data/output_proximity_filtration/psi_operon_full/",file), sep = "\t") %>% filter(!is.na(Query_label)) %>% select(Target_label, Query_label, Psiblast, operon, midas4_tax) %>% mutate(Psiblast = str_sub(file, end = -5))
 
 database_all <- database_all %>% bind_rows(temp) %>% filter(!is.na(Psiblast))

}

database_all <- database_all %>% mutate(x = paste0(Psiblast, operon))


data_upset_all <- database_all %>%
  mutate(label = substring(Target_label, 1, nchar(Target_label) - 6)) %>%
  mutate(
    Psiblast = str_replace(Psiblast, "alginate", "Alginate"),
    Psiblast = str_replace(Psiblast, "cellulose_All", "Cellulose"),
    Psiblast = str_replace(Psiblast, "HA_Pasteurella", "HA (pmHAS)"),
    Psiblast = str_replace(Psiblast, "HA_streptococcus", "HA (has)"),
    Psiblast = str_replace(Psiblast, "NulO_merged", "NulO"),
    Psiblast = str_replace(Psiblast, "pel_merged", "Pel"),
    Psiblast = str_replace(Psiblast, "pnag_pga", "PNAG (pga)"),
    Psiblast = str_replace(Psiblast, "B_subtilis_EPS", "B. subtilis EPS"),
    Psiblast = str_replace(Psiblast, "pnag_ica", "PNAG (ica)"),
    Psiblast = str_replace(Psiblast, "xanthan", "Xanthan"),
    Psiblast = str_replace(Psiblast, "psl", "Psl"),
    Psiblast = str_replace(Psiblast, "curdlan", "Curdlan"),
    Psiblast = str_replace(Psiblast, "diutan", "Diutan"),
    Psiblast = str_replace(Psiblast, "succinoglycan", "Succinoglycan"),
    Psiblast = str_replace(Psiblast, "gellan2", "Gellan2"),
    Psiblast = str_replace(Psiblast, "burkholderia_eps", "Burkholderia EPS"),
    Psiblast = str_replace(Psiblast, "amylovoran", "Amylovoran"),
    Psiblast = str_replace(Psiblast, "ColA", "Colanic Acid"),
    Psiblast = str_replace(Psiblast, "salecan", "Salecan"),
    Psiblast = str_replace(Psiblast, "stewartan", "Stewartan"),
    Psiblast = str_replace(Psiblast, "vps", "Vibrio EPS"),
    Psiblast = str_replace(Psiblast, "rhizobium_eps", "Rhizobium EPS"),
    Psiblast = str_replace(Psiblast, "gellan1", "Gellan1"),
    Psiblast = str_replace(Psiblast, "acetan", "Acetan"),
    Psiblast = str_replace(Psiblast, "s88", "s88"),
    Psiblast = str_replace(Psiblast, "levan", "Levan"),
    Psiblast = str_replace(Psiblast, "methanolan", "Methanolan"),
    Psiblast = str_replace(Psiblast, "synechan", "Synechan"),
    Psiblast = str_replace(Psiblast, "galactoglucan", "Galactoglucan"),
    Psiblast = str_replace(Psiblast, "succinoglycan", "Succinoglycan"),
    Psiblast = str_replace(Psiblast, "B_fragilis_PS_A", "B. fragilis PS A"),
    Psiblast = str_replace(Psiblast, "B_fragilis_PS_B", "B. fragilis PS B"),
    Psiblast = str_replace(Psiblast, "B_pseudomallei_EPS", "B. pseudomallei PS"),
    Psiblast = str_replace(Psiblast, "cepacian", "Cepacian"),
    Psiblast = str_replace(Psiblast, "E_faecalis_PS", "E. faecalis PS"),
    Psiblast = str_replace(Psiblast, "emulsan", "Emulsan"),
    Psiblast = str_replace(Psiblast, "EPS273", "EPS273"),
    Psiblast = str_replace(Psiblast, "GG", "GG"),
    Psiblast = str_replace(Psiblast, "glucorhamnan", "Glucorhamnan"),
    Psiblast = str_replace(Psiblast, "L_johnsonii_ATCC_33200_EPS_A", "L. johnsonii PS A"),
    Psiblast = str_replace(Psiblast, "L_johnsonii_ATCC_11506_EPS_B", "L. johnsonii PS B"),
    Psiblast = str_replace(Psiblast, "L_johnsonii_ATCC_2767_EPS_C", "L. johnsonii PS C"),
    Psiblast = str_replace(Psiblast, "L_lactis_EPS", "L. lactis PS"),
    Psiblast = str_replace(Psiblast, "L_plantarum_HePS", "L. plantarum PS"),
    Psiblast = str_replace(Psiblast, "phosphonoglycan", "Phosphonoglycan")
    )


##Upset plot

data_upset_all_unique <- data_upset_all %>% group_by(Psiblast) %>% distinct(Target_label) %>% ungroup() %>%
  group_by(Target_label) %>%
  summarise(queries = list(unique(Psiblast))) %>% ungroup()

hit_list <- unique(data_upset_all$Psiblast)


plot_upset <- function(d, sets) {
  ggplot(d, aes(x = queries)) +
    geom_bar() +
    scale_x_upset(n_intersections = 200, sets=sets) +
    ylab("") +
    xlab("") +
    # scale_y_continuous(limits = c(0, 10), breaks = c(1:10)) +
    theme(
      axis.text.y = element_text(size = 12)
    )
} 

upset_all <- plot_upset(data_upset_all_unique, sets = hit_list) + 
  ggtitle("PSI-BLAST hits shared between EPS queries") 


ggsave("./figures/upset_ProkkaNO.pdf", upset_all, limitsize = FALSE, width = 18, height = 12, dpi = 300)


# Which systems are found in which systems? 
systems_in_genomes <- data_upset_all %>% filter(Psiblast != "NulO") %>% group_by(Psiblast) %>% distinct(midas4_tax) %>% ungroup() %>%
  group_by(midas4_tax) %>%
  summarise(systems = list(unique(Psiblast))) 



# Creating a new dataframe with separate columns for each attribute
systems_in_genomes1 <- systems_in_genomes %>%
  unnest(systems) %>%
  mutate(has_attr = TRUE) %>%
  pivot_wider(names_from = systems, values_from = has_attr, values_fill = FALSE)

systems_in_genomes2 <- systems_in_genomes %>% group_by(midas4_tax) %>%
  mutate(systems = paste( unlist(systems), collapse=', '), systems)


work_book <- createWorkbook()



addWorksheet(work_book, sheetName =  "heatmap")
addWorksheet(work_book, sheetName =  "list")

writeData(work_book, "heatmap", systems_in_genomes1)
writeData(work_book, "list", systems_in_genomes2)


saveWorkbook(work_book,
             file= "./figures/systems_overview.xlsx",
             overwrite = TRUE)

