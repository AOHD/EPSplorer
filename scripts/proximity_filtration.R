
##---------------------------------------------------------------
##                         Main function                         
##---------------------------------------------------------------
proximity_filtration <- function(filename_psiblast,
                                 same_database = FALSE, #Used for cellulose because we filter the same psiblast results for 4 (5) different types
                                 perc_id = 20,
                                 max_dist_prok = 8,
                                 max_dist_gene = 2000,
                                 min_genes = 2,
                                 flanking_genes = 2,
                                 magstats = magstats.,
                                 query_metadata = query_metadata.,
                                 gff = gff.,
                                 exclude_gene = "",
                                 essential_genes = NA,
                                 PEgenes_filter = TRUE
                                 ){
  
  if (same_database == FALSE) {
  filename_psiblast_col <- paste(filename_psiblast, collapse = "_")
  } else {
    filename_psiblast_col <- same_database
  }
  
  
  message(glue("Percent identity filtration ({filename_psiblast_col})"))
    # Remove old results
    unlink(glue("./data/output_proximity_filtration/psi_operon_full/{filename_psiblast_col}.tsv"))
    unlink(glue("./data/output_proximity_filtration/psi_proxi_filt/{filename_psiblast_col}.tsv"))
    unlink(glue("./data/output_proximity_filtration/psi_percID_filt/{filename_psiblast_col}.tsv"))
    unlink(glue("./data/output_proximity_filtration/fasta_output/{filename_psiblast_col}/*"))
  # Loading raw psiblast results
  test1 <- filename_psiblast %>%
    paste0("./data/psiblast_results/concatenated/", ., ".faa") %>%
    lapply(
      fread,
      header = FALSE, sep = "\t", fill = TRUE, 
      col.names = c(
        "Query_label", "Target_label", "Percent_identity", "align_length",
        "mismatch", "gap_opens", "start_query", "end_query",
        "start_target", "end_target", "E_value", "Bit score")
    ) %>% bind_rows()
    
    if (ncol(test1) == 0) {
      message(glue("No psiblast results for {filename_psiblast_col}"))
      return()
    } 

    test1 <- test1 %>%
    bind_rows() %>%
    filter(!(Query_label %in% exclude_gene)) %>% 
  # Filtering
    subset(!is.na(E_value)) %>%
    filter(Percent_identity > perc_id) %>%
    group_by(Target_label) %>%
    filter(`Bit score` == max(`Bit score`)) %>%
    filter(!duplicated(Target_label)) %>%
  # Cleaning
    separate(Target_label, c("ID", "ProkkaNO"),
             sep = "_(?!.*_)", remove = FALSE, extra = "merge") %>%
    mutate(
      Target_label = paste(ID, ProkkaNO, sep = "_"),
      ProkkaNO = as.numeric(ProkkaNO)
      ) %>%
  # Merging gff, MAG statistics and query information
    inner_join(gff, by = c("Target_label", "ID", "ProkkaNO") ) %>%
    left_join(magstats, by = "ID", keep = FALSE) %>%
    left_join(query_metadata %>% filter(Psiblast %in% filename_psiblast) %>% select(Genename, Function, Psiblast), 
              by = c("Query_label" = "Genename"), keep=FALSE)  %>%
    arrange(ID, ProkkaNO) %>%
  # Operon Grouping
    group_by(seqname, ID) %>%
    # define gene and prokka distance to prior psiblast hit
    mutate(prio_prok = replace_na(ProkkaNO - shift(ProkkaNO) < max_dist_prok, FALSE),
           prio_gene = replace_na(start - shift(end, 1) < max_dist_gene, FALSE)) %>%
    ungroup() %>% 
    # If a gene don't satisfy distance to prior gene, it's start of operon
    mutate(operon_place = ifelse(!(prio_prok | prio_gene), "start", NA),
           operon = ifelse(operon_place == "start", row_number(), NA)) %>%
    fill(operon, .direction = "down") %>%
    assign(x = "psi_perc_ID_filt", value = ., pos = 1)
  
  # Write percent identity filtered results
  dir.create("./data/output_proximity_filtration/psi_percID_filt", showWarnings = FALSE)
  write.table(psi_perc_ID_filt, 
              file = glue("./data/output_proximity_filtration/psi_percID_filt/{filename_psiblast_col}.tsv"),
              quote = F, sep = "\t", row.names = F)
  
  ##---------------------------------------------------------------
  message(glue("Proximity filtration ({filename_psiblast_col})"))
  # Operon number of genes filtering
  
  if (PEgenes_filter == TRUE) {
    psi_perc_ID_filt <- psi_perc_ID_filt %>%
      group_by(operon) %>%
      filter(any(c("POLTRANS", "TRANS", "Export", "PE", "SY", "Polymerisation") %in% Function))
  }
  
  test2 <- psi_perc_ID_filt %>% group_by(operon) %>%
    filter(length(unique(Query_label)) >= min_genes & (any(essential_genes %in% Query_label) | any(is.na(essential_genes)))) %>%
#    {if(nrow(.) == 0) stop("All results were filtrated away, try less strict filtration.") else .} %>% 
    select(-prio_prok, -prio_gene) %>% 
    ungroup() %>% 
    nest(cols = !c(ID, operon)) %>% 
    group_by(ID) %>% 
    mutate(letter = toupper(letters)[seq_along(unique(operon))],
           ID2 = case_when(
             length(unique(operon)) > 1 ~ paste(ID, letter, sep = "_"),
             TRUE ~ ID)
           ) %>% 
    select(-letter) %>% unnest(cols = c(cols)) %>% 
    assign(x = "psi_proxi_filt", value = ., pos = 1)

  if (nrow(psi_proxi_filt) == 0) {
      message(glue("No {filename_psiblast_col} gene clusters"))
      return()
    } 
  
  # Write proximity filtered results
  dir.create("./data/output_proximity_filtration/psi_proxi_filt", showWarnings = FALSE)
  write.table(psi_proxi_filt, 
              file = glue("./data/output_proximity_filtration/psi_proxi_filt/{filename_psiblast_col}.tsv"), 
              quote = F, sep = "\t", row.names = F)
  
  message(glue("Operon expansion with intermediate and flanking genes ({filename_psiblast_col})"))
  ##---------------------------------------------------------------
  ## Expansion of filtered psiblast hits with surrounding genes
  # Identify surrounding genes
  if (ncol(psi_proxi_filt)*nrow(psi_proxi_filt) != 0 ) {
  test3 <- psi_proxi_filt %>% group_by(operon) %>%
    fill(ID2, .direction = "updown") %>% 
    summarise(
      ID = unique(ID),
      tig = unique(seqname),
      max = max(ProkkaNO),
      min = min(ProkkaNO)
      ) %>%
    full_join(gff, by = "ID") %>%
    group_by(operon) %>%
    filter(ProkkaNO <= (max + flanking_genes) &
           ProkkaNO >= (min - flanking_genes) & 
           seqname  == tig) %>%
    # Merging surrounding and psiblast genes
    full_join(psi_proxi_filt, by = c("operon", "ID", "seqname", "start", "end", "strand", "ProkkaNO", "Target_label")) %>%
    # Prepare for plotting in gggenes - see plot_operon.R
    group_by(operon) %>%
    mutate(
      # gggenes require 1 and -1 as strand direction
      strand = case_when(strand == "+" ~ 1, strand == "-" ~ -1, TRUE ~ 0),
      # Relative distances (start always minimum value of gene position)
      end = end - min(start),
      start = start - min(start),
      # converting AA numbers to nucleotide numbers
      end_target_plot =   ifelse(is.na(end_target),   end,   end_target   * 3 + start),
      start_target_plot = ifelse(is.na(start_target), start, start_target * 3 + start),
      # Assigning artificial percent identity to domains
      Percent_identity = ifelse(is.na(Percent_identity), 40, Percent_identity),
      # Order operons from highest bit score to lowest
      operon = factor(operon, levels = unique(arrange(., -`Bit score`)$operon))
      ) %>%
    fill(midas4_tax, .direction = "updown") %>%
    fill(ID2, .direction = "updown") %>% 
    # ->ASSIGNING "psi_operon_full"<- 
    assign(x = "psi_operon_full", value = ., pos = 1)
  # Write operons with surrounding and intermediate genes
  dir.create("./data/output_proximity_filtration/psi_operon_full", showWarnings = FALSE)
  write.table(psi_operon_full, file = glue("./data/output_proximity_filtration/psi_operon_full/{filename_psiblast_col}.tsv"), 
              quote = F, sep = "\t", row.names = F)
  
  ##------------------------------------------------------------------
  ##  Export fasta files of identified genes, with surrounding genes  
  ##------------------------------------------------------------------
  message(glue("Saving FASTA files of genes in expanded operons ({filename_psiblast_col})"))
  for (f in unique(psi_operon_full$ID2)){
    dir.create(glue("./data/output_proximity_filtration/fasta_output/{filename_psiblast_col}"), showWarnings = FALSE, recursive = TRUE)
    g <- unique(filter(psi_operon_full, ID2 == f)$ID)
    read.fasta(file = glue("./data/prokka/{g}.faa"), 
                            seqtype="AA", 
                            as.string=TRUE, 
                            set.attributes=FALSE) %>% 
      `[`(names(.) %in% psi_operon_full$Target_label) %>% 
      write.fasta(names=names(.), file.out=glue("./data/output_proximity_filtration/fasta_output/{filename_psiblast_col}/{f}.faa"))
  }
}

                                 }
