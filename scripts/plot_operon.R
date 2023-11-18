plot_operon <-  function(filename_psiblast,
                         same_database = FALSE, #Used for cellulose because we filter the same psiblast results for 4 (5) different types
                         name_addon = "",
                         query_title = filename_psiblast,
                         article_plot_domain = TRUE,
                         mags = "all",
                         width = 2,
                         textsize = 2,
                         expression = FALSE, 
                         limits = c(0.1, NA),
                         exclude_MAGs = "",
                         Processing_tank = "Anoxic stage"){
  ##---------------------------------------------------------------
  ##      Parameter definition
  ##---------------------------------------------------------------
  # The name of the .txt files
  if (same_database == FALSE) {
    filename_psiblast_col <- paste(filename_psiblast, collapse = "_")
  } else {
    filename_psiblast_col <- same_database
  }
  # Loading results from proximity filtration
  genes <- fread(glue("./data/output_proximity_filtration/psi_operon_full/{filename_psiblast_col}.tsv"))
  if(any("all" != mags)) {
    genes <- filter(genes, ID %in%  mags)
  }
  
  genes <- genes %>% filter(ID2 %ni% exclude_MAGs)
  
  # File with metadata on the query genes, e.g. function
  query_metadata. <- excel_sheets("./Query_figur.xlsx") %>%
    sapply(function(X) read_xlsx("./Query_figur.xlsx", sheet = X, skip = 1), USE.NAMES = T) %>% 
    lapply(as.data.frame) %>% 
    `[`(!(names(.) %in% c("Abbreveations", "HA_S_pyogenes"))) %>%
    rbindlist(fill = TRUE)  %>% 
    mutate(
      Function = str_replace(Function, "MOD", "Modification"),
      Function = str_replace(Function, "PE", "Polymerization & Export"),
      Function = str_replace(Function, "GT", "Glycosyl Transferase"),
      Function = str_replace(Function, "ABC", "ABC transporter"),
      Function = str_replace(Function, "SY", "Synthase"),
      Function = str_replace(Function, "PS", "Precursor Synthesis"),
      Function = str_replace(Function, "REG", "Regulation"),
      Function = str_replace(Function, "DEG", "Degradation"),
      Function = str_replace(Function, "TRANS", "Transport"),
      Function = str_replace(Function, "POL", "Polymerisation"),
      Function = str_replace(Function, "POLTRANS", "Polymerisation and Transport"),
      Function = str_replace(Function, "PRE", "Precursor synthesis"),
      Function = str_replace(Function, "NA", "Unknown function"),
    ) %>% mutate(Polysaccheride = ifelse(Polysaccheride == "S88","s88",Polysaccheride),
                 Psiblast = ifelse(Psiblast == "S88","s88",Psiblast))
  
  # colors in the "Function" legend 
  function_colors <- c("#c4a78b", "#944c2e", "#fbbb81",
                       "#305b46", "#fbc669", "#db6447",
                       "#e5e4b8", "#b3b887", "#ffe892", 
                       "#e4fbe2", "#c3dbce", "#89b5b0", 
                       "#425985", "#403675", "#f18009",
                       "#18a2bf", "#18a2bf", "#18a2bf")
  # print query_metadata.$Function %>% unique() to see order
  names(function_colors) <- query_metadata.$Function %>% unique()
  
  
  ##---------------------------------------------------------------
  ##  Loading interproscan data and merging with psiblast operons  
  ##---------------------------------------------------------------
  if (article_plot_domain == TRUE) {
  clean_gff3 <- function(df){
    # Extraction of relevant annotation information
    df %>% filter(V3 != "polypeptide") %>%
      mutate(
        domain = str_extract(V9, "signature_desc=[^;]*;"),
        domain = str_sub(domain, 1, -2),
        domain = gsub("signature_desc=", "", x = domain)) %>%
      subset(select = -c(V2, V3, V7, V8, V9)) %>%
      setNames(c("Target_label", "start1", "end1", "e_value", "Domain")) %>%
      filter(!is.na(Domain)) %>%
      distinct() %>%
      # Formating of domain names (a bit of confusing regular expressions)
      mutate(
        Domain = str_replace(Domain, "[gG]lycosyl.*transferase.*[fF]amily ", "GT family "),
        Domain = str_replace(Domain, "[gG]lycosyl.*transferase.*[gG]roup ", "GT group "),
        Domain = str_replace(Domain, "[gG]lycosyl.*transferase.[lL]ike.[fF]amily", "GT like family "),
        Domain = str_replace(Domain, "[gG]lycosyl.*transferase", "GT "),
        Domain = str_replace(Domain, "[gG]lycosyl [hH]ydrolase.*[fF]amily " , "GH"),
        Domain = str_replace(Domain, "[gG]lycosyl [hH]ydrolase" , "GH"),
        Domain = str_replace(Domain, ".*[cC]ellulose.*synth.*protein[^' ']*" , "CS "),
        Domain = str_replace(Domain, ".*[cC]ellulose.*synth[^' ']*", "CS "),
        Domain = str_replace(Domain, " N.terminal domain" , " N terminal"),
        Domain = str_replace(Domain, " C.terminal domain" , " C terminal"),
        Domain = str_replace(Domain, ".*BCSC_C.*", "BcsC"),
        Domain = str_replace(Domain, ".*GIL.*", "BcsE"),
        Domain = str_replace(Domain, ".*subunit D.*", "BcsD"),
        Domain = str_replace(Domain, ".*complementing protein A.*", "ccpA"),
        Domain = str_replace(Domain, "[iI]nitation [fF]actor" , "IF"),
        Domain = str_replace(Domain, "[eE]longation [fF]actor" , "EF"),
        Domain = str_replace(Domain, ".*Tetratrico.*|.*TPR.*", "Tetratrico"),
        Domain = str_replace(Domain, ".*[dD]omain of unknown function.*", "NA"),
        Domain = str_replace(Domain, "Phosphoglucomutase/phosphomannomutase, alpha/beta/alpha domain II", "PGM_PMM_II"),
        Domain = str_replace(Domain, "Phosphoglucomutase/phosphomannomutase, alpha/beta/alpha domain I", "PGM_PMM_I"),
        Domain = str_replace(Domain, "Phosphoglucomutase/phosphomannomutase, alpha/beta/alpha domain III", "PGM_PMM_III"),
        Domain = str_replace(Domain, "Phosphoglucomutase/phosphomannomutase,_C", "PGM_PMM_IV"),
        Domain = str_replace(Domain, "RTX calcium-binding nonapeptide repeat.*", "Hemolysn_Ca-bd"),
        Domain = str_replace(Domain, "UDP-glucose/GDP-mannose dehydrogenase family, UDP binding domain", "UDPG_MGDP_dh_C"),
        Domain = str_replace(Domain, "UDP-glucose/GDP-mannose dehydrogenase family, central domain", "UDP-Glc/GDP-Man_DH_dimer"),
        Domain = str_replace(Domain, "UDP-glucose/GDP-mannose dehydrogenase family, NAD binding domain", "UDPG_MGDP_dh_N"),
        Domain = str_replace(Domain, "SGNH hydrolase-like domain, acetyltransferase AlgX", "ALGX/ALGJ_SGNH-like"),
        Domain = str_replace(Domain, "MBOAT, membrane-bound O-acyltransferase family", "MBOAT_fam"),
        Domain = str_replace(Domain, "Periplasmic copper-binding protein.*", "NosD_dom"),
        Domain = str_replace(Domain, "UDP-N-acetylglucosamine 2-epimerase", "UDP-GlcNAc_Epase")
        
      ) 
  }
  ## Loading
  if (article_plot_domain == TRUE) {
  domains <- filename_psiblast_col %>% 
    lapply(function(query) {
      # Use this function on each polysaccharide name (e.g. "cellulose1)
      list.files(paste0("./data/interproscan_results/", query, "_fasta_removed/")) %>%
        lapply(function(mag) {
          # Use this function on each mag ID (e.g. "Ega_18-Q3-R5-49_MAXAC.199")
          mag_path = paste0("./data/interproscan_results/", query, "_fasta_removed/", mag)
          if (file.info(mag_path)$size > 0) {
            read.table(
              mag_path,
              sep = "\t",
              fill = TRUE,
              comment.char = "#",
              quote = "\""
            )
          }
        }) %>%
        bind_rows()
    }) %>%
    bind_rows() %>%
    clean_gff3() %>% 
    ## Combining information from genes to domains
    full_join(genes) %>%
    mutate(start2 = start + as.numeric(start1) * 3,
           end2 = start + as.numeric(end1) * 3) %>%
    subset(select = c(
      "start", "end", "start2", "end2", "Domain", "operon", "ID",
      "ID2", "Function", "strand", "midas4_tax")) %>%
    mutate(Percent_identity = 50)
  }
  }
  ##----------------------------------------------------------------
  ##                    Reversing strand direction                    
  ##----------------------------------------------------------------
  operon_reversed <- genes %>%
    group_by(operon) %>% 
    summarize(direction = sum(strand)) %>% 
    filter(direction < 0) %>% 
    pull(operon)
  
  
  genes <- genes %>% 
    group_by(operon) %>% 
    mutate(
      reverse = operon %in% operon_reversed,
      operon_middle = (min(c(start, end)) + max(c(start, end)))/2,
      start = ifelse(reverse, 2*operon_middle - start, start),
      end = ifelse(reverse, 2*operon_middle - end, end),
      min_operon = min(c(start, end)),
      start_target_plot = ifelse(reverse, 2*operon_middle - start_target_plot, start_target_plot) - min_operon,
      end_target_plot = ifelse(reverse, 2*operon_middle - end_target_plot, end_target_plot) - min_operon,
      end = end - min_operon,
      start = start - min_operon
    ) %>% 
    select(-reverse, -operon_middle)
  
  if (article_plot_domain == TRUE) {
  domains <- domains %>% 
    group_by(operon) %>% 
    mutate(
      reverse = operon %in% operon_reversed,
      # Reverse direction of matched gene, same as above
      operon_middle = (min(c(start, end)) + max(c(start, end)))/2,
      start = ifelse(reverse, 2*operon_middle - start, start),
      end = ifelse(reverse, 2*operon_middle - end, end),
      # Defining minimum in operon to ensure all start at 0
      min_operon = min(c(start, end)),
      gene_middle = (start + end)/2,
      # Reverse direction in operon
      start2 = ifelse(reverse, 2*operon_middle - start2, start2) - min_operon,
      end2 = ifelse(reverse, 2*operon_middle - end2, end2) - min_operon,
      # Reverse direction in gene
      start2 = ifelse(reverse, 2*gene_middle - start2, start2) - min_operon,
      end2 = ifelse(reverse, 2*gene_middle - end2, end2) - min_operon
      
      #strand = ifelse(reverse, strand, strand)
    ) %>% 
    select(-reverse, -operon_middle)
  }
  
  
  ##----------------------------------------------------------------
  ##                    Adding query information                    
  ##----------------------------------------------------------------
  genes <- query_metadata. %>% 
    mutate(
      start = Start,
      start_target_plot = Start,
      end = End,
      end_target_plot = End,
      title = glue("**{query_title}**"),
      ID2 = "Query",
      mi_phylum = "a",
      operon = 0,
      strand = ifelse(Strand == "+", 1, -1),
      Query_label = Genename,
    ) %>% 
    filter(Psiblast %in% filename_psiblast) %>% 
    full_join(genes) %>% 
    mutate(
      Function = str_replace(Function, "MOD", "Modification"),
      Function = str_replace(Function, "PE", "Polymerization & Export"),
      Function = str_replace(Function, "GT", "Glycosyl Transferase"),
      Function = str_replace(Function, "ABC", "ABC transporter"),
      Function = str_replace(Function, "SY", "Synthase"),
      Function = str_replace(Function, "PS", "Precursor Synthesis"),
      Function = str_replace(Function, "REG", "Regulation"),
      Function = str_replace(Function, "DEG", "Degradation"),
      Function = str_replace(Function, "TRANS", "Transport"),
      Function = str_replace(Function, "POL", "Polymerisation"),
      Function = str_replace(Function, "POLTRANS", "Polymerisation and Transport"),
      Function = str_replace(Function, "PRE", "Precursor synthesis"),
      Function = str_replace(Function, "UNKNOWN", "Unknown function"),
    )
  
  if (article_plot_domain == TRUE) {
  # Extraction of relevant annotation information
  domains <- read.table( 
    glue("./queries_ips/{filename_psiblast}.gff3"),
    sep = "\t", fill = TRUE, comment.char = "#", quote = "\""
  ) %>%
    clean_gff3() %>% 
    mutate(Genename = Target_label) %>% 
    # Adding gene information about query genes
    left_join(query_metadata.) %>% 
    mutate(
      start2 = Start + as.numeric(start1) * 3,
      end2 = Start + as.numeric(end1) * 3,
      #ID = "Query",
      ID2 = "Query",
      title = glue("**{query_title}**"),
      mi_phylum = "a",
      operon = 0,
      strand = ifelse(Strand == "+", 1, -1),
      start = Start,
      end = End
    ) %>% 
    select(c("start2", "end2", "Domain", "title",
             "ID2", "Function", "strand", "operon", "start", "end", "mi_phylum")) %>%
    mutate(Percent_identity = 50) %>%
    filter(!is.na(Domain)) %>% 
    full_join(domains) %>% 
    mutate(
      Function = str_replace(Function, "MOD", "Modification"),
      Function = str_replace(Function, "PE", "Polymerization & Export"),
      Function = str_replace(Function, "GT", "Glycosyl Transferase"),
      Function = str_replace(Function, "ABC", "ABC transporter"),
      Function = str_replace(Function, "SY", "Synthase"),
      Function = str_replace(Function, "PS", "Precursor Synthesis"),
      Function = str_replace(Function, "REG", "Regulation")
    )

  
  
  ##---------------------------------------------------------------
  ##       Make Domain colours
  ##---------------------------------------------------------------
  
  # colors in the "Domains" legend 
  # 30 distinct colors generated by https://mokole.com/palette.html
  domain_colors <- c("#696969", "#556b2f", "#228b22", "#7f0000", "#483d8b",
                     "#008b8b", "#4682b4", "#000080", "#d2691e", "#daa520",
                     "#8fbc8f", "#800080", "#b03060", "#ff4500", "#0000cd",
                     "#00ff00", "#00ff7f", "#dc143c", "#ffdead", "#00ffff",
                     "#f08080", "#adff2f", "#d8bfd8", "#ff00ff", "#1e90ff",
                     "#ffff54", "#90ee90", "#ff1493", "#7b68ee", "#ee82ee",
                     "#696969", "#556b2f", "#228b22", "#7f0000", "#483d8b",
                     "#008b8b", "#4682b4", "#000080", "#d2691e", "#daa520",
                     "#8fbc8f", "#800080", "#b03060", "#ff4500", "#0000cd",
                     "#00ff00", "#00ff7f", "#dc143c", "#ffdead", "#00ffff",
                     "#f08080", "#adff2f", "#d8bfd8", "#ff00ff", "#1e90ff",
                     "#ffff54", "#90ee90", "#ff1493", "#7b68ee", "#ee82ee")
  
  ##---------------------------------------------------------------
  ##       Filtering relevant domains
  ##---------------------------------------------------------------
  if (article_plot_domain == TRUE) {
  domains_filtered <- domains %>% filter(
    grepl("GT|[sS]accharide|[sS]ugar|[cC]arbohydrate|[eE][pP][sS]|ABC|[eE]pimerase|[wW]z|[Aa]cetyltransferase|[bB]eta[-| ]barrel|GH[0-9]|
          [sS]ps|[Aa]tr|[rR]hs|[pP]ga|[Ii]ca[A|B|C|D]|[gG]el|[rR]ml|[dD]ps|[hH]as[A|B|C]|pmHAS|[pP]sl|[sS]le|[bB]ep|
          [eE]xo|[gG]um|[bB]cs|[nN]eu|[cC]rd|[pP]el[A|B|C|D|E|F|G]|[aA]lg|[cC]ps|[aA]ms|[wW]ca|[pP]ss|
          [aA]ce[A|B|C|D|E|G|H|I|R|Q|P|M|F]|[tT]etratrico|[mM]annos|[rR]hamnos",
          Domain))
  
  domains_filtered <- domains_filtered %>%
    mutate(Domain = replace(Domain, str_detect(Domain, ".*GT.*4.*"), "GT family 4"),
           Domain = replace(Domain, str_detect(Domain, ".*GT.*2.*"), "GT family 2"),
           Domain = replace(Domain, str_detect(Domain, ".*GT.*1.*"), "GT family 1"),
           Domain = replace(Domain, str_detect(Domain, ".*GT.*9.*"), "GT family 9"),
           Domain = replace(Domain, str_detect(Domain, ".*ABC.*[tT]ransporter.*"), "ABC transporter"),
           Domain = replace(Domain, str_detect(Domain, ".*[eE]xosortase*"), "Exosortase"),
           Domain = replace(Domain, str_detect(Domain, ".*[eE]pimerase.*"), "Epimerase"),
           Domain = replace(Domain, str_detect(Domain, ".*[iI]somerase.*"), "Isomerase"),
           Domain = replace(Domain, str_detect(Domain, ".*[dD]ehydratase.*"), "Dehydratase"),
           Domain = replace(Domain, str_detect(Domain, ".*[aA]cetyltransferase.*"), "Acetyltransferase"),
           Domain = replace(Domain, str_detect(Domain, ".*[mM]annosyltransferase.*"), "Mannosyltransferase"),
           Domain = replace(Domain, str_detect(Domain, ".*[wW]za.*"), "wza"),
           Domain = replace(Domain, str_detect(Domain, ".*[bB]eta[-| ]barrel.*"), "Beta-barrel"),
           Domain = replace(Domain, str_detect(Domain, "^((?![0-9]).)*GT((?![0-9]).)*$"), "GT"),
           Domain = replace(Domain, str_detect(Domain, "[^((?![0-9]).)*pP]olysaccharide((?![0-9]).)*$"), "PS"),
           Domain = replace(Domain, str_detect(Domain, ".*Carbohydrate-binding.*module.*48.*(Isoamylase.*N.*terminal.*)"), "NA"),
           Domain = replace(Domain, str_detect(Domain, ".*CS.*BcsQ.*"), "BcsQ"),
           Domain = replace(Domain, str_detect(Domain, ".*CS.*BcsN,*"), "BcsN"),
           Domain = replace(Domain, str_detect(Domain, ".*CS.*BcsG.*"), "BcsG"))
  
  names(domain_colors) <- domains_filtered$Domain %>% unique()
  }
  
  ##---------------------------------------------------------------
  ##            Plotting operons with domain annotation            
  ##---------------------------------------------------------------
  assign("genes", genes, pos = 1)
  if (article_plot_domain == TRUE) {
  assign("domains", domains_filtered, pos = 1)
  }
  }

  ##---------------------------------------------------------------
  ##            Operon plotting            
  ##---------------------------------------------------------------
  
  genes <- genes %>% mutate(
    Psiblast = str_replace(Psiblast, "alginate", "Alginate"),
    Psiblast = ifelse(Psiblast == "cellulose" & same_database == "celluloseNA", "Unclassified cellulose", Psiblast),
    Psiblast = ifelse(Psiblast == "cellulose" & same_database == "celluloseII", "Cellulose II", Psiblast),
    Psiblast = str_replace(Psiblast, "celluloseII", "Cellulose II"),
    Psiblast = str_replace(Psiblast, "celluloseIII", "Cellulose III"),
    Psiblast = str_replace(Psiblast, "cellulose_Ac", "Acetylated cellulose"),
    Psiblast = str_replace(Psiblast, "cellulose_NA", "Unclassified cellulose"),
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
  
  Psiblast <- genes$Psiblast
  Psiblast <- unique(Psiblast[!is.na(Psiblast)])
  
  gene_height <- 5
  n_operon <- length(unique(genes$ID2))
  
  if (article_plot_domain == TRUE) {
  domains <- domains %>% filter(ID2 %ni% exclude_MAGs)
  domains_filtered <- domains_filtered %>% filter(ID2 %ni% exclude_MAGs) %>% mutate(mi_phylum = ifelse(ID2 == "Query", "a", mi_phylum))
  }
  genes <- genes %>% mutate(mi_phylum = ifelse(ID2 == "Query", "a", mi_phylum))
  
  operon_plot <- ggplot(
    genes, aes(xmin = start, xmax = end, y = paste(mi_phylum), forward = strand)
  ) +
    # Empty gene arrows
    geom_gene_arrow(
      arrowhead_height = unit(gene_height, "mm"),
      arrow_body_height = unit(gene_height, "mm"),
      arrowhead_width = unit(5, "mm")
    ) +
    # Colored gene arrows (match in psiblast)
    geom_subgene_arrow(
      data = genes,
      mapping = aes(xmin = start, xmax = end, y = paste(mi_phylum),
                    xsubmin = start_target_plot,
                    xsubmax = end_target_plot,
                    fill = Function,
                    forward = strand
      ),
      arrowhead_height = unit(gene_height, "mm"),
      arrow_body_height = unit(gene_height, "mm"),
      arrowhead_width = unit(5, "mm"),
      # position = position_nudge(y = 0.3)
    ) +
    facet_wrap(
      ~ mi_phylum + ID2, 
      scales = "free_y",
      ncol = 1
    ) +
    theme_genes() +
    theme(
      #legend.position = "right",
      legend.direction = "vertical",
      legend.spacing.x = unit(6, "mm"),
      legend.text = element_text(margin = margin(b  = -15)),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")
    ) + scale_fill_manual(
      values = function_colors[names(function_colors) %in% unique(genes$Function)], 
      na.value = "transparent") +
    guides(
      fill = guide_legend(
        title = "Function of Matched Query Gene", 
        title.position = "top", 
        title.hjust = 0.5, 
        label.vjust = 1,
        label.theme = element_text(size = 6),
        keyheight = 0.9,
        label.position = "bottom")
    ) +
    scale_x_continuous(expand = rep(max(genes$end, genes$start) * 0.0000003, 2))
  if (article_plot_domain == TRUE) {
    ##---------------------------------------------------------------
    ##            Operon plotting with domain            
    ##---------------------------------------------------------------
    operon_plot2 <- operon_plot + 
      # Domains boxes
      geom_text(
        data = genes %>% mutate(start = (start + end)/2),
        mapping = aes(x = start, label = Query_label),
        size = textsize,
        nudge_y = -0.25, 
        angle = -45,
        hjust = 0
      ) +
      geom_richtext(
        data = genes %>% group_by(operon) %>% slice(1),
        mapping = aes(x = 0, label = midas4_tax),
        size = 4,
        nudge_y = 0.45,
        hjust = 0,
        fill = NA, label.color = NA
      ) +
      new_scale_fill() +
      geom_gene_arrow(
        data = domains_filtered, #%>% filter(ID2 != "Query"),
        mapping = aes(
          fill = Domain,
          xmin = start2,
          xmax = end2,
          forward = strand
        ),
        arrowhead_height = unit(gene_height - 2.5, "mm"),
        arrow_body_height = unit(gene_height - 2.5, "mm"),
        arrowhead_width = unit(0, "mm"),
        position = position_nudge(y = +0.22)
      ) +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.spacing.x = unit(6, "mm"),
        legend.text = element_text(margin = margin(b  = -15)),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), plot.margin = unit(c(0,0,0,0), "cm"),
        legend.background = element_rect(linetype = 1, linewidth = 0.5, colour = 1)
      ) +
      scale_fill_manual(
        values = domain_colors[names(domain_colors) %in% unique(domains_filtered$Domain)], 
        na.value = "transparent") +
      guides(
        fill = guide_legend(
          title = "IPS Domain Annotations",
          title.position = "top", 
          title.hjust = 0.5, 
          label.vjust = 0.75,
          label.theme = element_text(size = 6),
          keyheight = 0.9,
          label.position = "bottom"),
      ) +
      ggtitle(Psiblast)
    
    
    ggsave(
      plot = operon_plot2,
      glue("./figures/operons/operon_", paste(filename_psiblast_col, collapse = "_"), "{name_addon}.pdf"), 
      width = unit(13, "mm"),
      height = unit(1.4 * n_operon + width, "mm"),
      limitsize = FALSE)
    
    
  } else{
    operon_plot <- operon_plot + 
      geom_text(
        data = genes %>% mutate(start = (start + end)/2),
        mapping = aes(x = start, label = Query_label),
        size = textsize,
       nudge_y = -0.25, 
        angle = -45,
        hjust = 0
      ) +
      geom_richtext(
        data = genes %>% group_by(operon) %>% slice(1),
        mapping = aes(x = 0, label = midas4_tax),
        size = 4,
        nudge_y = 0.45,
        hjust = 0,
        fill = NA, label.color = NA
      ) +
      ggtitle(Psiblast)

    ggsave(
      plot = operon_plot,
      glue("./figures/operons/operon_", paste(filename_psiblast_col, collapse = "_"), "{name_addon}.pdf"), 
      width = unit(13, "mm"),
      height = unit(1.4 * n_operon + 1, "mm"),
      limitsize = FALSE)
  }
}
