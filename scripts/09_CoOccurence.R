# ============================================================
# Script: 09_CoOccurence.R
# Description: NetCoMi SPIEC-EASI microbiome network analysis
#              for Litopenaeus vannamei production pond microbiome
#              data. Networks constructed per compartment
#              (Intestine, Water, Sediment) at genus level.
# Input: step4_phyloseq_object.rds
# Note: Plot aesthetics have been simplified for public release.
#       Final publication figures used customized formatting
#       not included in this repository.
#       Requires external packages: NetCoMi, SpiecEasi.
#       Install SpiecEasi from GitHub for Matrix >= 1.5.0 compatibility:
#       remotes::install_github("zdk123/SpiecEasi")
# ============================================================

# =========================================================================
# ENVIRONMENT CLEANUP AND INITIALIZATION
# =========================================================================

cleanup_environment <- function() {
  base_packages <- c("stats", "graphics", "grDevices", "utils", "datasets",
                     "methods", "base")
  loaded_packages    <- search()
  packages_to_detach <- loaded_packages[grepl("^package:", loaded_packages)]
  packages_to_detach <- packages_to_detach[
    !grepl(paste(base_packages, collapse = "|"), packages_to_detach)
  ]
  for (pkg in packages_to_detach) {
    tryCatch(
      detach(pkg, character.only = TRUE, unload = TRUE, force = TRUE),
      error = function(e) invisible(NULL)
    )
  }
  gc(verbose = FALSE)
}

cleanup_environment()

# =========================================================================
# PACKAGE INSTALLATION AND LOADING
# =========================================================================

required_packages <- c(
  "phyloseq", "NetCoMi", "SpiecEasi", "circlize", "tidyverse",
  "janitor", "parallel", "stringr", "readr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("phyloseq", "limma")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    library(pkg, character.only = TRUE)
  }
}

suppressWarnings({
  library(phyloseq)
  library(tidytree)
})

# =========================================================================
# SESSION INFO
# =========================================================================

cat("R Version:", R.version.string, "\n")
cat("Key Package Versions:\n")
for (pkg in c("phyloseq", "NetCoMi", "SpiecEasi", "circlize")) {
  cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
}
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

set.seed(54)

# File paths — update for your environment
physeq_path <- "step4_phyloseq_object.rds"
output_dir  <- "output/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Parallel processing
cores <- max(1, parallel::detectCores() - 1)

# Analysis parameters
filter_trial         <- "1"
prevalence_threshold <- 0.2
max_taxa_display     <- 15   # Max taxa in chord diagrams

# SPIEC-EASI parameters
spieceasi_params <- list(
  method        = "mb",
  nlambda       = 30,
  pulsar.params = list(rep.num = 100)
)

# Analysis scope — genus level only
tax_levels   <- c("genus")
sample_types <- c("Intestine", "Water", "Sediment")

# Colors — kept from publication version
type_colors <- c("Intestine" = "#E84646", "Sediment" = "#E69F00", "Water" = "#0072B2")

bright_palette <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
  "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
  "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02",
  "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
  "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
  "#A50026", "#D73027", "#F46D43", "#808080"
)

# Chord diagram parameters
radial_labels      <- FALSE
label_cex          <- 1.2
gap_degree         <- 7
track_height       <- 0.05
chord_transparency <- 0.4
canvas_width       <- 12
canvas_height      <- 12

# =========================================================================
# UTILITY FUNCTIONS
# =========================================================================

clean_taxonomy <- function(physeq_obj) {
  tax_df <- as.data.frame(tax_table(physeq_obj))
  for (col in colnames(tax_df)) {
    if (col %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) {
      tax_df[[col]] <- gsub("[a-z]__",              "", tax_df[[col]])
      tax_df[[col]] <- gsub("_.*",                  "", tax_df[[col]])
      tax_df[[col]] <- gsub("^(Bacteria|Archaea)_", "", tax_df[[col]])
      tax_df[[col]] <- trimws(tax_df[[col]])
      tax_df[[col]] <- ifelse(
        is.na(tax_df[[col]]) | tax_df[[col]] == "", "Unidentified", tax_df[[col]]
      )
    }
  }
  rownames(tax_df) <- rownames(tax_table(physeq_obj))
  tax_table(physeq_obj) <- as.matrix(tax_df)
  return(physeq_obj)
}

safe_sample_data <- function(physeq_obj) {
  tryCatch(
    data.frame(sample_data(physeq_obj)),
    error = function(e) data.frame(row.names = sample_names(physeq_obj))
  )
}

is_technical_name <- function(name) {
  if (is.null(name) || length(name) == 0 || is.na(name) ||
      name == "" || name == "Unidentified") return(TRUE)
  name <- as.character(name)
  technical_patterns <- c(
    "^UBA\\d+", "^ASV_", "^OTU_", "^[A-Z0-9]+-\\d+", "^[A-Z]{2,4}\\d+$",
    "^GCA-", "^SAR\\d+", "^MBA\\d+", "uncultured", "bacterium",
    "Bacteria$", "Archaea$"
  )
  any(sapply(technical_patterns, function(p) grepl(p, name, ignore.case = TRUE)))
}

get_hierarchical_labels <- function(taxa_ids, physeq_obj) {
  if (is.null(physeq_obj) || is.null(tax_table(physeq_obj))) {
    return(setNames(paste0("Taxon_", seq_along(taxa_ids)), taxa_ids))
  }
  
  tax_df     <- as.data.frame(tax_table(physeq_obj))
  taxa_names <- setNames(character(length(taxa_ids)), taxa_ids)
  
  clean_val <- function(x) {
    if (is.na(x)) return(NA)
    x <- gsub("^[a-z]__",              "", x)
    x <- gsub("_.*",                   "", x)
    x <- gsub("^(Bacteria|Archaea)_",  "", x)
    x <- trimws(x)
    if (x == "") return(NA)
    return(x)
  }
  
  for (i in seq_along(taxa_ids)) {
    id <- taxa_ids[i]
    if (!id %in% rownames(tax_df)) {
      taxa_names[i] <- paste0("Taxon_", substr(id, 1, 8))
      next
    }
    row    <- tax_df[id, ]
    genus  <- clean_val(if ("Genus"  %in% names(row)) as.character(row[["Genus"]])  else NA)
    family <- clean_val(if ("Family" %in% names(row)) as.character(row[["Family"]]) else NA)
    order  <- clean_val(if ("Order"  %in% names(row)) as.character(row[["Order"]])  else NA)
    class  <- clean_val(if ("Class"  %in% names(row)) as.character(row[["Class"]])  else NA)
    phylum <- clean_val(if ("Phylum" %in% names(row)) as.character(row[["Phylum"]]) else NA)
    
    assigned <- NA
    for (val in list(genus, family, order, class, phylum)) {
      if (!is.na(val) && !is_technical_name(val)) { assigned <- val; break }
    }
    taxa_names[i] <- if (is.na(assigned)) paste0("Taxon_", substr(id, 1, 8)) else assigned
  }
  
  if (length(unique(taxa_names)) < length(taxa_names)) {
    orig_names <- names(taxa_names)
    taxa_names <- make.unique(taxa_names, sep = "_")
    names(taxa_names) <- orig_names
  }
  return(taxa_names)
}

# =========================================================================
# STEP 1: LOAD AND PREPROCESS DATA
# =========================================================================

physeq <- readRDS(physeq_path)
physeq <- clean_taxonomy(physeq)

if (!is.null(sample_data(physeq))) {
  sample_df <- data.frame(sample_data(physeq))
  if ("type" %in% colnames(sample_df)) {
    sample_df$type <- ifelse(sample_df$type == "Gut", "Intestine",
                             as.character(sample_df$type))
  }
  if (!is.null(filter_trial) && "trial" %in% colnames(sample_df)) {
    sample_df <- sample_df[sample_df$trial == filter_trial, , drop = FALSE]
    physeq    <- prune_samples(rownames(sample_df), physeq)
  }
  sample_data(physeq) <- sample_data(sample_df)
}

physeq_filtered <- prune_samples(sample_sums(physeq) > 0, physeq)
physeq_filtered <- prune_taxa(taxa_sums(physeq_filtered) > 0, physeq_filtered)

# =========================================================================
# STEP 2: NETWORK CONSTRUCTION (SPIEC-EASI, GENUS LEVEL)
# =========================================================================

tax_level_results   <- list()
saved_analysis_keys <- character()

physeq_agg <- tryCatch(
  tax_glom(physeq_filtered, taxrank = "Genus"),
  error = function(e) {
    stop("tax_glom failed at Genus level: ", e$message)
  }
)

sample_results <- list()

for (sample_type in sample_types) {
  sample_df_sub <- safe_sample_data(physeq_agg)
  
  type_samples <- if ("type" %in% colnames(sample_df_sub)) {
    rownames(sample_df_sub)[sample_df_sub$type == sample_type]
  } else {
    rownames(sample_df_sub)
  }
  
  if (length(type_samples) < 10) next
  
  physeq_sub <- prune_samples(type_samples, physeq_agg)
  physeq_sub <- prune_taxa(taxa_sums(physeq_sub) > 0, physeq_sub)
  
  if (ntaxa(physeq_sub) < 10) next
  
  net_result <- tryCatch({
    NetCoMi::netConstruct(
      physeq_sub,
      filtTax     = "numbSamp",
      filtTaxPar  = list(numbSamp = prevalence_threshold),
      measure     = "spieceasi",
      measurePar  = spieceasi_params,
      normMethod  = "none",
      zeroMethod  = "none",
      sparsMethod = "none",
      verbose     = 0,
      cores       = cores,
      seed        = 54
    )
  }, error = function(e) {
    warning("netConstruct failed for ", sample_type, ": ", e$message)
    NULL
  })
  if (is.null(net_result)) next
  
  props_result <- tryCatch({
    NetCoMi::netAnalyze(
      net_result,
      centrLCC        = TRUE,
      avDissIgnoreInf = TRUE,
      sPathNorm       = TRUE,
      clustMethod     = "cluster_fast_greedy",
      hubPar          = "degree",
      hubQuant        = 0.9,
      normDeg         = TRUE,
      normBetw        = TRUE,
      normClose       = TRUE,
      normEigen       = TRUE,
      weightDeg       = FALSE,
      gcmHeat         = FALSE
    )
  }, error = function(e) {
    warning("netAnalyze failed for ", sample_type, ": ", e$message)
    NULL
  })
  if (is.null(props_result)) next
  
  key         <- paste(sample_type, "genus", sep = "_")
  result_file <- file.path(output_dir,
                           paste0("Network_Result_", key, "_SPIEC-EASI.rds"))
  saveRDS(list(network    = net_result,
               properties = props_result,
               physeq     = physeq_sub), result_file)
  
  sample_results[[sample_type]] <- list(
    network    = net_result,
    properties = props_result,
    physeq     = physeq_sub
  )
  saved_analysis_keys <- c(saved_analysis_keys, key)
}

tax_level_results[["genus"]] <- sample_results

# =========================================================================
# STEP 3: CHORD DIAGRAMS
# =========================================================================

chord_results <- list()

for (sample_type in sample_types) {
  key <- paste(sample_type, "genus", sep = "_")
  if (!key %in% saved_analysis_keys) next
  
  network_result <- tax_level_results[["genus"]][[sample_type]]
  if (is.null(network_result)) next
  
  network_obj   <- network_result$network
  physeq_subset <- network_result$physeq
  
  if (is.null(network_obj$edgelist1) || nrow(network_obj$edgelist1) == 0) next
  
  sig_edges <- network_obj$edgelist1 %>%
    dplyr::mutate(weight = if ("adja" %in% colnames(.)) adja else abs(asso)) %>%
    dplyr::filter(weight > 0) %>%
    dplyr::arrange(dplyr::desc(weight))
  
  if (nrow(sig_edges) == 0) next
  
  all_taxa <- unique(c(sig_edges$v1, sig_edges$v2))
  
  if (length(all_taxa) > max_taxa_display) {
    top_taxa  <- sig_edges %>%
      tidyr::pivot_longer(cols = c(v1, v2), values_to = "taxon") %>%
      dplyr::count(taxon, sort = TRUE) %>%
      dplyr::slice_head(n = max_taxa_display) %>%
      dplyr::pull(taxon)
    sig_edges <- dplyr::filter(sig_edges, v1 %in% top_taxa & v2 %in% top_taxa)
    all_taxa  <- unique(c(sig_edges$v1, sig_edges$v2))
  }
  
  taxa_labels <- get_hierarchical_labels(all_taxa, physeq_subset)
  valid_idx   <- !is.na(taxa_labels) & taxa_labels != "" &
    !grepl("^NA$", taxa_labels)
  all_taxa    <- all_taxa[valid_idx]
  taxa_labels <- taxa_labels[valid_idx]
  
  if (length(all_taxa) < 3) next
  
  sig_edges <- dplyr::filter(sig_edges, v1 %in% all_taxa & v2 %in% all_taxa)
  if (nrow(sig_edges) == 0) next
  
  label_lookup <- taxa_labels
  if (length(unique(label_lookup)) < length(label_lookup)) {
    orig         <- names(label_lookup)
    label_lookup <- make.unique(as.character(label_lookup), sep = "_")
    names(label_lookup) <- orig
  }
  
  edge_df <- data.frame(
    source = label_lookup[sig_edges$v1],
    target = label_lookup[sig_edges$v2],
    weight = sig_edges$weight * 100,
    stringsAsFactors = FALSE
  )
  edge_df <- edge_df[!is.na(edge_df$source) & !is.na(edge_df$target), ]
  if (nrow(edge_df) == 0) next
  
  unique_labels <- sort(unique(c(edge_df$source, edge_df$target)))
  tax_matrix    <- matrix(0,
                          nrow = length(unique_labels),
                          ncol = length(unique_labels),
                          dimnames = list(unique_labels, unique_labels))
  for (i in seq_len(nrow(edge_df))) {
    s <- edge_df$source[i]; t <- edge_df$target[i]; w <- edge_df$weight[i]
    if (s %in% rownames(tax_matrix) && t %in% colnames(tax_matrix)) {
      tax_matrix[s, t] <- w
      tax_matrix[t, s] <- w
    }
  }
  
  taxa_colors <- if (length(unique_labels) <= length(bright_palette)) {
    bright_palette[seq_along(unique_labels)]
  } else {
    colorRampPalette(bright_palette)(length(unique_labels))
  }
  names(taxa_colors) <- unique_labels
  
  png_filename <- file.path(output_dir,
                            paste0("Chord_", sample_type, "_Genus_SPIEC-EASI.png"))
  
  tryCatch({
    png(png_filename, width = canvas_width, height = canvas_height,
        units = "in", res = 300)
    circos.clear()
    circos.par(gap.after    = rep(gap_degree, ncol(tax_matrix)),
               track.height = track_height)
    
    chordDiagram(
      tax_matrix,
      directional       = FALSE,
      grid.col          = taxa_colors,
      transparency      = chord_transparency,
      annotationTrack   = c("grid"),
      preAllocateTracks = 1
    )
    
    if (radial_labels) {
      circos.track(track.index = 1, panel.fun = function(x, y) {
        xlim        <- get.cell.meta.data("xlim")
        ylim        <- get.cell.meta.data("ylim")
        sector.name <- get.cell.meta.data("sector.index")
        theta       <- circlize::circlize(mean(xlim), 1.3)[1, 1] %% 360
        text_angle  <- if (theta >= 90 && theta < 270) theta + 180 else theta
        hjust_val   <- if (theta >= 90 && theta < 270) 1 else 0
        circos.text(mean(xlim), ylim[1] + 0.15, sector.name,
                    adj = c(hjust_val, 0.5), cex = label_cex, srt = text_angle)
      }, bg.border = NA)
    } else {
      circos.track(track.index = 1, panel.fun = function(x, y) {
        xlim        <- get.cell.meta.data("xlim")
        ylim        <- get.cell.meta.data("ylim")
        sector.name <- get.cell.meta.data("sector.index")
        circos.text(mean(xlim), ylim[1], sector.name,
                    facing = "bending.outside", niceFacing = TRUE,
                    cex = label_cex)
      }, bg.border = NA)
    }
    
    title(main    = paste0(sample_type, " — Genus associations (SPIEC-EASI)"),
          cex.main = 1.2, font.main = 2)
    dev.off()
    
    chord_results[[key]] <- list(matrix = tax_matrix,
                                 labels = unique_labels,
                                 colors = taxa_colors)
  }, error = function(e) {
    warning("Chord diagram failed for ", key, ": ", e$message)
    if (dev.cur() > 1) dev.off()
    circos.clear()
  })
}

# =========================================================================
# STEP 4: NETWORK PROPERTIES CSV EXPORT
# =========================================================================

summary_stats <- tibble::tibble(
  sample_type  = character(), method       = character(),
  n_taxa       = integer(),   n_edges      = integer(),
  n_hubs       = integer(),   edge_density = numeric(),
  clustering   = numeric(),   modularity   = numeric(),
  top_keystone = character()
)

for (sample_type in sample_types) {
  key            <- paste(sample_type, "genus", sep = "_")
  network_result <- tax_level_results[["genus"]][[sample_type]]
  if (is.null(network_result)) next
  
  net           <- network_result$network
  props         <- network_result$properties
  physeq_subset <- network_result$physeq
  
  n_taxa  <- if (!is.null(net$adjaMat1))  nrow(net$adjaMat1) else 0
  n_edges <- if (!is.null(net$edgelist1)) sum(net$edgelist1$adja > 0) else 0
  n_hubs  <- if (!is.null(props$hubs1))   length(props$hubs1) else 0
  
  top_taxon_id <- if (!is.null(props$centralities) &&
                      nrow(props$centralities) > 0) {
    rownames(props$centralities)[which.max(props$centralities$degree)]
  } else "None"
  
  top_taxon <- if (top_taxon_id != "None") {
    lbl <- get_hierarchical_labels(top_taxon_id, physeq_subset)
    lbl[top_taxon_id]
  } else "None"
  
  gp           <- props$globalProps
  edge_density <- if (!is.null(gp)) {
    idx <- grep("density|Density", names(gp), value = TRUE)
    if (length(idx)) gp[[idx[1]]] else NA_real_
  } else NA_real_
  clustering <- if (!is.null(gp)) {
    idx <- grep("cluster|Cluster", names(gp), value = TRUE)
    if (length(idx)) gp[[idx[1]]] else NA_real_
  } else NA_real_
  modularity <- if (!is.null(gp)) {
    idx <- grep("modular|Modular", names(gp), value = TRUE)
    if (length(idx)) gp[[idx[1]]] else NA_real_
  } else NA_real_
  
  summary_stats <- dplyr::bind_rows(summary_stats, tibble::tibble(
    sample_type  = sample_type,
    method       = "SPIEC-EASI",
    n_taxa       = n_taxa,
    n_edges      = n_edges,
    n_hubs       = n_hubs,
    edge_density = edge_density,
    clustering   = clustering,
    modularity   = modularity,
    top_keystone = top_taxon
  ))
}

readr::write_csv(summary_stats, file.path(output_dir, "Analysis_Summary.csv"))

cat("Analysis complete.\n")
cat("Results saved to:", output_dir, "\n")
cat("Successful networks:", length(saved_analysis_keys), "\n")
cat("Chord diagrams:", length(chord_results), "\n")
