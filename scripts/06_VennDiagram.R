# ============================================================
# Script: 06_VennDiagram.R
# Description: Venn diagram analysis of shared and unique
#              genera across Intestine, Sediment, and Water
#              compartments in Litopenaeus vannamei production
#              ponds. Analysis uses iterative
#              rarefaction with consensus threshold.
# Input: phyloseq object (step4_phyloseq_object.rds)
# Note: Plot aesthetics have been simplified for public release.
#       Final publication figures used customized formatting
#       not included in this repository.
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
  "phyloseq", "vegan", "ggplot2", "patchwork", "ggtext",
  "ggVennDiagram", "janitor", "future.apply", "parallel", "future", "tidyverse"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("phyloseq", quietly = TRUE)) {
  BiocManager::install("phyloseq", update = FALSE, ask = FALSE)
  library(phyloseq)
}

# =========================================================================
# SESSION INFO
# =========================================================================

cat("R Version:", R.version.string, "\n")
cat("Key Package Versions:\n")
key_packages <- c("phyloseq", "vegan", "ggplot2", "ggVennDiagram", "dplyr")
for (pkg in key_packages) cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

# File paths
phyloseq_path <- "step4_phyloseq_object.rds"
output_dir    <- "output/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
n_iterations        <- 100
rarefaction_depth   <- 57841
consensus_threshold <- 0.2   # Taxon retained if present in >= this fraction of iterations
set.seed(54)

# Sample filtering
filter_ponds      <- c("C3", "C6", "D1", "D4")
filter_timepoints <- c("Week 0", "Week 2", "Week 4", "Week 6", "Week 8")
filter_shrimp_max <- 5

# Taxonomic level
tax_level <- "Genus"

# Placeholder color gradients — publication colors not included in public release
color_gradients <- list(
  "Intestine" = c("#FDEEEE", "#F8BEBE", "#F28E8E", "#EC5E5E", "#E84646"),
  "Sediment"  = c("#FDF6E6", "#FAE5B3", "#F4CC66", "#F0BB33", "#E69F00"),
  "Water"     = c("#E6F2FF", "#B3D9FF", "#66B2FF", "#3399FF", "#0072B2")
)

# =========================================================================
# UTILITY FUNCTIONS
# =========================================================================

clean_taxonomy <- function(tax_names) {
  tax_names <- gsub("^[a-z]__", "", tax_names)
  tax_names <- gsub("_.*", "", tax_names)
  tax_names[is.na(tax_names) | tax_names == ""] <- "Unidentified"
  return(tax_names)
}

process_sample_metadata <- function(physeq) {
  sample_metadata <- as_tibble(sample_data(physeq), rownames = "sample_id") %>%
    janitor::clean_names() %>%
    dplyr::filter(trial == "1")
  
  sample_metadata %>%
    dplyr::mutate(
      group_type = dplyr::case_when(
        grepl("gut|intestin", type, ignore.case = TRUE) ~ "Intestine",
        grepl("water",        type, ignore.case = TRUE) ~ "Water",
        grepl("soil|sediment",type, ignore.case = TRUE) ~ "Sediment",
        TRUE ~ as.character(type)
      ),
      group_time = dplyr::case_when(
        grepl("^T0$", timepoint) ~ "Week 0",
        grepl("^T1$", timepoint) ~ "Week 2",
        grepl("^T2$", timepoint) ~ "Week 4",
        grepl("^T3$", timepoint) ~ "Week 6",
        grepl("^T4$", timepoint) ~ "Week 8",
        TRUE ~ as.character(timepoint)
      )
    ) %>%
    dplyr::filter(
      pond %in% filter_ponds,
      group_time %in% filter_timepoints,
      is.na(shrimp) | shrimp <= filter_shrimp_max
    )
}

get_taxa_for_samples <- function(sample_ids, rarefied, tax_data) {
  if (length(sample_ids) == 0) return(character(0))
  
  if (length(sample_ids) == 1) {
    present_taxa <- rarefied[sample_ids, ] > 0
  } else {
    present_taxa <- colSums(rarefied[sample_ids, , drop = FALSE] > 0) > 0
  }
  
  present_taxa_ids  <- names(present_taxa)[present_taxa]
  matching_taxa     <- tax_data$taxa_id %in% present_taxa_ids
  tax_names         <- tax_data[[tax_level]][matching_taxa]
  tax_names         <- tax_names[!is.na(tax_names)]
  return(unique(tax_names))
}

# =========================================================================
# RAREFACTION FUNCTION
# =========================================================================

rarefaction_for_venn <- function(iteration, otu_matrix, depth, physeq_obj,
                                 sample_types, time_points) {
  set.seed(iteration * 100)
  rarefied <- vegan::rrarefy(otu_matrix, depth)
  
  tax_table_data         <- as.data.frame(tax_table(physeq_obj))
  tax_table_data$taxa_id <- rownames(tax_table_data)
  for (col in setdiff(names(tax_table_data), "taxa_id")) {
    tax_table_data[[col]] <- clean_taxonomy(tax_table_data[[col]])
  }
  
  sample_data_df            <- as.data.frame(sample_data(physeq_obj))
  sample_data_df$sample_id  <- rownames(sample_data_df)
  
  # Taxa by type at each timepoint
  results_by_timepoint <- list()
  for (tp in time_points) {
    taxa_by_type <- list()
    for (st in sample_types) {
      type_samples <- sample_data_df$sample_id[
        sample_data_df$group_time == tp & sample_data_df$group_type == st
      ]
      type_samples <- intersect(type_samples, rownames(rarefied))
      if (length(type_samples) > 0) {
        taxa_by_type[[st]] <- get_taxa_for_samples(type_samples, rarefied, tax_table_data)
      }
    }
    taxa_by_type <- taxa_by_type[lengths(taxa_by_type) > 0]
    if (length(taxa_by_type) >= 2) results_by_timepoint[[tp]] <- taxa_by_type
  }
  
  # Taxa by timepoint within each type
  results_by_sampletype <- list()
  for (st in sample_types) {
    taxa_by_time <- list()
    for (tp in time_points) {
      time_samples <- sample_data_df$sample_id[
        sample_data_df$group_type == st & sample_data_df$group_time == tp
      ]
      time_samples <- intersect(time_samples, rownames(rarefied))
      if (length(time_samples) > 0) {
        taxa_by_time[[tp]] <- get_taxa_for_samples(time_samples, rarefied, tax_table_data)
      }
    }
    taxa_by_time <- taxa_by_time[lengths(taxa_by_time) > 0]
    if (length(taxa_by_time) >= 2) results_by_sampletype[[st]] <- taxa_by_time
  }
  
  list(by_timepoint = results_by_timepoint, by_sampletype = results_by_sampletype)
}

# =========================================================================
# COMBINE RAREFACTION RESULTS
# =========================================================================

combine_rarefaction_results <- function(results_list, type = "by_timepoint") {
  all_results <- lapply(results_list, function(x) x[[type]])
  all_keys    <- unique(unlist(lapply(all_results, names)))
  
  combined_results <- list()
  for (key in all_keys) {
    key_results <- Filter(Negate(is.null), lapply(all_results, function(x) x[[key]]))
    if (length(key_results) == 0) next
    
    all_subkeys    <- unique(unlist(lapply(key_results, names)))
    consensus_list <- list()
    
    for (subkey in all_subkeys) {
      all_taxa      <- unlist(lapply(key_results, function(x) x[[subkey]]))
      taxon_counts  <- table(all_taxa)
      consensus_taxa <- names(taxon_counts)[taxon_counts >= (n_iterations * consensus_threshold)]
      if (length(consensus_taxa) > 0) consensus_list[[subkey]] <- consensus_taxa
    }
    
    if (length(consensus_list) >= 2) combined_results[[key]] <- consensus_list
  }
  
  return(combined_results)
}

# =========================================================================
# PARALLEL PROCESSING SETUP
# =========================================================================

setup_parallel_processing <- function() {
  tryCatch({
    n_cores <- max(1, parallel::detectCores() - 1)
    future::plan(future::multisession, workers = n_cores)
    return(TRUE)
  }, error = function(e) {
    future::plan(future::sequential)
    return(FALSE)
  })
}

# =========================================================================
# PLOTTING FUNCTIONS
# =========================================================================

create_venn_plot <- function(data_list, title, colors = NULL,
                             fill_limits = NULL, plot_type = "timepoint") {
  if (length(data_list) < 2) return(NULL)
  
  if (is.null(colors)) {
    colors <- c("#EDF8E9", "#BAE4B3", "#74C476", "#238B45")
  }
  
  label_size <- if (plot_type == "timepoint") 4 else 3
  
  p <- ggVennDiagram(
    data_list,
    label               = "both",
    label_alpha         = 1,
    label_color         = "black",
    label_geom          = "text",
    label_percent_digit = 2,
    label_size          = label_size,
    set_size            = 4,
    edge_lty            = "solid",
    edge_size           = 0.3,
    show_percentage     = TRUE,
    show_quantity       = TRUE
  ) +
    scale_fill_gradientn(colors = colors, limits = fill_limits) +
    labs(title = title) +
    theme_minimal() +
    theme(
      panel.grid      = element_blank(),
      axis.ticks      = element_blank(),
      plot.title      = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "bottom",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 8),
      axis.title      = element_blank(),
      axis.text       = element_blank()
    )
  
  return(p)
}

get_region_counts <- function(results_list) {
  all_counts <- unlist(lapply(results_list, function(group) {
    if (length(group) < 2) return(NULL)
    venn   <- ggVennDiagram::process_data(ggVennDiagram::Venn(group))
    ggVennDiagram::venn_region(venn)$count
  }))
  if (length(all_counts) == 0) return(NULL)
  c(min(all_counts, na.rm = TRUE), max(all_counts, na.rm = TRUE))
}

create_all_plots <- function(combined_results, sample_types, time_points) {
  limits_timepoint  <- get_region_counts(combined_results$by_timepoint)
  limits_sampletype <- get_region_counts(combined_results$by_sampletype)
  
  plots_timepoint  <- list()
  plots_sampletype <- list()
  
  for (tp in names(combined_results$by_timepoint)) {
    if (length(combined_results$by_timepoint[[tp]]) < 2) next
    plot <- create_venn_plot(
      combined_results$by_timepoint[[tp]],
      title       = tp,
      colors      = NULL,
      fill_limits = limits_timepoint,
      plot_type   = "timepoint"
    )
    if (!is.null(plot)) plots_timepoint[[tp]] <- plot
  }
  
  for (st in names(combined_results$by_sampletype)) {
    if (length(combined_results$by_sampletype[[st]]) < 2) next
    colors <- if (st %in% names(color_gradients)) color_gradients[[st]] else
      c("#D1D5DB", "#9CA3AF", "#6B7280", "#374151")
    plot <- create_venn_plot(
      combined_results$by_sampletype[[st]],
      title       = st,
      colors      = colors,
      fill_limits = limits_sampletype,
      plot_type   = "sampletype"
    )
    if (!is.null(plot)) plots_sampletype[[st]] <- plot
  }
  
  list(timepoint = plots_timepoint, sampletype = plots_sampletype)
}

# =========================================================================
# STEP 1: LOAD AND PREPARE DATA
# =========================================================================

physeq <- readRDS(phyloseq_path)

sample_metadata_clean <- process_sample_metadata(physeq)

sample_types <- sort(unique(sample_metadata_clean$group_type[!is.na(sample_metadata_clean$group_type)]))
time_points  <- sort(unique(sample_metadata_clean$group_time[!is.na(sample_metadata_clean$group_time)]))

sample_data_updated <- sample_metadata_clean %>%
  tibble::column_to_rownames("sample_id")
sample_data(physeq) <- sample_data(sample_data_updated)

# =========================================================================
# STEP 2: RAREFACTION SETUP
# =========================================================================

parallel_available <- setup_parallel_processing()

otu_mat <- otu_table(physeq) %>% as("matrix")
if (taxa_are_rows(physeq)) otu_mat <- t(otu_mat)

samples_sufficient <- rownames(otu_mat)[rowSums(otu_mat) >= rarefaction_depth]
otu_mat_filtered   <- otu_mat[samples_sufficient, ]

# =========================================================================
# STEP 3: RUN RAREFACTION ANALYSIS
# =========================================================================

if (parallel_available) {
  rarefaction_results <- tryCatch({
    future.apply::future_lapply(1:n_iterations, function(i) {
      rarefaction_for_venn(i, otu_mat_filtered, rarefaction_depth,
                           physeq, sample_types, time_points)
    }, future.seed = TRUE)
  }, error = function(e) {
    future::plan(future::sequential)
    lapply(1:n_iterations, function(i) {
      rarefaction_for_venn(i, otu_mat_filtered, rarefaction_depth,
                           physeq, sample_types, time_points)
    })
  })
} else {
  rarefaction_results <- lapply(1:n_iterations, function(i) {
    rarefaction_for_venn(i, otu_mat_filtered, rarefaction_depth,
                         physeq, sample_types, time_points)
  })
}

future::plan(future::sequential)

combined_by_timepoint  <- combine_rarefaction_results(rarefaction_results, "by_timepoint")
combined_by_sampletype <- combine_rarefaction_results(rarefaction_results, "by_sampletype")

combined_results <- list(
  by_timepoint  = combined_by_timepoint,
  by_sampletype = combined_by_sampletype
)

# =========================================================================
# STEP 4: CREATE PLOTS
# =========================================================================

plots <- create_all_plots(combined_results, sample_types, time_points)

# =========================================================================
# STEP 5: COMBINED FIGURES
# =========================================================================

# Figure 1: Sample types across timepoints
if (length(plots$timepoint) >= 5) {
  ordered_tp <- plots$timepoint[intersect(time_points, names(plots$timepoint))]
  
  fig_timepoint <- wrap_plots(ordered_tp, ncol = 5) +
    plot_annotation(
      title    = "Shared and unique genera between compartments at each timepoint",
      subtitle = paste0(
        "Rarefied to ", scales::comma(rarefaction_depth), " reads, ",
        n_iterations, " iterations, ",
        consensus_threshold * 100, "% consensus threshold. ",
        "Panels A-E: Week 0, 2, 4, 6, 8."
      ),
      tag_levels = "A"
    ) &
    theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9,  hjust = 0.5),
      plot.tag      = element_text(size = 10, face = "bold")
    )
  
  ggsave(
    filename = file.path(output_dir, "VennDiagram_Genus_TypesAcrossTime.png"),
    plot     = fig_timepoint,
    width    = 24,
    height   = 8,
    units    = "in",
    dpi      = 300
  )
}

# Figure 2: Temporal dynamics within each sample type
if (length(plots$sampletype) >= 3) {
  ordered_st <- plots$sampletype[intersect(sample_types, names(plots$sampletype))]
  
  fig_sampletype <- wrap_plots(ordered_st, ncol = 3) +
    plot_annotation(
      title    = "Temporal dynamics of shared and unique genera within each compartment",
      subtitle = paste0(
        "Rarefied to ", scales::comma(rarefaction_depth), " reads, ",
        n_iterations, " iterations, ",
        consensus_threshold * 100, "% consensus threshold. ",
        "Panels A-C: Intestine, Sediment, Water."
      ),
      tag_levels = "A"
    ) &
    theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9,  hjust = 0.5),
      plot.tag      = element_text(size = 10, face = "bold")
    )
  
  ggsave(
    filename = file.path(output_dir, "VennDiagram_Genus_TimeWithinTypes.png"),
    plot     = fig_sampletype,
    width    = 18,
    height   = 8,
    units    = "in",
    dpi      = 300
  )
}
