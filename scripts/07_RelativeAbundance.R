# ============================================================
# Script: RelativeAbundance.R
# Description: Alluvial diagram showing genus-level temporal
#              dynamics across Intestine, Sediment, and Water
#              compartments in Litopenaeus vannamei production
#              ponds.
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
  "phyloseq", "vegan", "ggplot2", "ggalluvial", "patchwork", "ggtext",
  "janitor", "tidyverse", "scales"
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
key_packages <- c("phyloseq", "ggplot2", "ggalluvial", "dplyr", "tidyr")
for (pkg in key_packages) cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

# File paths
phyloseq_path <- "step4_phyloseq_object.rds"
output_dir    <- "output/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Sample filtering
filter_ponds      <- c("C3", "C6", "D1", "D4")
filter_timepoints <- c("Week 0", "Week 2", "Week 4", "Week 6", "Week 8")
filter_shrimp_max <- 5

# Taxonomy
min_abundance <- 0.02   # Taxa below this threshold grouped as "Other"

# Placeholder palette — publication colors not included in public release
genus_palette <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
  "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
  "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02",
  "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
  "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
  "#A50026", "#D73027", "#F46D43", "#808080"
)

set.seed(54)

# =========================================================================
# UTILITY FUNCTIONS
# =========================================================================

clean_taxonomy <- function(tax_names) {
  tax_names <- gsub("^[a-z]__", "", tax_names)
  tax_names <- gsub("_.*", "", tax_names)
  tax_names[is.na(tax_names) | tax_names == ""] <- "Unidentified"
  return(tax_names)
}

process_metadata <- function(metadata_raw) {
  metadata_raw %>%
    dplyr::filter(trial == 1) %>%
    mutate(
      type = case_when(
        type == "Gut"   ~ "Intestine",
        type == "Water" ~ "Water",
        type == "Soil"  ~ "Sediment",
        TRUE ~ as.character(type)
      ),
      timepoint = case_when(
        timepoint == "T0" ~ "Week 0",
        timepoint == "T1" ~ "Week 2",
        timepoint == "T2" ~ "Week 4",
        timepoint == "T3" ~ "Week 6",
        timepoint == "T4" ~ "Week 8",
        TRUE ~ as.character(timepoint)
      ),
      time_numeric = as.numeric(gsub("Week ", "", timepoint)),
      timepoint    = factor(timepoint, levels = c("Week 0", "Week 2", "Week 4", "Week 6", "Week 8")),
      type         = factor(type, levels = c("Intestine", "Sediment", "Water")),
      pond         = as.factor(pond),
      time         = factor(time_numeric, levels = c(0, 2, 4, 6, 8))
    )
}

assign_colors <- function(taxa_vector) {
  n_taxa <- length(taxa_vector)
  colors <- genus_palette[1:min(n_taxa, length(genus_palette))]
  if (n_taxa > length(genus_palette)) {
    colors <- rep(genus_palette, ceiling(n_taxa / length(genus_palette)))[1:n_taxa]
  }
  names(colors) <- taxa_vector
  if ("Other" %in% taxa_vector) colors["Other"] <- "#808080"
  return(colors)
}

# =========================================================================
# RELATIVE ABUNDANCE FUNCTIONS
# =========================================================================

calculate_relative_abundance <- function(otu_matrix, tax_df, metadata_df) {
  taxa_counts <- list()
  
  for (i in 1:nrow(otu_matrix)) {
    sample_id     <- rownames(otu_matrix)[i]
    sample_counts <- otu_matrix[i, ]
    nonzero_idx   <- which(sample_counts > 0)
    
    for (j in nonzero_idx) {
      otu_id      <- colnames(otu_matrix)[j]
      count_value <- sample_counts[j]
      
      if (otu_id %in% rownames(tax_df)) {
        taxon <- tax_df[otu_id, "Genus"]
        taxa_counts[[length(taxa_counts) + 1]] <- list(
          sample = sample_id, taxon = taxon, count = count_value
        )
      }
    }
  }
  
  if (length(taxa_counts) == 0) stop("No taxa counts extracted.")
  
  taxa_df <- do.call(rbind, lapply(taxa_counts, function(x) {
    data.frame(sample = x$sample, genus = x$taxon, count = x$count,
               stringsAsFactors = FALSE)
  }))
  
  agg_counts <- aggregate(count ~ sample + genus, data = taxa_df, sum)
  result     <- left_join(agg_counts, metadata_df, by = "sample")
  return(result)
}

process_relative_abundance <- function(raw_counts) {
  raw_counts %>%
    group_by(sample) %>%
    mutate(
      total_abundance = sum(count, na.rm = TRUE),
      abundance       = ifelse(total_abundance > 0, count / total_abundance, 0)
    ) %>%
    ungroup() %>%
    mutate(
      timepoint = factor(timepoint, levels = c("Week 0", "Week 2", "Week 4", "Week 6", "Week 8")),
      type      = factor(type,      levels = c("Intestine", "Sediment", "Water")),
      time      = factor(time,      levels = c(0, 2, 4, 6, 8))
    )
}

# =========================================================================
# ALLUVIAL PLOT FUNCTION
# =========================================================================

create_alluvial_plot <- function(data, sample_type, min_abund = 0.02) {
  
  filtered <- data %>% dplyr::filter(type == sample_type, abundance > 0)
  if (nrow(filtered) == 0) {
    warning(paste("No data for", sample_type))
    return(NULL)
  }
  
  # Summarize to timepoint x genus proportions
  plot_data <- filtered %>%
    group_by(timepoint, genus) %>%
    summarise(total = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(timepoint) %>%
    mutate(proportion = total / sum(total, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(genus) %>%
    mutate(max_prop = max(proportion, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(taxon = ifelse(max_prop >= min_abund, genus, "Other")) %>%
    group_by(timepoint, taxon) %>%
    summarise(proportion = sum(proportion, na.rm = TRUE), .groups = "drop")
  
  taxa_order <- plot_data %>%
    group_by(taxon) %>%
    summarise(mean_prop = mean(proportion), .groups = "drop") %>%
    arrange(desc(mean_prop)) %>%
    pull(taxon)
  
  plot_data <- plot_data %>%
    mutate(taxon = factor(taxon, levels = taxa_order))
  
  colors <- assign_colors(levels(plot_data$taxon))
  
  ggplot(plot_data, aes(x = timepoint, y = proportion,
                        alluvium = taxon, stratum = taxon, fill = taxon)) +
    geom_alluvium(alpha = 0.7, width = 0.5, color = NA, knot.pos = 0.45) +
    geom_stratum(width = 0.6, alpha = 1.0, color = NA) +
    geom_text(stat = "stratum", aes(label = taxon), size = 3, fontface = "plain") +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0.01, 0.01)) +
    scale_x_discrete(expand = c(0.05, 0.05)) +
    labs(
      title = sample_type,
      x     = "Sampling Timepoint",
      y     = "Relative Abundance",
      fill  = "Genus"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5),
      axis.title.x    = element_blank(),
      axis.title.y    = element_text(face = "bold"),
      axis.text.x     = element_text(color = "black", face = "bold"),
      axis.text.y     = element_text(color = "black"),
      legend.title    = element_text(face = "bold"),
      panel.grid      = element_blank(),
      legend.position = "right"
    )
}

# =========================================================================
# STEP 1: LOAD AND PREPARE DATA
# =========================================================================

physeq <- readRDS(phyloseq_path)

sample_metadata_raw   <- phyloseq::sample_data(physeq) %>%
  as_tibble(rownames = "sample") %>%
  janitor::clean_names()

sample_metadata_clean <- process_metadata(sample_metadata_raw)

filtered_samples <- sample_metadata_clean %>%
  dplyr::filter(
    pond %in% filter_ponds,
    timepoint %in% filter_timepoints,
    is.na(shrimp) | shrimp <= filter_shrimp_max
  ) %>%
  pull(sample)

physeq_filtered <- phyloseq::prune_samples(filtered_samples, physeq)

sample_metadata <- process_metadata(
  phyloseq::sample_data(physeq_filtered) %>%
    as_tibble(rownames = "sample") %>%
    janitor::clean_names()
)

# =========================================================================
# STEP 2: CALCULATE RELATIVE ABUNDANCE
# =========================================================================

otu_mat <- phyloseq::otu_table(physeq_filtered) %>% as("matrix")
if (phyloseq::taxa_are_rows(physeq_filtered)) otu_mat <- t(otu_mat)

tax_df <- as.data.frame(phyloseq::tax_table(physeq_filtered))
for (col in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) {
  if (col %in% colnames(tax_df)) tax_df[[col]] <- clean_taxonomy(tax_df[[col]])
}

sample_metadata_filtered <- sample_metadata %>%
  dplyr::filter(sample %in% rownames(otu_mat))

raw_counts <- calculate_relative_abundance(otu_mat, tax_df, sample_metadata_filtered)
genus_data <- process_relative_abundance(raw_counts)

# =========================================================================
# STEP 3: CREATE ALLUVIAL PLOTS
# =========================================================================

plot_intestine <- create_alluvial_plot(genus_data, "Intestine", min_abund = min_abundance)
plot_sediment  <- create_alluvial_plot(genus_data, "Sediment",  min_abund = min_abundance)
plot_water     <- create_alluvial_plot(genus_data, "Water",     min_abund = min_abundance)

# =========================================================================
# STEP 4: COMBINED FIGURE
# =========================================================================

combined_plot <- plot_intestine / plot_sediment / plot_water +
  plot_annotation(
    title    = "Genus-level temporal dynamics in Litopenaeus vannamei production ponds",
    subtitle = paste0(
      "Relative abundance of bacterial genera across Intestine, Sediment, and Water compartments ",
      "over an 8-week production cycle. ",
      "Genera representing less than ", min_abundance * 100, "% relative abundance grouped as 'Other'."
    ),
    tag_levels = "A"
  ) &
  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    plot.tag      = element_text(size = 12, face = "bold")
  )

ggsave(
  filename = file.path(output_dir, "Alluvial_Genus.png"),
  plot     = combined_plot,
  width    = 14,
  height   = 18,
  units    = "in",
  dpi      = 300
)
