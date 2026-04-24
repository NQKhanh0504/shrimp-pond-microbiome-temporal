# ============================================================
# Script: 08_DifferentialAnalysis.R
# Description: ANCOM-BC2 differential abundance analysis for
#              Litopenaeus vannamei production pond microbiome
#              data. Intestine as reference group,
#              two pairwise comparisons: Water vs Intestine
#              and Sediment vs Intestine.
# Input: step4_phyloseq_object.rds
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

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_packages <- c("phyloseq", "ANCOMBC", "microbiome")
cran_packages <- c("tidyverse", "janitor", "ggrepel", "patchwork",
                   "ggtext", "readr", "scales", "parallel")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  library(pkg, character.only = TRUE)
}

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# =========================================================================
# SESSION INFO
# =========================================================================

cat("R Version:", R.version.string, "\n")
cat("Key Package Versions:\n")
for (pkg in c("phyloseq", "ANCOMBC", "tidyverse", "ggrepel")) {
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
results_dir <- file.path(output_dir, "Results")
plots_dir   <- file.path(output_dir, "Volcano_Plots")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)

# Statistical parameters
alpha_threshold      <- 0.05
prevalence_cutoff_high <- 0.1   # For Phylum and Class levels
prevalence_cutoff_low  <- 0.2   # For Order, Family, Genus levels
taxonomic_levels     <- c("Genus")

# Plotting parameters
max_labels_per_group <- 8   # Top taxa to label per group in volcano plot

# Colors and shapes — kept from publication version
type_colors <- c("Intestine" = "#E84646", "Sediment" = "#E69F00", "Water" = "#0072B2")
type_shapes <- c("Intestine" = 21,        "Sediment" = 24,        "Water" = 22)

# Parallel processing
n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)

# =========================================================================
# HELPER FUNCTIONS
# =========================================================================

clean_taxonomy_hierarchical <- function(tax_row) {
  clean_string <- function(s) {
    s <- gsub("^[a-z]__", "", s)
    s <- gsub("_.*", "", s)
    s <- gsub("^(Bacteria|Archaea)_", "", s)
    trimws(s)
  }
  
  genus  <- clean_string(tax_row["Genus"])
  family <- clean_string(tax_row["Family"])
  order  <- clean_string(tax_row["Order"])
  class  <- clean_string(tax_row["Class"])
  
  if (!is.na(genus) && genus != "" && !genus %in% c("Bacteria", "Archaea")) {
    return(genus)
  } else if (!is.na(family) && family != "" && !family %in% c("Bacteria", "Archaea")) {
    return(paste0("Unclassified_", family))
  } else if (!is.na(order) && order != "" && !order %in% c("Bacteria", "Archaea")) {
    return(paste0("Unclassified_", order))
  } else if (!is.na(class) && class != "" && !class %in% c("Bacteria", "Archaea")) {
    return(paste0("Unclassified_", class))
  } else {
    return("Unidentified")
  }
}

extract_significant_results <- function(ancombc_result, alpha = 0.05) {
  if (is.null(ancombc_result) || is.null(ancombc_result$res) ||
      nrow(ancombc_result$res) == 0) return(NULL)
  
  res_df <- ancombc_result$res
  
  res_df$NaN_Flag <- apply(
    res_df[, c("lfc_type_clean1", "lfc_type_clean2",
               "q_type_clean1",   "q_type_clean2")],
    1, function(row) any(is.nan(row) | is.infinite(row) | is.na(row))
  )
  res_df$Sensitivity_Pass <- res_df$passed_ss_type_clean1 & res_df$passed_ss_type_clean2
  
  # lfc_type_clean1 = log2(Water/Intestine)
  # lfc_type_clean2 = log2(Sediment/Intestine)
  comparisons <- list(
    "Water vs Intestine"   = list(lfc_col = "lfc_type_clean1", q_col = "q_type_clean1"),
    "Sediment vs Intestine" = list(lfc_col = "lfc_type_clean2", q_col = "q_type_clean2")
  )
  
  results_long <- do.call(rbind, lapply(names(comparisons), function(comp) {
    pair <- comparisons[[comp]]
    res_df %>%
      dplyr::select(taxon, NaN_Flag, Sensitivity_Pass,
                    all_of(c(pair$lfc_col, pair$q_col))) %>%
      dplyr::rename(lfc   = all_of(pair$lfc_col),
                    q_val = all_of(pair$q_col)) %>%
      dplyr::mutate(comparison = comp)
  }))
  
  sig_results <- results_long %>%
    dplyr::filter(q_val < alpha & !NaN_Flag & Sensitivity_Pass) %>%
    dplyr::arrange(q_val)
  
  list(full_results = results_long, significant_results = sig_results)
}

export_comprehensive_results <- function(ancombc_result, tax_level, output_dir) {
  if (is.null(ancombc_result) || is.null(ancombc_result$res)) return(NULL)
  
  res_df <- ancombc_result$res
  
  bias_corrected <- if (!is.null(ancombc_result$bias_correct_log_table)) {
    bc <- as.data.frame(ancombc_result$bias_correct_log_table)
    bc$taxon <- rownames(bc)
    bc
  } else NULL
  
  build_comparison_output <- function(lfc_col, se_col, w_col, p_col, q_col,
                                      diff_col, ss_col, comp_label, direction_pos) {
    out <- res_df %>%
      dplyr::select(taxon,
                    all_of(c(lfc_col, se_col,
                             grep(paste0("^", w_col, "$"), colnames(res_df), value = TRUE),
                             p_col, q_col, diff_col, ss_col))) %>%
      dplyr::rename(
        LFC                  = all_of(lfc_col),
        SE                   = all_of(se_col),
        W_statistic          = dplyr::matches(paste0("^", w_col, "$")),
        P_value              = all_of(p_col),
        Adjusted_P_value     = all_of(q_col),
        Differentially_Abundant = all_of(diff_col),
        Sensitivity_Pass     = all_of(ss_col)
      ) %>%
      dplyr::mutate(
        Comparison  = comp_label,
        Significant = Adjusted_P_value < 0.05 & Sensitivity_Pass,
        Direction   = dplyr::case_when(
          LFC > 0 ~ paste0(direction_pos, " enriched"),
          LFC < 0 ~ "Intestine enriched",
          TRUE    ~ "No change"
        )
      )
    if (!is.null(bias_corrected)) out <- dplyr::left_join(out, bias_corrected, by = "taxon")
    out
  }
  
  water_out <- build_comparison_output(
    "lfc_type_clean1", "se_type_clean1", "W_type_clean1",
    "p_type_clean1",   "q_type_clean1",  "diff_type_clean1", "passed_ss_type_clean1",
    "Water vs Intestine", "Water"
  )
  
  sediment_out <- build_comparison_output(
    "lfc_type_clean2", "se_type_clean2", "W_type_clean2",
    "p_type_clean2",   "q_type_clean2",  "diff_type_clean2", "passed_ss_type_clean2",
    "Sediment vs Intestine", "Sediment"
  )
  
  combined <- dplyr::bind_rows(water_out, sediment_out) %>%
    dplyr::mutate(Clean_Name = taxon)
  
  first_cols <- c("taxon", "Clean_Name", "Comparison", "LFC", "SE", "W_statistic",
                  "P_value", "Adjusted_P_value", "Significant", "Direction",
                  "Differentially_Abundant", "Sensitivity_Pass")
  combined <- combined %>%
    dplyr::select(all_of(c(first_cols, setdiff(names(combined), first_cols))))
  
  readr::write_csv(combined,
                   file.path(output_dir,
                             paste0("ANCOMBC_", tax_level, "_Comprehensive_Results.csv")))
  
  sig_only <- combined %>%
    dplyr::filter(Significant) %>%
    dplyr::arrange(Comparison, Adjusted_P_value)
  if (nrow(sig_only) > 0) {
    readr::write_csv(sig_only,
                     file.path(output_dir,
                               paste0("ANCOMBC_", tax_level, "_Significant_Only.csv")))
  }
  
  for (comp in c("Water vs Intestine", "Sediment vs Intestine")) {
    comp_data     <- dplyr::filter(combined, Comparison == comp)
    comp_filename <- gsub(" ", "_", comp)
    readr::write_csv(comp_data,
                     file.path(output_dir,
                               paste0("ANCOMBC_", tax_level, "_", comp_filename, ".csv")))
  }
  
  return(combined)
}

prepare_volcano_data <- function(data, comparison, sig_threshold = 0.05, max_labels = 5) {
  parts <- strsplit(comparison, " vs ")[[1]]
  type1 <- parts[1]
  type2 <- parts[2]
  
  comp_data <- data %>% dplyr::filter(comparison == !!comparison)
  
  # Negate Water vs Intestine so positive LFC = Water enriched
  if (comparison == "Water vs Intestine") {
    comp_data <- comp_data %>% dplyr::mutate(lfc = -lfc)
  }
  
  taxa_to_label_relaxed <- comp_data %>%
    dplyr::filter(q_val < sig_threshold, !NaN_Flag,
                  is.finite(lfc), is.finite(q_val), !is.na(Clean_Name)) %>%
    dplyr::mutate(
      importance = -log10(pmax(q_val, 1e-10)) * abs(lfc),
      group      = dplyr::if_else(lfc > 0, type1, type2)
    ) %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(dplyr::desc(importance), .by_group = TRUE) %>%
    dplyr::slice_head(n = max_labels) %>%
    dplyr::ungroup()
  
  taxa_to_label_strict <- dplyr::filter(taxa_to_label_relaxed, Sensitivity_Pass)
  
  taxa_to_label <- if (nrow(taxa_to_label_strict) >= 4) {
    taxa_to_label_strict$taxon
  } else if (nrow(taxa_to_label_relaxed) > 0) {
    taxa_to_label_relaxed$taxon
  } else {
    character(0)
  }
  
  plot_data <- comp_data %>%
    dplyr::filter(is.finite(lfc) & is.finite(q_val)) %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        q_val < sig_threshold & Sensitivity_Pass ~ "Significant",
        TRUE ~ "Not Significant"
      ),
      neg_log_q  = -log10(pmax(q_val, 1e-10)),
      label      = ifelse(taxon %in% taxa_to_label, Clean_Name, ""),
      sample_type = dplyr::if_else(lfc > 0, type1, type2)
    )
  
  list(plot_data = plot_data, type1 = type1, type2 = type2)
}

build_volcano_plot <- function(plot_list, sig_threshold = 0.05) {
  pd <- plot_list$plot_data
  
  ggplot(pd, aes(x = lfc, y = neg_log_q)) +
    geom_hline(yintercept = -log10(sig_threshold),
               linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_point(aes(color = sample_type, shape = sample_type,
                   fill  = ifelse(significance == "Significant", sample_type, NA)),
               size = 3, alpha = 0.6, stroke = 0.8) +
    ggrepel::geom_text_repel(
      aes(label = label),
      size              = 3,
      max.overlaps      = Inf,
      box.padding       = 0.4,
      point.padding     = 0.2,
      force             = 5,
      force_pull        = 0.1,
      min.segment.length = 0,
      segment.size      = 0.4,
      segment.color     = "black",
      segment.alpha     = 0.6,
      segment.linetype  = 1,
      segment.curvature = 0
    ) +
    scale_color_manual(values = type_colors, name = "Enriched in") +
    scale_fill_manual( values = type_colors, name = "Sample Type",
                       na.value = NA, guide = "none") +
    scale_shape_manual(values = type_shapes, name = "Enriched in") +
    scale_x_continuous(breaks = seq(-8, 8, by = 1), limits = c(-8, 8)) +
    scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
    labs(
      x = "Log Fold Change",
      y = "-log10(FDR-adjusted p-value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position  = "top",
      legend.title     = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# =========================================================================
# STEP 1: LOAD AND PREPARE DATA
# =========================================================================

physeq <- readRDS(physeq_path)

# Apply hierarchical taxonomic naming
tax_df <- as.data.frame(tax_table(physeq))
tax_df$Genus_Clean <- apply(tax_df, 1, clean_taxonomy_hierarchical)

new_tax_table <- tax_table(physeq)
for (i in seq_len(ncol(new_tax_table))) {
  new_tax_table[, i] <- gsub("^[a-z]__", "", new_tax_table[, i])
  new_tax_table[, i] <- gsub("_.*",       "", new_tax_table[, i])
  new_tax_table[, i][is.na(new_tax_table[, i]) | new_tax_table[, i] == ""] <- "Unidentified"
}
new_tax_table[, "Genus"] <- tax_df$Genus_Clean
tax_table(physeq) <- new_tax_table
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# Prepare metadata
meta_df <- as(sample_data(physeq), "data.frame") %>%
  janitor::clean_names() %>%
  dplyr::filter(trial == "1") %>%
  dplyr::mutate(
    type_clean = dplyr::case_when(
      stringr::str_detect(tolower(type), "gut|intestin") ~ "Intestine",
      stringr::str_detect(tolower(type), "water")        ~ "Water",
      stringr::str_detect(tolower(type), "soil|sediment") ~ "Sediment",
      TRUE ~ as.character(type)
    ),
    timepoint_clean = dplyr::case_when(
      timepoint == "T0" ~ "Week_0",
      timepoint == "T1" ~ "Week_2",
      timepoint == "T2" ~ "Week_4",
      timepoint == "T3" ~ "Week_6",
      timepoint == "T4" ~ "Week_8",
      TRUE ~ as.character(timepoint)
    ),
    pond_clean  = toupper(trimws(as.character(pond))),
    type_clean  = factor(type_clean,
                         levels = c("Intestine", "Water", "Sediment")),
    timepoint_clean = factor(timepoint_clean,
                             levels = c("Week_0", "Week_2", "Week_4",
                                        "Week_6", "Week_8"))
  )

sample_data(physeq)                    <- sample_data(meta_df)
rownames(sample_data(physeq))          <- rownames(meta_df)

physeq_final <- subset_samples(
  physeq,
  timepoint_clean %in% c("Week_0", "Week_2", "Week_4", "Week_6", "Week_8")
)
physeq_final <- prune_taxa(taxa_sums(physeq_final) > 0, physeq_final)

# =========================================================================
# STEP 2: ANCOM-BC2 ANALYSIS
# =========================================================================

sample_type_results  <- list()
ancombc_full_results <- list()

for (tax_level in taxonomic_levels) {
  prv_cutoff <- if (tax_level %in% c("Phylum", "Class")) {
    prevalence_cutoff_high
  } else {
    prevalence_cutoff_low
  }
  
  ancombc_result <- tryCatch({
    ANCOMBC::ancombc2(
      data        = physeq_final,
      tax_level   = tax_level,
      fix_formula = "type_clean + timepoint_clean",
      rand_formula = "(1 | pond_clean)",
      p_adj_method = "fdr",
      prv_cut     = prv_cutoff,
      lib_cut     = 0,
      s0_perc     = 0.5,
      group       = "type_clean",
      struc_zero  = TRUE,
      neg_lb      = TRUE,
      alpha       = alpha_threshold,
      n_cl        = n_cores,
      verbose     = FALSE,
      global      = FALSE,
      pairwise    = TRUE,
      dunnet      = FALSE,
      trend       = FALSE,
      pseudo_sens = TRUE
    )
  }, error = function(e) {
    warning("ANCOM-BC2 failed for ", tax_level, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(ancombc_result) && !is.null(ancombc_result$res)) {
    ancombc_full_results[[tax_level]] <- ancombc_result
    
    results <- extract_significant_results(ancombc_result, alpha_threshold)
    
    if (!is.null(results)) {
      readr::write_csv(
        results$full_results,
        file.path(results_dir, paste0("ANCOMBC_Type_", tax_level, "_Full.csv"))
      )
      if (nrow(results$significant_results) > 0) {
        readr::write_csv(
          results$significant_results,
          file.path(results_dir, paste0("ANCOMBC_Type_", tax_level, "_Significant.csv"))
        )
      }
      sample_type_results[[tax_level]] <- results
    }
  }
  gc()
}

# =========================================================================
# STEP 3: ADD CLEAN NAMES AND EXPORT CSV
# =========================================================================

for (lvl in names(sample_type_results)) {
  if (!is.null(sample_type_results[[lvl]]$full_results)) {
    sample_type_results[[lvl]]$full_results$Clean_Name <-
      sample_type_results[[lvl]]$full_results$taxon
  }
}

for (tax_level in names(ancombc_full_results)) {
  export_comprehensive_results(ancombc_full_results[[tax_level]], tax_level, results_dir)
}

# =========================================================================
# STEP 4: VOLCANO PLOTS
# =========================================================================

sig_threshold <- 0.05
max_labels    <- max_labels_per_group

water_intestine_data <- prepare_volcano_data(
  sample_type_results$Genus$full_results,
  "Water vs Intestine",
  max_labels = max_labels
)

soil_intestine_data <- prepare_volcano_data(
  sample_type_results$Genus$full_results,
  "Sediment vs Intestine",
  max_labels = max_labels
)

volcano_water_intestine <- build_volcano_plot(water_intestine_data, sig_threshold)
volcano_soil_intestine  <- build_volcano_plot(soil_intestine_data,  sig_threshold)

# Save individual plots
ggsave(file.path(plots_dir, "Volcano_Genus_Water_vs_Intestine.png"),
       volcano_water_intestine, width = 8, height = 7, dpi = 300)

ggsave(file.path(plots_dir, "Volcano_Genus_Sediment_vs_Intestine.png"),
       volcano_soil_intestine, width = 8, height = 7, dpi = 300)

# Combined figure
plot_main_title <- paste0(
  "Volcano plot of genus-level differential abundance in ",
  "*Litopenaeus vannamei* production ponds"
)
plot_subtitle <- paste0(
  "ANCOM-BC2 analysis comparing Intestine, Sediment, and Water compartments ",
  "during an 8-week production cycle. FDR-adjusted p-value threshold: 0.05. ",
  "Opaque points passed both significance and sensitivity analysis. ",
  "Labels show top genera per compartment. ",
  "Panel (**A**) Water vs Intestine; Panel (**B**) Sediment vs Intestine."
)

combined_volcano <- (volcano_water_intestine + volcano_soil_intestine) +
  patchwork::plot_annotation(
    title      = plot_main_title,
    subtitle   = plot_subtitle,
    tag_levels = "A"
  ) &
  theme(
    plot.title    = ggtext::element_markdown(face = "bold", size = 12, hjust = 0),
    plot.subtitle = ggtext::element_markdown(size = 9, hjust = 0, lineheight = 1.3,
                                             margin = margin(t = 5, b = 10)),
    plot.tag      = element_text(size = 12, face = "bold")
  )

ggsave(file.path(plots_dir, "ANCOMBC_Volcano_Genus.png"),
       combined_volcano, width = 16, height = 8, dpi = 300)

cat("Analysis complete.\n")
cat("Results saved to:", results_dir, "\n")
cat("Plots saved to:",   plots_dir,   "\n")
