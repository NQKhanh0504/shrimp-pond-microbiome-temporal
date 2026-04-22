# ============================================================
# Script: 05_BetaDiversity.R
# Description: Beta diversity analysis using iterative
#              rarefaction and PCoA ordination across
#              Intestine, Sediment, and Water compartments
#              in Litopenaeus vannamei production ponds.
#              Distance metrics: Jaccard, Aitchison,
#              Weighted UniFrac, Unweighted UniFrac.
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
  "phyloseq", "vegan", "dplyr", "ggplot2", "patchwork", "ggtext",
  "stringr", "parallel", "future.apply", "tidyverse", "tidyr", "janitor",
  "ape", "readr", "multcompView", "scales"
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
key_packages <- c("phyloseq", "vegan", "ggplot2", "patchwork", "dplyr")
for (pkg in key_packages) cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

# File paths
phyloseq_path <- "step4_phyloseq_object.rds"
output_dir    <- "output/"
results_dir   <- file.path(output_dir, "Statistical_Results")
plots_dir     <- file.path(output_dir, "Plots")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,   recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
n_iterations      <- 100
rarefaction_depth <- 57841
correction_method <- "fdr"
set.seed(54)
options(contrasts = c("contr.sum", "contr.poly"))

# Colors and shapes — kept from publication version
type_colors <- c("Intestine" = "#E84646", "Sediment" = "#E69F00", "Water" = "#0072B2")
timepoint_colors <- c(
  "Week 0" = "#E84646", "Week 2" = "#f6bc41",
  "Week 4" = "#3CAEA3", "Week 6" = "#20639B", "Week 8" = "#173F5F"
)
type_shapes <- c("Intestine" = 21, "Sediment" = 24, "Water" = 22)

# =========================================================================
# UTILITY FUNCTIONS
# =========================================================================

format_pvalue <- function(p) {
  ifelse(is.na(p), "= NA",
         ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p)))
}

clean_sample_metadata <- function(metadata) {
  metadata %>%
    dplyr::filter(trial == "1") %>%
    dplyr::mutate(
      type_clean = dplyr::case_when(
        stringr::str_detect(tolower(type), "gut|intestin") ~ "Intestine",
        stringr::str_detect(tolower(type), "water")        ~ "Water",
        stringr::str_detect(tolower(type), "soil|sediment")~ "Sediment",
        TRUE ~ as.character(type)
      ),
      timepoint_clean = dplyr::case_when(
        stringr::str_detect(timepoint, "^T0$") ~ "Week 0",
        stringr::str_detect(timepoint, "^T1$") ~ "Week 2",
        stringr::str_detect(timepoint, "^T2$") ~ "Week 4",
        stringr::str_detect(timepoint, "^T3$") ~ "Week 6",
        stringr::str_detect(timepoint, "^T4$") ~ "Week 8",
        TRUE ~ as.character(timepoint)
      ),
      time_numeric = as.numeric(stringr::str_extract(timepoint_clean, "[0-9]+")),
      pond_clean   = as.factor(toupper(trimws(as.character(pond))))
    ) %>%
    dplyr::mutate(
      type_clean      = factor(type_clean,      levels = c("Intestine", "Water", "Sediment")),
      timepoint_clean = factor(timepoint_clean, levels = c("Week 0", "Week 2", "Week 4", "Week 6", "Week 8")),
      time_factor     = as.factor(time_numeric)
    )
}

create_pvalue_annotation <- function(method, permanova_results, betadisp_results) {
  if (!is.null(permanova_results[[method]])) {
    perm_res      <- permanova_results[[method]]
    type_p        <- perm_res$`Pr(>F)`[1]
    time_p        <- perm_res$`Pr(>F)`[2]
    interaction_p <- perm_res$`Pr(>F)`[3]
  } else {
    type_p <- time_p <- interaction_p <- NA
  }
  disp_type_p <- tryCatch(anova(betadisp_results[[method]]$type)$`Pr(>F)`[1], error = function(e) NA)
  disp_time_p <- tryCatch(anova(betadisp_results[[method]]$time)$`Pr(>F)`[1], error = function(e) NA)
  
  list(
    type_annotation = stringr::str_glue(
      "Compartment: P {format_pvalue(type_p)}\n",
      "Time: P {format_pvalue(time_p)}\n",
      "Interaction: P {format_pvalue(interaction_p)}\n",
      "Dispersion: P {format_pvalue(disp_type_p)}"
    ),
    time_annotation = stringr::str_glue("Dispersion: P {format_pvalue(disp_time_p)}")
  )
}

# =========================================================================
# STEP 1: LOAD AND PREPARE DATA
# =========================================================================

physeq <- readRDS(phyloseq_path)

sample_metadata <- phyloseq::sample_data(physeq) %>%
  as("data.frame") %>%
  tibble::as_tibble(rownames = "sample_id") %>%
  janitor::clean_names()

sample_metadata_clean <- clean_sample_metadata(sample_metadata)

phyloseq::sample_data(physeq) <- phyloseq::sample_data(
  sample_metadata_clean %>%
    dplyr::select(-sample_id) %>%
    as.data.frame() %>%
    `rownames<-`(sample_metadata_clean$sample_id)
)

# =========================================================================
# STEP 2: RAREFACTION SETUP
# =========================================================================

otu_mat <- phyloseq::otu_table(physeq) %>% as("matrix")
if (phyloseq::taxa_are_rows(physeq)) otu_mat <- t(otu_mat)

samples_sufficient <- names(rowSums(otu_mat))[rowSums(otu_mat) >= rarefaction_depth]
otu_mat_filtered   <- otu_mat[samples_sufficient, , drop = FALSE]
meta_filtered      <- sample_metadata_clean %>% dplyr::filter(sample_id %in% samples_sufficient)
tree               <- phyloseq::phy_tree(physeq, errorIfNULL = FALSE)

# =========================================================================
# STEP 3: ITERATIVE RAREFACTION
# =========================================================================

run_rarefaction_iteration <- function(i, otu_matrix, depth, tree_obj, meta_data, physeq_obj) {
  if (!require("vegan",    quietly = TRUE)) stop("vegan required")
  if (!require("phyloseq", quietly = TRUE)) stop("phyloseq required")
  if (!require("ape",      quietly = TRUE)) stop("ape required")
  if (!require("dplyr",    quietly = TRUE)) stop("dplyr required")
  
  set.seed(i * 100)
  rarefied        <- vegan::rrarefy(otu_matrix, depth)
  rarefied_pseudo <- rarefied + 1
  clr_mat         <- t(apply(rarefied_pseudo, 1, function(x) log(x / exp(mean(log(x))))))
  
  distances <- list(
    jaccard   = vegan::vegdist(rarefied > 0, method = "jaccard"),
    aitchison = vegan::vegdist(clr_mat,      method = "euclidean")
  )
  
  if (!is.null(tree_obj)) {
    tryCatch({
      temp_otu    <- phyloseq::otu_table(rarefied, taxa_are_rows = FALSE)
      temp_tax    <- phyloseq::tax_table(physeq_obj)[colnames(rarefied), ]
      temp_sample <- phyloseq::sample_data(
        meta_data %>%
          dplyr::select(-sample_id) %>%
          as.data.frame() %>%
          `rownames<-`(meta_data$sample_id)
      )
      present_taxa      <- colnames(rarefied)[colSums(rarefied) > 0]
      tree_tips         <- intersect(present_taxa, tree_obj$tip.label)
      if (length(tree_tips) > 2) {
        tree_pruned   <- ape::keep.tip(tree_obj, tree_tips)
        temp_phyloseq <- phyloseq::phyloseq(
          temp_otu[, tree_tips], temp_tax[tree_tips, ], temp_sample, tree_pruned
        )
        distances$wunifrac <- phyloseq::distance(temp_phyloseq, method = "wunifrac")
        distances$uunifrac <- phyloseq::distance(temp_phyloseq, method = "unifrac")
      }
    }, error = function(e) invisible(NULL))
  }
  
  list(iteration = i, distances = distances, sample_names = rownames(rarefied))
}

future::plan(future::multisession)
distance_iterations <- future.apply::future_lapply(1:n_iterations, function(i) {
  run_rarefaction_iteration(i, otu_mat_filtered, rarefaction_depth, tree, meta_filtered, physeq)
}, future.seed = TRUE)
future::plan(future::sequential)

# =========================================================================
# STEP 4: AVERAGE DISTANCE MATRICES
# =========================================================================

distance_methods <- c("jaccard", "aitchison", "wunifrac", "uunifrac")
sample_names     <- rownames(otu_mat_filtered)
n_samples        <- length(sample_names)

averaged_distances <- list()
for (method in distance_methods) {
  sum_matrix       <- matrix(0, nrow = n_samples, ncol = n_samples,
                             dimnames = list(sample_names, sample_names))
  valid_iterations <- 0
  for (iter_result in distance_iterations) {
    if (method %in% names(iter_result$distances)) {
      dist_matrix <- as.matrix(iter_result$distances[[method]])
      if (all(rownames(dist_matrix) %in% sample_names) &&
          all(colnames(dist_matrix) %in% sample_names)) {
        sum_matrix       <- sum_matrix + dist_matrix[sample_names, sample_names]
        valid_iterations <- valid_iterations + 1
      }
    }
  }
  if (valid_iterations > 0) {
    averaged_distances[[method]] <- as.dist(sum_matrix / valid_iterations)
  }
}
averaged_distances <- averaged_distances[!sapply(averaged_distances, is.null)]

# =========================================================================
# STEP 5: AVERAGE INTESTINE SAMPLES TO POND LEVEL
# =========================================================================
# Intestine has up to 5 shrimp replicates per pond x timepoint.
# Pond is the true experimental unit; averaging collapses pseudoreplicates
# before ordination and statistical testing.
# =========================================================================

average_intestine_samples <- function(physeq_obj) {
  otu_mat_orig <- as(phyloseq::otu_table(physeq_obj), "matrix")
  if (phyloseq::taxa_are_rows(physeq_obj)) otu_mat_orig <- t(otu_mat_orig)
  
  meta_df         <- as(phyloseq::sample_data(physeq_obj), "data.frame")
  intestine_idx   <- meta_df$type_clean == "Intestine"
  intestine_meta  <- meta_df[intestine_idx, ]
  intestine_otu   <- otu_mat_orig[intestine_idx, ]
  
  intestine_meta$group_id <- paste(intestine_meta$pond_clean,
                                   intestine_meta$timepoint_clean, sep = "_")
  averaged_otu  <- list()
  averaged_meta <- list()
  
  for (group in unique(intestine_meta$group_id)) {
    grp       <- intestine_meta$group_id == group
    avg_abund <- if (sum(grp) == 1) intestine_otu[grp, ] else colMeans(intestine_otu[grp, ])
    grp_meta  <- intestine_meta[grp, ][1, ]
    new_id    <- paste("Avg_Intestine", group, sep = "_")
    averaged_otu[[new_id]]  <- avg_abund
    averaged_meta[[new_id]] <- grp_meta
  }
  
  avg_otu_mat <- do.call(rbind, averaged_otu)
  other_otu   <- otu_mat_orig[!intestine_idx, ]
  other_meta  <- meta_df[!intestine_idx, ]
  avg_meta_df <- do.call(rbind, averaged_meta)
  all_cols    <- union(names(avg_meta_df), names(other_meta))
  
  final_otu  <- rbind(avg_otu_mat, other_otu)
  final_meta <- dplyr::bind_rows(
    avg_meta_df[, intersect(all_cols, names(avg_meta_df))],
    other_meta[,  intersect(all_cols, names(other_meta))]
  )
  rownames(final_meta) <- rownames(final_otu)
  
  phyloseq::phyloseq(
    phyloseq::otu_table(final_otu,    taxa_are_rows = FALSE),
    phyloseq::sample_data(final_meta),
    phyloseq::tax_table(physeq_obj),
    phyloseq::phy_tree(physeq_obj)
  )
}

temp_physeq <- phyloseq::phyloseq(
  phyloseq::otu_table(otu_mat_filtered, taxa_are_rows = FALSE),
  phyloseq::sample_data(
    meta_filtered %>%
      dplyr::select(-sample_id) %>%
      as.data.frame() %>%
      `rownames<-`(meta_filtered$sample_id)
  ),
  phyloseq::tax_table(physeq)[colnames(otu_mat_filtered), ],
  phyloseq::phy_tree(physeq)
)

physeq_averaged <- average_intestine_samples(temp_physeq)

# =========================================================================
# STEP 6: FINAL DISTANCE MATRICES ON AVERAGED DATA
# =========================================================================

otu_matrix_final <- as(phyloseq::otu_table(physeq_averaged), "matrix")
if (phyloseq::taxa_are_rows(physeq_averaged)) otu_matrix_final <- t(otu_matrix_final)

clr_final <- t(apply(otu_matrix_final + 1, 1, function(x) log(x / exp(mean(log(x))))))

final_distances <- list(
  jaccard   = vegan::vegdist(otu_matrix_final > 0, method = "jaccard"),
  aitchison = vegan::vegdist(clr_final,             method = "euclidean")
)

if (!is.null(phyloseq::phy_tree(physeq_averaged, errorIfNULL = FALSE))) {
  tryCatch({
    final_distances$wunifrac <- phyloseq::distance(physeq_averaged, method = "wunifrac")
  }, error = function(e) invisible(NULL))
  tryCatch({
    final_distances$uunifrac <- phyloseq::distance(physeq_averaged, method = "unifrac")
  }, error = function(e) invisible(NULL))
}
final_distances <- final_distances[!sapply(final_distances, is.null)]

meta_plot_final <- as(phyloseq::sample_data(physeq_averaged), "data.frame") %>%
  dplyr::mutate(
    type_clean      = factor(type_clean,      levels = c("Intestine", "Water", "Sediment")),
    timepoint_clean = factor(timepoint_clean, levels = c("Week 0", "Week 2", "Week 4", "Week 6", "Week 8"))
  )

# =========================================================================
# STEP 7: PCoA ORDINATION
# =========================================================================

pcoa_results  <- list()
pcoa_variance <- list()
for (method in names(final_distances)) {
  pcoa_results[[method]]  <- cmdscale(final_distances[[method]], k = 2, eig = TRUE)
  pcoa_variance[[method]] <- pcoa_results[[method]]$eig / sum(pcoa_results[[method]]$eig)
}

# =========================================================================
# STEP 8: STATISTICAL ANALYSIS
# =========================================================================

permanova_results <- list()
betadisp_results  <- list()
betadisp_summary  <- tibble::tibble()

for (method in names(final_distances)) {
  dist_mat          <- final_distances[[method]]
  dist_sample_names <- rownames(as.matrix(dist_mat))
  meta_permanova    <- meta_plot_final[match(dist_sample_names, rownames(meta_plot_final)), ]
  
  permanova_results[[method]] <- tryCatch({
    if ("pond_clean" %in% names(meta_permanova)) {
      vegan::adonis2(dist_mat ~ type_clean * timepoint_clean,
                     data = meta_permanova, permutations = 9999,
                     by = "terms", strata = meta_permanova$pond_clean)
    } else {
      vegan::adonis2(dist_mat ~ type_clean * timepoint_clean,
                     data = meta_permanova, permutations = 9999, by = "terms")
    }
  }, error = function(e) NULL)
  
  betadisp_type <- vegan::betadisper(dist_mat, meta_permanova$type_clean)
  betadisp_time <- vegan::betadisper(dist_mat, meta_permanova$timepoint_clean)
  betadisp_results[[method]] <- list(type = betadisp_type, time = betadisp_time)
  
  anova_type <- anova(betadisp_type)
  anova_time <- anova(betadisp_time)
  
  betadisp_summary <- dplyr::bind_rows(betadisp_summary, tibble::tibble(
    Test           = paste0(method, "_", c("Type", "Time")),
    F_value        = c(anova_type$`F value`[1], anova_time$`F value`[1]),
    p_value        = c(anova_type$`Pr(>F)`[1],  anova_time$`Pr(>F)`[1]),
    Interpretation = dplyr::case_when(
      p_value < 0.05 ~ "Unequal dispersion",
      TRUE           ~ "Equal dispersion"
    )
  ))
}

# =========================================================================
# STEP 9: SAVE STATISTICAL RESULTS
# =========================================================================

permanova_summary <- tibble::tibble(Effect = c("Sample Type", "Timepoint", "Type x Time"))
for (method in names(permanova_results)) {
  if (!is.null(permanova_results[[method]])) {
    res <- permanova_results[[method]]
    permanova_summary[[paste0(method, "_R2")]] <- round(res$R2[1:3] * 100, 1)
    permanova_summary[[paste0(method, "_p")]]  <- res$`Pr(>F)`[1:3]
  }
}

readr::write_csv(permanova_summary, file.path(results_dir, "PERMANOVA_Results.csv"))
readr::write_csv(betadisp_summary,  file.path(results_dir, "Beta_Dispersion_Results.csv"))
readr::write_csv(
  tibble::tibble(
    Method                 = names(final_distances),
    PCoA_PC1_PC2_Variance  = sapply(pcoa_variance, function(x) round(sum(x[1:2]) * 100, 1)),
    Rarefaction_Depth      = rarefaction_depth,
    Rarefaction_Iterations = n_iterations,
    FDR_Correction         = correction_method
  ),
  file.path(results_dir, "Ordination_Quality_Summary.csv")
)
readr::write_csv(meta_plot_final, file.path(results_dir, "Final_Metadata.csv"))

# =========================================================================
# STEP 10: CREATE PCoA PLOTS
# =========================================================================

method_titles <- list(
  jaccard   = "Jaccard presence-absence dissimilarity",
  aitchison = "Aitchison compositional dissimilarity",
  wunifrac  = "Weighted UniFrac phylogenetic dissimilarity",
  uunifrac  = "Unweighted UniFrac phylogenetic dissimilarity"
)

plot_list <- list()

for (method in names(pcoa_results)) {
  pc1_var   <- round(pcoa_variance[[method]][1] * 100, 1)
  pc2_var   <- round(pcoa_variance[[method]][2] * 100, 1)
  scores    <- pcoa_results[[method]]$points
  pcoa_data <- cbind(PC1 = scores[, 1], PC2 = scores[, 2], meta_plot_final)
  pval_ann  <- create_pvalue_annotation(method, permanova_results, betadisp_results)
  
  # Type plot
  plot_list[[paste0(method, "_type")]] <- ggplot(pcoa_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = type_clean, fill = type_clean, shape = type_clean),
               size = 4, alpha = 0.6, stroke = 1) +
    stat_ellipse(aes(fill = type_clean), level = 0.95, type = "norm",
                 alpha = 0.15, linewidth = 0.6, geom = "polygon") +
    annotate("text", x = Inf, y = -Inf, label = pval_ann$type_annotation,
             hjust = 1.05, vjust = -0.2, size = 3, fontface = "italic", lineheight = 0.9) +
    scale_color_manual(values = type_colors, name = "Compartment") +
    scale_fill_manual(values  = type_colors, name = "Compartment") +
    scale_shape_manual(values = type_shapes, name = "Compartment") +
    labs(
      title = method_titles[[method]],
      x     = paste0("PC1 (", pc1_var, "%)"),
      y     = paste0("PC2 (", pc2_var, "%)")
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0, size = 10),
      axis.title       = element_text(face = "bold"),
      axis.text        = element_text(color = "black"),
      legend.position  = "right",
      legend.title     = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  # Time plot
  plot_list[[paste0(method, "_time")]] <- ggplot(pcoa_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = timepoint_clean, fill = timepoint_clean, shape = type_clean),
               size = 4, alpha = 0.6, stroke = 1) +
    annotate("text", x = Inf, y = -Inf, label = pval_ann$time_annotation,
             hjust = 1.05, vjust = -0.5, size = 3, fontface = "italic", lineheight = 0.9) +
    scale_color_manual(values = timepoint_colors, name = "Timepoint") +
    scale_fill_manual(values  = timepoint_colors, guide = "none") +
    scale_shape_manual(values = type_shapes,       name = "Compartment") +
    labs(
      title = paste(method_titles[[method]], "- Temporal"),
      x     = paste0("PC1 (", pc1_var, "%)"),
      y     = paste0("PC2 (", pc2_var, "%)")
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0, size = 10),
      axis.title       = element_text(face = "bold"),
      axis.text        = element_text(color = "black"),
      legend.position  = "right",
      legend.title     = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# =========================================================================
# STEP 11: COMPOSITE FIGURES
# =========================================================================

# Figure 1: Aitchison + Jaccard
if (all(c("aitchison_type", "aitchison_time", "jaccard_type", "jaccard_time") %in% names(plot_list))) {
  fig_aj <- (plot_list[["aitchison_type"]] + plot_list[["aitchison_time"]]) /
    (plot_list[["jaccard_type"]]            + plot_list[["jaccard_time"]]) +
    plot_annotation(
      title    = "PCoA of microbial community composition - Aitchison and Jaccard dissimilarity",
      subtitle = paste0(
        "Intestine, Sediment, and Water compartments over an 8-week production cycle. ",
        "Rarefied to ", scales::comma(rarefaction_depth), " reads, ", n_iterations, " iterations. ",
        "Panels A-B: Aitchison (CLR-transformed), C-D: Jaccard (presence-absence)."
      ),
      tag_levels = "A"
    ) &
    theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 9,  hjust = 0),
      plot.tag      = element_text(size = 11, face = "bold")
    )
  
  ggsave(
    filename = file.path(plots_dir, "PCoA_Aitchison_Jaccard.png"),
    plot     = fig_aj,
    width    = 16, height = 12, units = "in", dpi = 300
  )
}

# Figure 2: Weighted + Unweighted UniFrac
if (all(c("wunifrac_type", "wunifrac_time", "uunifrac_type", "uunifrac_time") %in% names(plot_list))) {
  fig_uf <- (plot_list[["wunifrac_type"]] + plot_list[["wunifrac_time"]]) /
    (plot_list[["uunifrac_type"]]          + plot_list[["uunifrac_time"]]) +
    plot_annotation(
      title    = "PCoA of phylogenetic diversity - Weighted and Unweighted UniFrac",
      subtitle = paste0(
        "Intestine, Sediment, and Water compartments over an 8-week production cycle. ",
        "Rarefied to ", scales::comma(rarefaction_depth), " reads, ", n_iterations, " iterations. ",
        "Panels A-B: Weighted UniFrac (abundance-based), C-D: Unweighted UniFrac (presence-absence)."
      ),
      tag_levels = "A"
    ) &
    theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 9,  hjust = 0),
      plot.tag      = element_text(size = 11, face = "bold")
    )
  
  ggsave(
    filename = file.path(plots_dir, "PCoA_UniFrac.png"),
    plot     = fig_uf,
    width    = 16, height = 12, units = "in", dpi = 300
  )
}
