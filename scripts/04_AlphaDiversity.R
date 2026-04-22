# ============================================================
# Script: 04_AlphaDiversity.R
# Description: Alpha diversity analysis for shrimp microbiome
#              temporal dynamics experiment.
#              Produces: Observed Richness, Shannon Entropy,
#              Berger-Parker Dominance, and Faith's PD raincloud
#              plots with LME statistics and CLD annotations.
# Input: phyloseq object (step4_phyloseq_object.rds)
# Note: Plot aesthetics have been simplified for public release.
#       Final publication figures used customized formatting
#       not included in this repository.
# ============================================================

# =========================================================================
# ENVIRONMENT CLEANUP AND INITIALIZATION
# =========================================================================

cleanup_environment <- function() {
  tryCatch({
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(.Random.seed, envir = .GlobalEnv)
    }
  }, error = function(e) invisible(NULL))
  
  tryCatch({
    if ("future" %in% loadedNamespaces()) {
      future::plan(future::sequential)
    }
  }, error = function(e) invisible(NULL))
  
  base_packages <- c("stats", "graphics", "grDevices", "utils", "datasets",
                     "methods", "base", "tools", "parallel")
  
  loaded_ns    <- loadedNamespaces()
  ns_to_unload <- setdiff(loaded_ns, base_packages)
  
  if (length(ns_to_unload) == 0) {
    gc(verbose = FALSE)
    return(invisible(NULL))
  }
  
  get_importers <- function(pkg) {
    importers <- character(0)
    for (ns in loadedNamespaces()) {
      if (ns == pkg) next
      imports <- tryCatch(getNamespaceImports(ns), error = function(e) list())
      if (pkg %in% unlist(imports)) importers <- c(importers, ns)
    }
    return(importers)
  }
  
  unloaded       <- character(0)
  max_iterations <- length(ns_to_unload) * 2
  iteration      <- 0
  
  while (length(ns_to_unload) > 0 && iteration < max_iterations) {
    iteration     <- iteration + 1
    made_progress <- FALSE
    
    for (pkg in ns_to_unload) {
      importers        <- get_importers(pkg)
      importers_loaded <- intersect(importers, loadedNamespaces())
      
      if (length(importers_loaded) == 0) {
        success <- tryCatch({
          pkg_attached <- paste0("package:", pkg)
          if (pkg_attached %in% search()) {
            detach(pkg_attached, character.only = TRUE, unload = FALSE, force = TRUE)
          }
          unloadNamespace(pkg)
          TRUE
        }, error = function(e) FALSE)
        
        if (success) {
          unloaded      <- c(unloaded, pkg)
          ns_to_unload  <- setdiff(ns_to_unload, pkg)
          made_progress <- TRUE
        }
      }
    }
    
    if (!made_progress) break
  }
  
  gc(verbose = FALSE)
}

cleanup_environment()

# =========================================================================
# PACKAGE INSTALLATION AND LOADING
# =========================================================================

cran_packages <- c(
  "tidyverse", "vegan", "nlme", "emmeans", "janitor", "picante", "ape",
  "parallel", "future.apply", "multcomp", "multcompView", "ggpubr",
  "patchwork", "broom", "ggtext", "ggrepel", "gghalves", "scales"
)

bioc_packages <- c("phyloseq")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  library(pkg, character.only = TRUE)
}

# =========================================================================
# SESSION INFO
# =========================================================================

cat("R Version:", R.version.string, "\n")
cat("Key Package Versions:\n")
key_packages <- c("phyloseq", "vegan", "nlme", "emmeans", "multcomp", "ggplot2", "dplyr", "tidyr")
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
random_seed       <- 54
n_iterations      <- 100
rarefaction_depth <- 57841
alpha_level       <- 0.05

# Parallel processing configuration
# Options: "future" (cross-platform, RStudio-safe), "mclapply" (Unix/Mac), "parLapply" (Windows)
parallel_method <- "future"

# Sample filtering
filter_ponds      <- c("C3", "C6", "D1", "D4")
filter_timepoints <- c("Week 0", "Week 2", "Week 4", "Week 6", "Week 8")
filter_shrimp_max <- 5

# Diversity metrics
primary_metrics <- c("observed_richness", "shannon", "berger_parker", "faith_pd")

display_names <- c(
  "observed_richness" = "Observed Richness",
  "shannon"           = "Shannon Entropy",
  "faith_pd"          = "Faith's Phylogenetic Diversity",
  "berger_parker"     = "Berger-Parker Dominance"
)

# Placeholder palettes — publication colors not included in public release
color_palette <- c("Intestine" = "#E41A1C", "Sediment" = "#E69F00", "Water" = "#377EB8")
type_shapes   <- c("Intestine" = 21, "Sediment" = 24, "Water" = 22)

set.seed(random_seed)
options(contrasts = c("contr.sum", "contr.poly"))

# =========================================================================
# UTILITY FUNCTIONS
# =========================================================================

format_pvalue <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))
}

process_metadata <- function(metadata_raw) {
  metadata_raw %>%
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
      time_numeric = as.numeric(stringr::str_remove(timepoint, "Week ")),
      timepoint    = factor(timepoint, levels = c("Week 0", "Week 2", "Week 4", "Week 6", "Week 8")),
      type         = factor(type, levels = c("Intestine", "Sediment", "Water")),
      pond         = as.factor(pond),
      time_factor  = factor(time_numeric, levels = c(0, 2, 4, 6, 8))
    )
}

filter_samples <- function(metadata_clean, ponds = NULL, timepoints = NULL, shrimp_max = NULL) {
  filtered <- metadata_clean
  if (!is.null(ponds))      filtered <- dplyr::filter(filtered, pond %in% ponds)
  if (!is.null(timepoints)) filtered <- dplyr::filter(filtered, timepoint %in% timepoints)
  if (!is.null(shrimp_max)) filtered <- dplyr::filter(filtered, is.na(shrimp) | shrimp <= shrimp_max)
  return(filtered)
}

calculate_alpha_diversity_one_iter <- function(i, otu_matrix, depth, tree) {
  set.seed(i * 100)
  rarefied <- vegan::rrarefy(otu_matrix, depth)
  
  base_metrics <- data.frame(
    sample_id         = rownames(rarefied),
    iteration         = i,
    observed_richness = vegan::specnumber(rarefied),
    shannon           = vegan::diversity(rarefied, index = "shannon")
  )
  
  pd_result <- picante::pd(rarefied, tree, include.root = FALSE)
  faith_pd  <- data.frame(sample_id = rownames(pd_result), faith_pd = pd_result$PD)
  
  additional_metrics <- do.call(rbind, lapply(1:nrow(rarefied), function(x) {
    row_data    <- rarefied[x, ]
    total_reads <- sum(row_data)
    data.frame(
      sample_id     = rownames(rarefied)[x],
      berger_parker = ifelse(total_reads > 0, max(row_data) / total_reads, NA_real_)
    )
  }))
  
  result <- merge(base_metrics, faith_pd,           by = "sample_id")
  result <- merge(result,       additional_metrics, by = "sample_id")
  return(result)
}

# =========================================================================
# PLOT FUNCTION
# =========================================================================

create_type_box_plot <- function(data, metric, cld_data, type_p, time_p, interaction_p,
                                 letter_offset = 0.15) {
  
  fmt_p <- function(p) {
    if (is.na(p))  return("P-value = NA")
    if (p < 0.001) return("P-value < 0.001")
    stringr::str_glue("P-value = {sprintf('%.3f', p)}")
  }
  
  p_annotation <- stringr::str_glue(
    "Compartment: {fmt_p(type_p)}\n",
    "Time: {fmt_p(time_p)}\n",
    "Interaction: {fmt_p(interaction_p)}"
  )
  
  show_cld <- !is.na(type_p) && type_p < alpha_level
  
  y_max_group <- data %>%
    group_by(type) %>%
    summarise(
      max_val = max(.data[[metric]], na.rm = TRUE),
      min_val = min(.data[[metric]], na.rm = TRUE),
      n       = dplyr::n(),
      .groups = "drop"
    )
  
  global_max   <- max(data[[metric]], na.rm = TRUE)
  global_min   <- min(data[[metric]], na.rm = TRUE)
  global_range <- global_max - global_min
  offset_val   <- global_max * letter_offset
  
  p <- ggplot(data, aes(x = type, y = .data[[metric]], fill = type)) +
    geom_half_violin(side = "r", trim = FALSE, alpha = 0.5, color = NA) +
    geom_boxplot(aes(color = type, fill = after_scale(alpha(color, 0.2))),
                 width = 0.4, outlier.shape = NA, lwd = 3) +
    stat_summary(aes(color = type), fun = mean, geom = "point",
                 shape = 23, size = 13, fill = "white", stroke = 4) +
    geom_jitter(aes(shape = type, color = type), position = position_jitter(width = 0.3),
                size = 13, alpha = 0.4, stroke = 4) +
    annotate("text", x = Inf, y = Inf, label = p_annotation,
             hjust = 1, vjust = 1.1, size = 4, fontface = "italic", lineheight = 1) +
    scale_fill_manual(values = color_palette, name = "Compartment") +
    scale_color_manual(values = color_palette, name = "Compartment") +
    scale_shape_manual(values = type_shapes, name = "Compartment") +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
    labs(title = display_names[metric], y = display_names[metric], x = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", hjust = 0.5),
      axis.ticks.x       = element_line(color = "black"),
      axis.text.x        = element_text(color = "black", face = "bold"),
      axis.text.y        = element_text(color = "black"),
      axis.title         = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "gray92"),
      legend.position    = "none"
    )
  
  if (show_cld) {
    cld_plot_data <- cld_data %>%
      mutate(Letter = stringr::str_trim(.group)) %>%
      dplyr::select(type, Letter) %>%
      left_join(y_max_group, by = "type") %>%
      mutate(y_pos = max_val + offset_val)
    
    p <- p + geom_text(data = cld_plot_data,
                       aes(x = type, y = y_pos, label = Letter),
                       size = 4, fontface = "bold", color = "black")
  }
  
  n_label_data <- y_max_group %>%
    mutate(
      label = stringr::str_glue("n = {n}"),
      y_pos = min_val - global_range * (letter_offset * 0.6)
    )
  
  p <- p + geom_text(
    data        = n_label_data,
    aes(x = type, y = y_pos, label = label, color = type),
    size        = 3,
    fontface    = "plain",
    inherit.aes = FALSE,
    show.legend = FALSE
  )
  
  return(p)
}

# =========================================================================
# STEP 1: LOAD AND PREPARE DATA
# =========================================================================

physeq <- readRDS(phyloseq_path)

sample_metadata <- phyloseq::sample_data(physeq) %>%
  as_tibble(rownames = "sample") %>%
  janitor::clean_names() %>%
  dplyr::filter(!is.na(trial) & trial == "1") %>%
  process_metadata() %>%
  filter_samples(ponds = filter_ponds, timepoints = filter_timepoints, shrimp_max = filter_shrimp_max)

physeq_filtered <- phyloseq::prune_samples(sample_metadata$sample, physeq)

sample_metadata <- phyloseq::sample_data(physeq_filtered) %>%
  as_tibble(rownames = "sample") %>%
  janitor::clean_names() %>%
  process_metadata()

# =========================================================================
# STEP 2: CALCULATE ALPHA DIVERSITY
# =========================================================================

otu_mat <- phyloseq::otu_table(physeq_filtered) %>% as("matrix")
if (phyloseq::taxa_are_rows(physeq_filtered)) otu_mat <- t(otu_mat)

otu_mat_filtered <- otu_mat[rowSums(otu_mat) >= rarefaction_depth, ]
tree             <- phyloseq::phy_tree(physeq_filtered)

test_result <- tryCatch({
  calculate_alpha_diversity_one_iter(1, otu_mat_filtered, rarefaction_depth, tree)
}, error = function(e) {
  stop(stringr::str_glue("Test iteration failed: {e$message}\nParallel processing cannot proceed."))
})

n_cores <- parallel::detectCores() - 1

if (parallel_method == "mclapply") {
  alpha_div_iterations <- parallel::mclapply(1:n_iterations, function(i) {
    tryCatch({
      calculate_alpha_diversity_one_iter(i, otu_mat_filtered, rarefaction_depth, tree)
    }, error = function(e) list(error = TRUE, message = as.character(e), iteration = i))
  }, mc.cores = n_cores, mc.set.seed = TRUE)
  
} else if (parallel_method == "parLapply") {
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(cl, c("calculate_alpha_diversity_one_iter", "otu_mat_filtered",
                                "rarefaction_depth", "tree"), envir = environment())
  parallel::clusterEvalQ(cl, { library(vegan); library(picante); library(dplyr) })
  alpha_div_iterations <- parallel::parLapply(cl, 1:n_iterations, function(i) {
    tryCatch({
      calculate_alpha_diversity_one_iter(i, otu_mat_filtered, rarefaction_depth, tree)
    }, error = function(e) list(error = TRUE, message = as.character(e), iteration = i))
  })
  parallel::stopCluster(cl)
  
} else {
  future::plan(future::multisession, workers = n_cores)
  alpha_div_iterations <- future.apply::future_lapply(1:n_iterations, function(i) {
    tryCatch({
      calculate_alpha_diversity_one_iter(i, otu_mat_filtered, rarefaction_depth, tree)
    }, error = function(e) list(error = TRUE, message = as.character(e), iteration = i))
  }, future.seed = TRUE)
  future::plan(future::sequential)
}

null_results  <- sapply(alpha_div_iterations, is.null)
error_results <- sapply(alpha_div_iterations, function(x) {
  if (is.list(x) && !is.null(x$error) && x$error) TRUE else FALSE
})

valid_results <- !null_results & !error_results
if (sum(valid_results) == 0) {
  alpha_div_iterations <- lapply(1:n_iterations, function(i) {
    calculate_alpha_diversity_one_iter(i, otu_mat_filtered, rarefaction_depth, tree)
  })
} else {
  alpha_div_iterations <- alpha_div_iterations[valid_results]
}

all_iterations <- do.call(rbind, alpha_div_iterations)

alpha_diversity <- all_iterations %>%
  group_by(sample_id) %>%
  summarise(across(all_of(primary_metrics), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  left_join(sample_metadata, by = c("sample_id" = "sample")) %>%
  dplyr::filter(complete.cases(type, timepoint, pond)) %>%
  mutate(across(c(type, timepoint, pond, time_factor), as.factor))

# =========================================================================
# STEP 2B: AVERAGE INTESTINE REPLICATES TO POND LEVEL
# =========================================================================
# Intestine has up to 5 shrimp replicates per pond x timepoint.
# Pond is the true experimental unit; averaging collapses pseudoreplicates
# to give 1 observation per pond x type x timepoint before modeling.
# Water and Sediment are unaffected (already 1 sample per pond x timepoint).
# =========================================================================

alpha_diversity <- alpha_diversity |>
  dplyr::group_by(pond, type, timepoint, time_numeric, time_factor) |>
  dplyr::summarise(
    across(all_of(primary_metrics), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

check <- alpha_diversity |>
  dplyr::count(pond, type, timepoint) |>
  dplyr::filter(n > 1)

if (nrow(check) > 0) {
  warning("Duplicate rows remain after averaging intestine replicates")
  print(check)
}

# =========================================================================
# STEP 3: FIT MODELS AND EXTRACT RESULTS
# =========================================================================

# optim optimizer more robust than default nlminb for this dataset
lme_control <- nlme::lmeControl(msMaxIter = 200, opt = "optim")

models_list            <- list()
anova_list             <- list()
type_only_models_list  <- list()
type_only_anova_list   <- list()
type_only_emmeans_list <- list()

for (metric in primary_metrics) {
  
  # Full model with interaction — used for Type and Time p-values on plots
  model <- nlme::lme(
    as.formula(stringr::str_glue("{metric} ~ type * time_factor")),
    random      = ~1 | pond,
    correlation = nlme::corAR1(form = ~time_numeric | pond / type),
    data        = alpha_diversity,
    na.action   = na.omit,
    method      = "REML",
    control     = lme_control
  )
  
  models_list[[metric]] <- model
  anova_list[[metric]]  <- anova(model, type = "marginal")
  
  # Compartment-only model — used for CLD letters on plots
  type_only_model <- nlme::lme(
    as.formula(stringr::str_glue("{metric} ~ type")),
    random      = ~1 | pond,
    correlation = nlme::corAR1(form = ~time_numeric | pond / type),
    data        = alpha_diversity,
    na.action   = na.omit,
    method      = "REML",
    control     = lme_control
  )
  
  type_only_models_list[[metric]] <- type_only_model
  type_only_anova_list[[metric]]  <- anova(type_only_model, type = "marginal")
  
  type_only_emm <- emmeans::emmeans(type_only_model, ~ type)
  type_only_cld <- multcomp::cld(type_only_emm, Letters = letters, adjust = "tukey")
  type_only_emmeans_list[[metric]] <- as_tibble(type_only_cld) |>
    dplyr::mutate(.group = stringr::str_trim(.group))
}

# =========================================================================
# STEP 4: CREATE PLOTS
# =========================================================================

type_box_plots <- list()

for (metric in primary_metrics) {
  type_box_plots[[metric]] <- create_type_box_plot(
    alpha_diversity, metric,
    type_only_emmeans_list[[metric]],
    type_only_anova_list[[metric]]["type", "p-value"],
    anova_list[[metric]]["time_factor", "p-value"],
    anova_list[[metric]]["type:time_factor", "p-value"]
  )
}

# =========================================================================
# STEP 5: COMBINED FIGURE
# =========================================================================

combined_alpha_plot <- (
  type_box_plots[["observed_richness"]] +
    type_box_plots[["shannon"]] +
    type_box_plots[["berger_parker"]] +
    type_box_plots[["faith_pd"]]
) + plot_layout(ncol = 4)

ggsave(
  filename = file.path(output_dir, "AlphaDiversity.png"),
  plot     = combined_alpha_plot,
  width    = 12,
  height   = 8,
  units    = "in",
  dpi      = 300
)
