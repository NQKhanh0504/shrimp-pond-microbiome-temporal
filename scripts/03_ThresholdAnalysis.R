# ============================================================
# Script: 03_ThresholdAnalysis.R
# Description: Abundance and prevalence filtering threshold
#              analysis for Litopenaeus vannamei production
#              pond microbiome data.
#              Tests combined abundance + prevalence filtering
#              and monitors critical taxa impacts.
# Input: phyloseq object (seqtab_nochim.rds)
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

required_packages <- c("phyloseq", "stringr", "dplyr", "tidyr", "readr")

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
key_packages <- c("phyloseq", "stringr", "dplyr", "readr")
for (pkg in key_packages) cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

# File paths
phyloseq_file <- "seqtab_nochim.rds"
output_dir    <- "output/"
csv_dir       <- file.path(output_dir, "csv_exports")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

# Filtering parameters
abundance_thresholds  <- seq(0, 100, by = 5)
prevalence_thresholds <- seq(0, 20,  by = 1)   # % of samples within group

# Analysis scope
analyze_global_prevalence <- FALSE  # Test prevalence across all samples
analyze_group_prevalence  <- TRUE   # Test prevalence within sample types

# Taxonomic analysis settings
analyze_levels  <- c("Phylum", "Class", "Order", "Family", "Genus")
primary_level   <- "Genus"
secondary_level <- "Family"

# Phylogenetic annotation handling
aggregate_phylogenetic_variants <- TRUE
phylo_suffix_pattern <- "(_[A-Z0-9]+|_UCG-[0-9]+|_NK[0-9A-Z]+|_unclassified|_incertae_sedis)$"

# Critical taxa for monitoring
critical_taxa <- list(
  Genus = list(
    list(name = "Vibrio",        patterns = c("Vibrio"),        exact_match = FALSE),
    list(name = "Aeromonas",     patterns = c("Aeromonas"),     exact_match = FALSE),
    list(name = "Pseudomonas",   patterns = c("Pseudomonas"),   exact_match = FALSE),
    list(name = "Lactobacillus", patterns = c("Lactobacillus"), exact_match = FALSE),
    list(name = "Bacillus",      patterns = c("Bacillus"),      exact_match = FALSE),
    list(name = "Enterococcus",  patterns = c("Enterococcus"),  exact_match = FALSE),
    list(name = "Shewanella",    patterns = c("Shewanella"),    exact_match = FALSE),
    list(name = "Flavobacterium",patterns = c("Flavobacterium"),exact_match = FALSE)
  ),
  Family = list(
    list(name = "Vibrionaceae",      patterns = c("Vibrionaceae"),      exact_match = FALSE),
    list(name = "Lactobacillaceae",  patterns = c("Lactobacillaceae"),  exact_match = FALSE),
    list(name = "Bacillaceae",       patterns = c("Bacillaceae"),       exact_match = FALSE),
    list(name = "Rhodobacteraceae",  patterns = c("Rhodobacteraceae"),  exact_match = FALSE),
    list(name = "Flavobacteriaceae", patterns = c("Flavobacteriaceae"), exact_match = FALSE)
  ),
  Order = list(
    list(name = "Vibrionales",    patterns = c("Vibrionales"),    exact_match = FALSE),
    list(name = "Lactobacillales",patterns = c("Lactobacillales"),exact_match = FALSE),
    list(name = "Bacillales",     patterns = c("Bacillales"),     exact_match = FALSE),
    list(name = "Rhodobacterales",patterns = c("Rhodobacterales"),exact_match = FALSE),
    list(name = "Pseudomonadales",patterns = c("Pseudomonadales"),exact_match = FALSE)
  )
)

# Critical taxa filtering thresholds
critical_abundance_threshold  <- 10
critical_prevalence_threshold <- 0.02
use_critical_filters          <- TRUE

# Taxonomic cleaning settings
remove_prefixes      <- TRUE
handle_unclassified  <- TRUE

set.seed(54)

# =========================================================================
# UTILITY FUNCTIONS
# =========================================================================

get_sample_type_column <- function(metadata) {
  potential_cols <- c("type", "sample_type", "sampletype", "group", "category")
  for (col in potential_cols) {
    if (col %in% colnames(metadata)) return(col)
  }
  for (col in colnames(metadata)) {
    if (is.character(metadata[[col]]) || is.factor(metadata[[col]])) {
      n_levels <- length(unique(metadata[[col]]))
      if (n_levels > 1 && n_levels < nrow(metadata) / 2) return(col)
    }
  }
  return(NULL)
}

clean_taxonomic_names <- function(taxa_vector, level_name, preserve_length = FALSE) {
  if (length(taxa_vector) == 0) return(character(0))
  
  taxa_clean  <- taxa_vector
  na_positions <- is.na(taxa_clean)
  
  if (remove_prefixes) {
    level_prefixes <- list(
      Kingdom = "^[dk]__", Phylum = "^p__", Class  = "^c__",
      Order   = "^o__",    Family = "^f__", Genus  = "^g__", Species = "^s__"
    )
    level_key <- names(level_prefixes)[
      tolower(names(level_prefixes)) == tolower(level_name)
    ]
    if (length(level_key) > 0) {
      taxa_clean <- stringr::str_remove(taxa_clean, level_prefixes[[level_key[1]]])
    }
  }
  
  unclassified_pos <- taxa_clean == "" | is.na(taxa_clean) | na_positions
  if (handle_unclassified || preserve_length) {
    taxa_clean[unclassified_pos] <- paste0("Unclassified_", level_name)
  } else if (!preserve_length) {
    taxa_clean <- taxa_clean[!unclassified_pos]
  }
  
  return(taxa_clean)
}

aggregate_phylogenetic_names <- function(taxa_vector) {
  if (!aggregate_phylogenetic_variants || length(taxa_vector) == 0) return(taxa_vector)
  stringr::str_remove(taxa_vector, phylo_suffix_pattern)
}

match_critical_taxon <- function(taxa_names, pattern, exact_match = FALSE) {
  if (exact_match) {
    which(stringr::str_detect(taxa_names, paste0("^", pattern, "$")))
  } else if (aggregate_phylogenetic_variants) {
    taxa_agg <- aggregate_phylogenetic_names(taxa_names)
    which(stringr::str_detect(taxa_agg, paste0("^", pattern, "$")))
  } else {
    which(stringr::str_detect(taxa_names, paste0("^", pattern, "($|_)")))
  }
}

check_taxonomy_columns <- function(taxa_matrix, levels_to_check = analyze_levels) {
  available <- character()
  for (level in levels_to_check) {
    if      (level            %in% colnames(taxa_matrix)) available <- c(available, level)
    else if (tolower(level)   %in% colnames(taxa_matrix)) available <- c(available, tolower(level))
  }
  if (length(available) == 0) stop("No taxonomic levels found. Expected: ",
                                   paste(levels_to_check, collapse = ", "))
  return(available)
}

apply_abundance_filter <- function(seqtab, threshold) {
  asv_counts   <- colSums(seqtab)
  keep_idx     <- which(asv_counts >= threshold)
  list(
    filtered_seqtab = seqtab[, keep_idx, drop = FALSE],
    removed_indices = which(asv_counts < threshold),
    kept_indices    = keep_idx,
    n_removed       = sum(asv_counts < threshold),
    n_kept          = length(keep_idx)
  )
}

apply_prevalence_filter <- function(seqtab, threshold_percent) {
  n_samples       <- nrow(seqtab)
  threshold_count <- ceiling(n_samples * threshold_percent / 100)
  asv_prev        <- colSums(seqtab > 0)
  keep_idx        <- which(asv_prev >= threshold_count)
  list(
    filtered_seqtab = seqtab[, keep_idx, drop = FALSE],
    removed_indices = which(asv_prev < threshold_count),
    kept_indices    = keep_idx,
    n_removed       = sum(asv_prev < threshold_count),
    n_kept          = length(keep_idx),
    threshold_count = threshold_count
  )
}

apply_combined_filter <- function(seqtab, abundance_thresh, prevalence_thresh_percent) {
  abund_result <- apply_abundance_filter(seqtab, abundance_thresh)
  
  if (prevalence_thresh_percent == 0) {
    return(list(
      filtered_seqtab      = abund_result$filtered_seqtab,
      removed_indices      = abund_result$removed_indices,
      kept_indices         = abund_result$kept_indices,
      n_removed            = abund_result$n_removed,
      n_kept               = abund_result$n_kept,
      removed_by_abundance = abund_result$n_removed,
      removed_by_prevalence = 0
    ))
  }
  
  prev_result   <- apply_prevalence_filter(abund_result$filtered_seqtab, prevalence_thresh_percent)
  original_kept <- abund_result$kept_indices
  final_kept    <- original_kept[prev_result$kept_indices]
  final_removed <- setdiff(seq_len(ncol(seqtab)), final_kept)
  
  list(
    filtered_seqtab       = prev_result$filtered_seqtab,
    removed_indices       = final_removed,
    kept_indices          = final_kept,
    n_removed             = length(final_removed),
    n_kept                = length(final_kept),
    removed_by_abundance  = abund_result$n_removed,
    removed_by_prevalence = prev_result$n_removed
  )
}

analyze_taxonomic_impact <- function(taxa_matrix, removed_indices, taxonomic_level,
                                     critical_taxa_list = NULL, sample_metadata = NULL,
                                     seqtab = NULL, abundance_threshold = NULL,
                                     prevalence_threshold = NULL) {
  empty_result <- list(
    level = taxonomic_level, total_taxa_removed = 0,
    taxa_counts = table(character(0)), critical_impacts = numeric(0),
    sample_type_impacts = data.frame()
  )
  
  if (length(removed_indices) == 0) return(empty_result)
  
  level_col <- if (taxonomic_level %in% colnames(taxa_matrix)) taxonomic_level
  else tolower(taxonomic_level)
  if (!level_col %in% colnames(taxa_matrix)) return(empty_result)
  
  taxa_removed <- taxa_matrix[removed_indices, level_col]
  taxa_clean   <- clean_taxonomic_names(taxa_removed, taxonomic_level, preserve_length = FALSE)
  if (aggregate_phylogenetic_variants && length(taxa_clean) > 0) {
    taxa_clean <- aggregate_phylogenetic_names(taxa_clean)
  }
  if (length(taxa_clean) == 0) return(empty_result)
  
  taxa_counts      <- table(taxa_clean)
  critical_impacts <- numeric(0)
  sample_type_impacts <- data.frame()
  
  if (!is.null(critical_taxa_list) && taxonomic_level %in% names(critical_taxa_list)) {
    critical_list  <- critical_taxa_list[[taxonomic_level]]
    all_taxa_names <- clean_taxonomic_names(taxa_matrix[, level_col], taxonomic_level,
                                            preserve_length = FALSE)
    
    for (critical_item in critical_list) {
      taxon_name <- critical_item$name
      all_taxon_idx <- unique(unlist(lapply(critical_item$patterns, function(pat) {
        match_critical_taxon(all_taxa_names, pat, critical_item$exact_match)
      })))
      
      if (length(all_taxon_idx) == 0) next
      
      if (use_critical_filters && !is.null(seqtab)) {
        total_abund  <- sum(colSums(seqtab[, all_taxon_idx, drop = FALSE]))
        prevalence   <- sum(rowSums(seqtab[, all_taxon_idx, drop = FALSE] > 0) > 0) / nrow(seqtab)
        abund_pass   <- total_abund >= critical_abundance_threshold
        prev_pass    <- prevalence  >= critical_prevalence_threshold
        if (!abund_pass || !prev_pass) next
      }
      
      removed_count <- length(intersect(all_taxon_idx, removed_indices))
      if (removed_count > 0) critical_impacts[taxon_name] <- removed_count
    }
    
    # Sample-type specific impacts
    if (length(critical_impacts) > 0 && !is.null(sample_metadata) && !is.null(seqtab)) {
      type_column <- get_sample_type_column(sample_metadata)
      if (!is.null(type_column)) {
        for (taxon_name in names(critical_impacts)) {
          critical_item <- NULL
          for (item in critical_list) {
            if (item$name == taxon_name) { critical_item <- item; break }
          }
          if (is.null(critical_item)) next
          
          all_taxon_idx <- unique(unlist(lapply(critical_item$patterns, function(pat) {
            match_critical_taxon(all_taxa_names, pat, critical_item$exact_match)
          })))
          
          for (type in unique(sample_metadata[[type_column]])) {
            removed_in_type <- intersect(all_taxon_idx, removed_indices)
            if (length(removed_in_type) > 0) {
              sample_type_impacts <- rbind(sample_type_impacts, data.frame(
                Taxon       = taxon_name,
                Level       = taxonomic_level,
                Sample_Type = type,
                ASVs_Lost   = length(removed_in_type),
                stringsAsFactors = FALSE
              ))
            }
          }
        }
      }
    }
  }
  
  list(
    level               = taxonomic_level,
    total_taxa_removed  = length(taxa_counts),
    taxa_counts         = taxa_counts,
    critical_impacts    = critical_impacts,
    sample_type_impacts = sample_type_impacts
  )
}

# =========================================================================
# STEP 1: LOAD AND PREPARE DATA
# =========================================================================

ps_object <- readRDS(phyloseq_file)

if (!inherits(ps_object, "phyloseq")) {
  stop("Loaded object is not a phyloseq object. Found: ", class(ps_object)[1])
}

seqtab.nochim <- as.matrix(phyloseq::otu_table(ps_object))
if (phyloseq::taxa_are_rows(ps_object)) seqtab.nochim <- t(seqtab.nochim)

if (!is.null(phyloseq::sample_data(ps_object, errorIfNULL = FALSE))) {
  sample_metadata <- data.frame(phyloseq::sample_data(ps_object))
} else {
  sample_metadata <- data.frame(
    Sample_ID = rownames(seqtab.nochim),
    group     = "All_Samples",
    stringsAsFactors = FALSE
  )
  rownames(sample_metadata) <- rownames(seqtab.nochim)
}

if (!is.null(phyloseq::tax_table(ps_object, errorIfNULL = FALSE))) {
  taxa_matrix <- as.matrix(phyloseq::tax_table(ps_object))
} else {
  stop("No taxonomy found in phyloseq object.")
}

available_levels <- check_taxonomy_columns(taxa_matrix, analyze_levels)

if (!primary_level %in% available_levels && !tolower(primary_level) %in% available_levels) {
  primary_level <- available_levels[1]
}
if (!secondary_level %in% available_levels && !tolower(secondary_level) %in% available_levels) {
  secondary_level <- if (length(available_levels) > 1) available_levels[2] else primary_level
}

type_column <- get_sample_type_column(sample_metadata)
if (is.null(type_column)) {
  sample_metadata$group <- "All_Samples"
  type_column           <- "group"
  analyze_group_prevalence <- FALSE
}

sample_groups <- unique(sample_metadata[[type_column]])
total_asvs    <- ncol(seqtab.nochim)

# =========================================================================
# STEP 2: BASELINE ABUNDANCE AND PREVALENCE ANALYSIS
# =========================================================================

asv_counts     <- colSums(seqtab.nochim)
asv_prevalence <- colSums(seqtab.nochim > 0)
n_samples      <- nrow(seqtab.nochim)

# Sample-type-specific analysis
type_asv_lists <- list()
type_summaries <- list()

for (group in sample_groups) {
  group_samples  <- sample_metadata[[type_column]] == group
  group_seqtab   <- seqtab.nochim[group_samples, , drop = FALSE]
  group_asv_cnt  <- colSums(group_seqtab)
  group_seqtab   <- group_seqtab[, group_asv_cnt > 0, drop = FALSE]
  
  n_group        <- nrow(group_seqtab)
  group_asvs     <- colnames(group_seqtab)
  group_prev     <- colSums(group_seqtab > 0)
  group_abund    <- colSums(group_seqtab)
  
  type_asv_lists[[group]] <- group_asvs
  type_summaries[[group]] <- list(
    n_samples        = n_group,
    n_asvs           = length(group_asvs),
    singleton_pct    = round(sum(group_prev == 1) / length(group_asvs) * 100, 1),
    core_pct         = round(sum(group_prev > n_group * 0.5) / length(group_asvs) * 100, 1),
    median_prevalence = median(group_prev),
    median_abundance  = median(group_abund)
  )
}

# Export baseline comparison
comparison_df <- data.frame(
  Sample_Type       = names(type_summaries),
  N_Samples         = sapply(type_summaries, function(x) x$n_samples),
  N_ASVs            = sapply(type_summaries, function(x) x$n_asvs),
  Singleton_Pct     = sapply(type_summaries, function(x) x$singleton_pct),
  Core_Pct          = sapply(type_summaries, function(x) x$core_pct),
  Median_Prevalence = sapply(type_summaries, function(x) x$median_prevalence),
  Median_Abundance  = sapply(type_summaries, function(x) x$median_abundance)
)

readr::write_csv(comparison_df, file.path(csv_dir, "sample_type_comparison.csv"))

# ASV sharing analysis (3 groups)
if (length(sample_groups) == 3) {
  g1 <- sample_groups[1]; g2 <- sample_groups[2]; g3 <- sample_groups[3]
  
  sharing_df <- data.frame(
    Category       = c(
      paste0(g1, "_only"), paste0(g2, "_only"), paste0(g3, "_only"),
      paste0(g1, "_", g2, "_only"), paste0(g1, "_", g3, "_only"),
      paste0(g2, "_", g3, "_only"), "shared_all"
    ),
    N_ASVs         = c(
      length(setdiff(type_asv_lists[[g1]], union(type_asv_lists[[g2]], type_asv_lists[[g3]]))),
      length(setdiff(type_asv_lists[[g2]], union(type_asv_lists[[g1]], type_asv_lists[[g3]]))),
      length(setdiff(type_asv_lists[[g3]], union(type_asv_lists[[g1]], type_asv_lists[[g2]]))),
      length(setdiff(intersect(type_asv_lists[[g1]], type_asv_lists[[g2]]), type_asv_lists[[g3]])),
      length(setdiff(intersect(type_asv_lists[[g1]], type_asv_lists[[g3]]), type_asv_lists[[g2]])),
      length(setdiff(intersect(type_asv_lists[[g2]], type_asv_lists[[g3]]), type_asv_lists[[g1]])),
      length(Reduce(intersect, list(type_asv_lists[[g1]], type_asv_lists[[g2]], type_asv_lists[[g3]])))
    )
  )
  sharing_df$Pct_of_Total <- round(sharing_df$N_ASVs / total_asvs * 100, 2)
  readr::write_csv(sharing_df, file.path(csv_dir, "asv_sharing_analysis.csv"))
}

# =========================================================================
# STEP 3: GLOBAL PREVALENCE ANALYSIS
# =========================================================================

if (analyze_global_prevalence) {
  global_critical_summary <- data.frame()
  filtering_matrix        <- data.frame()
  
  for (abundance_thresh in abundance_thresholds) {
    for (prevalence_thresh in prevalence_thresholds) {
      filter_result <- apply_combined_filter(seqtab.nochim, abundance_thresh, prevalence_thresh)
      
      filtering_matrix <- rbind(filtering_matrix, data.frame(
        Abundance_Threshold   = abundance_thresh,
        Prevalence_Threshold  = prevalence_thresh,
        ASVs_Remaining        = filter_result$n_kept,
        ASVs_Removed          = filter_result$n_removed,
        Removed_by_Abundance  = filter_result$removed_by_abundance,
        Removed_by_Prevalence = filter_result$removed_by_prevalence,
        Percent_Remaining     = round(filter_result$n_kept / total_asvs * 100, 2)
      ))
      
      for (level in available_levels) {
        taxa_impact <- analyze_taxonomic_impact(
          taxa_matrix       = taxa_matrix,
          removed_indices   = filter_result$removed_indices,
          taxonomic_level   = level,
          critical_taxa_list = critical_taxa,
          sample_metadata   = sample_metadata,
          seqtab            = seqtab.nochim,
          abundance_threshold  = abundance_thresh,
          prevalence_threshold = prevalence_thresh
        )
        
        if (length(taxa_impact$critical_impacts) > 0) {
          for (taxon in names(taxa_impact$critical_impacts)) {
            global_critical_summary <- rbind(global_critical_summary, data.frame(
              Scope                = "Global",
              Sample_Group         = "All",
              Abundance_Threshold  = abundance_thresh,
              Prevalence_Threshold = prevalence_thresh,
              Taxonomic_Level      = level,
              Taxon                = taxon,
              ASVs_Lost            = taxa_impact$critical_impacts[taxon]
            ))
          }
        }
      }
    }
  }
  
  readr::write_csv(filtering_matrix,
                   file.path(csv_dir, "filtering_impact_matrix_global.csv"))
  
  if (nrow(global_critical_summary) > 0) {
    readr::write_csv(global_critical_summary,
                     file.path(csv_dir, "critical_taxa_global.csv"))
  }
  
  # Safe threshold identification
  safe_thresholds <- data.frame()
  for (i in seq_len(nrow(filtering_matrix))) {
    row        <- filtering_matrix[i, ]
    abund_t    <- row$Abundance_Threshold
    prev_t     <- row$Prevalence_Threshold
    
    has_critical <- if (nrow(global_critical_summary) > 0) {
      nrow(global_critical_summary[
        global_critical_summary$Abundance_Threshold  == abund_t &
          global_critical_summary$Prevalence_Threshold == prev_t, ]) > 0
    } else FALSE
    
    if (!has_critical) {
      safe_thresholds <- rbind(safe_thresholds, data.frame(
        Abundance         = abund_t,
        Prevalence        = prev_t,
        ASVs_Remaining    = row$ASVs_Remaining,
        Percent_Remaining = row$Percent_Remaining
      ))
    }
  }
  
  if (nrow(safe_thresholds) > 0) {
    safe_thresholds <- safe_thresholds[order(safe_thresholds$ASVs_Remaining), ]
    readr::write_csv(safe_thresholds,
                     file.path(csv_dir, "threshold_combinations.csv"))
  }
}

# =========================================================================
# STEP 4: GROUP-SPECIFIC PREVALENCE ANALYSIS
# =========================================================================

if (analyze_group_prevalence) {
  group_critical_summary <- data.frame()
  
  for (group in sample_groups) {
    group_samples    <- rownames(sample_metadata)[sample_metadata[[type_column]] == group]
    group_seqtab     <- seqtab.nochim[group_samples, , drop = FALSE]
    group_asv_counts <- colSums(group_seqtab)
    group_seqtab     <- group_seqtab[, group_asv_counts > 0, drop = FALSE]
    group_n_asvs     <- ncol(group_seqtab)
    
    group_matrix <- data.frame()
    
    for (abundance_thresh in abundance_thresholds) {
      for (prevalence_thresh in prevalence_thresholds) {
        filter_result <- apply_combined_filter(group_seqtab, abundance_thresh, prevalence_thresh)
        
        removed_asv_names    <- colnames(group_seqtab)[filter_result$removed_indices]
        removed_global_idx   <- which(colnames(seqtab.nochim) %in% removed_asv_names)
        
        group_matrix <- rbind(group_matrix, data.frame(
          Sample_Group          = group,
          Abundance_Threshold   = abundance_thresh,
          Prevalence_Threshold  = prevalence_thresh,
          ASVs_Remaining        = filter_result$n_kept,
          ASVs_Removed          = filter_result$n_removed,
          Removed_by_Abundance  = filter_result$removed_by_abundance,
          Removed_by_Prevalence = filter_result$removed_by_prevalence,
          Percent_Remaining     = round(filter_result$n_kept / group_n_asvs * 100, 2)
        ))
        
        for (level in available_levels) {
          taxa_impact <- analyze_taxonomic_impact(
            taxa_matrix        = taxa_matrix,
            removed_indices    = removed_global_idx,
            taxonomic_level    = level,
            critical_taxa_list = critical_taxa,
            sample_metadata    = sample_metadata[group_samples, , drop = FALSE],
            seqtab             = seqtab.nochim,
            abundance_threshold  = abundance_thresh,
            prevalence_threshold = prevalence_thresh
          )
          
          if (length(taxa_impact$critical_impacts) > 0) {
            for (taxon in names(taxa_impact$critical_impacts)) {
              group_critical_summary <- rbind(group_critical_summary, data.frame(
                Scope                = "Group",
                Sample_Group         = group,
                Abundance_Threshold  = abundance_thresh,
                Prevalence_Threshold = prevalence_thresh,
                Taxonomic_Level      = level,
                Taxon                = taxon,
                ASVs_Lost            = taxa_impact$critical_impacts[taxon]
              ))
            }
          }
        }
      }
    }
    
    clean_group <- stringr::str_replace_all(group, "[^A-Za-z0-9]", "_")
    readr::write_csv(group_matrix,
                     file.path(csv_dir, paste0("filtering_impact_matrix_", clean_group, ".csv")))
  }
  
  if (nrow(group_critical_summary) > 0) {
    readr::write_csv(group_critical_summary,
                     file.path(csv_dir, "critical_taxa_by_group.csv"))
  }
}

# =========================================================================
# STEP 5: TAXONOMIC LOSS SUMMARIES BY LEVEL
# =========================================================================

if (analyze_global_prevalence && exists("global_critical_summary") &&
    nrow(global_critical_summary) > 0) {
  for (level in available_levels) {
    level_data <- global_critical_summary[global_critical_summary$Taxonomic_Level == level, ]
    if (nrow(level_data) > 0) {
      readr::write_csv(
        level_data,
        file.path(csv_dir, paste0("taxonomic_loss_", tolower(level), ".csv"))
      )
    }
  }
}
