# ============================================================
# Script: 00_DADA2.R
# Description: DADA2 16S rRNA V3-V4 denoising pipeline for
#              Litopenaeus vannamei production pond microbiome
#              data. Processes paired-end FASTQ files
#              through ASV inference and chimera removal to
#              generate sequence table and phyloseq object.
# Input: Paired-end FASTQ files (demultiplexed)
# Output: seqtab_nochim.rds (phyloseq object)
# Note: Trimming parameters (trim_F, trim_R) derived from
#       Figaro analysis of V3-V4 amplicon data.
# ============================================================

# =========================================================================
# PACKAGE INSTALLATION AND LOADING
# =========================================================================

required_packages <- c("dada2", "ShortRead", "Biostrings", "phyloseq", "tidyverse")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg, ask = FALSE)
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
}

# =========================================================================
# SESSION INFO
# =========================================================================

cat("R Version:", R.version.string, "\n")
cat("Key Package Versions:\n")
for (pkg in c("dada2", "ShortRead", "phyloseq")) {
  cat(sprintf("  %s: %s\n", pkg, packageVersion(pkg)))
}
cat("\n")

# =========================================================================
# CONFIGURATION
# =========================================================================

# File paths — update for your environment
path          <- "fastq/"
save_dir      <- "output/"
metadata_file <- file.path(save_dir, "metadata.txt")

dir.create(save_dir,                    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(path, "filtered"), recursive = TRUE, showWarnings = FALSE)

# Primer parameters (V3-V4 region)
# Forward: CCTACGGGNGGCWGCAG (16 bp)
# Reverse: GACTACHVGGGTATCTAATCC (24 bp)
primer_len_F <- 16
primer_len_R <- 24

# Trimming parameters — derived from Figaro analysis
# Figaro input: expected amplicon length 460 bp, V3-V4 primers
trim_F <- 304
trim_R <- 176

# Quality filtering parameters
max_expected_errors_F <- 2   # Forward reads (more lenient)
max_expected_errors_R <- 1   # Reverse reads (stricter due to quality degradation)
truncQ <- 2                  # Truncate at first quality score <= truncQ
maxN   <- 0                  # No ambiguous bases allowed
rm_phix <- TRUE              # Remove PhiX contamination

# Denoising parameters
pool_samples   <- "pseudo"      # Recommended for datasets >100 samples
chimera_method <- "consensus"

set.seed(54)

# =========================================================================
# STEP 1: FILE PREPARATION
# =========================================================================

fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz$", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)

if (length(fnFs) != length(fnRs)) {
  stop("Mismatch between forward and reverse file counts.")
}
if (length(fnFs) <= 2) {
  stop("Fewer than 2 samples detected. Ensure data is demultiplexed.")
}

# =========================================================================
# STEP 2: FILTER AND TRIM
# =========================================================================

filtFs <- file.path(path, "filtered", str_c(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", str_c(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen  = c(trim_F, trim_R),
  trimLeft  = c(primer_len_F, primer_len_R),
  maxN      = maxN,
  maxEE     = c(max_expected_errors_F, max_expected_errors_R),
  truncQ    = truncQ,
  rm.phix   = rm_phix,
  compress  = TRUE,
  multithread = TRUE
)

filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]

# =========================================================================
# STEP 3: LEARN ERROR RATES
# =========================================================================

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# =========================================================================
# STEP 4: SAMPLE INFERENCE (DENOISING)
# =========================================================================

dadaFs <- dada(filtFs, err = errF, multithread = TRUE, pool = pool_samples)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE, pool = pool_samples)

# =========================================================================
# STEP 5: MERGE PAIRED READS
# =========================================================================

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = FALSE)

# =========================================================================
# STEP 6: CONSTRUCT SEQUENCE TABLE
# =========================================================================

seqtab <- makeSequenceTable(mergers)

# =========================================================================
# STEP 7: REMOVE CHIMERAS
# =========================================================================

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method      = chimera_method,
  multithread = TRUE,
  verbose     = FALSE
)

# =========================================================================
# STEP 8: TRACK READS THROUGH PIPELINE
# =========================================================================

getN <- function(x) sum(getUniques(x))

track <- cbind(
  out,
  sapply(dadaFs,  getN),
  sapply(dadaRs,  getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.csv(track, file.path(save_dir, "read_tracking.csv"))

# =========================================================================
# STEP 9: CREATE PHYLOSEQ OBJECT WITH METADATA
# =========================================================================

otu <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)

if (file.exists(metadata_file)) {
  sam_data <- read.table(metadata_file, header = TRUE, row.names = 1, sep = "\t")
  sam      <- sample_data(sam_data)
  ps       <- phyloseq(otu, sam)
} else {
  ps <- phyloseq(otu)
}

# =========================================================================
# STEP 10: SAVE OUTPUT
# =========================================================================

saveRDS(ps, file.path(save_dir, "seqtab_nochim.rds"))

cat("Pipeline complete.\n")
cat("Samples:", nrow(seqtab.nochim), "\n")
cat("ASVs:", ncol(seqtab.nochim), "\n")
cat("Chimera fraction removed:", round((1 - sum(seqtab.nochim) / sum(seqtab)) * 100, 2), "%\n")
cat("Output saved to:", file.path(save_dir, "seqtab_nochim.rds"), "\n")
