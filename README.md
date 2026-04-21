# Temporal dynamics across the shrimp intestine, culture water, and sediment in production ponds.
16S rRNA amplicon sequencing analysis of microbial community temporal dynamics across intestine, sediment, and water compartments in Litopenaeus vannamei production ponds.

## Overview

This repository contains R scripts for reproducing the microbiome analyses reported in "Temporal dynamics across the shrimp intestine, culture water, and sediment in production ponds". The study characterizes temporal dynamics of microbial communities across three compartments, intestine, pond sediment, and pond water, in Pacific white shrimp (*Litopenaeus vannamei*) production ponds over an 8-week production cycle (Week 0, 2, 4, 6, 8) across four ponds.

Raw 16S rRNA amplicon sequencing data are deposited in the NCBI Sequence Read Archive under BioProject accession **PRJNA1444495**. The processed phyloseq object used as input for all analysis scripts is available at [Zenodo DOI — to be added upon acceptance].

---

## Repository Structure

```
.
├── README.md
├── metadata/
│   └── metadata.txt        # Sample metadata (see Metadata section)
└── scripts/
    ├── 00_DADA2.R                  # DADA2 denoising: raw FASTQ -> seqtab_nochim.rds
    ├── 01_16S_Pipeline.R           # 16S processing: seqtab_nochim.rds -> step4_phyloseq_object.rds
    ├── 02_Rarefaction.R            # Rarefaction curve analysis and depth selection
    ├── 03_ThresholdAnalysis.R      # Abundance and prevalence filtering threshold analysis
    ├── 04_AlphaDiversity.R         # Alpha diversity: richness, Shannon, Berger-Parker, Faith's PD
    ├── 05_BetaDiversity.R          # Beta diversity: PCoA, PERMANOVA, BETADISPER
    ├── 06_VennDiagram.R            # ASV sharing across compartments and timepoints
    ├── 07_Alluvial.R               # Relative abundance alluvial/stacked bar plots
    ├── 08_ANCOMBC.R                # ANCOM-BC2 differential abundance analysis
    └── 09_NetCoMi.R                # SPIEC-EASI co-occurrence network analysis
```

---

## Data Availability

Raw FASTQ reads at NCBI SRA: BioProject [PRJNA1444495](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1444495) |
Processed phyloseq object at Zenodo: [DOI to be added upon acceptance] |
Sample metadata | This repository: `metadata/Exp2433_Metadata.txt` |

---

## Metadata

The sample metadata file (`metadata/metadata.txt`) is a tab-delimited text file with the following columns:

| Column | Description |
|---|---|
| `sample` | Unique sample identifier |
| `fq1` | Forward FASTQ filename |
| `fq2` | Reverse FASTQ filename |
| `timepoint` | Sampling timepoint (T0–T4 corresponding to Week 0, 2, 4, 6, 8) |
| `type` | Sample compartment (Gut, Soil, Water) |
| `shrimp` | Shrimp individual within pond |
| `pond` | Pond identifier (C3, C6, D1, D4) |
| `code` | Sampling code |
| `om_percent` | Organic matter percentage |
| `nitrite` | Nitrite concentration |
| `ph` | pH |
| `ammonia` | Ammonia concentration |
| `feedinput` | Feed input |
| `temperature` | Water temperature |
| `do` | Dissolved oxygen |
| `trial` | Trial number |
| `shared` | Shared sample indicator |
| `trt` | Treatment group |

---

## Requirements

### R Environment

| Software | Version |
|---|---|
| R | 4.5.2 (2025-10-31) |
| RStudio | 2026.1.0.392 (Apple Blossom) |
| Platform | macOS Tahoe 26.3.1, aarch64-apple-darwin20 |

### R Packages

All packages are mentioned in each script.
> **Note:** Install `SpiecEasi` from GitHub for compatibility with Matrix >= 1.5.0:
> ```r
> remotes::install_github("zdk123/SpiecEasi")
> ```


### External Tools

| Tool | Version |
|---|---|
| VSEARCH | 2.30.0 |
| VeryFastTree | 4.0.5 |
| Figaro | 1.1.2 | 

>**Figaro** was used externally to determine optimal DADA2 trimming parameters prior to running `00_DADA2.R`. The resulting parameters (`trim_F = 304`, `trim_R = 176`) are hardcoded in the script. See [Figaro documentation](https://github.com/Zymo-Research/figaro) for installation and usage.
---

## How to Run

Scripts must be run in numbered order. Each script reads from and writes to an `output/` directory which should be created before running. Update file paths in the configuration section at the top of each script before running.

### Input and output summary

| Script | Input | Output |
|---|---|---|
| `00_DADA2.R` | Raw FASTQ files, `Exp2433_Metadata.txt` | `seqtab_nochim.rds`, `read_tracking.csv` |
| `01_16S_Pipeline.R` | `seqtab_nochim.rds` | `step4_phyloseq_object.rds`, QC summary CSVs |
| `02_Rarefaction.R` | `step4_phyloseq_object.rds` | Rarefaction plots (PNG), rarefaction CSVs |
| `03_ThresholdAnalysis.R` | `seqtab_nochim.rds` | Filtering impact matrix CSVs |
| `04_AlphaDiversity.R` | `step4_phyloseq_object.rds` | Alpha diversity plots (PNG), statistics CSVs |
| `05_BetaDiversity.R` | `step4_phyloseq_object.rds` | PCoA plots (PNG), PERMANOVA/BETADISPER CSVs |
| `06_VennDiagram.R` | `step4_phyloseq_object.rds` | Venn diagram plots (PNG) |
| `07_Alluvial.R` | `step4_phyloseq_object.rds` | Alluvial/stacked bar plots (PNG) |
| `08_ANCOMBC.R` | `step4_phyloseq_object.rds` | Volcano plots (PNG), differential abundance CSVs |
| `09_NetCoMi.R` | `step4_phyloseq_object.rds` | Chord diagrams (PNG), network summary CSVs |

> Scripts `02_Rarefaction.R` and `03_ThresholdAnalysis.R` can be run independently of each other after `01_16S_Pipeline.R`. Scripts `04` through `09` all take `step4_phyloseq_object.rds` as input and can be run in any order after `01_16S_Pipeline.R`. Alternatively, download the processed phyloseq object directly from Zenodo [DOI to be added] to skip scripts `00` and `01`.

---
## Contact

For questions regarding the analysis scripts, please open an issue on this repository or contact khanhnguyen@auburn.edu or qkhanh0504@gmail.com.
