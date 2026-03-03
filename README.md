# MetaPop - Python Edition

A pipeline for macro- and micro-diversity analyses of metagenomic-derived populations.

**This is a modernized fork** of the [original MetaPop](https://github.com/metaGmetapop/metapop) with significant improvements for ease of use, particularly on Google Colab.

---

## What's New in This Fork

This fork replaces R dependencies with pure Python implementations, making the pipeline:

- **Easier to install** - No R packages required
- **Cloud-ready** - Three deployment options: Google Colab (interactive), Modal (serverless), local (command-line)
- **Resilient to failures** - Multi-stage Modal pipeline with checkpointing and recovery support
- **More reliable** - Fixed path handling issues for cross-platform compatibility
- **Faster setup** - Fewer dependencies, simpler installation

### Key Changes from Original

| Feature | Original | This Fork |
|---------|----------|-----------|
| Microdiversity calculation | R scripts | Python (pandas/numpy) |
| Macrodiversity calculation | R scripts | Python (pandas/scipy) |
| FST calculation | R scripts | Python |
| Codon bias analysis | R scripts | Python |
| Path handling | Relative (R setwd) | Absolute paths |
| Colab support | Limited | Full support with interactive notebook |

---

## Description

MetaPop analyzes metagenomic data to assess diversity at both the sample level (macrodiversity) and within-population level (microdiversity).

### Pipeline Modules

1. **Preprocessing** - Filter reads by quality, length, and percent identity; calculate coverage metrics
2. **Microdiversity** - Identify SNPs, calculate nucleotide diversity (pi, theta), Tajima's D, pN/pS ratios
3. **Macrodiversity** - Compute normalized abundances, alpha diversity indices, beta diversity distances
4. **Visualization** - Generate plots and summary figures

### Diversity Metrics Calculated

**Alpha Diversity:**
- Richness, Chao index, ACE index
- Shannon's H, Simpson's index, Inverse Simpson's
- Fisher alpha diversity, Pielou's J evenness

**Beta Diversity:**
- Bray-Curtis dissimilarity
- Jaccard distances
- Center Log Ratio (CLR) Euclidean distances

**Microdiversity:**
- Pi (nucleotide diversity)
- Watterson's Theta
- Tajima's D
- pN/pS ratio (selection pressure)
- Fixation index (Fst)

---

## Choosing Your Deployment Option

| Option | Best For | Setup Time | Cost | Data Size |
|--------|----------|------------|------|-----------|
| **Google Colab** | Quick testing, <50 GB data, interactive use | 2 min | Free (with Pro: $10/mo) | Up to 50 GB (Pro: 100 GB) |
| **Modal** | Large datasets (>100 GB), automated runs, production | 10 min | $0.50/GPU-hour | Unlimited (volume-based) |
| **Local** | Full control, GPU available, small data | 5 min | Depends on hardware | As much as your machine can hold |

**Quick recommendation:**
- First time? **Colab** (easiest, no setup)
- Large dataset? **Modal** (recoverable stages, scales to 1TB+)
- Got a beefy machine? **Local** (full control, fast iteration)

---

## Quick Start: Google Colab (Recommended for First-Time Users)

The easiest way to run MetaPop is through our interactive Colab notebook.

### Step 1: Open the Notebook

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/espickle1/metapop/blob/main/notebooks/MetaPop_Colab.ipynb)

Or manually:
1. Go to [Google Colab](https://colab.research.google.com)
2. File → Open notebook → GitHub
3. Enter: `https://github.com/espickle1/metapop`
4. Select `notebooks/MetaPop_Colab.ipynb`

### Step 2: Prepare Your Data

Upload your data to Google Drive with this structure:

```
My Drive/
└── metapop_data/
    ├── bams/                    # Your BAM files
    │   ├── sample1.bam
    │   ├── sample2.bam
    │   └── ...
    └── references/              # Reference FASTA files
        └── reference_genomes.fna
```

**Requirements:**
- **BAM files**: Aligned reads (sorted/unsorted - MetaPop will sort them)
- **Reference FASTA**: The genome(s) your reads were aligned to (`.fa`, `.fasta`, or `.fna`)

### Step 3: Run the Notebook

1. **Run Cell 1** - Install dependencies (~2-3 minutes)
2. **Run Cell 2** - Mount Google Drive (authorize access)
3. **Run Cell 3** - Import libraries
4. **Configure paths** - Update the directory paths to match your data
5. **Adjust parameters** - Use the interactive sliders
6. **Click "Run MetaPop"** - Start the pipeline

### Step 4: View Results

Results are saved to your specified output directory:

```
metapop_results/
└── MetaPop/
    ├── 00.Log_and_Parameters/     # Run logs and settings
    ├── 01.Genomes_and_Genes/      # Processed reference files
    ├── 02.Filtered_Samples/       # Quality-filtered BAM files
    ├── 03.Breadth_and_Depth/      # Coverage statistics
    ├── 05.Variant_Calls/          # Raw variant calls
    ├── 07.Cleaned_SNPs/           # Filtered SNP data
    ├── 08.Codon_Bias/             # Codon usage analysis
    ├── 09.Linked_SNPs/            # Multi-SNP codon analysis
    ├── 10.Microdiversity/         # Pi, theta, Tajima's D, pN/pS
    ├── 11.Macrodiversity/         # Abundance tables, diversity indices
    └── 12.Visualizations/         # Generated plots (PDF)
```

---

## Local Installation

### Requirements

- Python 3.8+
- samtools
- bcftools
- prodigal

### Install via pip

```bash
# Install system dependencies (Ubuntu/Debian)
sudo apt-get install samtools bcftools prodigal

# Install MetaPop
pip install git+https://github.com/espickle1/metapop.git
```

### Install via conda

```bash
# Create environment
conda create -n metapop python=3.10
conda activate metapop

# Install dependencies
conda install -c bioconda samtools bcftools prodigal
pip install git+https://github.com/espickle1/metapop.git
```

---

## Command Line Usage

```bash
metapop --input_samples /path/to/bams \
        --reference /path/to/references \
        --output /path/to/output \
        [OPTIONS]
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--input_samples` | Directory containing BAM files |
| `--reference` | Directory containing reference FASTA file(s) |
| `--output` | Output directory (will be created if needed) |

### Optional Arguments

**Module Control:**

| Argument | Description |
|----------|-------------|
| `--no_micro` | Skip microdiversity analysis |
| `--no_macro` | Skip macrodiversity analysis |
| `--no_viz` | Skip visualization generation |
| `--preprocess_only` | Only run preprocessing |
| `--micro_only` | Only run microdiversity (requires prior preprocessing) |
| `--macro_only` | Only run macrodiversity (requires prior preprocessing) |

**Preprocessing Parameters:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--id_min` | 95 | Minimum percent identity for read filtering |
| `--min_len` | 30 | Minimum alignment length (bp) |
| `--min_cov` | 20 | Minimum breadth of coverage (%) |
| `--min_dep` | 10 | Minimum truncated average depth |
| `--trunc` | 10 | Truncation percentile for TAD calculation |

**Variant Calling Parameters:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--min_obs` | 4 | Minimum observations to call a SNP |
| `--min_pct` | 1 | Minimum population % for SNP call |
| `--min_qual` | 20 | Minimum PHRED quality score |
| `--ref_sample` | (none) | Sample to use as reference for SNPs |

**Other Options:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--threads` | auto | Number of parallel threads |
| `--genes` | (auto) | Pre-computed gene predictions (Prodigal format) |
| `--norm` | (auto) | Library normalization file |
| `--subsample_size` | 10 | Subsampling depth for diversity calculations |

### Example Commands

**Full analysis:**
```bash
metapop --input_samples ./bams --reference ./refs --output ./results --threads 8
```

**Preprocessing only:**
```bash
metapop --input_samples ./bams --reference ./refs --output ./results --preprocess_only
```

**Skip macrodiversity:**
```bash
metapop --input_samples ./bams --reference ./refs --output ./results --no_macro
```

**Stringent filtering:**
```bash
metapop --input_samples ./bams --reference ./refs --output ./results \
        --id_min 98 --min_cov 50 --min_dep 20
```

---

## Running on Modal (Cloud)

MetaPop can run on [Modal](https://modal.com) for large datasets that exceed local or Colab resources. Modal provides serverless cloud compute with persistent storage, ideal for processing multi-gigabyte BAM files without managing infrastructure.

### Modal Deployment Options

Two application files are available for different use cases:

| File | Approach | Stages | Recovery | Best For |
|------|----------|--------|----------|----------|
| **`modal_app.py`** | Multi-stage orchestration | 6 independent stages | ✅ Full recovery support | Large datasets, long pipelines |
| **`modal_app_phase_1.py`** | Monolithic execution | All stages in one function | ❌ No partial recovery | Small-to-medium datasets, simple runs |

**Recommendation:** Use `modal_app.py` for datasets >100 GB. It splits work into checkpointed stages, so a crash at stage 5 won't restart from stage 1.

### Prerequisites & Setup

1. **Install Modal:**
   ```bash
   pip install modal
   ```

2. **Authenticate with Modal:**
   ```bash
   modal setup
   ```
   This opens a browser window to sign in (free account available at https://modal.com).

3. **Create Modal Volume (one-time):**
   ```bash
   modal volume create metapop-data
   ```
   This creates persistent cloud storage where your BAM files and results live.

### Data Preparation: Upload to Modal Volume

Modal volumes are persistent cloud storage mounted at `/mnt/` inside functions. Upload your input data:

```bash
# Create the directory structure in the volume
modal volume put metapop-data /path/to/local/bams /bams
modal volume put metapop-data /path/to/local/refs /refs
```

**Local directory structure (before upload):**
```
~/my_data/
├── bams/
│   ├── sample1.bam
│   ├── sample2.bam
│   └── ...
└── refs/
    └── reference.fasta
```

**Modal volume structure (after upload):**
```
/mnt/metapop-data/
├── bams/
│   ├── sample1.bam
│   ├── sample2.bam
│   └── ...
└── refs/
    └── reference.fasta
```

### Running the Pipeline

Both `modal_app.py` and `modal_app_phase_1.py` use the same command-line interface. The only difference is whether stages are separate (recoverable) or combined.

**Launch the multi-stage pipeline (recommended):**
```bash
modal run modal_app.py \
    --input-samples /mnt/metapop-data/bams \
    --reference /mnt/metapop-data/refs \
    --output /mnt/metapop-data/results \
    --threads 8
```

**Launch the monolithic pipeline (simpler but no recovery):**
```bash
modal run modal_app_phase_1.py \
    --input-samples /mnt/metapop-data/bams \
    --reference /mnt/metapop-data/refs \
    --output /mnt/metapop-data/results \
    --threads 8
```

All standard MetaPop CLI flags are supported (e.g., `--id-min 98`, `--min-cov 50`, `--no-micro`). **Note:** Modal uses hyphens (`--input-samples`) rather than underscores (unlike local CLI).

**Full example with parameters:**
```bash
modal run modal_app.py \
    --input-samples /mnt/metapop-data/bams \
    --reference /mnt/metapop-data/refs \
    --output /mnt/metapop-data/results \
    --id-min 96 \
    --min-cov 30 \
    --min-dep 15 \
    --threads 8
```

### How the Multi-stage Pipeline Works

`modal_app.py` executes the pipeline as 6 independent stages, with checkpointing between each:

1. **Setup** — Initialize directories, combine reference FASTAs, run gene prediction (Prodigal)
2. **Preprocessing** — Filter BAM files by quality, length, coverage, and depth metrics
3. **Variant Calling** — Identify SNPs via samtools, correct consensus bases, analyze codon bias
4. **Microdiversity** — Calculate pi, Tajima's D, pN/pS ratios, FST; detect linked SNPs
5. **Macrodiversity** — Compute normalized abundances, alpha/beta diversity indices
6. **Visualizations** — Generate summary plots and statistical figures

**Key benefit:** Each stage commits results to the Modal volume via `vol.commit()`. If a stage fails, you can resume from the previous successful stage instead of restarting the entire pipeline.

### Recovery from Partial Runs

If a stage crashes (e.g., due to worker preemption or timeout), resume by skipping already-completed stages:

**Stage completion markers** (checked by the pipeline):
- Stage 1 always re-runs (cheap, regenerates paths)
- Stage 2: Check for `04.Depth_per_Pos/` directory
- Stage 3: Check for `05.Variant_Calls/` directory
- Stage 4: Check for `10.Microdiversity/` directory
- Stage 5: Check for `11.Macrodiversity/` directory

**Recovery flags:**

| Flag | Skips | Use When |
|------|-------|----------|
| `--skip-preproc` | Stage 2 | Preprocessing completed but later stage failed |
| `--skip-snp-calling` | Stage 3 | SNP calling completed but microdiversity failed |
| `--no-micro` | Stages 3+4 | Only want macrodiversity results |
| `--no-macro` | Stage 5 | Only want microdiversity results |
| `--no-viz` | Stage 6 | Skip visualization generation (fastest) |
| `--preprocess-only` | Stages 3+ | Only run preprocessing and stop |
| `--viz-only` | Stages 2-5 | Only generate visualizations (requires prior run) |

**Example: Resume after stage 3 (variant calling) crashes:**
```bash
modal run modal_app.py \
    --input-samples /mnt/metapop-data/bams \
    --reference /mnt/metapop-data/refs \
    --output /mnt/metapop-data/results \
    --skip-preproc \
    --threads 8
```
This skips preprocessing and variant calling, resuming at microdiversity.

**Monolithic pipeline note:** `modal_app_phase_1.py` has no partial recovery because it runs all stages in a single Modal function and only commits at the very end. If it crashes midway, you must restart entirely.

### Download Results

Once the pipeline completes, download results from the Modal volume to your local machine:

```bash
# Download the entire results directory
modal volume get metapop-data /results ./local_results

# Or download just specific files
modal volume get metapop-data /results/MetaPop/11.Macrodiversity ./local_macro_results
```

Results are organized identically to local runs (see [Output Files](#output-files) section below).

### Pipeline Stages Reference (Multi-stage)

| Stage | Name | CPU | Memory | Timeout | Description |
|-------|------|-----|--------|---------|-------------|
| 1 | Setup | 2 cores | 4 GB | 30 min | Directory structure, combine FASTAs, Prodigal gene calling |
| 2 | Preprocessing | 8 cores | 16 GB | 4 hrs | Filter reads by identity, length, coverage, depth metrics |
| 3 | Variant Calling | 8 cores | 32 GB | 2 hrs | SNP identification via samtools, consensus correction, codon bias |
| 4 | Microdiversity | 8 cores | 8 GB | 2 hrs | Linked SNPs, Fisher's exact test, pi/theta/Tajima's D/FST |
| 5 | Macrodiversity | 4 cores | 4 GB | 1 hr | Normalized abundances, alpha/beta diversity indices |
| 6 | Visualizations | 2 cores | 4 GB | 30 min | Summary plots, heatmaps, PCoA figures |

### Troubleshooting Modal-Specific Issues

**"Worker interrupted due to preemption"**
- This happens when Modal reclaims a spot instance (temporary cloud server)
- **Solution:** Use the recovery flags above (e.g., `--skip-preproc`) to resume from the last completed stage
- **Mitigation:** The multi-stage design makes preemption less painful — only the current stage restarts

**"Volume not found: metapop-data"**
- The Modal volume doesn't exist in your Modal workspace
- **Solution:** Create it with `modal volume create metapop-data`

**"No such file or directory: /mnt/metapop-data/bams"**
- Data wasn't uploaded to the volume or paths are incorrect
- **Solution:** Verify with `modal volume ls metapop-data` and upload with `modal volume put`

**Stage times out (4 hrs exceeded for preprocessing, etc.)**
- Large BAM files or too-strict filtering parameters
- **Solution:** Try looser parameters or reduce dataset size; monitor via `modal app logs modal_app.py`

**"Out of memory" during variant calling**
- Stage 3 processes all BAM files in parallel; large files + tight parameters = high memory usage
- **Solution:** Reduce `--threads` or use `--preprocess-only` to identify problematic samples first

---

## Input File Formats

### BAM Files

Standard BAM format from any aligner (bowtie2, BWA, minimap2, etc.).

- Files do NOT need to be sorted or indexed (MetaPop handles this)
- File names become sample identifiers (e.g., `sample1.bam` → "sample1")

### Reference FASTA

Standard FASTA format with contig/genome sequences.

```
>contig_name_1
ATCGATCGATCG...
>contig_name_2
GCTAGCTAGCTA...
```

**Important:** Sequence names must match the reference used for alignment.

### Normalization File (Optional)

Tab-separated file with sample names and library sizes:

```
sample1	52000000
sample2	43000000
sample3	73000000
```

If not provided, MetaPop auto-generates this from BAM read counts.

### Gene Predictions (Optional)

Prodigal-format FASTA file. If not provided, MetaPop runs Prodigal automatically.

```
>contig1_1 # 2 # 181 # 1 # ID=1_1;partial=10;start_type=Edge;...
ATGCGATCGATCG...
```

---

## Output Files

### Key Results Files

| File | Description |
|------|-------------|
| `11.Macrodiversity/Alpha_diversity_stats.tsv` | Alpha diversity metrics per sample |
| `11.Macrodiversity/normalized_abundances_table.tsv` | Normalized abundance matrix |
| `11.Macrodiversity/Beta_diversity_*.tsv` | Beta diversity distance matrices |
| `10.Microdiversity/global_gene_microdiversity.tsv` | Gene-level diversity metrics |
| `10.Microdiversity/global_contig_microdiversity.tsv` | Contig-level diversity metrics |
| `07.Cleaned_SNPs/genic_snps.tsv` | All identified SNPs in genes |
| `03.Breadth_and_Depth/*_breadth_and_depth.tsv` | Coverage statistics per sample |

### Visualization PDFs

| File | Contents |
|------|----------|
| `12.Visualizations/preprocessing_summaries.pdf` | Read filtering statistics |
| `12.Visualizations/alpha_diversity_scatterplots.pdf` | Alpha diversity plots |
| `12.Visualizations/PCoA_*.pdf` | Principal coordinate analysis |
| `12.Visualizations/codon_bias_plots.pdf` | Codon usage analysis |

---

## Troubleshooting

### Common Issues

**"FileNotFoundError: MetaPop/..."**
- This indicates a path issue. Make sure you're using the latest version:
  ```bash
  pip install --upgrade --force-reinstall git+https://github.com/espickle1/metapop.git
  ```

**"No BAM files found"**
- Verify your BAM directory path is correct
- Check that files have `.bam` extension (case-sensitive on Linux)

**"Reference directory not found"**
- The `--reference` argument should be a DIRECTORY containing FASTA files, not a file path

**Pipeline hangs or runs very slowly**
- Reduce thread count if running on limited resources
- For large datasets, consider running on a machine with more RAM

**Memory errors**
- Large BAM files require significant RAM
- Try processing fewer samples at once
- Use a machine with more memory (Colab Pro offers 25GB+)

### Getting Help

- **Issues:** [GitHub Issues](https://github.com/espickle1/metapop/issues)
- **Email:** James.Chang@bcm.edu

---

## Citation

If you use MetaPop in your research, please cite the original publication:

> [Original MetaPop citation - check the original repository]

And mention this fork if you used the Python/Colab version:

> MetaPop Python Edition: https://github.com/espickle1/metapop

---

## License

This project is licensed under the same terms as the original MetaPop repository.

---

## Acknowledgments

- Original MetaPop developers at Ohio State University
- Contributors to this Python port
