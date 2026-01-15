# MetaPop - Python Edition

A pipeline for macro- and micro-diversity analyses of metagenomic-derived populations.

**This is a modernized fork** of the [original MetaPop](https://github.com/metaGmetapop/metapop) with significant improvements for ease of use, particularly on Google Colab.

---

## What's New in This Fork

This fork replaces R dependencies with pure Python implementations, making the pipeline:

- **Easier to install** - No R packages required
- **Cloud-ready** - Optimized for Google Colab with interactive notebook
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

## Quick Start: Google Colab (Recommended)

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
