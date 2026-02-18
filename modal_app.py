"""
MetaPop on Modal — Multi-Stage Pipeline

Runs the MetaPop metagenomics pipeline on Modal (modal.com) as a series of
right-sized cloud functions. No existing MetaPop code is modified.

Usage:
    # Upload BAM directory and reference FASTA to Modal Volume ("metapop-data").
    # Paths after the volume name are relative to the volume root.
    # The volume is mounted at /mnt/metapop-data in the container (matching modal shell).
    modal volume put metapop-data /path/to/bam_dir /bams
    modal volume put metapop-data /path/to/ref_fastas /refs

    # Run full pipeline (all MetaPop flags are passed through).
    # Use /mnt/metapop-data/... paths because that is where the volume is mounted.
    modal run modal_app.py --input-samples /mnt/metapop-data/bams --reference /mnt/metapop-data/refs --output /mnt/metapop-data/results --threads 8

    # Run with preprocessing only
    modal run modal_app.py --input-samples /mnt/metapop-data/bams --reference /mnt/metapop-data/refs --output /mnt/metapop-data/results --preprocess-only

    # Download results (volume-relative path, not container path)
    modal volume get metapop-data /results ./local_results
"""

from __future__ import annotations

import os
import sys
import modal

app = modal.App("metapop")

# ---------------------------------------------------------------------------
# Image: install system tools + MetaPop from local source
# ---------------------------------------------------------------------------
metapop_image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("samtools", "bcftools", "prodigal")
    .pip_install(
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "pysam>=0.16.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "scikit-learn>=0.24.0",
    )
    .add_local_dir(".", remote_path="/root/metapop_src", copy=True, ignore=[".git", "toy_dataset", "notebooks", "__pycache__"])
    .run_commands("cd /root/metapop_src && pip install .")
)

vol = modal.Volume.from_name("metapop-data", create_if_missing=True)
VOLUME_MOUNT = "/mnt/metapop-data"

# ---------------------------------------------------------------------------
# Stage 1: Setup — directory prep, combine FASTAs, gene calling, MAG log
# ---------------------------------------------------------------------------
@app.function(image=metapop_image, volumes={VOLUME_MOUNT: vol}, cpu=2, memory=4096, timeout=1800)
def setup(config: dict) -> dict:
    """Create output dirs, combine references, call genes, build MAG log."""
    vol.reload()

    import metapop.metapop_helper_functions as hf

    out = config["output"]
    ref = config["reference"]
    genes_arg = config["genes"]
    no_mac = config["no_macro"]
    norm_arg = config["norm"]
    bams = config["input_samples"]
    threads = config["threads"]

    hf.main_dir_prep(out)

    # Normalization file (needed unless --no_macro)
    norm_file = norm_arg
    if not no_mac:
        if norm_file == "":
            norm_file = hf.produce_default_normalization_file(out, bams, threads)
        else:
            norm_file = os.path.abspath(norm_file)

    joined_fastas = hf.multi_to_single_fasta(out, ref)

    if genes_arg == "":
        reference_genes = hf.gene_calls(joined_fastas, out)
    else:
        reference_genes = os.path.abspath(genes_arg)
        print("Using", reference_genes, "as the genes file.")

    treat_as_mag = False
    mag_contig_dict, mag_length_dict = hf.create_mag_log(
        ref, out + "/MetaPop/00.Log_and_Parameters/mags.txt", treat_as_mag
    )

    # Write parameter log
    _write_run_settings(config, out, joined_fastas, reference_genes, norm_file)

    vol.commit()

    return {
        "joined_fastas": joined_fastas,
        "reference_genes": reference_genes,
        "mag_contig_dict": mag_contig_dict,
        "mag_length_dict": mag_length_dict,
        "norm_file": norm_file,
    }


def _write_run_settings(config, out, joined_fastas, reference_genes, norm_file):
    """Write the run_settings.tsv parameter log (mirrors main.py logic)."""
    path = out + "/MetaPop/00.Log_and_Parameters/run_settings.tsv"
    with open(path, "w") as param:
        print("parameter", "setting", sep="\t", file=param)
        print("Directory", os.path.normpath(out + "/MetaPop/"), sep="\t", file=param)
        print("Samtools Location", "", sep="\t", file=param)
        print("BCFTools Location", "", sep="\t", file=param)
        print("Library Location", "", sep="\t", file=param)
        print("Assembly", joined_fastas, sep="\t", file=param)
        print("Genes", reference_genes, sep="\t", file=param)
        print("Normalization File", norm_file, sep="\t", file=param)
        print("ID Cutoff", config["min_id"], sep="\t", file=param)
        print("Min. Read Length", config["min_len"], sep="\t", file=param)
        print("Coverage", config["min_cov"], sep="\t", file=param)
        print("Depth", config["min_dep"], sep="\t", file=param)
        print("Truncation", config["trunc"], sep="\t", file=param)
        print("Variant Base Call Cutoff", config["min_qual"], sep="\t", file=param)
        print("Subsample Size", config["subsample_size"], sep="\t", file=param)
        print("BP for Detect", config["min_bp_len"], sep="\t", file=param)
        print("Macrodiversity with whole genomes", config["whole_genome"], sep="\t", file=param)
        print("Macrodiversity percent coverage minimum", config["genome_cov_cutoff"], sep="\t", file=param)
        print("All Genomes Plotted", config["plot_all"], sep="\t", file=param)
        print("SNP Scale", config["snp_scale"], sep="\t", file=param)
        print("Threads", config["threads"], sep="\t", file=param)


# ---------------------------------------------------------------------------
# Stage 2: Preprocessing — filter reads by identity, length, coverage, depth
# ---------------------------------------------------------------------------
@app.function(image=metapop_image, volumes={VOLUME_MOUNT: vol}, cpu=8, memory=16384, timeout=14400)
def preprocess(config: dict, setup_result: dict):
    """Filter BAMs by percent identity, read length, breadth, and depth."""
    vol.reload()

    import metapop.metapop_filter as mf

    out = config["output"]
    bams = config["input_samples"]
    threads = config["threads"]

    filter_command_base = [
        "input_bam",
        out,
        int(config["min_len"]),
        float(config["min_id"]),
        config["is_global"],
    ]
    bdb = [
        int(config["min_cov"]),
        int(config["min_dep"]),
        False,  # treat_as_mag is always False currently
        float(config["trunc"]),
    ]

    mf.filt(
        bams,
        filter_command_base,
        bdb,
        threads,
        setup_result["mag_contig_dict"],
        setup_result["mag_length_dict"],
        setup_result["joined_fastas"],
    )

    vol.commit()


# ---------------------------------------------------------------------------
# Stage 3: Variant calling — identify SNPs, correct consensus, codon bias
# ---------------------------------------------------------------------------
@app.function(image=metapop_image, volumes={VOLUME_MOUNT: vol}, cpu=8, memory=16384, timeout=7200)
def call_variants(config: dict, setup_result: dict) -> dict:
    """Call variant positions, correct reference bases, compute codon bias."""
    vol.reload()

    import metapop.metapop_snp_call as snp

    out = config["output"]
    threads = config["threads"]

    joined_fastas, reference_genes = snp.call_variant_positions(
        out,
        setup_result["joined_fastas"],
        int(config["min_obs"]),
        int(config["min_qual"]),
        float(config["min_pct"]),
        threads,
        config["ref_sample"],
        setup_result["reference_genes"],
    )

    vol.commit()

    return {"joined_fastas": joined_fastas, "reference_genes": reference_genes}


# ---------------------------------------------------------------------------
# Stage 4: Microdiversity — linked SNPs, Fisher test, pi/theta, FST
# ---------------------------------------------------------------------------
@app.function(image=metapop_image, volumes={VOLUME_MOUNT: vol}, cpu=8, memory=8192, timeout=7200)
def run_microdiversity(config: dict, variant_result: dict):
    """Mine linked SNPs, run Fisher test, calculate microdiversity and FST."""
    vol.reload()

    import metapop.metapop_helper_functions as hf
    import metapop.metapop_mine_reads as mr
    import metapop.metapop_fisher as fisher
    import metapop.metapop_microdiversity as microdiv
    import metapop.metapop_fst as fst

    out = config["output"]
    threads = config["threads"]
    joined_fastas = variant_result["joined_fastas"]
    reference_genes = variant_result["reference_genes"]

    hf.micro_prep(out)

    linked_file = mr.do_mine_reads(out, threads)
    fisher.process_linked_snps(linked_file)

    microdiv.run_microdiversity(
        out,
        joined_fastas,
        reference_genes,
        int(config["min_cov"]),
        int(config["min_dep"]),
        int(config["subsample_size"]),
        threads,
    )

    microdiv_file = os.path.normpath(
        out + "/MetaPop/10.Microdiversity/global_raw_microdiversity_data_snp_loci_only.tsv"
    )
    lengths_file = os.path.normpath(
        out + "/MetaPop/10.Microdiversity/global_contig_microdiversity.tsv"
    )
    fst.perform_fst(microdiv_file, lengths_file, out, threads)

    vol.commit()


# ---------------------------------------------------------------------------
# Stage 5: Macrodiversity — abundance, alpha/beta diversity, ordination
# ---------------------------------------------------------------------------
@app.function(image=metapop_image, volumes={VOLUME_MOUNT: vol}, cpu=4, memory=4096, timeout=3600)
def run_macrodiversity(config: dict, setup_result: dict, variant_result: dict | None):
    """Calculate macro-level diversity metrics and generate ordination plots."""
    vol.reload()

    import metapop.metapop_helper_functions as hf
    import metapop.metapop_macrodiversity as macrodiv

    out = config["output"]
    threads = config["threads"]

    # Use base-corrected genomes/genes if variant calling ran, otherwise originals
    if variant_result is not None:
        joined_fastas = variant_result["joined_fastas"]
        reference_genes = variant_result["reference_genes"]
    else:
        joined_fastas = os.path.abspath(setup_result["joined_fastas"])
        reference_genes = os.path.abspath(setup_result["reference_genes"])

    hf.macro_prep(out)

    macro_results = macrodiv.run_macrodiversity(
        out,
        joined_fastas,
        reference_genes,
        setup_result["norm_file"],
        int(config["min_cov"]),
        int(config["min_dep"]),
        config["whole_genome"] == "1",
        int(config["genome_cov_cutoff"]),
        int(config["min_bp_len"]),
        threads,
    )

    viz_path = os.path.normpath(out + "/MetaPop/12.Visualizations")
    if not os.path.exists(viz_path):
        os.makedirs(viz_path, exist_ok=True)
    macrodiv.generate_macrodiversity_visualizations(macro_results, viz_path)

    vol.commit()


# ---------------------------------------------------------------------------
# Stage 6: Visualizations — preprocessing, microdiversity, codon bias plots
# ---------------------------------------------------------------------------
@app.function(image=metapop_image, volumes={VOLUME_MOUNT: vol}, cpu=2, memory=4096, timeout=1800)
def run_visualizations(config: dict):
    """Generate summary visualizations for preprocessing and microdiversity."""
    vol.reload()

    import metapop.metapop_helper_functions as hf
    import metapop.metapop_visualizations as viz

    out = config["output"]
    bams = config["input_samples"]
    no_mic = config["no_micro"]

    hf.viz_prep(out)

    samples = os.listdir(bams)
    samples = [os.path.normpath(bams + "/" + s) for s in samples]
    names = hf.get_base_names(samples)

    viz.visualize_preprocessing(out, os.path.abspath(os.path.normpath(bams)), names)

    if not no_mic:
        viz.visualize_microdiversity(out, config["plot_all"] == "1", config["snp_scale"])
        viz.visualize_codon_bias(out)

    vol.commit()


# ---------------------------------------------------------------------------
# Local entrypoint — parse CLI args and orchestrate remote stages
# ---------------------------------------------------------------------------
@app.local_entrypoint()
def main(
    input_samples: str = "",
    reference: str = "",
    genes: str = "",
    norm: str = "",
    output: str = "/mnt/metapop-data/results",
    threads: int = 8,
    # Filter options
    id_min: float = 95,
    min_len: int = 50,
    min_cov: int = 20,
    min_dep: int = 10,
    trunc: float = 10,
    is_global: bool = False,
    # Variant calling options
    min_qual: int = 20,
    min_obs: int = 2,
    min_pct: float = 1,
    ref_sample: str = "",
    subsample_size: int = 10,
    # Macrodiversity options
    whole_genomes: bool = False,
    genome_detection_cutoff: int = 0,
    minimum_bases_for_detection: int = 5000,
    # Visualization options
    plot_all: bool = False,
    snp_scale: str = "local",
    # Flow-control flags
    preprocess_only: bool = False,
    no_micro: bool = False,
    no_macro: bool = False,
    no_viz: bool = False,
    viz_only: bool = False,
    skip_preproc: bool = False,
    skip_snp_calling: bool = False,
):
    if not input_samples:
        print("MetaPop needs aligned reads. Provide --input-samples.")
        sys.exit(1)
    if not reference:
        print("MetaPop needs reference genomes. Provide --reference.")
        sys.exit(1)

    # Normalize output path (strip trailing slash)
    if output.endswith("/"):
        output = output[:-1]

    # Convert bool flags to "0"/"1" strings to match main.py convention
    whole_genome_str = "1" if whole_genomes else "0"
    plot_all_str = "1" if plot_all else "0"

    config = {
        "input_samples": input_samples,
        "reference": reference,
        "genes": genes,
        "norm": norm,
        "output": output,
        "threads": threads,
        "min_id": id_min,
        "min_len": min_len,
        "min_cov": min_cov,
        "min_dep": min_dep,
        "trunc": trunc,
        "is_global": is_global,
        "min_qual": min_qual,
        "min_obs": min_obs,
        "min_pct": min_pct,
        "ref_sample": ref_sample,
        "subsample_size": subsample_size,
        "whole_genome": whole_genome_str,
        "genome_cov_cutoff": genome_detection_cutoff,
        "min_bp_len": minimum_bases_for_detection,
        "plot_all": plot_all_str,
        "snp_scale": snp_scale,
        "no_macro": no_macro,
        "no_micro": no_micro,
    }

    # ---- Stage 1: Setup (always runs) ----
    print("[1/6] Setup: preparing directories, combining FASTAs, calling genes...")
    setup_result = setup.remote(config)
    print("       Setup complete.")

    # ---- Early exit: viz_only skips all data generation ----
    if viz_only:
        if not no_viz:
            print("[6/6] Visualizations (viz-only mode)...")
            run_visualizations.remote(config)
            print("       Visualizations complete.")
        print("MetaPop (Modal) finished.")
        return

    # ---- Stage 2: Preprocessing ----
    if not skip_preproc:
        print("[2/6] Preprocessing: filtering reads by identity, length, coverage, depth...")
        preprocess.remote(config, setup_result)
        print("       Preprocessing complete.")

        if preprocess_only:
            print("MetaPop (Modal) finished (preprocess-only mode).")
            return
    else:
        print("[2/6] Preprocessing: SKIPPED (--skip-preproc)")

    # ---- Stage 3: Variant calling ----
    variant_result = None
    if not no_micro:
        if not skip_snp_calling:
            print("[3/6] Variant calling: identifying SNPs and correcting consensus...")
            variant_result = call_variants.remote(config, setup_result)
            print("       Variant calling complete.")
        else:
            print("[3/6] Variant calling: SKIPPED (--skip-snp-calling)")
            # Use original paths (as absolute) when skipping
            variant_result = {
                "joined_fastas": os.path.abspath(setup_result["joined_fastas"]),
                "reference_genes": os.path.abspath(setup_result["reference_genes"]),
            }

        # ---- Stage 4: Microdiversity ----
        print("[4/6] Microdiversity: linked SNPs, Fisher test, pi/theta/FST...")
        run_microdiversity.remote(config, variant_result)
        print("       Microdiversity complete.")
    else:
        print("[3/6] Variant calling: SKIPPED (--no-micro)")
        print("[4/6] Microdiversity: SKIPPED (--no-micro)")

    # ---- Stage 5: Macrodiversity ----
    if not no_macro:
        print("[5/6] Macrodiversity: abundance, alpha/beta diversity...")
        run_macrodiversity.remote(config, setup_result, variant_result)
        print("       Macrodiversity complete.")
    else:
        print("[5/6] Macrodiversity: SKIPPED (--no-macro)")

    # ---- Stage 6: Visualizations ----
    if not no_viz:
        print("[6/6] Visualizations: generating summary plots...")
        run_visualizations.remote(config)
        print("       Visualizations complete.")
    else:
        print("[6/6] Visualizations: SKIPPED (--no-viz)")

    print("MetaPop (Modal) finished.")
