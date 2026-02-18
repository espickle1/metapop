"""
MetaPop on Modal — Monolithic Wrapper (Phase 1 / Option A)

Runs the entire MetaPop pipeline as a single Modal function.
No existing MetaPop code is modified.

Usage:
    # Upload BAM directory and reference FASTA to Modal Volume ("metapop-data").
    # Paths after the volume name are relative to the volume root.
    # The volume is mounted at /mnt/metapop-data in the container (matching modal shell).
    modal volume put metapop-data /path/to/bam_dir /bams
    modal volume put metapop-data /path/to/ref_fastas /refs

    # Run full pipeline (all MetaPop flags are passed through).
    # Use /mnt/metapop-data/... paths because that is where the volume is mounted.
    modal run modal_app_phase_1.py --input-samples /mnt/metapop-data/bams --reference /mnt/metapop-data/refs --output /mnt/metapop-data/results --threads 8

    # Run with preprocessing only
    modal run modal_app_phase_1.py --input-samples /mnt/metapop-data/bams --reference /mnt/metapop-data/refs --output /mnt/metapop-data/results --preprocess-only

    # Download results (volume-relative path, not container path)
    modal volume get metapop-data /results ./local_results
"""

from __future__ import annotations

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
# Single Modal function: run the entire MetaPop pipeline
# ---------------------------------------------------------------------------
@app.function(image=metapop_image, volumes={VOLUME_MOUNT: vol}, cpu=8, memory=16384, timeout=14400)
def run_metapop(cli_args: list[str]):
    """Run the full MetaPop pipeline by invoking main() with the given CLI args."""
    vol.reload()

    import metapop.metapop_main as mp

    # Inject CLI args so argparse sees them
    sys.argv = ["metapop"] + cli_args

    mp.main()

    vol.commit()


# ---------------------------------------------------------------------------
# Local entrypoint: collect CLI args and forward them to the remote function
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
    library: str = "",
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
    just_codon_bias: bool = False,
):
    if not input_samples:
        print("MetaPop needs aligned reads. Provide --input-samples.")
        sys.exit(1)
    if not reference:
        print("MetaPop needs reference genomes. Provide --reference.")
        sys.exit(1)

    # Build the CLI argument list that metapop main() expects
    cli_args = [
        "--input_samples", input_samples,
        "--reference", reference,
        "--output", output,
        "--threads", str(threads),
        "--id_min", str(id_min),
        "--min_len", str(min_len),
        "--min_cov", str(min_cov),
        "--min_dep", str(min_dep),
        "--trunc", str(trunc),
        "--min_qual", str(min_qual),
        "--min_obs", str(min_obs),
        "--min_pct", str(min_pct),
        "--subsample_size", str(subsample_size),
        "--genome_detection_cutoff", str(genome_detection_cutoff),
        "--minimum_bases_for_detection", str(minimum_bases_for_detection),
        "--snp_scale", snp_scale,
    ]

    # Optional string args (only include if non-empty)
    if genes:
        cli_args += ["--genes", genes]
    if norm:
        cli_args += ["--norm", norm]
    if ref_sample:
        cli_args += ["--ref_sample", ref_sample]
    if library:
        cli_args += ["--library", library]

    # Boolean flags (only include if True)
    if is_global:
        cli_args.append("--global")
    if whole_genomes:
        cli_args.append("--whole_genomes")
    if plot_all:
        cli_args.append("--plot_all")
    if preprocess_only:
        cli_args.append("--preprocess_only")
    if no_micro:
        cli_args.append("--no_micro")
    if no_macro:
        cli_args.append("--no_macro")
    if no_viz:
        cli_args.append("--no_viz")
    if viz_only:
        cli_args.append("--viz_only")
    if skip_preproc:
        cli_args.append("--skip_preproc")
    if skip_snp_calling:
        cli_args.append("--skip_snp_calling")
    if just_codon_bias:
        cli_args.append("--just_codon_bias")

    print("Running MetaPop on Modal with args:")
    print("  metapop", " ".join(cli_args))
    print()

    run_metapop.remote(cli_args)

    print("MetaPop (Modal) finished.")
