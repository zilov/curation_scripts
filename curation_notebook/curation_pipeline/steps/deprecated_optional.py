"""Optional curation steps (DEPRECATED - use individual commands instead)."""

import glob
import subprocess
from pathlib import Path

from curation_pipeline.context import CurationContext
from curation_pipeline.modules import module_cmd
from curation_pipeline.output import (
    console,
    print_done,
    print_info,
    print_next_step,
    print_step_header,
    print_warning,
)

# ---------------------------------------------------------------------------
# Farm script paths
# ---------------------------------------------------------------------------

_HAP_BEDGRAPH_SCRIPT = "/nfs/users/nfs_d/dz11/hap_bedgraph.py"
_GET_NEAREST_COMPARATOR = "/software/grit/projects/vgp_curation_scripts/get_nearest_comparator.rb"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _find_pretext_map_in_workdir(ctx: CurationContext) -> Path:
    """
    Returns the HR pretext map that was copied to workdir.

    Raises FileNotFoundError if not found (unless print_only).
    """
    pattern = str(ctx.workdir / f"{ctx.tol_id}*hr.pretext")
    if ctx.print_only:
        return ctx.workdir / f"{ctx.tol_id}_hr.pretext"
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(
            f"No HR pretext map found in workdir: {pattern}\n"
            "Run copy_pretext_maps first."
        )
    return Path(sorted(matches)[-1])


def _run(cmd: str, print_only: bool = False) -> None:
    """Print cmd; execute it unless print_only is True."""
    console.print(f"\n[yellow]Command:[/yellow] [green]{cmd}[/green]")
    if print_only:
        return
    subprocess.run(cmd, shell=True, check=True)


# ---------------------------------------------------------------------------
# Public step functions
# ---------------------------------------------------------------------------


def add_gap_track(ctx: CurationContext) -> None:
    """
    Adds a gap track to the pretext map using hap_bedgraph.py + PretextGraph.

    Prerequisite: the pretext map has already been copied to ``ctx.workdir``.

    Steps:
        1. Locate the HR pretext map in workdir.
        2. Build and run::

               module purge && module load pretextgraph/0.0.7--h4ac6f70_0 && \\
               python3 /nfs/users/nfs_d/dz11/hap_bedgraph.py \\
                   {ctx.workdir}/original.fa | \\
               PretextGraph -i {pretext_map_path} -n gap

    Prints:
        Step header, command executed.
    Next step hint: ``add_telo_track(ctx)``
    """
    print_step_header(ctx.ticket_id, ctx.tol_id, "Add gap track")

    pretext_map = _find_pretext_map_in_workdir(ctx)
    original_fa = ctx.workdir / "original.fa"
    ml = module_cmd("PRETEXTGRAPH")

    cmd = (
        f"{ml} && "
        f"python3 {_HAP_BEDGRAPH_SCRIPT} {original_fa} | "
        f"PretextGraph -i {pretext_map} -n gap"
    )
    _run(cmd, ctx.print_only)
    print_done("Gap track added.")
    print_next_step("add_telo_track(ctx)")


def add_telo_track(ctx: CurationContext) -> None:
    """
    Adds a telomere track to the pretext map if telo_*.bed.gz is available.

    Steps:
        1. Glob for telomere BED file from TreeVAL output.
        2. If not found: print warning and return.
        3. Build and run::

               module purge && module load pretextgraph/0.0.7--h4ac6f70_0 && \\
               zcat {telo_bed_gz} | \\
               awk '{ print $1"\\t"$2"\\t"$3"\\t"($3-$2) }' | \\
               PretextGraph -i {pretext_map_path} -n telomere

    Prints:
        Step header, telo file found (or warning if absent), command executed.
    """
    print_step_header(ctx.ticket_id, ctx.tol_id, "Add telo track")

    telo_pattern = str(
        ctx.assembly_draft_dir / "treeval" / "*" / "tv_output1" / "treeval_upload" / "telo_*.bed.gz"
    )

    if ctx.print_only:
        print_info("Telo pattern", telo_pattern)
    else:
        telo_files = glob.glob(telo_pattern)
        if not telo_files:
            print_warning(f"No telo BED file found at: {telo_pattern} — skipping telo track.")
            return
        telo_bed_gz = Path(sorted(telo_files)[-1])
        print_info("Telo file", str(telo_bed_gz))

    pretext_map = _find_pretext_map_in_workdir(ctx)
    ml = module_cmd("PRETEXTGRAPH")

    telo_arg = telo_pattern if ctx.print_only else str(telo_bed_gz)
    cmd = (
        f"{ml} && "
        f"zcat {telo_arg} | "
        r"awk '{ print $1\"\t\"$2\"\t\"$3\"\t\"($3-$2) }' | "
        f"PretextGraph -i {pretext_map} -n telomere"
    )
    _run(cmd, ctx.print_only)
    print_done("Telo track added.")


def add_bedgraph_track(ctx: CurationContext, bedgraph_path: str) -> None:
    """
    Adds an arbitrary bedgraph track to the pretext map.

    Steps:
        1. Verify ``bedgraph_path`` exists (skipped in print_only mode).
        2. Derive a track name from the bedgraph filename stem.
        3. Build and run::

               module purge && module load pretextgraph/0.0.7--h4ac6f70_0 && \\
               cat {bedgraph_path} | \\
               PretextGraph -i {pretext_map_path} -n {track_name}

    Prints:
        Step header, track name, command executed.
    """
    print_step_header(ctx.ticket_id, ctx.tol_id, "Add bedgraph track")

    bg_path = Path(bedgraph_path)
    if not ctx.print_only and not bg_path.exists():
        raise FileNotFoundError(f"Bedgraph file not found: {bg_path}")

    track_name = bg_path.stem
    print_info("Bedgraph file", str(bg_path))
    print_info("Track name", track_name)

    pretext_map = _find_pretext_map_in_workdir(ctx)
    ml = module_cmd("PRETEXTGRAPH")

    cmd = (
        f"{ml} && "
        f"cat {bg_path} | "
        f"PretextGraph -i {pretext_map} -n {track_name}"
    )
    _run(cmd, ctx.print_only)
    print_done(f"Bedgraph track '{track_name}' added.")


# run_sex_matcher and _print_sex_summary live in sex_matcher.py
from curation_pipeline.steps.sex_matcher import _print_sex_summary, run_sex_matcher  # noqa: F401

# microchromosome steps live in microchromosome.py
from curation_pipeline.steps.microchromosome import (  # noqa: F401
    run_microchromosome_curation,
    run_microchromosome_post_curation,
)


def _clean_species_name(species: str) -> str:
    """
    Normalise a species name for use with get_nearest_comparator.rb.

    Rules:
    - Strip anything in parentheses (alternative names).
    - Take the first two words.
    - If the second word is "sp." or contains any digit, use only the first word.

    Examples:
        "Anopheles rufipes"                       -> "Anopheles rufipes"
        "Anopheles sp. 123"                        -> "Anopheles"
        "Heliconius melpomene (postman butterfly)" -> "Heliconius melpomene"
        "Genus sp. (some form)"                    -> "Genus"
    """
    import re
    # Remove parenthetical remarks
    cleaned = re.sub(r"\(.*?\)", "", species).strip()
    words = cleaned.split()
    if len(words) == 0:
        return species.strip()
    if len(words) == 1:
        return words[0]
    second = words[1]
    if second == "sp." or any(ch.isdigit() for ch in second):
        return words[0]
    return f"{words[0]} {second}"


def find_closest_reference(ctx: CurationContext, number: int = 1) -> None:
    """
    Finds (and downloads) the closest reference genome from NCBI for the
    species being curated.

    The reference FASTA is downloaded into ``{ctx.workdir}/reference/``.
    The script must be run from that directory, so we ``cd`` into it first.

    Command::

        mkdir -p {ctx.workdir}/reference && \\
        cd {ctx.workdir}/reference && \\
        /software/grit/projects/vgp_curation_scripts/get_nearest_comparator.rb \\
            -s "{ctx.species}" -d -n {number}

    Prints:
        Step header, command executed, path to reference directory.
    """
    print_step_header(ctx.ticket_id, ctx.tol_id, "Find closest reference")

    ref_dir = ctx.workdir / "reference"
    species_query = _clean_species_name(ctx.species)
    print_info("Reference dir", str(ref_dir))
    print_info("Species (raw)", ctx.species)
    print_info("Species (query)", species_query)

    cmd = (
        f"mkdir -p {ref_dir} && "
        f"cd {ref_dir} && "
        f"{_GET_NEAREST_COMPARATOR} -s \"{species_query}\" -d -n {number}"
    )
    _run(cmd, ctx.print_only)
    print_done(f"Reference downloaded to {ref_dir}")


def run_fastga(ctx: CurationContext, reference_path: str | None = None) -> None:
    """
    Runs a FastGA dot-plot comparison of the curated assembly vs. a reference.

    Notebook source: ``run_fastga()`` and ``scp_fastga()`` functions.

    Steps:
        1. Determine the primary curated hap1 FASTA in ``ctx.workdir``.
        2. If reference is gzipped and not yet decompressed/reheadered::

               ml grit && gunzip {ref} && \
               reheader {ref.replace('.gz', '')} > {ref_reheader}

        3. Build bsub command::

               bsub -G team135 -n 8 -e e_fastga -o o_fastga \
                   -M 24000 -R'...' \
                   /software/grit/projects/vgp_curation_scripts/FastGA_dot_dgenies.sh \
                   {ref_reheader} {hap1_fa} {run_prefix} {outdir}

        4. Print command and submit via subprocess.
        5. If a ``fastga/`` output directory already exists: print scp commands
           for the index and PAF files::

               scp {ctx.farm_host}:{fastga_outdir}/*.f*a.idx ~/curations/.../
               scp {ctx.farm_host}:{fastga_outdir}/*FastGA.paf ~/curations/.../

    Prints:
        Step header, bsub command, scp commands (if output exists).
    """
    print_step_header(ctx.ticket_id, ctx.tol_id, "Run FastGA")

    # --- find curated hap1 fa ---
    hap1_pattern = str(ctx.workdir / f"{ctx.tol_id}*{ctx.hap1_prefix}*.curated.fa")
    if ctx.print_only:
        hap1_fa = ctx.workdir / f"{ctx.tol_id}.{ctx.hap1_prefix}.primary.curated.fa"
    else:
        hap1_matches = glob.glob(hap1_pattern)
        if not hap1_matches:
            raise FileNotFoundError(f"No curated hap1 FASTA found: {hap1_pattern}")
        hap1_fa = Path(sorted(hap1_matches)[-1])
    print_info("Curated hap1 FASTA", str(hap1_fa))

    # --- find reference ---
    ref_dir = ctx.workdir / "reference"
    ref_patterns = [
        str(ctx.workdir / "GC*.fna.gz"),
        str(ctx.workdir / "GC*.fna"),
        str(ref_dir / "*.fna.gz"),
        str(ref_dir / "*.fna"),
    ]
    ref_path = None
    for pattern in ref_patterns:
        if ctx.print_only:
            ref_path = Path(pattern.replace("*", "example"))
            break
        matches = glob.glob(pattern)
        if matches:
            ref_path = Path(sorted(matches)[-1])
            break

    if ref_path is None or (not ctx.print_only and not ref_path.exists()):
        print_info("No reference found, downloading closest reference")
        find_closest_reference(ctx)
        # After download, find again
        ref_matches = glob.glob(str(ref_dir / "*.fna.gz")) + glob.glob(str(ref_dir / "*.fna"))
        if not ref_matches:
            raise FileNotFoundError(f"No reference downloaded to {ref_dir}")
        ref_path = Path(sorted(ref_matches)[-1])

    print_info("Reference FASTA", str(ref_path))

    # --- prepare reference (gunzip + reheader if needed) ---
    ref_prefix = ref_path.stem.split(".")[0]
    assembly_prefix = hap1_fa.stem.split(".")[0]
    run_prefix = f"{ref_prefix}_vs_{assembly_prefix}"
    outdir = ctx.workdir / "fastga"
    ref_reheader = ctx.workdir / f"{ref_prefix}_reheader.fna"

    if ref_path.suffix == ".gz":
        gunzip_cmd = f"gunzip {ref_path}"
        reheader_cmd = f"reheader {ref_path.with_suffix('')} > {ref_reheader}"
        prep_cmd = f"{gunzip_cmd} && {reheader_cmd}"
    else:
        prep_cmd = f"reheader {ref_path} > {ref_reheader}"

    ml_grit = module_cmd("GRIT")
    fastga_script = "/software/grit/projects/vgp_curation_scripts/FastGA_dot_dgenies.sh"
    bsub_cmd = (
        f"bsub -G team135 -n 8 -e e_fastga -o o_fastga "
        f"-M 24000 -R'select[mem>24000] rusage[mem=24000] span[hosts=1]' "
        f"{fastga_script} {ref_reheader} {hap1_fa} {run_prefix} {outdir}"
    )

    full_cmd = f"{ml_grit} && {prep_cmd} && {bsub_cmd}"
    _run(full_cmd, ctx.print_only)

    # --- if output exists, print scp commands ---
    if ctx.print_only or outdir.exists():
        scp_local_dir = f"~/curations/work/{ctx.tol_id}"
        idx_files = glob.glob(str(outdir / "*f*a.idx")) if not ctx.print_only else [str(outdir / "example.f*a.idx")]
        paf_files = glob.glob(str(outdir / "*FastGA.paf")) if not ctx.print_only else [str(outdir / "example.FastGA.paf")]
        files_to_scp = idx_files + paf_files
        if files_to_scp:
            scp_cmds = [f"scp {ctx.farm_host}:{f} {scp_local_dir}" for f in files_to_scp]
            print_info("Scp FastGA results to local", " && ".join(scp_cmds))

    print_done("FastGA submitted.")



