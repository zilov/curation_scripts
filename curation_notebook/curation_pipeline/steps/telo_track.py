"""Add telomere track to pretext map."""

import glob
from pathlib import Path

from curation_pipeline.context import CurationContext
from curation_pipeline.modules import module_cmd
from curation_pipeline.output import (
    print_done,
    print_info,
    print_step_header,
    print_warning,
)
from curation_pipeline.steps._helpers import _find_pretext_map_in_workdir, _run


# ---------------------------------------------------------------------------
# Public step functions
# ---------------------------------------------------------------------------


def add_telo_track(ctx: CurationContext) -> None:
    """
    Adds a telomere track to the pretext map if telo_*.bed.gz is available.

    Steps:
        1. Glob for telomere BED file from TreeVAL output.
        2. If not found: print warning and return.
        3. Build and run::

               module purge && module load pretextgraph/0.0.7--h4ac6f70_0 && \\
               zcat {telo_bed_gz} | \\
               awk '{ print $1\\t$2\\t$3\\t($3-$2) }' | \\
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