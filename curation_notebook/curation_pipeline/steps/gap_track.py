"""Add gap track to pretext map."""

from curation_pipeline.context import CurationContext
from curation_pipeline.modules import module_cmd
from curation_pipeline.output import (
    print_done,
    print_next_step,
    print_step_header,
)
from curation_pipeline.steps._helpers import _find_pretext_map_in_workdir, _run

# ---------------------------------------------------------------------------
# Farm script paths
# ---------------------------------------------------------------------------

_HAP_BEDGRAPH_SCRIPT = "/nfs/users/nfs_d/dz11/hap_bedgraph.py"


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