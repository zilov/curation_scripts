"""Add bedgraph track to pretext map."""

from pathlib import Path

from curation_pipeline.context import CurationContext
from curation_pipeline.modules import module_cmd
from curation_pipeline.output import (
    print_done,
    print_info,
    print_step_header,
)
from curation_pipeline.steps._helpers import _find_pretext_map_in_workdir, _run


# ---------------------------------------------------------------------------
# Public step functions
# ---------------------------------------------------------------------------


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