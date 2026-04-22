"""Step: submit QV and k-mer completeness analysis via bsub."""

from __future__ import annotations

from curation_pipeline.context import CurationContext
from curation_pipeline.output import print_done, print_next_step, print_step_header
from curation_pipeline.steps._helpers import _submit_bsub


def run_qv(ctx: CurationContext) -> None:
    """
    Submits QV and k-mer completeness analysis via bsub.

    Notebook source: ``pre_and_post_curation()`` — ``run_qv_analysis`` section.

    Steps:
        1. Build and submit::

               cd {ctx.workdir} && kmer_completeness.bash {ctx.tol_id} {ctx.release_version}

    Prints:
        Step header, bsub command, job ID.
    Next step hint: ``validate_curated_files(ctx)``
    """
    print_step_header(ctx.ticket_id, ctx.tol_id, "QV analysis")

    inner_cmd = (
        f"cd {ctx.workdir} && kmer_completeness.bash {ctx.tol_id} {ctx.release_version}"
    )
    bsub_opts = "-q normal -M 8000 -R\"select[mem>8000] rusage[mem=8000]\""
    _submit_bsub(inner_cmd, bsub_opts, ctx.print_only)

    print_done("QV job submitted")
    print_next_step("validate_curated_files(ctx)")
