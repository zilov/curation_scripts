"""
Command-line interface for the curation pipeline.

Usage:
    curate <ticket_id> summary
    curate <ticket_id> pre
    curate <ticket_id> post

Config is loaded from ~/.grit_curation_config.yaml automatically.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml


def _load_user_config() -> dict:
    config_path = Path.home() / ".grit_curation_config.yaml"
    if not config_path.exists():
        print(
            f"[error] User config not found: {config_path}\n"
            "Create it from the template — see description.md for the required fields.",
            file=sys.stderr,
        )
        sys.exit(1)
    with open(config_path) as f:
        return yaml.safe_load(f)


def _build(args: argparse.Namespace, user_config: dict):
    from curation_pipeline.context import CurationContext

    yaml_override = None
    if args.yaml:
        yaml_path = Path(args.yaml)
        if not yaml_path.exists():
            print(f"[error] YAML file not found: {yaml_path}", file=sys.stderr)
            sys.exit(1)
        with open(yaml_path) as f:
            yaml_override = yaml.safe_load(f)

    return CurationContext.from_ticket(
        args.ticket_id,
        user_config,
        yaml_override=yaml_override,
        print_only=args.print_only,
    )


def cmd_summary(args: argparse.Namespace, user_config: dict) -> None:
    from curation_pipeline.steps.pre_curation import print_curation_summary

    ctx = _build(args, user_config)
    print_curation_summary(ctx)


def cmd_pre(args: argparse.Namespace, user_config: dict) -> None:
    from curation_pipeline.steps.pre_curation import (
        copy_pretext_maps,
        print_curation_summary,
        setup_curation,
    )

    ctx = _build(args, user_config)
    print_curation_summary(ctx)
    setup_curation(ctx)
    copy_pretext_maps(ctx)


def cmd_post(args: argparse.Namespace, user_config: dict) -> None:
    ctx = _build(args, user_config)

    # Post-curation steps will be added in task 7
    from curation_pipeline.output import console
    console.print("[yellow]Post-curation steps not yet implemented (task 7).[/yellow]")
    _ = ctx


def cmd_optional(args: argparse.Namespace, user_config: dict) -> None:
    from curation_pipeline.steps.gap_track import add_gap_track
    from curation_pipeline.steps.telo_track import add_telo_track
    from curation_pipeline.steps.bedgraph_track import add_bedgraph_track

    ctx = _build(args, user_config)
    track = args.track

    if track == "gap":
        add_gap_track(ctx)
    elif track == "telo":
        add_telo_track(ctx)
    elif track == "bedgraph":
        if not args.file:
            print(
                "[error] --file / FILE is required for the 'bedgraph' track.",
                file=sys.stderr,
            )
            sys.exit(1)
        add_bedgraph_track(ctx, args.file)
    else:
        print(f"[error] Unknown track: {track}", file=sys.stderr)
        sys.exit(1)


def cmd_sex_matcher(args: argparse.Namespace, user_config: dict) -> None:
    from curation_pipeline.steps.sex_matcher import run_sex_matcher

    ctx = _build(args, user_config)
    run_sex_matcher(ctx)


def cmd_microchromosome(args: argparse.Namespace, user_config: dict) -> None:
    from curation_pipeline.steps.microchromosome import run_microchromosome_curation

    ctx = _build(args, user_config)
    run_microchromosome_curation(ctx)


def cmd_microchromosome_post(args: argparse.Namespace, user_config: dict) -> None:
    from curation_pipeline.steps.microchromosome import run_microchromosome_post_curation

    ctx = _build(args, user_config)
    run_microchromosome_post_curation(ctx)


def cmd_find_reference(args: argparse.Namespace, user_config: dict) -> None:
    from curation_pipeline.steps.find_reference import find_closest_reference

    ctx = _build(args, user_config)
    find_closest_reference(ctx, number=args.number)


def cmd_fastga(args: argparse.Namespace, user_config: dict) -> None:
    from curation_pipeline.steps.fastga import run_fastga

    ctx = _build(args, user_config)
    run_fastga(ctx, reference_path=None)  # reference_path not used, we find it


def cmd_blast_contaminants(args: argparse.Namespace, user_config: dict) -> None:
    from curation_pipeline.steps.blast_contaminants import run_blast_contaminants

    ctx = _build(args, user_config)
    run_blast_contaminants(ctx)


def cmd_busco_synteny(args: argparse.Namespace, user_config: dict) -> None:
    from curation_pipeline.output import console
    console.print("[yellow]busco-synteny not yet implemented.[/yellow]")


def cmd_busco_curated(args: argparse.Namespace, user_config: dict) -> None:
    from curation_pipeline.output import console
    console.print("[yellow]busco-curated not yet implemented.[/yellow]")


def cmd_rename_and_orient(args: argparse.Namespace, user_config: dict) -> None:
    from curation_pipeline.steps.rename_and_orient import run_rename_and_orient

    ctx = _build(args, user_config)
    run_rename_and_orient(ctx)


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="curate",
        description="Genome curation pipeline CLI",
    )
    parser.add_argument("ticket_id", help="Jira ticket ID, e.g. RC-1234 or GRIT-567")
    parser.add_argument(
        "command",
        choices=["summary", "pre", "post", "optional", "sex-matcher", "microchromosome", "microchromosome-post", "find-reference", "fastga", "blast-contaminants", "busco-synteny", "busco-curated", "rename-and-orient"],
        help=(
            "summary — print ticket summary; "
            "pre — run pre-curation steps; "
            "post — run post-curation steps; "
            "optional — add optional pretext tracks (gap, telo, bedgraph); "
            "sex-matcher — run sex_matcher and print top contigs with genes of interest; "
            "microchromosome — run second-shot microchromosome pre-curation (split + HiC remapping); "
            "microchromosome-post — run post-curation after micro map is curated (pretext-to-asm + merge); "
            "find-reference — find and download the closest reference genome from NCBI; "
            "fastga — run FastGA dot-plot comparison against closest reference; "
            "blast-contaminants — run blast contaminants search in shrapnel scaffolds; "
            "busco-synteny — run BUSCO synteny analysis against reference; "
            "busco-curated — run BUSCO analysis on curated genome; "
            "rename-and-orient — rename and orient chromosomes to match reference"
        ),
    )
    parser.add_argument(
        "track",
        nargs="?",
        choices=["gap", "telo", "bedgraph"],
        help="Track type for the 'optional' command: gap | telo | bedgraph",
    )
    parser.add_argument(
        "file",
        nargs="?",
        metavar="FILE",
        help="Bedgraph file path (required for 'optional bedgraph FILE')",
    )
    parser.add_argument(
        "--number",
        type=int,
        metavar="N",
        default=1,
        help="Number of references to download (for 'find-reference', default: 1)",
    )
    parser.add_argument(
        "--lineage",
        metavar="LINEAGE",
        help="BUSCO lineage dataset (for 'busco-synteny' and 'busco-curated', e.g., vertebrata_odb10)",
    )
    parser.add_argument(
        "--yaml",
        metavar="PATH",
        help="Path to a local YAML file instead of fetching from Jira (useful for testing)",
    )
    parser.add_argument(
        "--print-only",
        action="store_true",
        default=False,
        help="Print commands without executing them",
    )

    args = parser.parse_args()

    if args.command == "optional" and not args.track:
        parser.error("'optional' command requires a track: gap | telo | bedgraph [FILE]")

    user_config = _load_user_config()

    dispatch = {
        "summary": cmd_summary,
        "pre": cmd_pre,
        "post": cmd_post,
        "optional": cmd_optional,
        "sex-matcher": cmd_sex_matcher,
        "microchromosome": cmd_microchromosome,
        "microchromosome-post": cmd_microchromosome_post,
        "find-reference": cmd_find_reference,
        "fastga": cmd_fastga,
        "blast-contaminants": cmd_blast_contaminants,
        "busco-synteny": cmd_busco_synteny,
        "busco-curated": cmd_busco_curated,
        "rename-and-orient": cmd_rename_and_orient,
    }
    dispatch[args.command](args, user_config)


if __name__ == "__main__":
    main()
