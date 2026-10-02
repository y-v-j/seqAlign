#!/usr/bin/env python3
"""
run_gapdh_comparison.py
=======================
Runnable example: all-vs-all comparison of the oxidoreductase GAPDH
(glyceraldehyde-3-phosphate dehydrogenase, EC 1.2.1.12) across 15
organisms, for both the protein and the coding-DNA sequences.

It uses the FASTA files in ``tests/data/oxidoreductases/`` (downloaded by
``download_sequences.py`` in this folder; pass ``--download`` to fetch them
again) and, for each sequence type, writes:

    oxidoreductase_protein_heatmap.png   + results_protein.txt
    oxidoreductase_dna_heatmap.png       + results_dna.txt

Usage (from the repository root; seqAlign does not need to be installed):

    python tests/oxidoreductases/run_gapdh_comparison.py
    python tests/oxidoreductases/run_gapdh_comparison.py --outdir my_results --protein-theme classic
    python tests/oxidoreductases/run_gapdh_comparison.py --download   # re-download the sequences first

By default the heatmaps go to ``docs/images/`` (regenerating the images
shown in the README) and the text results to ``tests/data/oxidoreductases/``;
with ``--outdir DIR`` everything is written to DIR.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import os
import sys
from typing import List, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA_DIR = os.path.join(REPO_ROOT, "tests", "data", "oxidoreductases")
DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "docs", "images")

# Make the in-repo package importable without installing it.
sys.path.insert(0, REPO_ROOT)

from seqalign import cli  # noqa: E402
from seqalign.heatmap import THEMES  # noqa: E402


def _organisms(data_dir: str) -> List[str]:
    """Organism labels in manifest order (mammals first, bacteria last)."""
    with open(os.path.join(data_dir, "manifest.tsv"), encoding="utf-8") as fh:
        return [row["label"] for row in csv.DictReader(fh, delimiter="\t")]


def _run(kind: str, organisms: List[str], data_dir: str, image_dir: str, results_dir: str,
         matrix: str, theme: str, dpi: int, verbose: bool) -> int:
    files = [os.path.join(data_dir, kind, f"{org}.fasta") for org in organisms]
    heatmap = os.path.join(image_dir, f"oxidoreductase_{kind}_heatmap.png")
    results = os.path.join(results_dir, f"results_{kind}.txt")
    label = "protein" if kind == "protein" else "coding DNA"
    argv = files + [
        "-i", "-d",
        "-s", matrix,
        "--theme", theme,
        "--dpi", str(dpi),
        "--title", f"GAPDH (EC 1.2.1.12) {label} - {len(files)} organisms",
        "--heatmap-output", heatmap,
        "-O", results,
    ]
    print(f"[{kind}] aligning {len(files)} sequences "
          f"({len(files) * (len(files) - 1) // 2} pairs, matrix {matrix}, theme {theme}) ...")
    capture = io.StringIO()
    with contextlib.redirect_stdout(sys.stdout if verbose else capture):
        status = cli.main(argv)
    if status == 0:
        print(f"[{kind}] heatmap: {os.path.relpath(heatmap)}")
        print(f"[{kind}] results: {os.path.relpath(results)}")
    return status


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default=None,
                        help="write heatmaps and results here (default: heatmaps to docs/images, "
                             "results to tests/data/oxidoreductases)")
    parser.add_argument("--data-dir", default=DATA_DIR, help="downloaded sequences (default: tests/data/oxidoreductases)")
    parser.add_argument("--download", action="store_true", help="(re-)download the sequences before comparing")
    parser.add_argument("--only", choices=["protein", "dna"], help="run only one sequence type")
    parser.add_argument("--protein-matrix", default="BLOSUM50", help="protein scoring matrix (default: BLOSUM50)")
    parser.add_argument("--dna-matrix", default="NUC.4.4", help="DNA scoring matrix (default: NUC.4.4)")
    parser.add_argument("--protein-theme", default="sunset-ocean", choices=list(THEMES), help="protein heatmap theme")
    parser.add_argument("--dna-theme", default="forest-berry", choices=list(THEMES), help="DNA heatmap theme")
    parser.add_argument("--dpi", type=int, default=150, help="heatmap resolution (default: 150)")
    parser.add_argument("-v", "--verbose", action="store_true", help="also print every pairwise result")
    args = parser.parse_args(argv)

    if args.download or not os.path.isfile(os.path.join(args.data_dir, "manifest.tsv")):
        from tests.oxidoreductases import download_sequences

        print("Downloading GAPDH sequences (UniProt / ENA) ...")
        if download_sequences.main(["--outdir", args.data_dir]) != 0:
            print("error: download failed for some organisms", file=sys.stderr)
            return 1

    image_dir = args.outdir or DEFAULT_OUTDIR
    results_dir = args.outdir or args.data_dir
    os.makedirs(image_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    organisms = _organisms(args.data_dir)
    status = 0
    if args.only in (None, "protein"):
        status |= _run("protein", organisms, args.data_dir, image_dir, results_dir,
                       args.protein_matrix, args.protein_theme, args.dpi, args.verbose)
    if args.only in (None, "dna"):
        status |= _run("dna", organisms, args.data_dir, image_dir, results_dir,
                       args.dna_matrix, args.dna_theme, args.dpi, args.verbose)
    return status


if __name__ == "__main__":
    sys.exit(main())
