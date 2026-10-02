"""
cli.py
======
The command-line driver. This is the direct analogue of align.c's
``main()``/``initenv()``: parse options, read sequence files, run every
pairwise comparison (the original's nested `for(k...) for(kj=k+1...)`
all-vs-all loop, driven there by `-a`), print/save results, and — new in
this port — draw a heatmap when more than two sequences are given.

Flag compatibility with the original
-------------------------------------
    -f GAP_OPEN     gap opening penalty              (was gdelval)
    -g GAP_EXTEND    gap extension penalty            (was ggapval)
    -m MARKX         output format selector (0/3 supported; others fall
                     back to 0 with a warning — the original's markx
                     values 1/2/4/9/10 depended on library-search-specific
                     formatting this driver's all-vs-all mode doesn't use)
    -n               treat input as DNA/RNA           (was dnaseq flag)
    -O FILE          write results to FILE as well as stdout (was resfile)
    -s NAME|FILE     scoring matrix name or NCBI-format matrix file
                     (BLOSUM50 default for protein, simple +5/-4 "DNA"
                     scheme default for nucleotides; see --list-matrices)
    -w WIDTH         alignment display line width      (was llen)
    -x OFF0 OFF1     numbering offsets for the two sequences (was sq0off/
                     sq1off)
    -d               suppress the aligned-sequence display (was showAlignment=0)
    -i               print percent similarity as well as identity
    -a N             kept for backward compatibility; with modern argparse
                     you can just list N files directly, so -a is optional
                     and, if given, is only validated against the actual
                     file count rather than driving a manual off-by-one loop.

New flags (not in the original)
--------------------------------
    --engine {auto,affine,linear-space}   choose the alignment algorithm
    --dna / --protein                     force sequence type detection
    --heatmap-output PATH                 where to save the multi-sequence
                                           heatmap (default: heatmap.png)
    --no-heatmap                          skip heatmap generation
    --theme NAME                          heatmap color theme (--list-themes)
    --upper-cmap / --lower-cmap CMAP      override a theme's color maps
    --label-rotation DEG                  angle of the heatmap's x labels
    --title TEXT                          heatmap title
    --list-matrices / --list-themes       print the available choices and exit
"""

from __future__ import annotations

import argparse
import itertools
import sys
from dataclasses import dataclass
from typing import List, Optional

from .align import align_sequences
from .consensus import format_alignment_blocks, score_alignment
from . import __version__
from .heatmap import (
    DEFAULT_THEME,
    DEFAULT_TITLE,
    THEMES,
    build_pairwise_matrices,
    describe_themes,
    get_theme,
    plot_heatmap,
)
from .io_utils import SequenceFileError, SequenceRecord, read_sequence
from .matrices import ScoringScheme, describe_matrices, get_scoring_scheme

VERSTR = f"seqAlign {__version__} (Python reimplementation of align.c's driver logic)"
REFSTR = "Original algorithm reference: Myers and Miller, CABIOS (1988)."


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="seqAlign",
        description=(
            "Pairwise / all-vs-all global sequence alignment "
            "(identity + similarity), with an optional heatmap for 3+ sequences."
        ),
        epilog=REFSTR,
    )
    parser.add_argument("files", nargs="*", help="FASTA file(s); 2+ for pairwise, "
                                                   "3+ triggers all-vs-all + heatmap")
    parser.add_argument("-f", "--gap-open", type=float, default=None, metavar="PENALTY",
                         help="gap opening penalty (default: 12 protein / 16 DNA)")
    parser.add_argument("-g", "--gap-extend", type=float, default=None, metavar="PENALTY",
                         help="gap extension penalty (default: 2 protein / 4 DNA)")
    parser.add_argument("-m", "--markx", type=int, default=0, choices=[0, 3],
                         help="output format (0=blocks, 3=blocks+match line); default 0")
    parser.add_argument("-n", "--dna", dest="force_dna", action="store_true",
                         help="treat sequences as DNA/RNA")
    parser.add_argument("--protein", dest="force_protein", action="store_true",
                         help="treat sequences as protein (overrides auto-detection)")
    parser.add_argument("-O", "--output-file", type=str, default=None, metavar="FILE",
                         help="also write full results to FILE")
    parser.add_argument("-s", "--matrix", type=str, default=None, metavar="NAME|FILE",
                         help="scoring matrix name or NCBI-format matrix file (default: BLOSUM50 "
                              "for protein, DNA for nucleotides; see --list-matrices)")
    parser.add_argument("--list-matrices", action="store_true",
                         help="list the bundled protein and nucleotide scoring matrices and exit")
    parser.add_argument("-w", "--line-width", type=int, default=60, metavar="N",
                         help="line width for alignment display (default: 60)")
    parser.add_argument("-x", "--offsets", type=int, nargs=2, default=(1, 1),
                         metavar=("OFF0", "OFF1"), help="numbering offsets for seq1/seq2")
    parser.add_argument("-d", "--no-display", dest="show_alignment", action="store_false",
                         help="suppress the aligned-sequence display")
    parser.add_argument("-i", "--similarity", dest="with_similarity", action="store_true",
                         help="also report percent similarity")
    parser.add_argument("-a", "--num-sequences", type=int, default=None, metavar="N",
                         help="(compat) expected number of input files; validated, not required")
    parser.add_argument("--engine", choices=["auto", "affine", "linear-space"], default="auto",
                         help="alignment algorithm (default: auto, chooses by sequence length)")
    parser.add_argument("--heatmap-output", type=str, default="heatmap.png", metavar="PATH",
                         help="output path for the all-vs-all heatmap (default: heatmap.png)")
    parser.add_argument("--no-heatmap", action="store_true", help="skip heatmap generation")
    parser.add_argument("--dpi", type=int, default=150, help="heatmap resolution (default: 150)")
    parser.add_argument("--theme", type=str.lower, default=DEFAULT_THEME, choices=list(THEMES),
                         help=f"heatmap color theme (default: {DEFAULT_THEME}; see --list-themes)")
    parser.add_argument("--list-themes", action="store_true",
                         help="list the heatmap color themes and exit")
    parser.add_argument("--upper-cmap", type=str, default=None, metavar="CMAP",
                         help="matplotlib colormap for the upper (similarity) triangle; overrides --theme")
    parser.add_argument("--lower-cmap", type=str, default=None, metavar="CMAP",
                         help="matplotlib colormap for the lower (identity) triangle; overrides --theme")
    parser.add_argument("--label-rotation", type=float, default=None, metavar="DEG",
                         help="rotation of the heatmap's x-axis labels in degrees "
                              "(default: auto - 45 when names are wider than their column)")
    parser.add_argument("--title", type=str, default=DEFAULT_TITLE, help="heatmap title")
    return parser


@dataclass
class PairResult:
    name0: str
    name1: str
    identity: float
    similarity: float


def _load_sequences(paths: List[str], force_dna: Optional[bool]) -> List[SequenceRecord]:
    records = []
    for path in paths:
        try:
            records.append(read_sequence(path, force_dna=force_dna))
        except SequenceFileError as exc:
            print(f"error: {exc}", file=sys.stderr)
            sys.exit(1)
    return records


def _run_pair(rec0: SequenceRecord, rec1: SequenceRecord, scheme: ScoringScheme,
              args: argparse.Namespace, out_handle) -> PairResult:
    gap_open = args.gap_open if args.gap_open is not None else scheme.default_gaps.gap_open
    gap_extend = args.gap_extend if args.gap_extend is not None else scheme.default_gaps.gap_extend

    result = align_sequences(
        rec0.sequence, rec1.sequence, scheme.matrix, gap_open, gap_extend, engine=args.engine
    )
    stats = score_alignment(result.aligned0, result.aligned1, scheme.score_fn)

    header = (
        f"{rec0.title:<40} {rec0.length:>6} {'DNA' if rec0.is_dna else 'aa'} vs.\n"
        f"{rec1.title:<40} {rec1.length:>6} {'DNA' if rec1.is_dna else 'aa'}\n"
        f"scoring: {scheme.name}, gap penalties: {gap_open}/{gap_extend}, engine: {result.engine}\n"
    )
    summary = f"{stats.percent_identity:4.1f}% identity"
    if args.with_similarity:
        summary += f"; {stats.percent_similarity:4.1f}% similarity"
    summary += f"\t'GA' score: {result.score:.1f}\n"

    lines = [header, summary]
    if args.show_alignment:
        lines.append(
            format_alignment_blocks(
                result.aligned0, result.aligned1, rec0.identifier, rec1.identifier,
                scheme.score_fn, line_width=args.line_width,
                offset0=args.offsets[0], offset1=args.offsets[1],
            )
        )

    text = "\n".join(lines)
    print(text)
    if out_handle:
        out_handle.write(text + "\n")

    return PairResult(rec0.identifier, rec1.identifier, stats.percent_identity, stats.percent_similarity)


def _check_colormaps(theme) -> None:
    """Fail before aligning anything if a colormap name is invalid."""
    try:
        from matplotlib import colormaps
    except ImportError:
        return  # reported later as a warning when the heatmap is drawn
    for cmap in (theme.upper_cmap, theme.lower_cmap):
        if cmap not in colormaps:
            raise ValueError(f"Unknown matplotlib colormap '{cmap}'.")


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.list_matrices or args.list_themes:
        if args.list_matrices:
            print(describe_matrices())
        if args.list_themes:
            print(describe_themes())
        return 0
    if not args.files:
        parser.error("at least two sequence files are required")

    print(f"{VERSTR}\n{REFSTR}\n")

    if args.force_dna and args.force_protein:
        print("error: --dna and --protein are mutually exclusive", file=sys.stderr)
        return 1
    force_dna = True if args.force_dna else (False if args.force_protein else None)

    if args.num_sequences is not None and args.num_sequences != len(args.files):
        print(
            f"warning: -a {args.num_sequences} does not match the {len(args.files)} "
            f"file(s) actually given; proceeding with the files provided.",
            file=sys.stderr,
        )

    records = _load_sequences(args.files, force_dna)
    if len(records) < 2:
        print("error: at least two sequence files are required.", file=sys.stderr)
        return 1

    kinds = {r.is_dna for r in records}
    if len(kinds) > 1:
        dna = [r.source_path for r in records if r.is_dna]
        prot = [r.source_path for r in records if not r.is_dna]
        print(
            "error: a mix of nucleotide and protein sequences was detected "
            f"(DNA: {', '.join(dna)}; protein: {', '.join(prot)}). "
            "Pass --dna or --protein if the auto-detection is wrong.",
            file=sys.stderr,
        )
        return 1
    try:
        scheme = get_scoring_scheme(args.matrix, dnaseq=records[0].is_dna)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    make_heatmap = len(records) > 2 and not args.no_heatmap
    if make_heatmap:
        try:
            theme = get_theme(args.theme, args.upper_cmap, args.lower_cmap)
            _check_colormaps(theme)
        except (ValueError, ImportError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1

    out_handle = None
    if args.output_file:
        try:
            out_handle = open(args.output_file, "w", encoding="utf-8")
        except OSError as exc:
            print(f"error: cannot open -O output file: {exc}", file=sys.stderr)
            return 1

    n = len(records)
    identity_matrix = [[0.0] * n for _ in range(n)]
    similarity_matrix = [[0.0] * n for _ in range(n)]

    try:
        for i, j in itertools.combinations(range(n), 2):
            pair = _run_pair(records[i], records[j], scheme, args, out_handle)
            identity_matrix[i][j] = identity_matrix[j][i] = pair.identity
            similarity_matrix[i][j] = similarity_matrix[j][i] = pair.similarity
    finally:
        if out_handle:
            out_handle.close()

    if make_heatmap:
        labels = [r.identifier for r in records]
        combined = build_pairwise_matrices(labels, identity_matrix, similarity_matrix)
        try:
            path = plot_heatmap(
                labels, combined, args.heatmap_output, dpi=args.dpi, title=args.title,
                theme=theme, label_rotation=args.label_rotation,
            )
            print(f"Heatmap written to: {path}")
        except ImportError as exc:
            print(f"warning: {exc}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
