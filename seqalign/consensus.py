"""
consensus.py
============
Direct analogues of the original's ``calcons()`` and ``discons()``.

``calcons()`` in align.c walked the raw edit-script (``res[]``) alongside
aa0/aa1, built the gapped consensus strings ``seqc0``/``seqc1``, and counted
``nd`` (identical columns) out of ``nc`` (total aligned columns) to compute
percent identity. Here the aligned strings already come back directly from
the aligner, so ``score_alignment()`` plays the same role: count identical
columns for percent identity, and (since the original's
``CHECK_SIMILARITY`` was never provided) count positive-scoring columns for
percent similarity, using the standard bioinformatics convention (as seen
in EMBOSS ``needle``/``water``): a column is "similar" if it's an exact
match OR the substitution matrix gives it a positive score.

``discons()`` printed the aligned sequences in fixed-width blocks;
``format_alignment_blocks()`` reproduces that, plus a EMBOSS-style match
line (``|`` identical, ``:`` similar, blank otherwise) for readability.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

ScoreFn = Callable[[str, str], int]


@dataclass
class AlignmentStats:
    aligned_length: int
    identical: int
    similar: int
    percent_identity: float
    percent_similarity: float


def score_alignment(aligned0: str, aligned1: str, score_fn: ScoreFn) -> AlignmentStats:
    """Equivalent of calcons()'s nd/nc bookkeeping, plus a similarity count."""
    if len(aligned0) != len(aligned1):
        raise ValueError("Aligned sequences must be the same length.")

    nc = len(aligned0)
    identical = 0
    similar = 0
    for a, b in zip(aligned0, aligned1):
        if a == "-" or b == "-":
            continue
        if a == b:
            identical += 1
            similar += 1
        elif score_fn(a, b) > 0:
            similar += 1

    # Denominator matches the classic convention: aligned columns excluding
    # gaps counted against the shorter/ungapped comparison length, expressed
    # as a percentage of alignment length (mirrors `percent = nd*100/nc` in
    # the original, generalised with a similarity numerator too).
    percent_identity = (identical * 100.0 / nc) if nc else 0.0
    percent_similarity = (similar * 100.0 / nc) if nc else 0.0

    return AlignmentStats(
        aligned_length=nc,
        identical=identical,
        similar=similar,
        percent_identity=percent_identity,
        percent_similarity=percent_similarity,
    )


def format_alignment_blocks(
    aligned0: str,
    aligned1: str,
    name0: str,
    name1: str,
    score_fn: ScoreFn,
    line_width: int = 60,
    offset0: int = 1,
    offset1: int = 1,
) -> str:
    """Block-formatted pairwise alignment display (discons() equivalent)."""
    if line_width <= 0:
        line_width = 60
    lines = []
    n = len(aligned0)
    pos0, pos1 = offset0, offset1
    label_w = max(len(name0), len(name1), 1)

    for start in range(0, n, line_width):
        seg0 = aligned0[start : start + line_width]
        seg1 = aligned1[start : start + line_width]
        match_line = "".join(
            "|" if a == b and a != "-" else (":" if a != "-" and b != "-" and score_fn(a, b) > 0 else " ")
            for a, b in zip(seg0, seg1)
        )
        end0 = pos0 + sum(1 for c in seg0 if c != "-") - 1
        end1 = pos1 + sum(1 for c in seg1 if c != "-") - 1

        lines.append(f"{name0:<{label_w}} {pos0:>6} {seg0} {end0}")
        lines.append(f"{'':<{label_w}} {'':>6} {match_line}")
        lines.append(f"{name1:<{label_w}} {pos1:>6} {seg1} {end1}")
        lines.append("")

        pos0 = end0 + 1
        pos1 = end1 + 1

    return "\n".join(lines).rstrip() + "\n"
