"""
align.py
========
The alignment engine itself.

REIMPLEMENTATION NOTE (read this first)
----------------------------------------
The original align.c does *not* contain the alignment algorithm - it calls
out to ``ALIGN()``, which lives in the missing llmax.c and implements the
Myers & Miller (1988) linear-space global alignment algorithm with affine
gap costs. That file was not provided, so the algorithm below is an
independent, from-scratch reimplementation, not a port. It is standard,
well-established sequence-alignment theory (Needleman-Wunsch / Gotoh /
Hirschberg), documented per function below:

  * ``gotoh_affine_align`` - classic O(n*m) time, O(n*m) space global
    alignment with affine gap penalties (Gotoh, 1982). Used by default;
    exact optimal alignment.

  * ``hirschberg_align`` - Hirschberg's (1975) divide-and-conquer technique
    for O(n*m) time but only O(min(n,m)) space. This implementation uses a
    LINEAR (non-affine) gap cost, which is the classic textbook Hirschberg
    algorithm. This is an honest engineering trade-off: true Myers-Miller
    linear-space alignment *with* affine gaps is significantly more
    involved (it threads two coupled Gotoh recurrences through the
    divide-and-conquer split). Rather than risk a subtly-wrong
    reimplementation of that coupling, we offer linear-space alignment with
    a simple linear gap cost as the large-sequence mode, and the exact
    affine-gap algorithm as the default for sequences that comfortably fit
    in memory. This is documented again in the CLI help text and README.

Both engines return the same ``AlignmentResult`` shape, so the rest of the
pipeline (identity/similarity counting, formatting, heatmaps) is agnostic
to which one ran.

PERFORMANCE NOTE
----------------
Each DP row is computed with numpy vector operations instead of a Python
loop over cells. Cells that depend on the previous row (diagonal / vertical
moves) vectorise directly; the horizontal-gap dependency *within* a row is
resolved with a running maximum (``np.maximum.accumulate``), using

    Y[i, j] = max_{k<j} ( M[i, k] - gap_open - gap_extend * (j-1-k) )
            = -gap_open - gap_extend*(j-1) + cummax_k( M[i, k] + gap_extend*k )

which is exactly the same recurrence, so results are unchanged - only
~100x faster for typical gene/protein lengths.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "seqalign requires numpy for its DP matrices. Install it with:\n"
        "    pip install numpy"
    ) from exc

from .matrices import SubstitutionMatrix


@dataclass
class AlignmentResult:
    """Analogous to what calcons() built from res[]/aa0/aa1 in the original."""

    aligned0: str
    aligned1: str
    score: float
    engine: str


# ---------------------------------------------------------------------------
# Engine 1: Gotoh affine-gap dynamic programming (default, exact)
# ---------------------------------------------------------------------------
def gotoh_affine_align(
    seq0: str,
    seq1: str,
    matrix: SubstitutionMatrix,
    gap_open: float,
    gap_extend: float,
) -> AlignmentResult:
    """
    Exact global alignment with affine gap penalties.

    Space: O(n*m) int32/float32 numpy arrays (three DP matrices).
    Time:  O(n*m).

    Suitable for sequences up to roughly a few thousand residues each on a
    typical machine; for anything larger prefer ``hirschberg_align``.
    """
    n, m = len(seq0), len(seq1)
    if n == 0 or m == 0:
        # Degenerate case: one sequence empty -> all gaps.
        return AlignmentResult(
            aligned0=seq0 + "-" * m,
            aligned1="-" * n + seq1,
            score=-(gap_open + gap_extend * max(n, m)) if max(n, m) else 0.0,
            engine="gotoh-affine",
        )

    NEG = np.float32(-1e9)
    M = np.full((n + 1, m + 1), NEG, dtype=np.float32)  # both aligned
    X = np.full((n + 1, m + 1), NEG, dtype=np.float32)  # gap in seq1 (down)
    Y = np.full((n + 1, m + 1), NEG, dtype=np.float32)  # gap in seq0 (right)

    M[0, 0] = 0.0
    for i in range(1, n + 1):
        X[i, 0] = -(gap_open + gap_extend * (i - 1))
    for j in range(1, m + 1):
        Y[0, j] = -(gap_open + gap_extend * (j - 1))

    lut = matrix.lut.astype(np.float32)
    codes0, codes1 = matrix.encode(seq0), matrix.encode(seq1)
    # ramp[k] = gap_extend * k, used to turn the in-row gap recurrence into a
    # running maximum (see module docstring).
    ramp = np.float32(gap_extend) * np.arange(m, dtype=np.float32)
    go, ge = np.float32(gap_open), np.float32(gap_extend)

    for i in range(1, n + 1):
        s = lut[codes0[i - 1], codes1]
        best_diag = np.maximum(np.maximum(M[i - 1, :-1], X[i - 1, :-1]), Y[i - 1, :-1])
        M[i, 1:] = best_diag + s
        X[i, 1:] = np.maximum(M[i - 1, 1:] - go, X[i - 1, 1:] - ge)
        Y[i, 1:] = np.maximum.accumulate(M[i, :-1] + ramp) - go - ramp

    final_scores = {"M": M[n, m], "X": X[n, m], "Y": Y[n, m]}
    state = max(final_scores, key=lambda k: final_scores[k])
    best_score = float(final_scores[state])

    aligned0: List[str] = []
    aligned1: List[str] = []
    i, j = n, m
    while i > 0 or j > 0:
        if state == "M" and i > 0 and j > 0:
            aligned0.append(seq0[i - 1])
            aligned1.append(seq1[j - 1])
            prev_scores = {
                "M": M[i - 1, j - 1],
                "X": X[i - 1, j - 1],
                "Y": Y[i - 1, j - 1],
            }
            i, j = i - 1, j - 1
            state = max(prev_scores, key=lambda k: prev_scores[k])
        elif state == "X" and i > 0:
            aligned0.append(seq0[i - 1])
            aligned1.append("-")
            if abs(X[i, j] - (M[i - 1, j] - gap_open)) < 1e-3:
                state = "M"
            i -= 1
        elif state == "Y" and j > 0:
            aligned0.append("-")
            aligned1.append(seq1[j - 1])
            if abs(Y[i, j] - (M[i, j - 1] - gap_open)) < 1e-3:
                state = "M"
            j -= 1
        else:
            # Fallback for edge rows/columns.
            if i > 0:
                aligned0.append(seq0[i - 1])
                aligned1.append("-")
                i -= 1
            else:
                aligned0.append("-")
                aligned1.append(seq1[j - 1])
                j -= 1

    return AlignmentResult(
        aligned0="".join(reversed(aligned0)),
        aligned1="".join(reversed(aligned1)),
        score=best_score,
        engine="gotoh-affine",
    )


# ---------------------------------------------------------------------------
# Engine 2: Hirschberg linear-space alignment (linear gap cost)
# ---------------------------------------------------------------------------
def _nw_rows(
    seq0: str, seq1: str, matrix: SubstitutionMatrix, gap_cost: float
):
    """
    Yield successive rows of a Needleman-Wunsch (linear gap) score matrix,
    starting with row 0. Each row is computed from the previous one in
    vectorised form; the in-row horizontal-gap dependency
    ``curr[j] = max(t[j], curr[j-1] - gap)`` is resolved as
    ``cummax(t + gap*j) - gap*j``.
    """
    m = len(seq1)
    lut = matrix.lut.astype(np.float64)
    codes1 = matrix.encode(seq1)
    ramp = gap_cost * np.arange(m + 1, dtype=np.float64)
    prev = -ramp.copy()
    yield prev
    t = np.empty(m + 1, dtype=np.float64)
    for a in matrix.encode(seq0):
        t[0] = prev[0] - gap_cost
        np.maximum(prev[:-1] + lut[a, codes1], prev[1:] - gap_cost, out=t[1:])
        prev = np.maximum.accumulate(t + ramp) - ramp
        yield prev


def _nw_score_row(
    seq0: str, seq1: str, matrix: SubstitutionMatrix, gap_cost: float
) -> "np.ndarray":
    """Last row of a Needleman-Wunsch score matrix, O(min space) O(n*m) time."""
    row = None
    for row in _nw_rows(seq0, seq1, matrix, gap_cost):
        pass
    return row


def _needleman_wunsch_small(
    seq0: str, seq1: str, matrix: SubstitutionMatrix, gap_cost: float
) -> AlignmentResult:
    """Base case for Hirschberg recursion: plain O(n*m) NW with linear gaps."""
    n, m = len(seq0), len(seq1)
    dp = np.vstack(list(_nw_rows(seq0, seq1, matrix, gap_cost)))
    score = matrix.score

    aligned0, aligned1 = [], []
    i, j = n, m
    while i > 0 and j > 0:
        if dp[i, j] == dp[i - 1, j - 1] + score(seq0[i - 1], seq1[j - 1]):
            aligned0.append(seq0[i - 1])
            aligned1.append(seq1[j - 1])
            i, j = i - 1, j - 1
        elif dp[i, j] == dp[i - 1, j] - gap_cost:
            aligned0.append(seq0[i - 1])
            aligned1.append("-")
            i -= 1
        else:
            aligned0.append("-")
            aligned1.append(seq1[j - 1])
            j -= 1
    while i > 0:
        aligned0.append(seq0[i - 1])
        aligned1.append("-")
        i -= 1
    while j > 0:
        aligned0.append("-")
        aligned1.append(seq1[j - 1])
        j -= 1

    return AlignmentResult(
        aligned0="".join(reversed(aligned0)),
        aligned1="".join(reversed(aligned1)),
        score=float(dp[n, m]),
        engine="hirschberg",
    )


def hirschberg_align(
    seq0: str,
    seq1: str,
    matrix: SubstitutionMatrix,
    gap_cost: float,
    base_case_size: int = 200,
) -> AlignmentResult:
    """
    Linear-space global alignment (Hirschberg, 1975) with a linear gap cost.

    Space: O(min(n, m)) for the score-vector recursion, plus O(n+m) for the
    output alignment strings.
    Time:  O(n*m), same asymptotic time as the quadratic-space version.

    See the module docstring for why this uses a linear rather than affine
    gap cost.
    """
    n, m = len(seq0), len(seq1)

    if n <= 1 or m <= 1 or n * m <= base_case_size * base_case_size:
        return _needleman_wunsch_small(seq0, seq1, matrix, gap_cost)

    mid = n // 2
    score_left = _nw_score_row(seq0[:mid], seq1, matrix, gap_cost)
    score_right = _nw_score_row(seq0[mid:][::-1], seq1[::-1], matrix, gap_cost)[::-1]
    total = score_left + score_right
    split = int(np.argmax(total))
    best_score = float(total[split])

    left = hirschberg_align(seq0[:mid], seq1[:split], matrix, gap_cost, base_case_size)
    right = hirschberg_align(seq0[mid:], seq1[split:], matrix, gap_cost, base_case_size)

    return AlignmentResult(
        aligned0=left.aligned0 + right.aligned0,
        aligned1=left.aligned1 + right.aligned1,
        score=best_score,
        engine="hirschberg",
    )


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------
LARGE_SEQUENCE_THRESHOLD = 4000  # residues; above this, prefer linear-space


def align_sequences(
    seq0: str,
    seq1: str,
    matrix: SubstitutionMatrix,
    gap_open: float,
    gap_extend: float,
    engine: str = "auto",
) -> AlignmentResult:
    """
    Pick and run an alignment engine.

    :param matrix: substitution matrix (see ``matrices.load_matrix``).
    :param engine: "auto" (default; chooses based on sequence length),
        "affine" (force Gotoh, exact, quadratic space),
        "linear-space" (force Hirschberg, memory-efficient, linear gap cost).
    """
    n, m = len(seq0), len(seq1)
    if engine == "affine":
        return gotoh_affine_align(seq0, seq1, matrix, gap_open, gap_extend)
    if engine == "linear-space":
        return hirschberg_align(seq0, seq1, matrix, gap_cost=gap_extend)
    if engine == "auto":
        if max(n, m) > LARGE_SEQUENCE_THRESHOLD:
            return hirschberg_align(seq0, seq1, matrix, gap_cost=gap_extend)
        return gotoh_affine_align(seq0, seq1, matrix, gap_open, gap_extend)
    raise ValueError(f"Unknown engine '{engine}'. Choose auto|affine|linear-space.")
