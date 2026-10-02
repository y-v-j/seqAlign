import random

import pytest

from seqalign.align import align_sequences, gotoh_affine_align, hirschberg_align
from seqalign.consensus import score_alignment
from seqalign.matrices import load_matrix

PROTEIN = "ACDEFGHIKLMNPQRSTVWY"


def _reference_affine_score(s0, s1, score, go, ge):
    """Plain cell-by-cell Gotoh recurrence - the reference the vectorised engine must match."""
    neg = float("-inf")
    n, m = len(s0), len(s1)
    M = [[neg] * (m + 1) for _ in range(n + 1)]
    X = [[neg] * (m + 1) for _ in range(n + 1)]
    Y = [[neg] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0
    for i in range(1, n + 1):
        X[i][0] = -(go + ge * (i - 1))
    for j in range(1, m + 1):
        Y[0][j] = -(go + ge * (j - 1))
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            M[i][j] = max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]) + score(s0[i - 1], s1[j - 1])
            X[i][j] = max(M[i - 1][j] - go, X[i - 1][j] - ge)
            Y[i][j] = max(M[i][j - 1] - go, Y[i][j - 1] - ge)
    return max(M[n][m], X[n][m], Y[n][m])


def _reference_linear_score(s0, s1, score, g):
    n, m = len(s0), len(s1)
    dp = [[-g * (i + j) if i == 0 or j == 0 else 0 for j in range(m + 1)] for i in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            dp[i][j] = max(dp[i - 1][j - 1] + score(s0[i - 1], s1[j - 1]), dp[i - 1][j] - g, dp[i][j - 1] - g)
    return dp[n][m]


def _alignment_score_affine(a0, a1, score, go, ge):
    total, in_gap = 0, None
    for x, y in zip(a0, a1):
        if x == "-" or y == "-":
            kind = 0 if x == "-" else 1
            total -= ge if in_gap == kind else go
            in_gap = kind
        else:
            total += score(x, y)
            in_gap = None
    return total


def _mutate(rng, seq, alphabet):
    return "".join(c if rng.random() > 0.25 else rng.choice(alphabet) for c in seq if rng.random() > 0.1)


@pytest.mark.parametrize("seed", range(40))
def test_affine_engine_matches_reference(seed):
    rng = random.Random(seed)
    dna = seed % 2 == 0
    matrix = load_matrix("DNA" if dna else "BLOSUM62")
    alphabet = "ACGT" if dna else PROTEIN
    s0 = "".join(rng.choice(alphabet) for _ in range(rng.randint(1, 50)))
    s1 = _mutate(rng, s0, alphabet) or "A"
    go, ge = rng.choice([(12, 2), (16, 4), (10, 1)])
    result = gotoh_affine_align(s0, s1, matrix, go, ge)
    assert result.aligned0.replace("-", "") == s0
    assert result.aligned1.replace("-", "") == s1
    assert result.score == pytest.approx(_reference_affine_score(s0, s1, matrix.score, go, ge))
    # The returned alignment really achieves the reported score.
    assert _alignment_score_affine(result.aligned0, result.aligned1, matrix.score, go, ge) == pytest.approx(result.score)


@pytest.mark.parametrize("seed", range(20))
def test_linear_space_engine_matches_reference(seed):
    rng = random.Random(1000 + seed)
    matrix = load_matrix("BLOSUM50")
    s0 = "".join(rng.choice(PROTEIN) for _ in range(rng.randint(20, 80)))
    s1 = _mutate(rng, s0, PROTEIN) or "A"
    result = hirschberg_align(s0, s1, matrix, gap_cost=2, base_case_size=6)
    assert result.aligned0.replace("-", "") == s0
    assert result.aligned1.replace("-", "") == s1
    assert result.score == pytest.approx(_reference_linear_score(s0, s1, matrix.score, 2))


def test_identical_sequences_are_100_percent():
    matrix = load_matrix("BLOSUM50")
    seq = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF"
    for engine in ("affine", "linear-space"):
        r = align_sequences(seq, seq, matrix, 12, 2, engine=engine)
        stats = score_alignment(r.aligned0, r.aligned1, matrix.score)
        assert stats.percent_identity == 100.0 and stats.percent_similarity == 100.0


def test_unknown_engine():
    with pytest.raises(ValueError):
        align_sequences("AC", "AC", load_matrix("DNA"), 16, 4, engine="nope")
