"""
matrices.py
===========
Scoring matrices for the alignment engine.

REIMPLEMENTATION NOTE
----------------------
The original align.c reads its substitution matrix from a compiled-in PAM
table (`upam.gbl`, populated by pam.c) or an external ``*.mat`` file via
``initpam()``/``initpam2()``. Those files were not available, so this module
ships the standard, publicly-published NCBI matrices as plain-text files in
``seqalign/data/matrices/`` and loads them on demand:

Protein
    BLOSUM45, BLOSUM50 (default), BLOSUM62, BLOSUM80, BLOSUM90,
    PAM30, PAM70, PAM250

Nucleic acid (DNA/RNA)
    DNA (built-in +5/-4 match/mismatch; default, mirrors the original's
    simple ``npam`` table set by ``resetp()``), NUC.4.4 (a.k.a. EDNAFULL,
    with IUPAC ambiguity codes), BLASTN (NCBI BLASTN +2/-3 with IUPAC
    codes), TRANS (transition/transversion matrix)

Just like the original's ``-s`` option, a path to any NCBI-format matrix
file can be given instead of a bundled name.

Every matrix is expanded into a 128x128 ASCII lookup table, so scoring a
whole row of the DP matrix is a single numpy fancy-index operation.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional

import numpy as np

_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "matrices")

PROTEIN_MATRICES = ["BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90", "PAM30", "PAM70", "PAM250"]
NUCLEOTIDE_MATRICES = ["DNA", "NUC.4.4", "BLASTN", "TRANS"]

DEFAULT_PROTEIN_MATRIX = "BLOSUM50"
DEFAULT_NUCLEOTIDE_MATRIX = "DNA"

# Aliases accepted by -s (case-insensitive).
_ALIASES = {"EDNAFULL": "NUC.4.4", "NUC44": "NUC.4.4"}

_MATRIX_DESCRIPTIONS = {
    "BLOSUM45": "distantly related proteins",
    "BLOSUM50": "general purpose protein matrix (default)",
    "BLOSUM62": "NCBI BLASTP default; moderately related proteins",
    "BLOSUM80": "closely related proteins",
    "BLOSUM90": "very closely related proteins",
    "PAM30": "short, highly similar protein sequences",
    "PAM70": "similar protein sequences",
    "PAM250": "distantly related proteins (Dayhoff)",
    "DNA": "simple match +5 / mismatch -4 (default)",
    "NUC.4.4": "EDNAFULL; +5/-4 with IUPAC ambiguity codes",
    "BLASTN": "NCBI BLASTN +2/-3 with IUPAC ambiguity codes",
    "TRANS": "transition/transversion matrix (A<->G, C<->T score higher)",
}

# Simple DNA/RNA match-mismatch scheme (mirrors the original's simple `npam`
# nucleotide table: identity is rewarded, everything else penalised equally).
DNA_MATCH = 5
DNA_MISMATCH = -4

_NUCLEOTIDE_CHARS = set("ACGTUNRYSWKMBDHV")


@dataclass(frozen=True)
class SubstitutionMatrix:
    """A substitution matrix plus a precomputed ASCII lookup table."""

    name: str
    alphabet: str
    scores: "np.ndarray"
    is_nucleotide: bool
    lut: "np.ndarray" = field(repr=False, compare=False)

    def encode(self, sequence: str) -> "np.ndarray":
        """Sequence -> array of ASCII codes, usable as indices into ``lut``."""
        return np.frombuffer(sequence.upper().encode("ascii", "replace"), dtype=np.uint8) & 0x7F

    def score(self, a: str, b: str) -> int:
        return int(self.lut[ord(a) & 0x7F, ord(b) & 0x7F])

    @property
    def min_score(self) -> int:
        return int(self.scores.min())

    @property
    def max_score(self) -> int:
        return int(self.scores.max())


def _build_lut(alphabet: str, scores: "np.ndarray", is_nucleotide: bool) -> "np.ndarray":
    """
    Expand an alphabet x alphabet matrix into a 128x128 lookup table.

    Characters missing from the matrix are scored like the matrix's
    wildcard ('N' for nucleotides, 'X' for proteins) if it has one, U is
    scored like T for RNA, and anything else gets the matrix's lowest score.
    Lower-case letters score like their upper-case counterparts.
    """
    index = {c: i for i, c in enumerate(alphabet)}
    wildcard = "N" if is_nucleotide else "X"

    def resolve(c: str) -> Optional[int]:
        c = c.upper()
        if c in index:
            return index[c]
        if is_nucleotide and c == "U" and "T" in index:
            return index["T"]
        if c.isalpha() or c == "*":
            return index.get(wildcard)
        return None

    rows: List[Optional[int]] = [resolve(chr(code)) for code in range(128)]
    lut = np.full((128, 128), scores.min(), dtype=np.int16)
    valid = np.array([r is not None for r in rows])
    mapped = np.array([r if r is not None else 0 for r in rows])
    lut[np.ix_(valid, valid)] = scores[np.ix_(mapped[valid], mapped[valid])]
    return lut


def _make_matrix(name: str, alphabet: str, scores: "np.ndarray", is_nucleotide: bool) -> SubstitutionMatrix:
    return SubstitutionMatrix(
        name=name,
        alphabet=alphabet,
        scores=scores,
        is_nucleotide=is_nucleotide,
        lut=_build_lut(alphabet, scores, is_nucleotide),
    )


def parse_matrix_file(path: str, name: Optional[str] = None) -> SubstitutionMatrix:
    """
    Parse an NCBI-format matrix file: '#' comment lines, a header line of
    column letters, then one row per letter (row letter followed by scores).
    """
    try:
        with open(path, "r", encoding="utf-8") as fh:
            lines = [ln.split() for ln in fh if ln.strip() and not ln.lstrip().startswith("#")]
    except OSError as exc:
        raise ValueError(f"Cannot read scoring matrix file '{path}': {exc}") from exc

    if not lines:
        raise ValueError(f"Scoring matrix file '{path}' is empty.")
    columns = [c.upper() for c in lines[0]]
    rows: Dict[str, List[int]] = {}
    for parts in lines[1:]:
        try:
            rows[parts[0].upper()] = [int(float(v)) for v in parts[1:]]
        except (ValueError, IndexError) as exc:
            raise ValueError(f"Malformed row in scoring matrix '{path}': {' '.join(parts)}") from exc

    if set(rows) != set(columns) or any(len(r) != len(columns) for r in rows.values()):
        raise ValueError(f"Scoring matrix '{path}' is not square / row and column labels differ.")

    alphabet = "".join(columns)
    scores = np.array([rows[c] for c in columns], dtype=np.int16)
    is_nucleotide = set(alphabet) <= _NUCLEOTIDE_CHARS
    label = name or os.path.splitext(os.path.basename(path))[0]
    return _make_matrix(label, alphabet, scores, is_nucleotide)


def _simple_dna_matrix() -> SubstitutionMatrix:
    """Identical letters (including IUPAC codes, and T/U) score +5, all else -4."""
    alphabet = "ACGTURYSWKMBDHVN"
    scores = np.full((len(alphabet), len(alphabet)), DNA_MISMATCH, dtype=np.int16)
    np.fill_diagonal(scores, DNA_MATCH)
    t, u = alphabet.index("T"), alphabet.index("U")
    scores[t, u] = scores[u, t] = DNA_MATCH
    return _make_matrix("DNA", alphabet, scores, is_nucleotide=True)


_CACHE: Dict[str, SubstitutionMatrix] = {}


def canonical_matrix_name(name: str) -> str:
    upper = name.upper()
    upper = _ALIASES.get(upper, upper)
    for known in PROTEIN_MATRICES + NUCLEOTIDE_MATRICES:
        if known.upper() == upper:
            return known
    return name


def load_matrix(name: str) -> SubstitutionMatrix:
    """Load a bundled matrix by name (case-insensitive), or any matrix file by path."""
    canonical = canonical_matrix_name(name)
    if canonical in _CACHE:
        return _CACHE[canonical]

    if canonical == "DNA":
        matrix = _simple_dna_matrix()
    elif canonical in PROTEIN_MATRICES or canonical in NUCLEOTIDE_MATRICES:
        matrix = parse_matrix_file(os.path.join(_DATA_DIR, f"{canonical}.mat"), name=canonical)
    elif os.path.isfile(name):
        matrix = parse_matrix_file(name)
    else:
        raise ValueError(
            f"Unknown scoring matrix '{name}'. Bundled matrices: "
            f"protein: {', '.join(PROTEIN_MATRICES)}; "
            f"nucleotide: {', '.join(NUCLEOTIDE_MATRICES)}. "
            f"You can also pass the path to an NCBI-format matrix file."
        )
    _CACHE[canonical] = matrix
    return matrix


def describe_matrices() -> str:
    """Human-readable list of bundled matrices (for --list-matrices)."""
    lines = ["Protein matrices:"]
    lines += [f"  {n:<10} {_MATRIX_DESCRIPTIONS[n]}" for n in PROTEIN_MATRICES]
    lines.append("Nucleic acid matrices:")
    lines += [f"  {n:<10} {_MATRIX_DESCRIPTIONS[n]}" for n in NUCLEOTIDE_MATRICES]
    lines.append("Any NCBI-format matrix file can also be passed by path, e.g. -s my_matrix.mat")
    return "\n".join(lines)


@dataclass(frozen=True)
class GapPenalties:
    """Affine gap cost: gap_open is charged once, gap_extend per residue."""

    gap_open: float
    gap_extend: float


@dataclass(frozen=True)
class ScoringScheme:
    """Bundles a substitution matrix with default gap penalties."""

    name: str
    matrix: SubstitutionMatrix
    default_gaps: GapPenalties
    is_dna: bool

    @property
    def score_fn(self) -> Callable[[str, str], int]:
        return self.matrix.score


def get_scoring_scheme(name: Optional[str], dnaseq: bool) -> ScoringScheme:
    """
    Resolve a scoring scheme by name, mirroring initenv()'s `-s` option and
    resetp()'s DNA-vs-protein switch in the original C driver.

    :param name: matrix name or file path requested via -s; ``None`` selects
                 BLOSUM50 for protein / the simple DNA scheme for nucleotides.
    :param dnaseq: True if the sequences are nucleotides (-n or auto-detected).
    :raises ValueError: for an unknown matrix, or a protein matrix used with
                        nucleotide sequences (and vice versa).
    """
    matrix = load_matrix(name or (DEFAULT_NUCLEOTIDE_MATRIX if dnaseq else DEFAULT_PROTEIN_MATRIX))
    if matrix.is_nucleotide != dnaseq:
        raise ValueError(
            f"Scoring matrix '{matrix.name}' is a "
            f"{'nucleotide' if matrix.is_nucleotide else 'protein'} matrix but the sequences are "
            f"{'DNA/RNA' if dnaseq else 'protein'}. Choose a matching matrix (see --list-matrices), "
            f"or pass --dna/--protein if sequence type auto-detection was wrong."
        )

    if dnaseq:
        label = "DNA (match=+5 / mismatch=-4)" if matrix.name == "DNA" else matrix.name
        gaps = GapPenalties(gap_open=16, gap_extend=4)
    else:
        label = matrix.name
        gaps = GapPenalties(gap_open=12, gap_extend=2)
    return ScoringScheme(name=label, matrix=matrix, default_gaps=gaps, is_dna=dnaseq)
