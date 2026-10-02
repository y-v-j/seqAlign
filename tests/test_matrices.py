import pytest

from seqalign.matrices import (
    NUCLEOTIDE_MATRICES,
    PROTEIN_MATRICES,
    get_scoring_scheme,
    load_matrix,
    parse_matrix_file,
)


@pytest.mark.parametrize("name", PROTEIN_MATRICES)
def test_protein_matrices_load_and_are_symmetric(name):
    m = load_matrix(name)
    assert not m.is_nucleotide
    assert set("ACDEFGHIKLMNPQRSTVWY") <= set(m.alphabet)
    assert (m.scores == m.scores.T).all()


@pytest.mark.parametrize("name", NUCLEOTIDE_MATRICES)
def test_nucleotide_matrices_load(name):
    m = load_matrix(name)
    assert m.is_nucleotide
    assert set("ACGT") <= set(m.alphabet)
    assert all(m.score(c, c) > 0 for c in "ACGT")


@pytest.mark.parametrize(
    "name,a,b,expected",
    [
        ("BLOSUM62", "W", "W", 11), ("BLOSUM62", "A", "A", 4), ("BLOSUM62", "E", "D", 2),
        ("BLOSUM50", "W", "W", 15), ("BLOSUM50", "C", "C", 13), ("BLOSUM50", "V", "B", -4),
        ("PAM250", "W", "W", 17), ("PAM30", "A", "R", -7),
        ("NUC.4.4", "A", "A", 5), ("NUC.4.4", "A", "W", 1), ("NUC.4.4", "A", "T", -4),
        ("BLASTN", "A", "A", 2), ("BLASTN", "A", "C", -3),
        ("TRANS", "A", "G", 4), ("TRANS", "A", "T", 0),
        ("DNA", "A", "A", 5), ("DNA", "A", "C", -4),
    ],
)
def test_published_values(name, a, b, expected):
    assert load_matrix(name).score(a, b) == expected


def test_case_rna_and_unknown_residues():
    nuc = load_matrix("NUC.4.4")
    assert nuc.score("a", "A") == 5
    assert nuc.score("U", "T") == 5  # RNA scored like DNA
    prot = load_matrix("BLOSUM62")
    assert prot.score("J", "A") == prot.score("X", "A")  # unknown -> wildcard X
    assert load_matrix("DNA").score("U", "T") == 5


def test_aliases_and_case_insensitive_names():
    assert load_matrix("ednafull").name == "NUC.4.4"
    assert load_matrix("blosum62").name == "BLOSUM62"


def test_scheme_rejects_mismatched_matrix_type():
    with pytest.raises(ValueError, match="nucleotide matrix"):
        get_scoring_scheme("NUC.4.4", dnaseq=False)
    with pytest.raises(ValueError, match="protein matrix"):
        get_scoring_scheme("BLOSUM62", dnaseq=True)
    with pytest.raises(ValueError, match="Unknown scoring matrix"):
        get_scoring_scheme("NOT_A_MATRIX", dnaseq=False)


def test_default_schemes():
    assert get_scoring_scheme(None, dnaseq=False).matrix.name == "BLOSUM50"
    assert get_scoring_scheme(None, dnaseq=True).matrix.name == "DNA"


def test_custom_matrix_file(tmp_path):
    path = tmp_path / "tiny.mat"
    path.write_text("# custom\n   A  C  G  T\nA  3 -1 -1 -1\nC -1  3 -1 -1\nG -1 -1  3 -1\nT -1 -1 -1  3\n")
    m = parse_matrix_file(str(path))
    assert m.is_nucleotide and m.score("G", "G") == 3
    assert get_scoring_scheme(str(path), dnaseq=True).matrix.name == "tiny"
