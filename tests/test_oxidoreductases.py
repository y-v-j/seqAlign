"""
End-to-end test on real data: GAPDH (EC 1.2.1.12, an oxidoreductase)
protein and coding-DNA sequences from 15 organisms, as downloaded by
``tests/oxidoreductases/download_sequences.py``.
"""

import csv
import os

import pytest

from seqalign import cli

DATA = os.path.join(os.path.dirname(__file__), "data", "oxidoreductases")


def _manifest():
    path = os.path.join(DATA, "manifest.tsv")
    if not os.path.isfile(path):
        pytest.skip("oxidoreductase test data missing; run tests/oxidoreductases/download_sequences.py")
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _parse_results(text):
    """(name0, name1) -> (identity, similarity) from seqAlign's text output."""
    lines = text.splitlines()
    results = {}
    for i, line in enumerate(lines):
        if line.endswith(" vs."):
            name0 = line.split()[0]
            name1 = lines[i + 1].split()[0]
            summary = next(ln for ln in lines[i + 2:] if "% identity" in ln)
            ident = float(summary.split("% identity")[0])
            sim = float(summary.split("; ")[1].split("% similarity")[0])
            results[(name0, name1)] = results[(name1, name0)] = (ident, sim)
    return results


@pytest.mark.parametrize("kind,matrix", [("protein", None), ("dna", "NUC.4.4")])
def test_fifteen_organism_all_vs_all(kind, matrix, tmp_path, capsys):
    rows = _manifest()
    assert len(rows) == 15
    files = [os.path.join(DATA, kind, f"{r['label']}.fasta") for r in rows]
    heatmap = tmp_path / f"{kind}.png"
    argv = files + ["-i", "-d", "--heatmap-output", str(heatmap), "--dpi", "60"]
    if matrix:
        argv += ["-s", matrix]

    assert cli.main(argv) == 0
    results = _parse_results(capsys.readouterr().out)

    assert len(results) == 15 * 14  # 105 pairs, stored both ways
    assert heatmap.stat().st_size > 0
    for ident, sim in results.values():
        assert 0 < ident <= sim <= 100

    # Biology sanity checks: GAPDH tracks the species tree.
    ident = {k: v[0] for k, v in results.items()}
    assert ident[("Mus_musculus", "Rattus_norvegicus")] > ident[("Mus_musculus", "Danio_rerio")]
    assert ident[("Homo_sapiens", "Mus_musculus")] > ident[("Homo_sapiens", "Escherichia_coli")]
    assert ident[("Oryza_sativa", "Zea_mays")] > ident[("Oryza_sativa", "Arabidopsis_thaliana")]
    assert ident[("Oryza_sativa", "Arabidopsis_thaliana")] > ident[("Oryza_sativa", "Homo_sapiens")]


def test_downloaded_cds_translate_to_downloaded_proteins():
    from tests.oxidoreductases.download_sequences import translate
    from seqalign.io_utils import read_sequence

    for row in _manifest():
        protein = read_sequence(os.path.join(DATA, "protein", f"{row['label']}.fasta")).sequence
        cds = read_sequence(os.path.join(DATA, "dna", f"{row['label']}.fasta"), force_dna=True).sequence
        peptide = translate(cds).rstrip("*")
        assert peptide[1:] == protein[1:], row["label"]
