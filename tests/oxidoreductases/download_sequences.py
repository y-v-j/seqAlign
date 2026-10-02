#!/usr/bin/env python3
"""
download_sequences.py
=====================
Test-data helper (NOT part of the seqAlign package): downloads the protein
and coding-DNA sequence of one oxidoreductase (EC class 1) for 15 organisms
spanning mammals, birds, fish, amphibians, insects, nematodes, fungi,
plants and bacteria, and writes one FASTA file per organism:

    tests/data/oxidoreductases/protein/<Organism>.fasta
    tests/data/oxidoreductases/dna/<Organism>.fasta
    tests/data/oxidoreductases/manifest.tsv

The default enzyme is glyceraldehyde-3-phosphate dehydrogenase (GAPDH,
EC 1.2.1.12) - an NAD-dependent oxidoreductase that is present and highly
conserved in every domain of life, which makes it a good all-vs-all
identity/similarity test case. Any other EC number can be used via --ec.

Sources (public REST APIs, standard library only - no extra dependencies):
  * Protein: UniProtKB (reviewed / Swiss-Prot entries)  https://rest.uniprot.org
  * DNA:     the coding sequence (CDS) behind the UniProt entry's EMBL/ENA
             cross-reference, from ENA (https://www.ebi.ac.uk/ena), with NCBI
             E-utilities as a fallback. Every CDS is translated and checked
             against the UniProt protein before it is accepted.

Usage:
    python tests/oxidoreductases/download_sequences.py            # GAPDH
    python tests/oxidoreductases/download_sequences.py --ec 1.15.1.1 --outdir /tmp/sod
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Dict, List, Optional, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUTDIR = os.path.join(HERE, "..", "data", "oxidoreductases")
DEFAULT_EC = "1.2.1.12"

# (label used as the FASTA identifier / heatmap label, NCBI taxonomy id)
ORGANISMS: List[Tuple[str, int]] = [
    ("Homo_sapiens", 9606),
    ("Mus_musculus", 10090),
    ("Rattus_norvegicus", 10116),
    ("Bos_taurus", 9913),
    ("Gallus_gallus", 9031),
    ("Xenopus_laevis", 8355),
    ("Danio_rerio", 7955),
    ("Drosophila_melanogaster", 7227),
    ("Caenorhabditis_elegans", 6239),
    ("Saccharomyces_cerevisiae", 559292),
    ("Arabidopsis_thaliana", 3702),
    ("Oryza_sativa", 39947),
    ("Zea_mays", 4577),
    ("Escherichia_coli", 83333),
    ("Bacillus_subtilis", 224308),
]

# Words that indicate a tissue/organelle-specific paralog rather than the
# main (cytosolic, glycolytic) enzyme; such entries are only used if
# nothing else is available.
_DEPRIORITISE = ("testis", "spermat", "chloroplast", "plastid", "mitochondrial", "glycosomal")

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
ENA_FASTA = "https://www.ebi.ac.uk/ena/browser/api/fasta/{}"
NCBI_CDS = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            "?db=protein&id={}&rettype=fasta_cds_na&retmode=text")

_CODON_TABLE = {
    a + b + c: aa
    for (a, b, c), aa in zip(
        ((x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
}


def _get(url: str, retries: int = 3, timeout: int = 60) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "seqAlign-test-downloader/1.0"})
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as resp:
                return resp.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt == retries:
                raise RuntimeError(f"GET {url} failed: {exc}") from exc
            time.sleep(2 * attempt)
    return ""


def _parse_fasta(text: str) -> Tuple[str, str]:
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines or not lines[0].startswith(">"):
        return "", ""
    seq = []
    for ln in lines[1:]:
        if ln.startswith(">"):
            break  # first record only
        seq.append(ln)
    return lines[0][1:], "".join(seq).upper()


def translate(cds: str) -> str:
    """Standard-code translation (unknown codons -> X), stop codons kept as '*'."""
    return "".join(_CODON_TABLE.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))


def _cds_matches(cds: str, protein: str) -> bool:
    if len(cds) % 3 or not cds:
        return False
    peptide = translate(cds).rstrip("*")
    # Alternative start codons (e.g. GTG/TTG in bacteria) still encode Met.
    return peptide[1:] == protein[1:] and len(peptide) == len(protein)


def _entry_rank(entry: dict, target_len: Optional[int]) -> Tuple:
    name = json.dumps(entry.get("proteinDescription", {})).lower()
    special = any(word in name for word in _DEPRIORITISE)
    length = entry["sequence"]["length"]
    return (
        special,
        -entry.get("annotationScore", 0),
        abs(length - target_len) if target_len else 0,
        entry["primaryAccession"],
    )


def find_protein(ec: str, taxid: int, target_len: Optional[int]) -> Optional[dict]:
    query = f"ec:{ec} AND organism_id:{taxid} AND reviewed:true"
    url = f"{UNIPROT_SEARCH}?{urllib.parse.urlencode({'query': query, 'format': 'json', 'size': 25})}"
    results = json.loads(_get(url)).get("results", [])
    if not results:
        return None
    return sorted(results, key=lambda e: _entry_rank(e, target_len))[0]


def _protein_name(entry: dict) -> str:
    desc = entry.get("proteinDescription", {})
    rec = desc.get("recommendedName") or (desc.get("submissionNames") or [{}])[0]
    return rec.get("fullName", {}).get("value", "unknown protein")


def _gene_name(entry: dict) -> str:
    genes = entry.get("genes") or [{}]
    return genes[0].get("geneName", {}).get("value", "")


def find_cds(entry: dict) -> Optional[Tuple[str, str, str]]:
    """Return (protein_id, source, cds) for the first EMBL xref whose CDS translates to the protein."""
    protein = entry["sequence"]["value"]
    xrefs = [x for x in entry.get("uniProtKBCrossReferences", []) if x["database"] == "EMBL"]
    # mRNA-derived CDS first, then genomic.
    xrefs.sort(key=lambda x: {p["key"]: p["value"] for p in x.get("properties", [])}.get("MoleculeType") != "mRNA")
    for xref in xrefs[:12]:
        props = {p["key"]: p["value"] for p in xref.get("properties", [])}
        protein_id = props.get("ProteinId", "-")
        if protein_id in ("-", ""):
            continue
        for source, template in (("ENA", ENA_FASTA), ("NCBI", NCBI_CDS)):
            try:
                _, cds = _parse_fasta(_get(template.format(protein_id)))
            except RuntimeError:
                continue
            if _cds_matches(cds, protein):
                return protein_id, source, cds
            time.sleep(0.34)
    return None


def _write_fasta(path: str, header: str, sequence: str, width: int = 60) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(sequence), width):
            fh.write(sequence[i:i + width] + "\n")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ec", default=DEFAULT_EC, help=f"EC number of the oxidoreductase (default: {DEFAULT_EC}, GAPDH)")
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR, help="output directory (default: tests/data/oxidoreductases)")
    parser.add_argument("--target-length", type=int, default=None,
                        help="preferred protein length when an organism has several entries (default: human entry length)")
    args = parser.parse_args(argv)

    if not args.ec.startswith("1."):
        print(f"warning: EC {args.ec} is not an oxidoreductase (EC class 1).", file=sys.stderr)

    prot_dir = os.path.join(args.outdir, "protein")
    dna_dir = os.path.join(args.outdir, "dna")
    os.makedirs(prot_dir, exist_ok=True)
    os.makedirs(dna_dir, exist_ok=True)

    manifest: List[Dict[str, str]] = []
    target_len = args.target_length
    failures = []
    for label, taxid in ORGANISMS:
        entry = find_protein(args.ec, taxid, target_len)
        if entry is None:
            print(f"  {label:<26} no reviewed UniProt entry for EC {args.ec}", file=sys.stderr)
            failures.append(label)
            continue
        acc = entry["primaryAccession"]
        protein = entry["sequence"]["value"]
        target_len = target_len or len(protein)
        gene = _gene_name(entry) or acc
        name = _protein_name(entry)
        organism = entry["organism"]["scientificName"]

        cds_hit = find_cds(entry)
        if cds_hit is None:
            print(f"  {label:<26} {acc}: no CDS matching the protein found", file=sys.stderr)
            failures.append(label)
            continue
        protein_id, source, cds = cds_hit

        _write_fasta(os.path.join(prot_dir, f"{label}.fasta"),
                     f"{label} {gene} UniProt:{acc} {name} [EC {args.ec}] OS={organism}", protein)
        _write_fasta(os.path.join(dna_dir, f"{label}.fasta"),
                     f"{label} {gene} CDS {source}:{protein_id} (UniProt:{acc}) OS={organism}", cds)
        manifest.append({
            "label": label, "organism": organism, "taxon_id": str(taxid), "gene": gene,
            "uniprot": acc, "protein_name": name, "protein_length": str(len(protein)),
            "cds_protein_id": protein_id, "cds_source": source, "cds_length": str(len(cds)),
        })
        print(f"  {label:<26} {acc:<8} {gene:<8} {len(protein):>4} aa   CDS {protein_id} ({source}) {len(cds)} nt")
        time.sleep(0.34)

    with open(os.path.join(args.outdir, "manifest.tsv"), "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest[0]) if manifest else ["label"], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest)

    print(f"\n{len(manifest)} organisms written to {os.path.abspath(args.outdir)}")
    if failures:
        print(f"failed: {', '.join(failures)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
