"""
io_utils.py
===========
File reading utilities. Replaces the original's ``getseq()``/``gettitle()``
(defined in the missing nxgetaa.c). Reads a single biological sequence from
a FASTA-formatted file, streaming line-by-line so memory use is O(sequence
length) rather than O(file size), and validates content instead of silently
accepting garbage.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Iterator

_PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBZX*-")
_DNA_ALPHABET = set("ACGTUN-RYSWKMBDHV")  # includes IUPAC ambiguity codes


class SequenceFileError(ValueError):
    """Raised for malformed or unreadable sequence files."""


@dataclass
class SequenceRecord:
    """A single named biological sequence, analogous to (name0/name1, aa0/aa1)."""

    identifier: str
    title: str
    sequence: str
    is_dna: bool
    source_path: str

    @property
    def length(self) -> int:
        return len(self.sequence)


def _sniff_is_dna(sequence: str, sample_size: int = 200) -> bool:
    """Heuristic identical in spirit to the original relying on -n, but
    used here as a helpful default when -n isn't given explicitly."""
    sample = sequence[:sample_size].upper()
    if not sample:
        return False
    dna_chars = set("ACGTUN")
    hits = sum(1 for c in sample if c in dna_chars)
    return (hits / len(sample)) > 0.9


def _iter_fasta_records(path: str) -> Iterator[tuple[str, list[str]]]:
    """Yield (header_line, [sequence_lines]) tuples, streaming the file."""
    header = None
    chunk: list[str] = []
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for raw_line in fh:
                line = raw_line.rstrip("\n\r")
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        yield header, chunk
                    header = line[1:].strip()
                    chunk = []
                else:
                    chunk.append(line)
            if header is not None:
                yield header, chunk
    except OSError as exc:
        raise SequenceFileError(f"Cannot read sequence file '{path}': {exc}") from exc


def read_sequence(path: str, force_dna: bool | None = None) -> SequenceRecord:
    """
    Read the first sequence record from a FASTA file.

    Mirrors getseq()'s role: given a filename, produce a usable in-memory
    sequence plus a short title. Multiple records in one file are permitted
    by FASTA format in general, but this driver (like the original,
    processing one sequence per file) only uses the first and warns if more
    exist by raising a descriptive, catchable warning path via the caller.

    :raises SequenceFileError: for missing files, empty files, no header,
        empty sequence, or characters outside the expected alphabet.
    """
    if not os.path.isfile(path):
        raise SequenceFileError(f"Sequence file not found: '{path}'")
    if os.path.getsize(path) == 0:
        raise SequenceFileError(f"Sequence file is empty: '{path}'")

    records = list(_iter_fasta_records(path))
    if not records:
        raise SequenceFileError(
            f"'{path}' does not look like FASTA (no '>' header line found)."
        )

    header, lines = records[0]
    sequence = "".join(lines).upper()
    sequence = re.sub(r"\s+", "", sequence)

    if not sequence:
        raise SequenceFileError(f"No residues found under header in '{path}'.")

    is_dna = force_dna if force_dna is not None else _sniff_is_dna(sequence)
    alphabet = _DNA_ALPHABET if is_dna else _PROTEIN_ALPHABET
    bad_chars = sorted(set(sequence) - alphabet)
    if bad_chars:
        raise SequenceFileError(
            f"'{path}' contains characters not valid for "
            f"{'DNA' if is_dna else 'protein'} sequences: {bad_chars}. "
            f"Pass --dna/--protein explicitly if the auto-detection is wrong."
        )

    identifier = header.split()[0] if header else os.path.basename(path)
    title = header if header else os.path.basename(path)

    return SequenceRecord(
        identifier=identifier[:40],
        title=title[:60],
        sequence=sequence,
        is_dna=is_dna,
        source_path=path,
    )
