"""Core sequence-analysis functions for HELIX.

The implementation is based on the original university HELIX translator project,
refactored from a notebook/web prototype into importable and testable functions.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, asdict
from typing import Dict, Iterable, List

_BASE = {
    'TTT':'Phe','TTC':'Phe','TTA':'Leu','TTG':'Leu',
    'CTT':'Leu','CTC':'Leu','CTA':'Leu','CTG':'Leu',
    'ATT':'Ile','ATC':'Ile','ATA':'Ile','ATG':'Met',
    'GTT':'Val','GTC':'Val','GTA':'Val','GTG':'Val',
    'TCT':'Ser','TCC':'Ser','TCA':'Ser','TCG':'Ser',
    'CCT':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro',
    'ACT':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr',
    'GCT':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala',
    'TAT':'Tyr','TAC':'Tyr','TAA':'Stop','TAG':'Stop',
    'CAT':'His','CAC':'His','CAA':'Gln','CAG':'Gln',
    'AAT':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys',
    'GAT':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu',
    'TGT':'Cys','TGC':'Cys','TGA':'Stop','TGG':'Trp',
    'CGT':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg',
    'AGT':'Ser','AGC':'Ser','AGA':'Arg','AGG':'Arg',
    'GGT':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly',
}

def _patch(base: Dict[str,str], overrides: Dict[str,str]) -> Dict[str,str]:
    table = dict(base)
    table.update(overrides)
    return table

GENETIC_CODES = {
    '1':  ('Standard (NCBI Table 1)', _BASE),
    '2':  ('Vertebrate Mitochondrial', _patch(_BASE, {'TGA':'Trp','ATA':'Met','AGA':'Stop','AGG':'Stop'})),
    '3':  ('Yeast Mitochondrial', _patch(_BASE, {'TGA':'Trp','ATA':'Met','CTT':'Thr','CTC':'Thr','CTA':'Thr','CTG':'Thr'})),
    '4':  ('Mold/Protozoan Mitochondrial', _patch(_BASE, {'TGA':'Trp'})),
    '5':  ('Invertebrate Mitochondrial', _patch(_BASE, {'TGA':'Trp','ATA':'Met','AGA':'Ser','AGG':'Ser'})),
    '6':  ('Ciliate Nuclear', _patch(_BASE, {'TAA':'Gln','TAG':'Gln'})),
    '9':  ('Echinoderm Mitochondrial', _patch(_BASE, {'AAA':'Asn','AGA':'Ser'})),
    '10': ('Euplotid Nuclear', _patch(_BASE, {'TGA':'Cys'})),
    '11': ('Bacterial/Archaeal', _BASE),
    '12': ('Alternative Yeast Nuclear', _patch(_BASE, {'CTG':'Ser'})),
    '13': ('Ascidian Mitochondrial', _patch(_BASE, {'TGA':'Trp','ATA':'Met','AGA':'Gly','AGG':'Gly'})),
}

AA_1LETTER = {
    'Phe':'F','Leu':'L','Ile':'I','Met':'M','Val':'V','Ser':'S','Pro':'P',
    'Thr':'T','Ala':'A','Tyr':'Y','His':'H','Gln':'Q','Asn':'N','Lys':'K',
    'Asp':'D','Glu':'E','Cys':'C','Trp':'W','Arg':'R','Gly':'G','Stop':'*',
}

COMPLEMENT = str.maketrans(
    'ATGCatgcRYSWKMBDHVNryswkmbdhvn',
    'TACGtacgYRSWMKVHDBNyrswmkvhdbn',
)

@dataclass(frozen=True)
class CodonCall:
    codon: str
    aa3: str
    aa1: str
    pos: int

@dataclass(frozen=True)
class ORF:
    start_codon_idx: int
    stop_codon_idx: int
    length_nt: int


def clean_sequence(raw: str) -> str:
    """Normalize a FASTA/plain DNA or RNA sequence to uppercase DNA notation.

    FASTA headers, whitespace and digits are ignored. U is converted to T.
    IUPAC ambiguous nucleotide codes are accepted; codons containing ambiguity
    translate to '?' rather than being guessed.
    """
    lines = raw.strip().splitlines()
    sequence = ''.join(line for line in lines if not line.lstrip().startswith('>'))
    sequence = re.sub(r'[\s0-9]', '', sequence).upper().replace('U', 'T')
    invalid = set(sequence) - set('ATGCRYSWKMBDHVN')
    if invalid:
        raise ValueError(f"Invalid character(s): {', '.join(sorted(invalid))}")
    if len(sequence) < 3:
        raise ValueError('Sequence must contain at least 3 nucleotides.')
    return sequence


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def translate_frame(sequence: str, offset: int, table: Dict[str,str]) -> List[CodonCall]:
    calls: List[CodonCall] = []
    for i in range(offset, len(sequence)-2, 3):
        codon = sequence[i:i+3]
        aa3 = table.get(codon, '?')
        calls.append(CodonCall(codon, aa3, AA_1LETTER.get(aa3, '?'), i))
    return calls


def find_orfs(codons: List[CodonCall], min_length_nt: int = 60) -> List[ORF]:
    """Find simple ATG-to-stop ORFs in one translated frame.

    ORFs are non-overlapping in this educational implementation: after the first
    start codon, the next in-frame stop closes the ORF.
    """
    orfs: List[ORF] = []
    start = None
    for i, call in enumerate(codons):
        if call.aa3 == 'Met' and start is None:
            start = i
        elif call.aa3 == 'Stop' and start is not None:
            length_nt = (i - start + 1) * 3
            if length_nt >= min_length_nt:
                orfs.append(ORF(start, i, length_nt))
            start = None
    return orfs


def sequence_statistics(sequence: str) -> dict:
    length = len(sequence)
    canonical = {base: sequence.count(base) for base in 'ATGC'}
    gc = (canonical['G'] + canonical['C']) / length * 100 if length else 0
    composition = {base: round(canonical[base] / length * 100, 2) for base in 'ATGC'}
    ambiguous = length - sum(canonical.values())
    return {
        'length': length,
        'gc_percent': round(gc, 2),
        'composition_percent': composition,
        'ambiguous_bases': ambiguous,
    }


def analyze_six_frames(raw_sequence: str, genetic_code: str='1', min_orf_length_nt: int=60) -> dict:
    sequence = clean_sequence(raw_sequence)
    if genetic_code not in GENETIC_CODES:
        raise ValueError(f'Unsupported genetic-code table: {genetic_code}')
    table_name, table = GENETIC_CODES[genetic_code]
    rc = reverse_complement(sequence)
    frames = []
    for strand, seq in [('+', sequence), ('-', rc)]:
        for offset in range(3):
            codons = translate_frame(seq, offset, table)
            orfs = find_orfs(codons, min_orf_length_nt)
            frames.append({
                'label': f'{strand}{offset+1}',
                'strand': strand,
                'offset': offset,
                'protein': ''.join(c.aa1 for c in codons),
                'codons': [asdict(c) for c in codons],
                'orfs': [asdict(o) for o in orfs],
            })
    stats = sequence_statistics(sequence)
    return {
        **stats,
        'genetic_code': genetic_code,
        'genetic_code_name': table_name,
        'total_orfs': sum(len(frame['orfs']) for frame in frames),
        'frames': frames,
    }
