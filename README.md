# HELIX — DNA/RNA Six-Frame Translator

[🌐 Live Demo](https://helix-dna-rna-translator.vercel.app/) · [Source Code](https://github.com/mehreenfatimah/helix-dna-rna-translator)

HELIX is an educational bioinformatics web tool for validating nucleotide sequences, translating all six reading frames, detecting simple open reading frames (ORFs), and summarizing sequence composition.

## What it does

- accepts DNA, RNA, plain-text or FASTA input
- normalizes RNA (`U`) to DNA notation (`T`) for translation
- supports IUPAC ambiguous nucleotide symbols
- computes reverse complement, GC content and nucleotide composition
- translates the three forward and three reverse-complement reading frames
- supports several commonly used NCBI genetic-code tables included in the original project
- detects simple in-frame `ATG` → stop ORFs above a chosen minimum length
- can fetch nucleotide FASTA records from NCBI by accession when internet access is available
- exports/illustrates results through the browser interface

## Why this project matters

Six-frame translation is a foundational bioinformatics operation when the protein-coding frame of a nucleotide sequence is unknown. HELIX demonstrates sequence parsing, strand handling, codon translation, ORF detection, biological input validation and web application development.

## Repository structure

```text
app.py                         Flask web application
src/helix/core.py              sequence-analysis logic
templates/index.html            browser UI
tests/test_core.py              automated biological/unit tests
examples/input/                 sample FASTA
examples/output/                example exports from the original project
docs/                           project history and architecture image
archive/                        original final notebook preserved for provenance
```

## Run locally

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pip install -e .
python app.py
```

Open `http://127.0.0.1:3030`.

## Run tests

```bash
pip install -e .
pytest -q
```

The tests verify sequence normalization, reverse complementation, standard-code translation, six-frame output and handling of ambiguous codons.

## Example

Input:

```text
ATGGCTTAA
```

Frame `+1` translates to:

```text
MA*
```

where `*` represents a stop codon.

## Project background

HELIX originated as a university bioinformatics software project and was later refactored into this cleaner, tested public version. The original final notebook is preserved in the `archive/` directory, with additional development history documented in [`docs/PROJECT_HISTORY.md`](docs/PROJECT_HISTORY.md).

## Scientific limitations

- ORF detection is intentionally simple and uses `ATG` as the start codon; biologically valid alternative start codons are not exhaustively modeled.
- Ambiguous codons are reported as `?` rather than probabilistically resolved.
- The tool performs translation and descriptive sequence analysis; it is not a gene predictor or annotation pipeline.
- NCBI fetching depends on an internet connection and the public E-utilities service.

## Skills demonstrated

Python · Flask · sequence parsing · six-frame translation · reverse complement · ORF detection · FASTA · NCBI E-utilities · testing · refactoring
