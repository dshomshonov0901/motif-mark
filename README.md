# motif-mark

Object-oriented Python program that visualizes motifs from a FASTA file.  

Outputs a single PNG figure per FASTA file using `pycairo`.
---

## *Overview*

`motif-mark-oop.py`:

- Reads a FASTA file (≤ 10 sequences, ≤ 1000 bases each)
- Reads a motif file (≤ 5 motifs, ≤ 10 bases each, one per line)
- Finds motif matches (including ambiguous IUPAC nucleotide codes: https://www.promega.com/resources/guides/nucleic-acid-analysis/restriction-enzyme-resource/restriction-enzyme-resource-tables/iupac-ambiguity-codes-for-nucleotide-degeneracy/)
- Draws:
  - Introns as a baseline (Introns are defined by lowercase bases in the FASTA sequence)
  - Exons as filled rectangles (Exons are defined by uppercase bases in the FASTA sequence)
  - Motifs as colored bars (staggered if overlapping)
- Outputs a single PNG image to scale

---

## *Environment Setup*

Only `pycairo` is required. See documentation: https://pycairo.readthedocs.io/en/latest/

```bash
conda create -n my_pycairo pycairo
conda activate my_pycairo
```
## Usage

```bash
python motif-mark-oop.py -f input.fa -m motifs.txt
```

Arguments:

-f : FASTA file

-m : Motifs file

Output:

PNG file with same prefix as FASTA

Example: Figure_1.fa → Figure_1.png
