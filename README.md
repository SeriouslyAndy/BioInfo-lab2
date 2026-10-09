# Lab 2 – From gene to protein

Translate the coding region of a gene into an amino acid sequence using the
genetic code, keep only **valid** genes (AUG … stop), and count how many valid
genes each section of an unlimited set of sequences contains.

## Task

1. Convert the coding region of a gene into an amino acid sequence using the
   genetic code table.
2. Extra rules: **AUG** is the start of a molecule; **UAA**, **UAG** and
   **UGA** are the stops. Only valid sequences are extracted.
3. From an AUG, read on until a stop is met. Every AUG met on the way is
   disregarded as a new start until that stop is found.
4. A second module takes an unlimited number of sections (e.g. 100) and shows a
   chart with the number of genes found in each one and the maximum over the
   set. Inside one section there may be 3–4 possible starts, but some of them
   never reach a stop, so they are not genes.

---

## 1. Theory

### 1.1 Central dogma

$$
\text{DNA} \xrightarrow{\text{transcription}} \text{mRNA} \xrightarrow{\text{translation}} \text{protein}
$$

- **DNA** uses the alphabet $\Sigma_{DNA} = \{A, C, G, T\}$.
- **RNA** uses $\Sigma_{RNA} = \{A, C, G, U\}$: thymine (T) is replaced by uracil (U).
- In this lab, transcription is just $T \rightarrow U$ (the coding strand is
  given, so no complementing is needed). In code: `seq.upper().replace('T', 'U')`.

### 1.2 Codons

The ribosome reads mRNA in non-overlapping groups of **3 nucleotides** called
**codons**. Each codon maps to one amino acid or to a stop signal.

$$
\text{number of codons} = 4^3 = 64
$$

There are 61 sense codons for **20 amino acids** and 3 stop codons. Several
codons give the same amino acid, which is why the code is called **degenerate**
(redundant). For example, Leu has 6 codons, but Met and Trp have only one each.

### 1.3 The genetic code table

| 1st \ 2nd | U | C | A | G | 3rd |
|:-:|:-:|:-:|:-:|:-:|:-:|
| **U** | UUU Phe (F)<br>UUC Phe (F)<br>UUA Leu (L)<br>UUG Leu (L) | UCU Ser (S)<br>UCC Ser (S)<br>UCA Ser (S)<br>UCG Ser (S) | UAU Tyr (Y)<br>UAC Tyr (Y)<br>**UAA Stop**<br>**UAG Stop** | UGU Cys (C)<br>UGC Cys (C)<br>**UGA Stop**<br>UGG Trp (W) | U<br>C<br>A<br>G |
| **C** | CUU Leu (L)<br>CUC Leu (L)<br>CUA Leu (L)<br>CUG Leu (L) | CCU Pro (P)<br>CCC Pro (P)<br>CCA Pro (P)<br>CCG Pro (P) | CAU His (H)<br>CAC His (H)<br>CAA Gln (Q)<br>CAG Gln (Q) | CGU Arg (R)<br>CGC Arg (R)<br>CGA Arg (R)<br>CGG Arg (R) | U<br>C<br>A<br>G |
| **A** | AUU Ile (I)<br>AUC Ile (I)<br>AUA Ile (I)<br>**AUG Met (M) – Start** | ACU Thr (T)<br>ACC Thr (T)<br>ACA Thr (T)<br>ACG Thr (T) | AAU Asn (N)<br>AAC Asn (N)<br>AAA Lys (K)<br>AAG Lys (K) | AGU Ser (S)<br>AGC Ser (S)<br>AGA Arg (R)<br>AGG Arg (R) | U<br>C<br>A<br>G |
| **G** | GUU Val (V)<br>GUC Val (V)<br>GUA Val (V)<br>GUG Val (V) | GCU Ala (A)<br>GCC Ala (A)<br>GCA Ala (A)<br>GCG Ala (A) | GAU Asp (D)<br>GAC Asp (D)<br>GAA Glu (E)<br>GAG Glu (E) | GGU Gly (G)<br>GGC Gly (G)<br>GGA Gly (G)<br>GGG Gly (G) | U<br>C<br>A<br>G |

Amino acids are written either with **one letter** (`MAMGV`) or **three
letters** (`Met-Ala-Met-Gly-Val`). The program prints both.

### 1.4 Start and stop codons

| Signal | Codons | Meaning |
|---|---|---|
| Start | `AUG` | Translation begins here; it is also the amino acid Met (M) |
| Stop | `UAA`, `UAG`, `UGA` | Translation ends; no amino acid is added |

### 1.5 Reading frames

Since codons do not overlap, a sequence can be read in **3 reading frames**,
depending on where the first codon begins (position 1, 2 or 3):

```
RNA:      G G G C A U G G C C A U G G G C G U A U A A C C C
frame 1:  GGG CAU GGC CAU GGG CGU AUA ACC
frame 2:   GGC AUG GCC AUG GGC GUA UAA CCC      <- AUG ... UAA
frame 3:    GCA UGG CCA UGG GCG UAU AAC
```

The same letters give completely different proteins in each frame, so all
three are scanned independently. A frame shift of even one nucleotide changes
every codon after it.

### 1.6 Open reading frame (ORF) – what counts as a valid gene

A **valid gene** (ORF) here is a stretch of one reading frame that

1. starts with `AUG`,
2. continues codon by codon in the same frame,
3. ends at the **first** stop codon met in that frame.

$$
\text{gene} = \texttt{AUG}\;(c_1\,c_2\,\dots\,c_k)\;\text{STOP}, \qquad c_i \notin \{\texttt{UAA}, \texttt{UAG}, \texttt{UGA}\}
$$

Its length in nucleotides is a multiple of 3, and the protein has

$$
L_{aa} = \frac{L_{nt}}{3} - 1
$$

amino acids (the stop codon adds nothing).

**Rule "disregard all AUG until you find a valid stop":** once a gene has
started, any AUG met inside it is just the amino acid M, not a new gene.
Scanning for the next start resumes only **after** the stop codon.

**Rule "only valid sequences":** if an AUG never meets a stop before the end of
the sequence, it is a **possible** gene but **not a valid** one, so it is not
counted. In each frame there is at most one such start: the first AUG after
the last stop. Any later AUG in that frame lies inside its reading, so it is
disregarded.

Example `AUGAAAUGACCCAUGGGG`:

```
frame 1:  AUG AAA UGA | CCC AUG GGG        -> gene "MK" (valid), then AUG at 13 with no stop
frame 2:  UGA AAU GAC CCA UGG              -> no AUG
frame 3:  GAA AUG ACC CAU GGG              -> AUG at 6 with no stop
```

Result: **1 valid gene** and **2 starts with no stop** (not counted).

### 1.7 The algorithm

For each frame $f \in \{1, 2, 3\}$:

```
i = f - 1
while a whole codon fits at i:
    if codon(i) != AUG:            i += 3; continue
    j = i
    while a whole codon fits at j and codon(j) is not a stop:
        j += 3
    if no stop was found:          stop scanning this frame (invalid ORF)
    record gene  rna[i : j+3],  protein = translate(rna[i : j])
    i = j + 3                      # resume after the stop
```

Every nucleotide is visited a constant number of times per frame, so one
section takes $O(N)$ time, where $N$ is its length.

### 1.8 Counting genes over many sections

For a set of sections $S_1, S_2, \dots, S_n$ (with $n$ unlimited):

$$
g_k = \text{number of valid genes in } S_k, \qquad
g_{max} = \max_{1 \le k \le n} g_k
$$

The chart shows one bar per section with height $g_k$. The tallest bar or bars
are highlighted, and a dashed line marks $g_{max}$. Sections are read one at a
time, and only the counts are kept, so memory does not grow with the number of
sequences.

**What to expect from random DNA.** In random sequence, each codon is a stop
with probability $p = 3/64 \approx 4.7\%$. An ORF that has already started
therefore stops after on average

$$
E[\text{codons until stop}] = \frac{1}{p} = \frac{64}{3} \approx 21 \text{ codons}
$$

So random ORFs are short, about 20 aa. Real genes, such as influenza M1 at
252 aa or neuraminidase at 447 aa here, are much longer than chance would
give. Length is a simple signal for telling real genes apart from noise.

---

## 2. Data – Influenza A (`influenza_na.fasta`)

| Accession | Content | Valid genes | Starts with no stop | Longest gene |
|---|---|:-:|:-:|---|
| X52226 | A/FPV/Rostock/34 (H7N1) – neuraminidase | 22 | 0 | 447 aa, frame 3, nt 21–1364 (neuraminidase) |
| X53029 | A/USSR/90/1977 (H1N1) – matrix proteins 1 and 2 | 11 | 0 | 252 aa, frame 2, nt 26–784 (**M1** protein, 252 aa in reality too) |
| X59778 | A/NIB/4/1988 (H1N1) – polyprotein | 20 | 2 | 34 aa – the real coding region is a fragment with no stop, so it is **not valid** |

Maximum over the set: **22 genes** (X52226).

---

## 3. Files

| File | What it does |
|---|---|
| `lab2.py` | Genetic code table, translation (1- and 3-letter), `find_genes`, `find_open_starts`, streaming FASTA reader; prints every gene and the longest protein per section |
| `lab2_coding_region.py` | Translates one given coding region; valid only if it starts with AUG and has an in-frame stop |
| `lab2_sections.py` | Module 2: any number of FASTA files/sections (or 100 random ones) → gene count per section + chart `genes_per_section.png` |
| `main6.py` | Earlier version, influenza set: per-gene breakdown + chart |
| `main7.py` | Earlier version, general case: 100 synthetic sections + chart |
| `influenza_na.fasta` | The 3 influenza sequences |

All three versions apply the same rules and give the same counts.

## 4. Usage

```bash
pip install matplotlib

python lab2.py                               # influenza file, every gene + longest protein
python lab2.py my.fasta

python lab2_coding_region.py ATGGCCATGGGCGTATAA
# Protein (1-letter): MAMGV
# Protein (3-letter): Met-Ala-Met-Gly-Val

python lab2_sections.py                      # 100 random sections
python lab2_sections.py influenza_na.fasta other.fasta ...

python main6.py [file.fasta]
python main7.py [file.fasta]
```

## 5. Simplifications

- Only the given (forward) strand is scanned. Real gene finders also scan the
  3 frames of the reverse complement, which makes 6 frames in total.
- Only AUG is accepted as a start. Some organisms also use alternative starts
  such as GUG or UUG.
- There is no minimum gene length, so even `AUG UAA` counts as a valid gene
  (protein `M`, 1 aa).
- Introns, splicing and the real promoter and ribosome-binding signals are
  ignored.
