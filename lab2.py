import os
import sys

# genetic code table: RNA codon -> one-letter amino acid
CODON_TABLE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'C', 'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# one-letter -> three-letter amino acid names
THREE_LETTER = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    '*': 'Stop', '?': '???',
}

START = 'AUG'
STOPS = ('UAA', 'UAG', 'UGA')


def clean_sequence(text):
    # drop FASTA header lines and whitespace, convert DNA (T) to RNA (U)
    seq = ""
    for line in text.splitlines():
        if not line.startswith(">"):
            seq = seq + line.strip().upper()
    return seq.replace('T', 'U')


def translate(rna):
    # translate codon by codon (no start/stop rules here)
    protein = ""
    for i in range(0, len(rna) - 2, 3):
        aa = CODON_TABLE.get(rna[i : i + 3], '?')
        if aa == 'Stop':
            aa = '*'
        protein = protein + aa
    return protein


def to_three_letter(protein):
    # MAMGV -> Met-Ala-Met-Gly-Val
    names = []
    for aa in protein:
        names.append(THREE_LETTER[aa])
    return "-".join(names)


def translate_coding_region(sequence):
    # the input is already a coding region: it must start with AUG and
    # contain an in-frame stop. Returns the protein, or None if not valid.
    rna = clean_sequence(sequence)
    if not rna.startswith(START):
        return None
    for i in range(0, len(rna) - 2, 3):
        if rna[i : i + 3] in STOPS:
            return translate(rna[:i])
    return None  # no stop codon -> not a valid gene


def find_genes(rna):
    # a valid gene starts at AUG and ends at the first in-frame stop.
    # AUGs inside a gene are just M: they are ignored as new starts
    # until the stop is reached. An AUG with no stop after it is not a gene.
    genes = []
    for frame in range(3):
        i = frame
        while i + 3 <= len(rna):
            if rna[i : i + 3] != START:
                i = i + 3
                continue

            # walk from this AUG looking for a stop in the same frame
            j = i
            while j + 3 <= len(rna) and rna[j : j + 3] not in STOPS:
                j = j + 3

            if j + 3 > len(rna):
                # ran off the end without a stop -> invalid, and every later
                # AUG in this frame would hit the same end, so stop the frame
                break

            genes.append({
                "frame": frame + 1,
                "start": i + 1,          # 1-based position of the A in AUG
                "end": j + 3,            # 1-based position of the last stop base
                "stop": rna[j : j + 3],
                "protein": translate(rna[i:j]),
            })
            i = j + 3  # continue after the stop
    genes.sort(key=lambda g: g["start"])
    return genes


def longest_gene(genes):
    best = None
    for g in genes:
        if best is None or len(g["protein"]) > len(best["protein"]):
            best = g
    return best


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "influenza_na.fasta")
    with open(path) as f:
        rna = clean_sequence(f.read())

    genes = find_genes(rna)

    print("File:", path)
    print("Len:", len(rna), "nt")
    print("Valid genes (AUG ... stop):", len(genes))
    print()
    for n, g in enumerate(genes, 1):
        print(n, ": frame", g["frame"], "| pos", g["start"], "-", g["end"],
              "|", g["stop"], "|", len(g["protein"]), "aa")

    best = longest_gene(genes)
    if best is None:
        print("No valid gene found.")
    else:
        print()
        print("Longest gene: frame", best["frame"], "pos", best["start"], "-", best["end"],
              "(", len(best["protein"]), "aa )")
        print(best["protein"])
        print()
        print(to_three_letter(best["protein"]))
