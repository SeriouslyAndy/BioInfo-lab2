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
    parts = []
    for line in text.splitlines():
        if not line.startswith(">"):
            parts.append(line.strip().upper())
    return "".join(parts).replace('T', 'U')


def read_fasta(lines):
    # yield (name, rna) one section at a time from any iterable of lines
    # (an open file or text.splitlines()), so the number of sections is
    # unlimited: only one sequence is held in memory at a time
    name = None
    chunks = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                yield name, clean_sequence("\n".join(chunks))
            name = line[1:]
            chunks = []
        else:
            chunks.append(line)
    if name is not None:
        yield name, clean_sequence("\n".join(chunks))
    elif chunks:
        # plain sequence without a header
        yield "sequence", clean_sequence("\n".join(chunks))


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


def find_open_starts(rna):
    # AUGs that never reach a stop (at most one per frame: the first AUG
    # after the last stop of that frame; later AUGs are inside its reading
    # and are disregarded). These are possible genes that are NOT valid.
    starts = []
    for frame in range(3):
        last_stop = frame - 3
        for i in range(frame, len(rna) - 2, 3):
            if rna[i : i + 3] in STOPS:
                last_stop = i
        for i in range(last_stop + 3, len(rna) - 2, 3):
            if rna[i : i + 3] == START:
                starts.append({"frame": frame + 1, "start": i + 1})
                break
    return starts


def longest_gene(genes):
    best = None
    for g in genes:
        if best is None or len(g["protein"]) > len(best["protein"]):
            best = g
    return best


def report(name, rna):
    genes = find_genes(rna)
    open_starts = find_open_starts(rna)

    print("=" * 70)
    print(name)
    print("Len:", len(rna), "nt")
    print("Valid genes (AUG ... stop):", len(genes))
    print("AUG with no stop (not valid):", len(open_starts))
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
    print()


if __name__ == "__main__":
    # usage: python lab2.py [file.fasta]   (default: the influenza file)
    here = os.path.dirname(os.path.abspath(__file__))
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "influenza_na.fasta")
    with open(path) as f:
        for name, rna in read_fasta(f):
            report(name, rna)
