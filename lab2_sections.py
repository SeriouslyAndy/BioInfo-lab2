import os
import random
import sys

import matplotlib.pyplot as plt

from lab2 import clean_sequence, find_genes

ALPHABET = "ACGT"


def generate_dna(length):
    # build a random DNA sequence, one letter at a time
    sequence = ""
    for i in range(length):
        sequence = sequence + random.choice(ALPHABET)
    return sequence


def read_fasta(text):
    # split a FASTA text into (name, sequence) sections, any number of them
    sections = []
    name = None
    lines = []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                sections.append((name, clean_sequence("\n".join(lines))))
            name = line[1:].strip()
            lines = []
        else:
            lines.append(line)
    if name is not None:
        sections.append((name, clean_sequence("\n".join(lines))))
    return sections


def random_sections(count, min_len, max_len):
    sections = []
    for n in range(count):
        dna = generate_dna(random.randint(min_len, max_len))
        sections.append(("random_" + str(n + 1), clean_sequence(dna)))
    return sections


def plot_gene_counts(names, counts, out_file):
    best = max(counts)
    # one hue; the section(s) with the most genes get the darker step
    colors = ["#184f95" if c == best else "#86b6ef" for c in counts]

    fig, ax = plt.subplots(figsize=(max(8, len(counts) * 0.12), 4.5))
    x = range(1, len(counts) + 1)
    ax.bar(x, counts, color=colors, width=0.8)

    ax.axhline(best, color="#52514e", linewidth=1, linestyle="--")
    ax.text(len(counts) + 0.5, best, " max = " + str(best),
            va="center", ha="left", color="#0b0b0b", fontsize=9)

    ax.set_title("Valid genes (AUG ... stop) per section", loc="left")
    ax.set_xlabel("Section")
    ax.set_ylabel("Genes found")
    ax.set_xlim(0.3, len(counts) + 0.7)
    ax.yaxis.get_major_locator().set_params(integer=True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e5e4e0", linewidth=0.8)
    ax.set_axisbelow(True)

    fig.tight_layout()
    fig.savefig(out_file, dpi=150)
    print("Chart saved to", out_file)
    plt.show()


if __name__ == "__main__":
    # usage: python lab2_sections.py [file1.fasta file2.fasta ...]
    # with no files, 100 random sections are generated
    if len(sys.argv) > 1:
        sections = []
        for path in sys.argv[1:]:
            with open(path) as f:
                sections = sections + read_fasta(f.read())
    else:
        sections = random_sections(100, 300, 1500)

    names = []
    counts = []
    for name, rna in sections:
        genes = find_genes(rna)
        names.append(name)
        counts.append(len(genes))
        print(len(names), ":", name[:50], "|", len(rna), "nt |", len(genes), "genes")

    best = max(counts)
    print()
    print("Sections:", len(sections))
    print("Max genes in one section:", best)
    for name, c in zip(names, counts):
        if c == best:
            print("  ->", name)

    here = os.path.dirname(os.path.abspath(__file__))
    plot_gene_counts(names, counts, os.path.join(here, "genes_per_section.png"))
