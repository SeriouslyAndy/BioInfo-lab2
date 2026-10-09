import os
import matplotlib.pyplot as plt

# Mapping codons to single-letter amino acids
CODON_MAP = {
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


def stream_fasta(file_path_or_text):
    """
    MODIFIED: Replaced list-based parsing with a generator to process an
    unlimited number of sequences without loading all sequences into RAM at once.
    """
    lines = file_path_or_text.splitlines() if "\n" in file_path_or_text else open(file_path_or_text, 'r')
    header = None
    chunks = []

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
                chunks = []
            header = line[1:]
        else:
            chunks.append(line)

    if header is not None:
        yield header, "".join(chunks)

    if not isinstance(file_path_or_text, str) or "\n" not in file_path_or_text:
        lines.close()


def count_valid_genes(sequence):
    """
    MODIFIED: Extracted from anterior extract_valid_proteins.
    Counts valid ORFs obeying the professor's rule:
    - Initiates translation at AUG.
    - Disregards any internal AUG within the frame until a stop codon is met.
    - Increments count only if terminated by a valid stop codon.
    """
    clean_seq = sequence.upper().replace('T', 'U').replace(" ", "").replace("\n", "")
    seq_len = len(clean_seq)
    gene_count = 0
    i = 0

    while i < seq_len - 2:
        if clean_seq[i:i + 3] == 'AUG':
            curr = i
            found_stop = False

            while curr + 3 <= seq_len:
                codon = clean_seq[curr:curr + 3]
                amino_acid = CODON_MAP.get(codon, '?')

                if amino_acid == 'Stop':
                    found_stop = True
                    i = curr + 3  # Skip past the stop codon
                    break
                curr += 3

            if found_stop:
                gene_count += 1
            else:
                i += 1
        else:
            i += 1

    return gene_count


def analyze_dataset_and_plot(records_source, chart_filename="gene_distribution.png"):
    """
    NEW MODULE: Processes an arbitrary number of sequences, calculates the
    maximum number of genes in the dataset, and generates the requested chart.
    """
    labels = []
    gene_counts = []

    # Stream through sequences (supports 100+ records)
    for idx, (header, sequence) in enumerate(stream_fasta(records_source), start=1):
        count = count_valid_genes(sequence)
        # Use short identifier for plotting clarity
        short_id = header.split()[0] if header else f"Seq_{idx}"
        labels.append(short_id)
        gene_counts.append(count)

    if not gene_counts:
        print("No valid sequences found to analyze.")
        return

    # Find dataset metrics
    max_genes = max(gene_counts)
    max_indices = [i for i, count in enumerate(gene_counts) if count == max_genes]
    max_labels = [labels[i] for i in max_indices]

    print(f"Total sequences analyzed: {len(gene_counts)}")
    print(f"Maximum genes found in any sequence: {max_genes}")
    print(f"Sequence(s) reaching maximum: {', '.join(max_labels)}")

    # Build the chart
    fig, ax = plt.subplots(figsize=(14, 6))
    indices = range(len(gene_counts))

    # Base bars
    colors = ['#1f77b4' if count < max_genes else '#d62728' for count in gene_counts]
    bars = ax.bar(indices, gene_counts, color=colors, width=0.7, edgecolor='black', linewidth=0.5)

    # Highlight line for the maximum value
    ax.axhline(
        y=max_genes,
        color='#d62728',
        linestyle='--',
        linewidth=1.5,
        label=f'Max Genes = {max_genes} (held by {len(max_indices)} sequence(s))'
    )

    # Annotate peak sequence(s)
    for max_idx in max_indices:
        ax.annotate(
            f'Max: {max_genes}\n({labels[max_idx]})',
            xy=(max_idx, max_genes),
            xytext=(max_idx, max_genes + (max_genes * 0.08 if max_genes > 0 else 1)),
            ha='center',
            fontsize=8,
            fontweight='bold',
            arrowprops=dict(facecolor='#d62728', arrowstyle='->', lw=1)
        )

    ax.set_title("Gene Frequency per Sequence Across Dataset", fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel("Sequence Index / ID", fontsize=11, labelpad=8)
    ax.set_ylabel("Number of Valid Genes Detected", fontsize=11, labelpad=8)

    # Adjust X-axis labels depending on sample size
    if len(gene_counts) <= 40:
        ax.set_xticks(indices)
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
    else:
        # Avoid cluttered tick labels when analyzing 100+ sequences
        step = max(1, len(gene_counts) // 20)
        ax.set_xticks(list(indices)[::step])
        ax.set_xticklabels([labels[i] for i in list(indices)[::step]], rotation=45, ha='right', fontsize=8)

    ax.set_ylim(0, max_genes * 1.25 if max_genes > 0 else 5)
    ax.legend(loc='upper right', frameon=True)
    ax.grid(axis='y', linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig(chart_filename, dpi=300)
    print(f"Chart successfully saved to {chart_filename}")
    plt.show()


# --- Demonstration with a synthetic 100-sequence test set ---
if __name__ == "__main__":
    import random

    # Generate synthetic mock FASTA data with 100 sequences to test scaling
    mock_fasta_records = []
    bases = ['A', 'C', 'G', 'T']

    for seq_num in range(1, 101):
        # Varying lengths between 300 and 1500 nt
        seq_len = random.randint(300, 1500)
        random_dna = "".join(random.choices(bases, k=seq_len))
        mock_fasta_records.append(f">Sequence_{seq_num:03d} synthetic test entry\n{random_dna}")

    mock_dataset = "\n".join(mock_fasta_records)

    # Run analysis and render chart
    analyze_dataset_and_plot(mock_dataset, chart_filename="gene_distribution.png")
