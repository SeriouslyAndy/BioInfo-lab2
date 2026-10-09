import matplotlib.pyplot as plt

# Mapping codons to single-letter amino acid abbreviations
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


def stream_fasta(fasta_source):
    """
    Generator that parses FASTA input line-by-line, yielding (header, sequence)
    without loading all sequence contents into memory simultaneously.
    """
    lines = fasta_source.splitlines() if "\n" in fasta_source else open(fasta_source, 'r')
    header = None
    seq_chunks = []

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_chunks)
                seq_chunks = []
            header = line[1:]
        else:
            seq_chunks.append(line)

    if header is not None:
        yield header, "".join(seq_chunks)

    if not isinstance(fasta_source, str) or "\n" not in fasta_source:
        lines.close()


def find_genes_in_frame(clean_seq, frame_offset):
    """
    MODIFIED: Extracted single-frame ORF scanning into a dedicated helper.
    Scans a sequence along a specific reading frame (offset 0, 1, or 2):
    - Initiates only on AUG.
    - Disregards any downstream AUGs within this frame until an in-frame stop codon is reached.
    - Discards incomplete ORFs that do not encounter a stop codon.
    """
    genes = []
    seq_len = len(clean_seq)
    pos = frame_offset

    while pos + 3 <= seq_len:
        codon = clean_seq[pos:pos + 3]
        if codon == 'AUG':
            protein = []
            start_pos = pos
            curr = pos
            found_stop = False

            while curr + 3 <= seq_len:
                c = clean_seq[curr:curr + 3]
                aa = CODON_MAP.get(c, '?')

                if aa == 'Stop':
                    found_stop = True
                    # Set the outer frame pointer to start after this stop codon
                    pos = curr + 3
                    break

                protein.append(aa)
                curr += 3

            if found_stop:
                genes.append({
                    'start': start_pos,
                    'stop': curr + 3,
                    'length_aa': len(protein),
                    'frame': frame_offset + 1,
                    'protein': "".join(protein)
                })
            else:
                # Ran out of sequence without hitting stop: discard incomplete ORF
                pos += 3
        else:
            pos += 3

    return genes


def count_valid_genes(sequence, check_all_frames=True):
    """
    MODIFIED: Now accepts a check_all_frames flag.
    When True, scans all 3 forward reading frames (offsets 0, 1, 2) independently,
    preventing frame collisions and identifying all possible coding subsequences.
    """
    clean_seq = sequence.upper().replace('T', 'U').replace(" ", "").replace("\n", "")
    frames_to_check = [0, 1, 2] if check_all_frames else [0]
    all_genes = []

    for frame in frames_to_check:
        all_genes.extend(find_genes_in_frame(clean_seq, frame))

    return all_genes


def analyze_influenza_set(fasta_text, check_all_frames=True, chart_output="influenza_genes_chart.png"):
    """
    NEW: Processes the provided Influenza sequences, computes max gene count,
    prints detailed breakdown per sequence, and generates the chart.
    """
    labels = []
    full_headers = []
    gene_counts = []
    details_per_seq = []

    for idx, (header, seq) in enumerate(stream_fasta(fasta_text), start=1):
        genes = count_valid_genes(seq, check_all_frames=check_all_frames)
        # Extract accession or short identifier (e.g., gi|58576|gb|X52226)
        short_id = header.split()[0] if header else f"Seq_{idx}"

        labels.append(short_id)
        full_headers.append(header)
        gene_counts.append(len(genes))
        details_per_seq.append(genes)

    if not gene_counts:
        print("No valid sequences found.")
        return

    max_genes = max(gene_counts)
    max_indices = [i for i, c in enumerate(gene_counts) if c == max_genes]

    # Print analysis table to terminal
    print("=" * 80)
    print(f"INFLUENZA GENE DETECTION (Reading Frames Scanned: {'1, 2, 3' if check_all_frames else '1 only'})")
    print("=" * 80)
    for idx, (header, count, genes) in enumerate(zip(full_headers, gene_counts, details_per_seq), start=1):
        print(f"\n[{idx}] {header}")
        print(f"    Total Valid Genes Detected: {count}")
        for g_idx, g in enumerate(genes, start=1):
            print(f"      - Gene {g_idx}: Frame +{g['frame']}, Start: nt {g['start']}, "
                  f"Stop: nt {g['stop']}, Length: {g['length_aa']} aa")

    print("\n" + "-" * 80)
    print(f"Maximum genes found in any sequence: {max_genes}")
    print(f"Sequence(s) reaching maximum: {', '.join(labels[i] for i in max_indices)}")
    print("-" * 80)

    # MODIFIED: Build visualization highlighting max sequence(s)
    fig, ax = plt.subplots(figsize=(10, 6))
    indices = range(len(gene_counts))

    bar_colors = ['#d62728' if c == max_genes else '#1f77b4' for c in gene_counts]
    bars = ax.bar(indices, gene_counts, color=bar_colors, width=0.55, edgecolor='black', linewidth=0.8)

    # Draw maximum reference line
    ax.axhline(
        y=max_genes,
        color='#d62728',
        linestyle='--',
        linewidth=1.5,
        label=f'Max Genes = {max_genes}'
    )

    # Annotate peak bars
    for max_idx in max_indices:
        ax.annotate(
            f'Max: {max_genes}',
            xy=(max_idx, max_genes),
            xytext=(max_idx, max_genes + 0.35),
            ha='center',
            fontsize=10,
            fontweight='bold',
            arrowprops=dict(facecolor='#d62728', arrowstyle='->', lw=1.2)
        )

    # Configure axes
    ax.set_title("Influenza Dataset: Valid Gene Count per Sequence", fontsize=13, fontweight='bold', pad=15)
    ax.set_xlabel("Genomic Sequence ID", fontsize=11, labelpad=10)
    ax.set_ylabel("Valid Genes Count (AUG to Stop)", fontsize=11, labelpad=10)
    ax.set_xticks(indices)
    ax.set_xticklabels(labels, rotation=20, ha='right', fontsize=9)
    ax.set_ylim(0, max(max_genes + 1.5, 4))
    ax.legend(loc='upper right')
    ax.grid(axis='y', linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig(chart_output, dpi=300)
    print(f"\nChart saved to {chart_output}")
    plt.show()


# --- Execution with the 3 Influenza Sequences ---
if __name__ == "__main__":
    influenza_dataset = """
>gi|58576|gb|X52226|Influenza A virus (A/FPV/Rostock/34(H7N1)) gene for neuraminidase, genomic RNA
AGCAAAAGCAGGAGTTCAAAATGAATCCAAATCAGAAAATAATAACCATTGGGTCAATCTGTATGGGGAT
CGGAATAATCAGCCTAATATTACAAATTGGAAACATAATCTCAATGTGGGTTAGTCATTCAATTCAGACT
GAAAATCAAAATCACCATGAAGCATGCAACCCAAGCATTGCTGGACAGGATGCAGCTTCAGTGGCACTAG
CAGGCAATTCCTCTCTTTGTCCCATTAGTGGGTGGGCTATATACAGTAAAGACAATGGTATAAGAATTGG
CTCCAAAGGAGACGTATTTGTCATAAGAGAGCCATTTATTTCATGCTCTCACTTGGAATGCAGGACCTTT
TTTCTGACTCAGGGCGCCTTGTTGAATGACAAACACTCAAATGGAACCGTTAAAGACAGAAGCCCTTATA
GAACCTTGATGAGCTGTCCTGTTGGTGAAGCTCCTTCTCCCTACAATTCAAGGTTCGTGTCGGTTGCATG
GTCAGCAAGTGCCTGCCATGATGGCATGGGTTGGCTAACAATCGGAATTTCTGGTCCAGATAATGGAGCA
GTGGCTGTATTAAAATACAATGGTATAATAACAGACACCATCAAAAGTTGGAAAAATAACATATTGAGAA
CGCAAGAGTCTGAATGTGCCTGTATAAATGGTTCATGTTTCACTATAATGACCGATGGCCCAAGTAATGG
GCAGGCCTCGTACAAAATTTTCAAGATAGAAAAGGGGAAGGTAGTCAAATCAAGTGAATTAAATGCACCT
AATTACCACTACGAGGAATGTTCTTGTTATCCTGATGCAGGTGAAGTAATGTGTGTATGCAGGGATAATT
GGCATGGTTCGAATCGACCATGGGTGTCTTTCAATAAAAACCTTGATTATCAAATAGGTTACATCTGCAG
TGGGGTTTTTGGTGATAATCCACGACCCAATGATGGAACAGGCAGCTGTGGTCCAGTGTCTTCTAATGGA
GCATATGGGATAAAGGGGTTCTCATTTAAGTATGGTAATGGTGTTTGGATAGGGAGAACCAAAAGCACTA
GTTCCAGAAGCGGGTTTGAGATGATTTGGGATCCCAATGGATGGACAGAAACTGACAGCAGTTTCTCTGT
GAAGCAAGATATTGTAGCAATAACTGATTGGTCGGGATATAGCGGGAGTTTCGTCCAACACCCTGAATTA
ACAGGGTTGGACTGCATGAGGCCTTGTTTCTGGGTTGAACTGATCAGGGGACGGCCCAACCACAATACGA
TCTGGACTAGTGGGAGCAGCATTTCCTTCTGTGGTGTGAACAGCGATACTGTAGGTTGGTCTTGGCCAGA
CGGTGCTGAGTTGCCATTCACCATTGACAAGTAGTTTGTTCAAAAAAACTCCTTGTTTCTACT
>gi|59292|gb|X53029|Influenza A virus (A/USSR/90/1977(H1N1)) genes for matrix proteins 1 and 2, genomic RNA
AGCAAAAGCAGGTAGATGTTGAAAGATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTC
CCGTCAGGCCCCCTCAAAGCCGAGATCGCACAGAGACTTGAAGATGTCTTTGCTGGGAAGAACACCGATC
TTGAGGCTCTCATGGAATGGCTAAAGACAAGACCAATCCTGTCACCTCTGACTAAGGGGATTTTAGGATT
TGTGTTCACGCTCACCGTGCCCAGTGAGCGAGGACTGCAGCGTAGACGCTTTGTCCAAAATGCCCTTAAT
GGGAATGGGGATCCAAATAACATGGACAGAGCAGTTAAACTGTATAGAAAGCTTAAGAGGGAGATAACAT
TCCATGGGGCCAAAGAAATAGCACTCAGTTATTCTGCTGGTGCACTTGCCAGTTGTATGGGCCTCATATA
CAACAGGATGGGGGCTGTGACCACCGAAGCGGCATTTGGCCTGATATGCGCAACCTGTGAACAGATTGCT
GACTCCCAGCATAGGTCTCATAGGCAAATGGTGACAACAACCAATCCACTAATAAGACATGAGAACAGAA
TGGTTCTGGCCAGCACTACAGCTAAGGCTATGGAGCAAATGGCTGGATCGAGTGAGCAAGCAGCAGAGGC
CATGGAGGTTGCTAGTCAGGCCAGGCAAATGGTGCAGGCAATGAGAGCCATTGGGACTCATCCTAGCTCC
AGTGCTGGTCTGAAAAATGATCTTCTTGAAAATTTGCAGGCCTATCAGAAACGAATGGGGGTGCAGATGC
AACGATTCAAGTGATCCTCTTGTTGTTGCCGCAAGTATCATTGGGATTTTGCACTTGATATTGTGGATTC
TTGATCGTCTTTTTTTCAAATGCATTTATCGTCTCTTTAAACACGGTCTGAAAAGAGGGCCTTCTACGGA
AGGAGTACCAGAGTCTATGAGGGAAGAATATCGAAAGGAACAGCAGAATGCTGTGGATGCTGACGATAGT
CATTTTGTCAACATAGAGCTAGAGTAAAAAACTACCTTGTTTCTACT
>gi|59294|gb|X59778|Influenza A virus (A/NIB/4/1988(H1N1)) gene for polyprotein, genomic RNA
GGATCCTGTTATGTGCATTTACAGCTACAGATGCAGACACAATATGTATAGGCTACCATGCGAACAACTC
AACCGACACTGTTGACACAGTACTTGAGAAGAACGTGACAGTGACACACTCTGTCAACCTACTTGAGGAC
AGTCACAACAGAAAACTATGTCGACTAAAAGGAATAGCCCCACTACAATTGGGTAATTGCAGCATTGCCG
GATGGATCTTAGGAAACCCAGAATGCGAATCACTGTTTTCTAAGGAATCATGGTCCTACATTGCAGAAAC
ACCAAACTCCGAGAATGGAACATGTTACCCAGGGTATTTCGCCGACTATGAGGAACTGAGGGAGCAATTG
AGTTCAGTATCATCATTCGAGGGATTCGAAATATTCCCCGAAGAAAGCTCATGGCCCAACCACACCGTAA
CCAAGGGAGTAACGGCATCATGCTCCCATAATGGGAAAAGCAGTTTTTACAGAAATTTGCTATGGCTGAC
GGAGAAGAATGGCTTGTACCCAAATCTGAGCAAGTCCTATGTAAACAACAAAGAGAAAGAAGTCCTTGTA
CTATGGGGTGTTCATCACCCGTCTAACATAGGGGACCAAAGGGCCATCTATCATACAGAAAATGCTTATG
TCTCTGTAGTGTCTTCACATTATAGCAGGAGATTCACCCCAGAAATAGCAAAAAGACCCAAAGTAAGAGA
TCAAGAAGGAAGAATTAACTACTACTGGACTCTGCTGGAACCCGGGGACACAATAATATTTGAGGCAAAT
GGAAATCTAATAGCGCCATGGTATGCTTTCGCACTGAGTAGAGGCTTTGGGTCAGGAATCATCACCTCAA
ACGCATCAATGGATGAATGTGACGCGAAGTGTCAAACACCCCAGGGAGCTATAAACAGTAGTCTTCCTTT
CCAGAATGTACACCCAGTCACAATAGGAGAGTGTCCAAAGTATGTCAGGAGTACAAAATTAAGGATGGTT
ACAGGACTAAGGAACATCCCATCCATTCAATCCAGAGGTTTGTTTGGAGCCATTGCCGGTTTCATTGAAG
GGGGATGGACTGGAATTC
"""

    analyze_influenza_set(influenza_dataset, check_all_frames=True)
