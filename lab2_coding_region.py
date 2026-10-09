import sys

from lab2 import clean_sequence, to_three_letter, translate_coding_region

# translate one coding region (AUG ... stop), in one-letter and three-letter form
# usage: python lab2_coding_region.py ATGGCCATGGGCGTATAA

if __name__ == "__main__":
    seq = sys.argv[1] if len(sys.argv) > 1 else "ATGGCCATGGGCGTATAA"
    protein = translate_coding_region(seq)

    print("Seq:", seq)
    print("RNA:", clean_sequence(seq))
    if protein is None:
        print("Not a valid coding region (must start with AUG and have a stop).")
    else:
        print("Protein (1-letter):", protein)
        print("Protein (3-letter):", to_three_letter(protein))
