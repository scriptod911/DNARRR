import sys

def dna_to_rna(dna_sequence):
    rna_sequence = dna_sequence.replace("T", "U")
    return rna_sequence

if __name__ == "__main__":
    dna_sequence = sys.argv[1]
    rna_sequence = dna_to_rna(dna_sequence)
    print(rna_sequence)
