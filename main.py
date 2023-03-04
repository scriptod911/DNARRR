from dna_to_rna import dna_to_rna
from rna_to_protein import rna_to_protein

def main(dna_sequence, mode):
	rna_sequence = dna_to_rna(dna_sequence)
	protein_sequence = rna_to_protein(rna_sequence)

	if mode == "normal":
		print(f"Protein sequence:\n{protein_sequence}")

if __name__ == "__main__":
	dna_sequence = input("Enter a DNA sequence: ")
	mode = input("Enter 'normal' to just translate to protein sequence: ")
	main(dna_sequence, mode)
