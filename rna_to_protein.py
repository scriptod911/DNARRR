import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna
from Bio.Alphabet import generic_protein

def rna_to_protein(rna_sequence):
    # Create RNA sequence object
    rna = Seq(rna_sequence, generic_rna)

    # Translate RNA to protein using standard genetic code
    protein = rna.translate(to_stop=False)

    # Convert protein to string and return
    return str(protein)

if __name__ == "__main__":
    rna_sequence = sys.argv[1]
    protein_sequence = rna_to_protein(rna_sequence)
    print(protein_sequence)
