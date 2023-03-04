from Bio.Seq import Seq

def rna_to_protein(rna_sequence):
    seq_obj = Seq(rna_sequence)
    return seq_obj.translate()
