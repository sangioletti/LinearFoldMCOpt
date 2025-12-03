import numpy as np
from codons import *
from mRNA import *

verbose = False

# Download the luciferase amino acid sequence from the PDB FASTA, then convert to nucleotide sequence
from biotite.database import rcsb
from biotite.sequence import ProteinSequence
from biotite.sequence.io import fasta

# Example: use luciferase PDB id "1BA3"
pdb_id = "1BA3"

# Download the FASTA for the PDB entry
fasta_file = rcsb.fetch(pdb_id, "fasta", target_path=".")

# Read the sequence from the FASTA file
fasta_file_obj = fasta.FastaFile.read(fasta_file)
# Get the first sequence (FastaFile is dictionary-like)
if len(fasta_file_obj) == 0:
    raise ValueError("No sequences found in FASTA file")
# Get first header and corresponding sequence
first_header = list(fasta_file_obj.keys())[0]
aa_sequence = str(fasta_file_obj[first_header])

# Now, convert this amino acid sequence into a nucleotide sequence using one possible codon per amino acid
# We'll use the 'codon_table' from codons.py, which should provide a mapping for each amino acid

def aa_to_codon_sequence(aa_sequence):
    # Convert amino acid 1-letter code to nucleotide sequence
    # Use aa_1L_to_codons to get possible codons for each amino acid
    nucleotide_seq = ""
    for aa in aa_sequence:
        if aa == "*" or aa == "X":
            # Stop codon or unknown, use UAA
            codon = "UAA"
        else:
            # Get possible codons for this amino acid
            possible_codons = aa_1L_to_codons(aa)
            if not possible_codons:
                raise ValueError(f"No codon found for amino acid: {aa}")
            # Use the first codon possibility for simplicity
            codon = possible_codons[0]
        nucleotide_seq += codon
    return nucleotide_seq

# Convert amino acid sequence to nucleotide sequence
nucleotide_sequence = aa_to_codon_sequence(aa_sequence)



sequence = nucleotide_sequence
sequence = "AUUUGGUGGAGG" + "AUG" + sequence + "UAA"
weights = {
                #'mfe': 1.0, 
                'cai': 10.0, 
                'fe': 1.0, 
                'cpg': 5.0, 
                'stem': 10.0, 
                'utr_hybridisation': 10.0, 
                'initial_hybridisation': 10.0, 
                #'restriction_sites': 0.0 
                }
system = mRNA( 
                sequence = sequence, 
                species = 'human', 
                aa_to_codon_cai=human_aa_to_codon_cai,
                verbose = verbose, 
                loss_weights=weights, 
                modify_utr=False,
                initial_region_end_index = 30,
                T_K = 310,
                bpp_cutoff = 0.01,
                beamsize=100,
                )
print(f"Sequence length: {len(sequence)}")
print(f"McCaskill free-energy (*codon): {system.free_energy} kcal/mol per codon")
system.optimize_codon_usage( 
    T_opt=0.3, 
    nsteps=1000, 
    sample_frequency = 100, 
    verbose = False,
    output_filename = "opt_statistics.txt",
    average_stem_num_samples = 10,
    )
print(f"Optimized mRNA: {system.sequence}")
print(f"Optimized codon usage: {system._codons}")
print(f"Initial amino acid sequence: {system.initial_aminoacid}")
print(f"Optimized amino acid usage: {system.codons_to_amino_acids()}")
print(f"Optimized loss: {system.loss}")
print(f"Optimized mfe (*codon): {system.mfe} kcal/mol per codon")
print(f"Optimized mfe (total): {system.mfe*len(system.codons)} kcal/mol")
print(f"McCaskill free-energy (*codon): {system.free_energy} kcal/mol per codon")
print(f"Optimized cai: {np.exp(system.calculate_CAI_log)}")
print(f"Structure: {system._structure}")
if system._prob_matrix is not None:
    print(f"P(i,j) matrix:")
    for i in range(system._prob_matrix.shape[0]):
        for j in range(system._prob_matrix.shape[1]):
            if system._prob_matrix[i,j] > 10**(-2):  # Only print significant probabilities
                print(f"P({i},{j}) = {system._prob_matrix[i,j]:.4f}")
else:
    print("P(i,j) matrix not available (call free_energy first)")
print(f"Printing structure with minimum energy in structure.svg")
system.visualize_structure(filename="structure.svg", format="svg")
print(f"Structure visualized in structure.svg")

