# This code implements the codon table for the 20 standard amino acids.
# Complete codon table with all 64 codons (mRNA notation with U instead of T):

codon_table = {
    # Alanine (Ala)
    'GCU': 'Alanine',
    'GCC': 'Alanine',
    'GCA': 'Alanine',
    'GCG': 'Alanine',
    
    # Arginine (Arg)
    'CGU': 'Arginine',
    'CGC': 'Arginine',
    'CGA': 'Arginine',
    'CGG': 'Arginine',
    'AGA': 'Arginine',
    'AGG': 'Arginine',
    
    # Asparagine (Asn)
    'AAU': 'Asparagine',
    'AAC': 'Asparagine',
    
    # Aspartic acid (Asp)
    'GAU': 'Aspartic acid',
    'GAC': 'Aspartic acid',
    
    # Cysteine (Cys)
    'UGU': 'Cysteine',
    'UGC': 'Cysteine',
    
    # Glutamic acid (Glu)
    'GAA': 'Glutamic acid',
    'GAG': 'Glutamic acid',
    
    # Glutamine (Gln)
    'CAA': 'Glutamine',
    'CAG': 'Glutamine',
    
    # Glycine (Gly)
    'GGU': 'Glycine',
    'GGC': 'Glycine',
    'GGA': 'Glycine',
    'GGG': 'Glycine',
    
    # Histidine (His)
    'CAU': 'Histidine',
    'CAC': 'Histidine',
    
    # Isoleucine (Ile)
    'AUU': 'Isoleucine',
    'AUC': 'Isoleucine',
    'AUA': 'Isoleucine',
    
    # Leucine (Leu)
    'CUA': 'Leucine',
    'CUG': 'Leucine',
    'CUU': 'Leucine',
    'CUC': 'Leucine',
    'UUA': 'Leucine',
    'UUG': 'Leucine',
    
    # Lysine (Lys)
    'AAA': 'Lysine',
    'AAG': 'Lysine',
    
    # Methionine (Met) - also Start codon
    'AUG': 'Methionine',
    
    # Phenylalanine (Phe)
    'UUU': 'Phenylalanine',
    'UUC': 'Phenylalanine',
    
    # Proline (Pro)
    'CCU': 'Proline',
    'CCC': 'Proline',
    'CCA': 'Proline',
    'CCG': 'Proline',
    
    # Serine (Ser)
    'UCU': 'Serine',
    'UCC': 'Serine',
    'UCA': 'Serine',
    'UCG': 'Serine',
    'AGU': 'Serine',
    'AGC': 'Serine',
    
    # Threonine (Thr)
    'ACU': 'Threonine',
    'ACC': 'Threonine',
    'ACA': 'Threonine',
    'ACG': 'Threonine',
    
    # Tryptophan (Trp)
    'UGG': 'Tryptophan',
    
    # Tyrosine (Tyr)
    'UAU': 'Tyrosine',
    'UAC': 'Tyrosine',
    
    # Valine (Val)
    'GUU': 'Valine',
    'GUC': 'Valine',
    'GUA': 'Valine',
    'GUG': 'Valine',
    
    # Stop codons
    'UAA': 'Stop',
    'UAG': 'Stop',
    'UGA': 'Stop',
}

# Amino acid name → tuple of codons
aminoacid_to_codon_table = {
    'Alanine': ('GCU', 'GCC', 'GCA', 'GCG'),
    'Arginine': ('CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
    'Asparagine': ('AAU', 'AAC'),
    'Aspartic acid': ('GAU', 'GAC'),
    'Cysteine': ('UGU', 'UGC'),
    'Glutamic acid': ('GAA', 'GAG'),
    'Glutamine': ('CAA', 'CAG'),
    'Glycine': ('GGU', 'GGC', 'GGA', 'GGG'),
    'Histidine': ('CAU', 'CAC'),
    'Isoleucine': ('AUU', 'AUC', 'AUA'),
    'Leucine': ('CUA', 'CUG', 'CUU', 'CUC', 'UUA', 'UUG'),
    'Lysine': ('AAA', 'AAG'),
    'Methionine': ('AUG',),
    'Phenylalanine': ('UUU', 'UUC'),
    'Proline': ('CCU', 'CCC', 'CCA', 'CCG'),
    'Serine': ('UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'),
    'Threonine': ('ACU', 'ACC', 'ACA', 'ACG'),
    'Tryptophan': ('UGG',),
    'Tyrosine': ('UAU', 'UAC'),
    'Valine': ('GUU', 'GUC', 'GUA', 'GUG'),
    'Stop': ('UAA', 'UAG', 'UGA'),
}

# Full name → three-letter code
amino_acid_to_3_letter = {
    'Alanine': 'Ala',
    'Arginine': 'Arg',
    'Asparagine': 'Asn',
    'Aspartic acid': 'Asp',
    'Cysteine': 'Cys',
    'Glutamic acid': 'Glu',
    'Glutamine': 'Gln',
    'Glycine': 'Gly',
    'Histidine': 'His',
    'Isoleucine': 'Ile',
    'Leucine': 'Leu',
    'Lysine': 'Lys',
    'Methionine': 'Met',
    'Phenylalanine': 'Phe',
    'Proline': 'Pro',
    'Serine': 'Ser',
    'Threonine': 'Thr',
    'Tryptophan': 'Trp',
    'Tyrosine': 'Tyr',
    'Valine': 'Val',
    'Stop': None,
}

amino_acid_3L_to_full_name = {v: k for k, v in amino_acid_to_3_letter.items()}

# Full name → one-letter code (use None for Stop)
amino_acid_to_1_letter = {
    'Alanine': 'A',
    'Arginine': 'R',
    'Asparagine': 'N',
    'Aspartic acid': 'D',
    'Cysteine': 'C',
    'Glutamic acid': 'E',
    'Glutamine': 'Q',
    'Glycine': 'G',
    'Histidine': 'H',
    'Isoleucine': 'I',
    'Leucine': 'L',
    'Lysine': 'K',
    'Methionine': 'M',
    'Phenylalanine': 'F',
    'Proline': 'P',
    'Serine': 'S',
    'Threonine': 'T',
    'Tryptophan': 'W',
    'Tyrosine': 'Y',
    'Valine': 'V',
    'Stop': "*",   # or *
}

amino_acid_1L_to_full_name = {short: full for full, short in amino_acid_to_1_letter.items()}

# Codon Adaptation Index (CAI) for Humans
# Values normalized so the most preferred codon for each amino acid = 1.0
# Based on codon usage in highly expressed human genes
human_cai = {
    # Alanine
    'GCU': 0.34,
    'GCC': 1.00,
    'GCA': 0.36,
    'GCG': 0.20,
    
    # Arginine
    'CGU': 1.00,
    'CGC': 0.36,
    'CGA': 0.07,
    'CGG': 0.11,
    'AGA': 0.15,
    'AGG': 0.15,
    
    # Asparagine
    'AAU': 0.41,
    'AAC': 0.59,
    
    # Aspartic acid
    'GAU': 1.00,
    'GAC': 0.59,
    
    # Cysteine
    'UGU': 1.00,
    'UGC': 0.59,
    
    # Glutamic acid
    'GAA': 1.00,
    'GAG': 0.59,
    
    # Glutamine
    'CAA': 0.34,
    'CAG': 1.00,
    
    # Glycine
    'GGU': 0.34,
    'GGC': 1.00,
    'GGA': 0.36,
    'GGG': 0.20,
    
    # Histidine
    'CAU': 0.41,
    'CAC': 0.59,
    
    # Isoleucine
    'AUU': 0.47,
    'AUC': 1.00,
    'AUA': 0.08,
    
    # Leucine
    'CUA': 0.07,
    'CUG': 1.00,
    'CUU': 0.36,
    'CUC': 0.20,
    'UUA': 0.07,
    'UUG': 0.13,
    
    # Lysine
    'AAA': 0.33,
    'AAG': 1.00,
    
    # Methionine (Start)
    'AUG': 1.00,
    
    # Phenylalanine
    'UUU': 1.00,
    'UUC': 0.59,
    
    # Proline
    'CCU': 1.00,
    'CCC': 0.31,
    'CCA': 0.36,
    'CCG': 0.20,
    
    # Serine
    'UCU': 1.00,
    'UCC': 0.31,
    'UCA': 0.36,
    'UCG': 0.20,
    'AGU': 0.35,
    'AGC': 0.65,
    
    # Threonine
    'ACU': 0.34,
    'ACC': 1.00,
    'ACA': 0.36,
    'ACG': 0.20,
    
    # Tryptophan
    'UGG': 1.00,
    
    # Tyrosine
    'UAU': 1.00,
    'UAC': 0.59,
    
    # Valine
    'GUU': 0.46,
    'GUC': 0.47,
    'GUA': 0.25,
    'GUG': 1.00,
    
    # Stop codons
    'UAA': 0.30,
    'UAG': 0.24,
    'UGA': 0.47,
}

# Codon Adaptation Index (CAI) for E. coli
# Values normalized so the most preferred codon for each amino acid = 1.0
# Based on codon usage in highly expressed E. coli genes (Sharp & Li, 1987)
ecoli_cai = {
    # Alanine
    'GCU': 1.00,
    'GCC': 0.56,
    'GCA': 0.25,
    'GCG': 0.13,
    
    # Arginine
    'CGU': 0.36,
    'CGC': 1.00,
    'CGA': 0.04,
    'CGG': 0.04,
    'AGA': 0.02,
    'AGG': 0.02,
    
    # Asparagine
    'AAU': 1.00,
    'AAC': 0.51,
    
    # Aspartic acid
    'GAU': 1.00,
    'GAC': 0.45,
    
    # Cysteine
    'UGU': 1.00,
    'UGC': 0.45,
    
    # Glutamic acid
    'GAA': 1.00,
    'GAG': 0.45,
    
    # Glutamine
    'CAA': 0.30,
    'CAG': 1.00,
    
    # Glycine
    'GGU': 1.00,
    'GGC': 0.56,
    'GGA': 0.25,
    'GGG': 0.13,
    
    # Histidine
    'CAU': 1.00,
    'CAC': 0.45,
    
    # Isoleucine
    'AUU': 1.00,
    'AUC': 0.60,
    'AUA': 0.04,
    
    # Leucine
    'CUA': 0.04,
    'CUG': 1.00,
    'CUU': 0.25,
    'CUC': 0.13,
    'UUA': 0.04,
    'UUG': 0.13,
    
    # Lysine
    'AAA': 1.00,
    'AAG': 0.36,
    
    # Methionine (Start)
    'AUG': 1.00,
    
    # Phenylalanine
    'UUU': 1.00,
    'UUC': 0.45,
    
    # Proline
    'CCU': 0.30,
    'CCC': 0.19,
    'CCA': 0.25,
    'CCG': 1.00,
    
    # Serine
    'UCU': 1.00,
    'UCC': 0.30,
    'UCA': 0.25,
    'UCG': 0.13,
    'AGU': 1.00,
    'AGC': 0.38,
    
    # Threonine
    'ACU': 1.00,
    'ACC': 0.56,
    'ACA': 0.28,
    'ACG': 0.13,
    
    # Tryptophan
    'UGG': 1.00,
    
    # Tyrosine
    'UAU': 1.00,
    'UAC': 0.45,
    
    # Valine
    'GUU': 0.40,
    'GUC': 0.35,
    'GUA': 0.25,
    'GUG': 1.00,
    
    # Stop codons
    'UAA': 0.30,
    'UAG': 0.24,
    'UGA': 0.47,
}

# Codon Adaptation Index dictionary for different organisms
codon_adaptation_index = {
    'human': human_cai,
    'mouse': None,
    'rat': None,
    'yeast': None,
    'e_coli': ecoli_cai,
    'other': None,
}

# Amino acid name → {codon: CAI_value} for humans
human_aa_to_codon_cai = {
    amino_acid: {codon: human_cai[codon] for codon in codons}
    for amino_acid, codons in aminoacid_to_codon_table.items()
}

# Amino acid name → {codon: CAI_value} for E. coli
ecoli_aa_to_codon_cai = {
    amino_acid: {codon: ecoli_cai[codon] for codon in codons}
    for amino_acid, codons in aminoacid_to_codon_table.items()
}

restriction_sites = {'Sbf1':'CCTGCAGG','Nhe1':'GCTAGC','Age1':'ACCGGT'}
