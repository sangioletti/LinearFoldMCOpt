# LinearFold MCOpt Source Package
# Core modules for mRNA codon optimization

from .mRNA import mRNA
from .codons import (
    codon_table,
    codon_adaptation_index,
    human_aa_to_codon_cai,
    ecoli_aa_to_codon_cai,
    aminoacid_to_codon_table,
    amino_acid_to_1_letter,
)
from .utils import (
    codon_to_amino_acid_1L,
    codon_to_amino_acid_3L,
    aa_1L_to_codons,
    aa_3L_to_codons,
    aa_to_codon_sequence,
)
from .linearpartition_wrapper import LinearPartitionWrapper

__all__ = [
    'mRNA',
    'codon_table',
    'codon_adaptation_index',
    'human_aa_to_codon_cai',
    'ecoli_aa_to_codon_cai',
    'aminoacid_to_codon_table',
    'amino_acid_to_1_letter',
    'codon_to_amino_acid_1L',
    'codon_to_amino_acid_3L',
    'aa_1L_to_codons',
    'aa_3L_to_codons',
    'aa_to_codon_sequence',
    'LinearPartitionWrapper',
]
