# This code implements the class mRNA and its methods, used to optimize 
# the codon usage of an mRNA sequence. 

import numpy as np
import time
from codons import *
from linearpartition_wrapper import LinearPartitionWrapper
from utils import *
import copy

class mRNA:
    def __init__(
        self, 
        cds:str|None, 
        fivep_utr:str|None,
        threep_utr:str|None,
        species:str, 
        aa_to_codon_cai:dict, 
        loss_weights:dict = {'mfe': 1.0, 'cai': 1.0},
        T_K : float = 310, 
        modify_fivep_utr:bool = False,
        modify_threep_utr:bool = False,
        immutable_fivep_utr_range : tuple[int, int]|None = None,
        immutable_threep_utr_range : tuple[int, int]|None = None,
        verbose:bool = False,
        initial_region_end_index:int|None = None,
        beamsize:int = 100,
        bpp_cutoff : float = 0.0,
        start_from_optimal_cai : bool = True,
        ):
        # Retrieve codon adaptation index dictionary for the species
        if species not in ['human', 'mouse', 'rat', 'yeast', 'e_coli', 'other']:
            raise ValueError("Invalid species")

        # Validate the sequence
        self.validate_and_store_sequence(fivep_utr = fivep_utr, threep_utr = threep_utr, cds = cds)

        # Then set the codons. This is purely done using the cds sequence
        self.codons = self.cds_sequence
        self._initial_codons = self.codons.copy()  # Store initial codon sequence for statistics

        # Set the initial region end index
        if initial_region_end_index is not None:
            self.initial_region_end_index = initial_region_end_index

        self._initial_aminoacid = self.codons_to_amino_acids()

        # Initialize mutability array (all mutable by default)
        self.nucleotides_mutability = np.ones(len(self.sequence))
        

        # Initialize codon mutability array
        self.codon_mutability = np.ones(len(self.codons))
        
        # Always make Kozak immutable
        for i in range(self.kozak_index_range[0], self.kozak_index_range[1]+1):
            self.nucleotides_mutability[i] = 0

        if not modify_fivep_utr:
            for i in range(self.fivep_utr_index_range[1]+1):
                    self.nucleotides_mutability[i] = 0

        if not modify_threep_utr:
            first_index = self.threep_utr_index_range[0]
            last_index = self.threep_utr_index_range[1]
            for i in range(first_index, last_index+1):
                    self.nucleotides_mutability[i] = 0 
        
        # Apply immutable_range if specified
        if immutable_fivep_utr_range is not None:
            start, end = immutable_fivep_utr_range
            for i in range(start, min(end+1, len(self.sequence))):
                self.nucleotides_mutability[i] = 0
        
        # Apply immutable_range if specified
        if immutable_threep_utr_range is not None:
            start, end = immutable_threep_utr_range
            first_index = self.threep_utr_index_range[0]
            for i in range(first_index + start, first_index + end+1):
                self.nucleotides_mutability[i] = 0

        # Set these attributes before start_from_optimal_cai() is called
        self._species = species
        self.aa_to_codon_cai = aa_to_codon_cai

        if start_from_optimal_cai:
            self.start_from_optimal_cai()
            print(f"Initial (Optimal CAI) codon sequence: {self._initial_codons}")
        else:
            print(f"Initial (Random) codon sequence: {self._initial_codons}")

        self.verbose = verbose
        self.RT = 0.001987 * T_K  # kcal/(mol*K) * K
        self._minimum_folding_energy = None
        self._free_energy = None
        self._structure = None
        self._prob_matrix = None
        self._cpg_count = None
        self._cai_log = None
        self._loss_weights = loss_weights
        self.validate_weights()
        # Optimization: Sequence-based caching (initialize early)
        self._cached_seq_hash = None  # For free_energy caching
        self._cached_mfe_seq_hash = None  # For mfe caching
        self._cached_bpp_seq_hash = None  # For prob_matrix caching
        self._cached_structure_seq_hash = None  # For structure caching
        self._cached_codons_string = None  # For codons_string caching
        self._codons_changed = True  # Track if codons changed
        # Optimization: Initialize LinearPartition wrapper
        # Use a reasonable beam size - smaller is faster but less accurate
        # For optimization, we can use a smaller beam size for speed
        self._linearpartition = LinearPartitionWrapper(use_vienna=True, beamsize=beamsize, verbose=verbose)
        self._bpp_cutoff = bpp_cutoff  # BPP cutoff - can be increased to speed up BPP calculation
        self._pf_computed = False  # Track if partition function has been computed
        # Note: _fivep_utr and _threep_utr attributes already initialized above (before codon_mutability)
        self._largest_stem = None
        self.set_largest_stem() # This will initialize self.largest_stem and self.largest_stem_length
        self.modify_fivep_utr = modify_fivep_utr
        self.modify_threep_utr = modify_threep_utr
        self._skip_svg = False  # SVG format is now used (may cause segfaults in some ViennaRNA versions)
        self._use_average_stem_length = False  # Flag to use average stem length from ensemble vs longest stem from MFE
        self._cached_average_stem_length = None  # Cache for average stem length
        self._average_stem_num_samples = 100  # Number of samples for average stem calculation

    @property
    def codons_string(self) -> str:
        """Return the sequence string from codons."""
        if self._codons is None or len(self._codons) == 0:
            raise ValueError("Codons not initialized")
        
        # Cache the string to avoid repeated joins
        if self._cached_codons_string is None or self._codons_changed:
            self._cached_codons_string = ''.join(self.codons)
            self._codons_changed = False
        
        seq = self._cached_codons_string
        
        # Validate sequence contains only valid RNA nucleotides
        valid_bases = set('AUCG')
        if not all(base in valid_bases for base in seq):
            invalid = [base for base in seq if base not in valid_bases]
            raise ValueError(f"Sequence contains invalid nucleotides: {set(invalid)}")
        return seq

    def validate_and_store_sequence(self, fivep_utr:str, threep_utr:str, cds:str):
        # First convert T to U and capitalize
        fivep_utr = fivep_utr.upper().replace('T', 'U')
        threep_utr = threep_utr.upper().replace('T', 'U')
        cds = cds.upper().replace('T', 'U')
        print( "Sequence converted to uppercase and T replaced with U")
        
        # Check length of each part is multiple of 3
        if len(cds) % 3 != 0:
            raise ValueError(f"CDS length is not divisible by 3: {len(cds)}")

        # Check that the 5'-utr does not contain more than a start codon at the end
        if fivep_utr.count('AUG') > 1:
            print(f"WARNING: 5' UTR contains more than a start codon: number of AUG: {fivep_utr.count('AUG')}")

        # Check that fivep_utr ends with kozak sequence
        kozak = ['ACCAUG','GCCAUG']
        if fivep_utr[-6:] not in kozak:
            raise ValueError(f"5' UTR does not end with kozak sequence: {fivep_utr[-6:]}")
        strict_kozak = ['CACCAUG','CGCCAUG']
        if fivep_utr[-7:] not in strict_kozak:
            print(f"WARNING: 5' UTR does not end with strict kozak sequence: {fivep_utr[-7:]}")

        # If start codon is repeated in CDS, remove the second one
        if cds[:3] == 'AUG':
            cds = cds[3:]
            print(f"Removed start codon from CDS because already present in 5' UTR: {cds}")
        
        # Check that CDS ends with a stop codon, otherwise add one 
        if cds[-3:] not in ['UAA', 'UAG', 'UGA']:
            cds = cds + 'UAA'
            print(f"Added stop codon to CDS because it does not end with a stop codon: {cds}")

        # NOW set codons after all CDS modifications
        self.codons = cds
        print(f"Initial codons: {self._codons}")
        
        # Check that CDS contains only valid codons
        for codon in self._codons:
            if codon not in codon_table:
                raise ValueError(f"Invalid codon: {codon}")

        # Check the coding sequence has a single stop codon
        stop_codons = ['UAA', 'UAG', 'UGA']
        stop_count = sum(1 for c in self._codons if c in stop_codons)
        
        if stop_count > 1:
            raise ValueError(f"Too many stop codons in the CDS: {stop_count} found")

        # Define the final sequence
        # Store the sequences in a dictionary
        self.sequence = {'fivep_utr': fivep_utr, 'cds': cds, 'threep_utr': threep_utr}
        
        # set the indices of the fivep_utr, cds and threep_utr
        self._start_index = len(fivep_utr)
        self._stop_index = self._start_index + len(cds)
        self.cds_index_range = (self._start_index, self._stop_index)
        self.fivep_utr_index_range = (0, self._start_index - 1)
        self.threep_utr_index_range = ( self.cds_index_range[1] + 1, len(self.sequence) - 1 )

        
        # Check that all nucleotides are valid
        if not all(base in 'AUCG' for base in self.sequence):
            raise ValueError("Sequence contains invalid nucleotides")

        return self.sequence


    @property
    def sequence(self) -> str:
        return self._sequence['fivep_utr'] + self._sequence['cds'] + self._sequence['threep_utr']

    @sequence.setter
    def sequence(self, value:dict):
        self._sequence = value
        return
    
    @property
    def kozak_sequence(self) -> str:
        return self._sequence['fivep_utr'][-9:]
    
    @property
    def fivep_utr_sequence(self) -> str:
        return self._sequence['fivep_utr']
    
    @fivep_utr_sequence.setter
    def fivep_utr_sequence(self, value:str):
        self._sequence['fivep_utr'] = value
        return
    
    @property
    def cds_sequence(self) -> str:
        return self._sequence['cds']
    
    @cds_sequence.setter
    def cds_sequence(self, value:str):
        self._sequence['cds'] = value
        return
    
    @property
    def threep_utr_sequence(self) -> str:
        return self._sequence['threep_utr']
    
    @threep_utr_sequence.setter
    def threep_utr_sequence(self, value:str):
        self._sequence['threep_utr'] = value
        return

    @property
    def fivep_utr_index_range(self) -> tuple:
        return self._fivep_utr_index_range
    
    @fivep_utr_index_range.setter
    def fivep_utr_index_range(self, value:tuple):
        self._fivep_utr_index_range = value
        return
    
    @property
    def kozak_index_range(self) -> tuple:
        """Return Kozak sequence range (last 9 nucleotides of 5' UTR including start codon)."""
        _, last = self._fivep_utr_index_range
        return (last - 8, last)  # 9 positions: indices [last-8, last] inclusive

    @kozak_index_range.setter
    def kozak_index_range(self, value:tuple):
        self._kozak_index_range = value
        return
    
    @property
    def cds_index_range(self) -> tuple:
        return self._cds_index_range

    @cds_index_range.setter
    def cds_index_range(self, value:tuple):
        self._cds_index_range = value
        return
    
    @property
    def threep_utr_index_range(self) -> tuple:
        return self._threep_utr_index_range

    @threep_utr_index_range.setter
    def threep_utr_index_range(self, value:tuple):
        self._threep_utr_index_range = value
        return

    @property
    def nucleotides_mutability(self):
        return self._nucleotides_mutability

    @nucleotides_mutability.setter
    def nucleotides_mutability(self, value:np.ndarray):
        self._nucleotides_mutability = value
        return

    @property
    def codon_mutability(self):
        return self._codon_mutability

    @codon_mutability.setter
    def codon_mutability(self, value:np.ndarray):
        self._codon_mutability = value
        return  

    @property
    def n_nucleotides(self):
        return len(self.sequence)
    
    def validate_weights(self):
        all_keys = set(self._loss_weights.keys())
        valid_keys = ['mfe', 'fe', 'cai', 'cpg', 'stem', 'threep_utr_hybridisation', 
                      'fivep_utr_hybridisation', 'initial_hybridisation', 'restriction_sites', 
                      'kozak', 'codon_divergence']
        for k in all_keys:
            if k not in valid_keys: 
                raise ValueError(f"Invalid loss weight key: {k}. Valid keys are: {valid_keys}")
        
        if self._loss_weights.get('mfe',0.0) != 0.0 and self._loss_weights.get('fe',0.0) != 0.0:
            raise ValueError("Only one of mfe or fe can be non-zero")
        return

    def codons_to_amino_acids(self):
        codons = self.codons
        amino_acids = []
        for codon in codons:
            amino_acids.append(codon_to_amino_acid_1L(codon))
            print(f"Codon: {codon}, Amino acid: {amino_acids[-1]}")
        return amino_acids
        
    @property
    def p_5p(self):
        n_5p_utr = self.nucleotides_mutability[self.fivep_utr_index_range[0]:self.fivep_utr_index_range[1]].sum()
        n_cds = self.nucleotides_mutability[self.cds_index_range[0]:self.cds_index_range[1]].sum()
        n_3p_utr = self.nucleotides_mutability[self.threep_utr_index_range[0]:self.threep_utr_index_range[1]].sum()
        
        return n_5p_utr / (n_5p_utr + n_cds + n_3p_utr)

    @property
    def p_cds(self):
        n_5p_utr = self.nucleotides_mutability[self.fivep_utr_index_range[0]:self.fivep_utr_index_range[1]].sum()
        n_cds = self.nucleotides_mutability[self.cds_index_range[0]:self.cds_index_range[1]].sum()
        n_3p_utr = self.nucleotides_mutability[self.threep_utr_index_range[0]:self.threep_utr_index_range[1]].sum()
        
        return n_cds / (n_5p_utr + n_cds + n_3p_utr)

    @property
    def p_3p(self):
        n_5p_utr = self.nucleotides_mutability[self.fivep_utr_index_range[0]:self.fivep_utr_index_range[1]].sum()
        n_cds = self.nucleotides_mutability[self.cds_index_range[0]:self.cds_index_range[1]].sum()
        n_3p_utr = self.nucleotides_mutability[self.threep_utr_index_range[0]:self.threep_utr_index_range[1]].sum()
        
        return n_3p_utr / (n_5p_utr + n_cds + n_3p_utr)

    def make_mutation(self):
        region = np.random.choice(['5p_utr', 'cds', '3p_utr'], p = [self.p_5p, self.p_cds, self.p_3p])
        
        if region == '5p_utr':
            # Get mutable positions in 5' UTR (respects nucleotides_mutability)
            start, end = self.fivep_utr_index_range
            mutable_positions = np.where(self.nucleotides_mutability[start:end+1] > 0)[0]
            if len(mutable_positions) == 0:
                return  # No mutable positions
            index = np.random.choice(mutable_positions)
            seq = self._sequence['fivep_utr']
            nucleotide = seq[index]
            possible_nucleotides = ['A', 'C', 'G', 'U']
            possible_nucleotides.remove(nucleotide)
            new_nucleotide = np.random.choice(possible_nucleotides)
            if self.verbose:
                print(f"Mutating 5' UTR: {seq[index]} -> {new_nucleotide}")
            # Strings are immutable - rebuild the string
            self._sequence['fivep_utr'] = seq[:index] + new_nucleotide + seq[index+1:]
            
        elif region == 'cds':
            iii = 0
            while True:
                if self.verbose:
                    print(f"Mutating CDS, trial {iii}")
                current_codon, new_codon, index = self.propose_codon_mutation()
                if current_codon != new_codon:
                    break
                iii += 1
            self._codons[index] = new_codon
            self._sequence['cds'] = ''.join(self._codons)

        elif region == '3p_utr':
            # Get mutable positions in 3' UTR (respects nucleotides_mutability)
            start_nt, end_nt = self.threep_utr_index_range
            # Convert to local index within 3' UTR
            threep_utr_len = len(self._sequence['threep_utr'])
            mutable_positions = np.where(self.nucleotides_mutability[start_nt:start_nt+threep_utr_len] > 0)[0]
            if len(mutable_positions) == 0:
                return  # No mutable positions
            index = np.random.choice(mutable_positions)
            seq = self._sequence['threep_utr']
            nucleotide = seq[index]
            possible_nucleotides = ['A', 'C', 'G', 'U']
            possible_nucleotides.remove(nucleotide)
            new_nucleotide = np.random.choice(possible_nucleotides)
            if self.verbose:
                print(f"Mutating 3' UTR: {seq[index]} -> {new_nucleotide}")
            # Strings are immutable - rebuild the string
            self._sequence['threep_utr'] = seq[:index] + new_nucleotide + seq[index+1:]

        # Invalidate all caches
        self._loss = None  # Critical: invalidate cached loss
        self._cached_seq_hash = None
        self._cached_mfe_seq_hash = None
        self._cached_bpp_seq_hash = None
        self._cached_structure_seq_hash = None
        self._cached_codons_string = None
        self._codons_changed = True
        
        return 
        

    def propose_codon_mutation(self):
        index_modifiable_codons = np.where(self.codon_mutability)[0]
        index = np.random.choice(index_modifiable_codons)
        current_codon = self.codons[index]
        current_amino_acid = codon_table[current_codon]
        possible_codons = list(self.aa_to_codon_cai[current_amino_acid].keys())
        probabilities = np.array(list(self.aa_to_codon_cai[current_amino_acid].values()))
        if self.verbose:
            print(f"Possible codons: {possible_codons}")
            print(f"Probabilities: {probabilities}")
        probabilities = probabilities / np.sum(probabilities)
        new_codon = np.random.choice(possible_codons, p = probabilities)
        if self.verbose:
            print(f"Proposed mutation: {current_codon} to {new_codon} at index {index}")
        return current_codon, new_codon, index


    @property
    def codons(self):
        return self._codons
    
    @codons.setter
    def codons(self, cds_sequence:str):
        codons = []
        for i in range(0, len(cds_sequence), 3):
            codons.append(cds_sequence[i:i+3])
        self._codons = np.array(codons, dtype=str)
        self._codons_changed = True  # Mark codons as changed
        self._cached_codons_string = None  # Invalidate cached string
        # Invalidate cached results when codons change
        if hasattr(self, '_pf_computed'):
            self._pf_computed = False
        return

    def calculate_CAI(self, form = 'log', normalise = False):
        if form == 'log':
            if normalise:
                result = self.calculate_CAI_log()/len(self.codons)
            else:
                result = self.calculate_CAI_log()
        elif form == 'linear':
            if normalise:
                result = np.exp(self.calculate_CAI_log()/len(self.codons))
            else:
                result = np.exp(self.calculate_CAI_log())
        else:
            raise ValueError("Invalid form")
        return result

    def calculate_CAI_log(self):
        """
        Vectorized calculation of the average log CAI of the mRNA sequence.
        """
        if self._cai_log is None or self._codons_changed:
            caidict = codon_adaptation_index.get(self._species, None)
        
            codon_arr = self.codons
            # Get CAI values for each codon (use np.vectorize for clarity)
            cai_vec = np.vectorize(lambda codon: caidict.get(codon, 0.0))
            cai_values = cai_vec(codon_arr)

            valid = ( cai_values > 0 ).all()
            if not valid:
                raise ValueError(f"Some codons have a CAI of 0: {codon_arr[cai_values <= 0]}")
            self._cai_log = np.log(cai_values).sum()
            # Note: _codons_changed is reset in codons_string property

        return self._cai_log

    def get_observed_codon_distribution(self) -> dict:
        """
        Calculate the observed codon frequency distribution per amino acid.
        
        Returns:
            Dict[amino_acid: Dict[codon: frequency]]
            Frequencies are normalized to sum to 1 for each amino acid.
        """
        # Count codons per amino acid
        codon_counts = {}  # {amino_acid: {codon: count}}
        
        for codon in self.codons:
            if codon not in codon_table:
                continue
            aa = codon_table[codon]
            if aa not in codon_counts:
                codon_counts[aa] = {}
            codon_counts[aa][codon] = codon_counts[aa].get(codon, 0) + 1
        
        # Normalize to frequencies
        observed_dist = {}
        for aa, counts in codon_counts.items():
            total = sum(counts.values())
            if total > 0:
                observed_dist[aa] = {codon: count / total for codon, count in counts.items()}
        
        return observed_dist

    def get_expected_codon_distribution(self) -> dict:
        """
        Get the expected codon probability distribution from CAI values.
        CAI values are normalized per amino acid to form a probability distribution.
        
        Returns:
            Dict[amino_acid: Dict[codon: probability]]
        """
        expected_dist = {}
        
        for aa, codon_cai in self.aa_to_codon_cai.items():
            # Normalize CAI values to probabilities
            total_cai = sum(codon_cai.values())
            if total_cai > 0:
                expected_dist[aa] = {codon: cai / total_cai for codon, cai in codon_cai.items()}
        
        return expected_dist

    @staticmethod
    def _jensen_shannon_divergence(p: np.ndarray, q: np.ndarray) -> float:
        """
        Calculate Jensen-Shannon divergence between two probability distributions.
        JS(p, q) = 0.5 * KL(p || m) + 0.5 * KL(q || m) where m = 0.5 * (p + q)
        
        This is symmetric and always finite (no division by zero issues).
        Uses log base 2 so the result is bounded between 0 and 1.
        
        Args:
            p: First probability distribution (numpy array)
            q: Second probability distribution (numpy array)
            
        Returns:
            Jensen-Shannon divergence (bounded between 0 and 1)
        """
        # Ensure valid probability distributions
        p = np.asarray(p, dtype=float)
        q = np.asarray(q, dtype=float)
        
        # Compute mixture distribution
        m = 0.5 * (p + q)
        
        # Compute KL divergences using safe log2 (0 * log(0) = 0)
        def safe_kl(a, b):
            """KL(a || b) with safe handling of zeros, using log base 2"""
            result = 0.0
            for ai, bi in zip(a, b):
                if ai > 0 and bi > 0:
                    result += ai * np.log2(ai / bi)
            return result
        
        js = 0.5 * safe_kl(p, m) + 0.5 * safe_kl(q, m)
        return js

    @property
    def codon_distribution_divergence(self) -> float:
        """
        Calculate the divergence between observed and expected codon distributions.
        Uses Jensen-Shannon divergence which is symmetric and always finite.
        
        The divergence is summed over all amino acids, weighted by the number of
        occurrences of each amino acid in the sequence.
        
        Returns:
            Total weighted JS divergence (lower = more similar to expected)
        """
        observed = self.get_observed_codon_distribution()
        expected = self.get_expected_codon_distribution()
        
        total_divergence = 0.0
        total_codons = 0
        
        for aa in observed:
            # Skip single-codon amino acids (no choice to optimize)
            if aa in ['Methionine', 'Tryptophan']:
                continue
            
            if aa not in expected:
                continue
            
            # Get all possible codons for this amino acid
            all_codons = list(expected[aa].keys())
            
            # Build probability vectors (same order for both)
            p_obs = np.array([observed[aa].get(codon, 0.0) for codon in all_codons])
            p_exp = np.array([expected[aa].get(codon, 0.0) for codon in all_codons])
            
            # Count of this amino acid in sequence (for weighting)
            aa_count = sum(1 for c in self.codons if codon_table.get(c) == aa)
            
            # Calculate JS divergence for this amino acid
            js_div = self._jensen_shannon_divergence(p_obs, p_exp)
            
            # Weight by amino acid frequency
            total_divergence += js_div * aa_count
            total_codons += aa_count
        
        # Normalize by total codons to get per-codon divergence
        if total_codons > 0:
            return total_divergence / total_codons
        return 0.0

    @property
    def mfe(self) -> float:
        # Check sequence-based cache
        seq = self.sequence
        seq_hash = hash(seq)
        
        if self._minimum_folding_energy is not None and \
           hasattr(self, '_cached_mfe_seq_hash') and \
           self._cached_mfe_seq_hash == seq_hash:
            return self._minimum_folding_energy
        
        if self._minimum_folding_energy is None or \
           (hasattr(self, '_cached_mfe_seq_hash') and self._cached_mfe_seq_hash != seq_hash):
            try:
                if not seq or len(seq) == 0:
                    raise ValueError("Empty sequence for MFE calculation")
                
                # Use LinearPartition to get MEA structure and ensemble energy
                # For MFE, we use the ensemble free energy as an approximation
                # (LinearPartition doesn't directly provide MFE structure)
                mea_structure, ensemble_energy = self._linearpartition.calculate_mea_structure(seq)
                
                # Store structure (this will also cache it via the structure property)
                self._structure = mea_structure
                self._cached_structure_seq_hash = seq_hash
                
                # Use ensemble free energy as MFE approximation
                # Convert from kcal/mol to normalized units
                self._minimum_folding_energy = ensemble_energy / self.RT 
                
                # Validate structure length
                if len(self._structure) != len(seq):
                    raise ValueError(f"Structure length {len(self._structure)} doesn't match sequence length {len(seq)}")
                
                # Validate structure contains only valid characters
                valid_chars = set('().')
                if not all(c in valid_chars for c in self._structure):
                    invalid_chars = set(self._structure) - valid_chars
                    print(f"WARNING: Structure from LinearPartition contains invalid characters: {invalid_chars}")
                    # Clean the structure by replacing invalid characters with '.'
                    self._structure = ''.join(c if c in valid_chars else '.' for c in self._structure)
                
                # Cache sequence hash
                self._cached_mfe_seq_hash = seq_hash
                
            except Exception as e:
                print(f"Error computing MFE: {e}")
                import traceback
                traceback.print_exc()
                raise
        return self._minimum_folding_energy

    @property
    def free_energy(self) -> float:
        """
        Calculate the free energy of the mRNA sequence using LinearPartition.
        This is the ensemble free energy computed over all possible secondary structures.
        
        This method ONLY computes the partition function. Structure and prob_matrix
        are computed separately when needed.
        
        Returns:
            float: Free energy in kcal/mol (normalized per codon)
        """
        # Check sequence-based cache
        seq = self.sequence
        seq_hash = hash(seq)
        
        if self._free_energy is not None and hasattr(self, '_cached_seq_hash') and self._cached_seq_hash == seq_hash:
            return self._free_energy
        
        if self._free_energy is None or (hasattr(self, '_cached_seq_hash') and self._cached_seq_hash != seq_hash):
            try:
                # Validate sequence before calling LinearPartition
                if not seq or len(seq) == 0:
                    raise ValueError("Empty sequence")
                
                # Only compute partition function if not already computed
                if not self._pf_computed:
                    ensemble_energy = self._linearpartition.calculate_partition_function(seq)
                    # Convert from kcal/mol to normalized units
                    self._free_energy = ensemble_energy / self.RT 
                    self._pf_computed = True
                else:
                    # If pf was already computed, we need to recompute it (shouldn't happen often)
                    ensemble_energy = self._linearpartition.calculate_partition_function(seq)
                    self._free_energy = ensemble_energy / self.RT 
                
                # Cache sequence hash
                self._cached_seq_hash = seq_hash
                
            except Exception as e:
                print(f"Error computing free energy: {e}")
                import traceback
                traceback.print_exc()
                # Set to None to allow retry
                self._free_energy = None
                if hasattr(self, '_cached_seq_hash'):
                    self._cached_seq_hash = None
                self._pf_computed = False
                raise
        return self._free_energy
    
    @property
    def structure(self) -> str:
        """
        Get the RNA secondary structure in dot-bracket notation.
        Computes MEA structure using LinearPartition when needed.
        
        Returns:
            str: Structure in dot-bracket notation
        """
        # Check sequence-based cache
        seq = self.sequence
        seq_hash = hash(seq)
        
        if self._structure is not None and \
           hasattr(self, '_cached_structure_seq_hash') and \
           self._cached_structure_seq_hash == seq_hash:
            return self._structure
        
        if self._structure is None or \
           (hasattr(self, '_cached_structure_seq_hash') and self._cached_structure_seq_hash != seq_hash):
            try:
                if not seq or len(seq) == 0:
                    raise ValueError("Empty sequence for structure calculation")
                
                # Get structure from LinearPartition (MEA structure)
                self._structure, _ = self._linearpartition.calculate_mea_structure(seq)
                
                # Validate structure length matches sequence length
                if len(self._structure) != len(seq):
                    raise ValueError(f"Structure length {len(self._structure)} doesn't match sequence length {len(seq)}")
                
                # Validate structure contains only valid characters
                valid_chars = set('().')
                if not all(c in valid_chars for c in self._structure):
                    invalid_chars = set(self._structure) - valid_chars
                    print(f"WARNING: Structure from LinearPartition contains invalid characters: {invalid_chars}")
                    # Clean the structure by replacing invalid characters with '.'
                    self._structure = ''.join(c if c in valid_chars else '.' for c in self._structure)
                
                # Cache sequence hash
                self._cached_structure_seq_hash = seq_hash
                
            except Exception as e:
                print(f"Error computing structure: {e}")
                import traceback
                traceback.print_exc()
                raise
        
        return self._structure
    
    @property
    def prob_matrix(self) -> np.ndarray:
        """
        Get the base pair probability matrix.
        Computes BPP matrix using LinearPartition when needed.
        Requires partition function to be computed first.
        
        Returns:
            np.ndarray: Base pair probability matrix (n x n)
        """
        # Check sequence-based cache
        seq = self.sequence
        seq_hash = hash(seq)
        
        if self._prob_matrix is not None and \
           hasattr(self, '_cached_bpp_seq_hash') and \
           self._cached_bpp_seq_hash == seq_hash:
            return self._prob_matrix
        
        if self._prob_matrix is None or \
           (hasattr(self, '_cached_bpp_seq_hash') and self._cached_bpp_seq_hash != seq_hash):
            try:
                if not seq or len(seq) == 0:
                    raise ValueError("Empty sequence for prob_matrix calculation")

                if not self._pf_computed:
                    _ = self.free_energy  # This will compute partition function and set _pf_computed
                
                # Get base pair probability matrix using LinearPartition
                # Use cutoff to speed up calculation - only compute probabilities above threshold
                # Note: For hybridisation_penalty, we need all probabilities, so cutoff should be 0.0
                # But we can use a small cutoff to speed up without losing much accuracy
                self._prob_matrix = self._linearpartition.calculate_bpp_matrix(seq, cutoff=self._bpp_cutoff)
                
                # Validate matrix
                if self._prob_matrix is None:
                    raise ValueError("LinearPartition bpp calculation returned None")
                
                # Process the matrix
                if self._prob_matrix.ndim != 2:
                    print(f"Warning: Expected 2D probability matrix, got {self._prob_matrix.ndim}D")
                    self._prob_matrix = None
                else:
                    # LinearPartition returns n x n matrix (0-indexed)
                    seq_len = len(seq)
                    if self._prob_matrix.shape[0] != seq_len or self._prob_matrix.shape[1] != seq_len:
                        print(f"Warning: prob_matrix shape {self._prob_matrix.shape} doesn't match expected size (seq_len={seq_len})")
                
                # Cache sequence hash
                self._cached_bpp_seq_hash = seq_hash
                
            except Exception as e:
                print(f"Error getting base pair probability matrix: {e}")
                import traceback
                traceback.print_exc()
                self._prob_matrix = None
                raise
        
        return self._prob_matrix
    
    
    def reset(self, what=None):
        """
        Reset cached values. If 'what' is None or 'all', resets everything.
        Otherwise, 'what' can be a list of strings: ['energy', 'structure', 'bpp', 'cai', 'cpg', 'stem']
        """
        # Reset everything (original behavior)
        self._structure = None
        self._minimum_folding_energy = None
        self._free_energy = None
        self._prob_matrix = None
        self._cai_log = None
        self._largest_stem = None
        self._cpg_count = None
        self._loss = None  # Critical: invalidate cached loss
        self._cached_average_stem_length = None
        if hasattr(self, '_cached_seq_hash'):
            self._cached_seq_hash = None
        if hasattr(self, '_cached_mfe_seq_hash'):
            self._cached_mfe_seq_hash = None
        if hasattr(self, '_cached_bpp_seq_hash'):
            self._cached_bpp_seq_hash = None
        if hasattr(self, '_cached_structure_seq_hash'):
            self._cached_structure_seq_hash = None
        # Reset partition function computation flag
        self._pf_computed = False
        
        return

    def save_state(self) -> dict:
        """
        Save current mutable state for lightweight backup.
        Much faster than copy.deepcopy() for optimization loops.
        
        Returns:
            dict: State dictionary that can be passed to restore_state()
        """
        return {
            'codons': self._codons.copy(),
            'sequence': dict(self._sequence),  # shallow copy of dict with string values
            'loss': self._loss,
            'cai_log': self._cai_log,
            'cpg_count': self._cpg_count,
        }

    def restore_state(self, state: dict):
        """
        Restore mutable state from a saved state dictionary.
        
        Args:
            state: State dictionary from save_state()
        """
        self._codons = state['codons']
        self._sequence = state['sequence']
        self._loss = state['loss']
        self._cai_log = state['cai_log']
        self._cpg_count = state['cpg_count']
        self._codons_changed = True
        # Reset computed values that depend on sequence
        self.reset()

    @property
    def cpg_count(self) -> int:
        """
        Count the number of CpG dinucleotides in the codon sequence.
        CpG refers to a cytosine (C) followed by a guanine (G).
        
        Returns:
            int: Number of CpG dinucleotides in the sequence
        """
        if self._cpg_count is None or self._codons_changed:
            seq = self.sequence
            # Use optimized C string method instead of Python loop
            self._cpg_count = seq.count('CG')
        return self._cpg_count

    @property
    def loss(self) -> float:
        if self._loss is not None:
            return self._loss

        loss = 0.0
        w_cai = self._loss_weights.get('cai', 0.0)
        if w_cai != 0.0:
            loss += -w_cai * self.calculate_CAI_log()
        w_mfe = self._loss_weights.get('mfe', 0.0)
        if w_mfe != 0.0:
            loss += w_mfe * self.mfe
        w_fe = self._loss_weights.get('fe', 0.0)
        if w_fe != 0.0:
            loss += w_fe * self.free_energy
        w_cpg = self._loss_weights.get('cpg', 0.0)
        if w_cpg != 0.0:
            loss += w_cpg * self.cpg_count
        w_stem = self._loss_weights.get('stem', 0.0)
        if w_stem != 0.0:
            loss += w_stem * self.stem_penalty 
        w_fivep_utr = self._loss_weights.get('fivep_utr_hybridisation', 0.0)
        if w_fivep_utr != 0.0:
            loss += w_fivep_utr * self.fivep_utr_hybridisation_penalty
        w_threep_utr = self._loss_weights.get('threep_utr_hybridisation', 0.0)
        if w_threep_utr != 0.0:
            loss += w_threep_utr * self.threep_utr_hybridisation_penalty
        w_hairpin = self._loss_weights.get('initial_hybridisation', 0.0)
        if w_hairpin != 0.0:
            loss += w_hairpin * self.initial_hybridisation_penalty
        w_codon_div = self._loss_weights.get('codon_divergence', 0.0)
        if w_codon_div != 0.0:
            loss += w_codon_div * self.codon_distribution_divergence
        self._loss = loss
        return self._loss

    def visualize_structure(self, filename="structure.svg", format="svg"):
        print(f"Visualizing structure with minimum energy in {filename}")
        try:
            seq = self.sequence
            if not seq or len(seq) == 0:
                raise ValueError("Empty sequence for visualization")
            
            # Use structure property (will compute if needed)
            structure = self.structure
            
            # Validate structure length matches sequence length
            if len(structure) != len(seq):
                raise ValueError(f"Structure length {len(structure)} doesn't match sequence length {len(seq)}")
            
            # Note: LinearPartition doesn't provide visualization functions
            # This method is kept for compatibility but visualization would need
            # to be done using external tools (e.g., ViennaRNA, VARNA, etc.)
            if format == "svg" or format == "eps":
                print(f"Warning: Visualization not directly supported by LinearPartition.")
                print(f"Structure saved to {filename}.txt for use with external visualization tools.")
                with open(filename + ".txt", 'w') as f:
                    f.write(f"{seq}\n{structure}\n")
            else:
                raise ValueError("Format must be 'svg' or 'eps'")
        except Exception as e:
            print(f"Error visualizing structure: {e}")
            import traceback
            traceback.print_exc()
            raise

    @property
    def fivep_utr_hybridisation_penalty(self):
        return self.hybridisation_penalty(self._fivep_utr_index_range[0], self._fivep_utr_index_range[1])

    @property
    def threep_utr_hybridisation_penalty(self):
        return self.hybridisation_penalty(self._threep_utr_index_range[0], self._threep_utr_index_range[1])

    @property
    def initial_hybridisation_penalty(self):
        if self.initial_region_end_index is None:
            print(f"Initial region end index not set")
            raise ValueError("Initial region end index not set but weight not equal to 0")
        return self.hybridisation_penalty(0, self.initial_region_end_index)

    def hybridisation_penalty(self, start_index, end_index ):
        if start_index == end_index:
            return 0.0
        
        # Use prob_matrix property (will compute if needed)
        try:
            prob_matrix = self.prob_matrix
        except Exception as e:
            print(f"Warning: Could not compute prob_matrix for hybridisation_penalty: {e}")
            return 0.0
        
        # Safety check
        if prob_matrix is None or not hasattr(prob_matrix, 'shape'):
            return 0.0
        
        # Bounds checking
        seq_len = len(self.sequence)
        if start_index < 0 or end_index >= seq_len or start_index >= seq_len:
            print(f"Warning: Invalid indices for hybridisation_penalty: start={start_index}, end={end_index}, seq_len={seq_len}")
            return 0.0
        
        # Ensure indices are within prob_matrix bounds
        matrix_size = prob_matrix.shape[0]
        start_index = max(0, min(start_index, matrix_size - 1))
        end_index = max(0, min(end_index, matrix_size - 1))
        
        if start_index > end_index:
            start_index, end_index = end_index, start_index
        
        try:
            # Only extract the relevant region of the probability matrix
            # This is more memory efficient than slicing the entire matrix
            region_size = end_index - start_index + 1
            if region_size <= 0:
                return 0.0
            
            # Extract the region: rows from start_index to end_index, all columns
            prob_region = prob_matrix[start_index:end_index+1, :]
            
            # Penalize the hybridisation of the region
            # Sum all probabilities in the region and normalize by region size
            penalty = np.sum(prob_region) / max(region_size, 1)
            return penalty
        except Exception as e:
            print(f"Warning: Error computing hybridisation_penalty: {e}")
            return 0.0

    @property
    def largest_stem_length(self):
        if self._largest_stem is None:
            self.set_largest_stem()
        start, end, length, stem_struct, stem_seq = self._largest_stem
        return length

    @property
    def largest_stem(self):
        if self._largest_stem is None:
            self.set_largest_stem()
        return self._largest_stem

    def set_largest_stem(self):
        """
        Find the largest (contiguous) double-stranded region in dot-bracket notation.
        
        Args:
            structure: String in dot-bracket notation. If None, uses structure property to get structure.
        
        Returns:
            tuple: (start_pos, end_pos, length, stem_structure, stem_sequence)
                   where start_pos and end_pos are the positions of the stem
        """
        # Use structure property (will compute if needed)
        structure = self.structure
        n = len(structure)
        
        # Find all base pairs
        pairs = {}  # position -> paired position
        stack = []
        
        for i, char in enumerate(structure):
            if char == '(':
                stack.append(i)
            elif char == ')':
                if stack:
                    j = stack.pop()
                    pairs[j] = i
                    pairs[i] = j
        
        # Find contiguous stems
        # A stem is a contiguous region where every position is paired
        # and pairs are adjacent (i pairs with j, i+1 pairs with j-1, etc.)
        
        visited = set()
        best_stem = None
        best_length = 0
        
        for start in range(n):
            if start in visited or structure[start] != '(':
                continue
                
            # Try to extend a stem starting from this position
            i = start
            stem_pairs = []
            
            while i < n and i in pairs:
                j = pairs[i]
                if j < i:  # j should be to the right
                    break
                    
                # Check if this forms a contiguous stem
                # (i, j) should be adjacent to previous pair (i-1, j+1)
                if stem_pairs:
                    prev_i, prev_j = stem_pairs[-1]
                    if i != prev_i + 1 or j != prev_j - 1:
                        break
                
                stem_pairs.append((i, j))
                visited.add(i)
                visited.add(j)
                i += 1
            
            if len(stem_pairs) > best_length:
                best_length = len(stem_pairs)
                best_stem = stem_pairs
        
        if best_stem and best_length > 0:
            start_pos = best_stem[0][0]
            end_pos = best_stem[-1][1]
            
            # Create stem structure representation
            stem_struct = list('.' * n)
            for i, j in best_stem:
                stem_struct[i] = '('
                stem_struct[j] = ')'
            stem_struct = ''.join(stem_struct)
            
            # Get the sequence of the stem
            seq = self.sequence
            stem_seq_5prime = seq[start_pos:start_pos + best_length]
            stem_seq_3prime = seq[end_pos - best_length + 1:end_pos + 1]

            self._largest_stem = (start_pos, end_pos, best_length, stem_struct, (stem_seq_5prime, stem_seq_3prime))
            return self._largest_stem
        
        self._largest_stem = (0, 0, 0, '.' * n, ('', ''))
        return self._largest_stem

    @property
    def stem_penalty(self):
        """
        Penalty associated with the stem of the mRNA sequence.
        
        Can use either:
        - Largest stem length from MEA structure (default, faster)
        """
        length = self.largest_stem_length
        
        return abs( min( 30 - length, 0 ) ) # Penalize stems longer than 30 bp

    def restriction_sites_count(self):
        """
        Count the number of restriction sites in the mRNA sequence.
        
        Returns:
            int: Total number of restriction sites found
        """
        curr_sequence = get_sequence_as_dna(self.sequence)
        total_count = 0
    
        for enzyme_name, site_seq in restriction_sites.items():
            # Convert T to U for RNA sequences if needed
            site_seq = get_sequence_as_dna(site_seq)
            count = curr_sequence.count(site_seq)
            total_count += count
        
        return total_count

    def calculate_hybridised_contacts_percentage(self):
        """
        Calculate the percentage of hybridised contacts.
        This is the sum of all probabilities in prob_matrix divided by the number of nucleotides.
        
        Returns:
            float: Percentage of hybridised contacts (0-100)
        """
        # Use prob_matrix property (will compute if needed)
        try:
            prob_matrix = self.prob_matrix
        except Exception as e:
            print(f"Warning: Could not compute prob_matrix for hybridised contacts: {e}")
            return 0.0
        
        # Safety check
        if prob_matrix is None or not hasattr(prob_matrix, 'shape'):
            return 0.0
        
        try:
            # Sum all probabilities in the matrix
            total_prob = np.sum(prob_matrix)
            # Number of nucleotides = 3 * number of codons
            n_nucleotides = 3 * len(self.codons)
            
            if n_nucleotides == 0:
                return 0.0
            
            # Return as percentage
            return (total_prob / n_nucleotides) * 100.0
        except Exception as e:
            print(f"Warning: Error calculating hybridised contacts percentage: {e}")
            return 0.0

    def calculate_sequence_identity(self):
        """
        Calculate sequence identity (codon-wise) compared to initial codon sequence.
        
        Returns:
            float: Sequence identity as percentage (0-100)
        """
        current_codons = self.codons
        initial_codons = self._initial_codons
        
        if len(current_codons) != len(initial_codons):
            raise ValueError("Current and initial codon sequences have different lengths")
        
        # Count identical codons
        identical = np.sum(current_codons == initial_codons)
        total = len(current_codons)
        
        # Return as percentage
        return (identical / total) * 100.0

    def save_statistics(self, step, acceptance_rate, output_file, loss_components=None):
        """
        Save statistics to output file.
        
        Args:
            step: Current step number
            acceptance_rate: Acceptance rate over last n_sample steps
            output_file: File object or path to write statistics
            loss_components: Optional dict of pre-computed loss component values.
                           If None, will compute from cached properties (slower).
        """
        #hybridised_pct = self.calculate_hybridised_contacts_percentage()
        sequence_identity = self.calculate_sequence_identity()
        
        # Collect loss component values for non-zero weights
        loss_values = []
        
        # Use pre-computed values if provided, otherwise access cached properties
        if loss_components is None:
            loss_components = {}
        
        # Access cached properties - they should be fast since loss was just computed
        # Properties use sequence-based caching, so accessing them after loss computation is efficient
        w_cai = self._loss_weights.get('cai', 0.0)
        if w_cai != 0.0:
            # For CAI, we actually recalculate it to store the linear form
            cai_value = self.calculate_CAI(form="linear",normalise=True)
            loss_values.append(cai_value)
        
        w_mfe = self._loss_weights.get('mfe', 0.0)
        if w_mfe != 0.0:
            if 'mfe' in loss_components:
                mfe_value = loss_components['mfe'] / self.n_nucleotides
            else:
                try:
                    mfe_value = self.mfe / self.n_nucleotides
                except Exception:
                    mfe_value = 0.0
            loss_values.append(mfe_value)
        
        w_fe = self._loss_weights.get('fe', 0.0)
        if w_fe != 0.0:
            if 'fe' in loss_components:
                fe_value = loss_components['fe'] / self.n_nucleotides
            else:
                try:
                    fe_value = self.free_energy / self.n_nucleotides
                except Exception:
                    fe_value = 0.0
            loss_values.append(fe_value)
        
        w_cpg = self._loss_weights.get('cpg', 0.0)
        if w_cpg != 0.0:
            if 'cpg' in loss_components:
                cpg_value = loss_components['cpg']
            else:
                cpg_value = self.cpg_count
            loss_values.append(cpg_value)
        
        w_stem = self._loss_weights.get('stem', 0.0)
        if w_stem != 0.0:
            if 'stem' in loss_components:
                stem_value = loss_components['stem']
            else:
                try:
                    stem_value = self.stem_penalty
                except Exception:
                    stem_value = 0.0
            loss_values.append(stem_value)
        
        w_5p_utr = self._loss_weights.get('fivep_utr_hybridisation', 0.0)
        if w_5p_utr != 0.0:
            if 'fivep_utr_hybridisation' in loss_components:
                fivep_utr_value = loss_components['fivep_utr_hybridisation']
            else:
                try:
                    fivep_utr_value = self.fivep_utr_hybridisation_penalty
                except Exception:
                    fivep_utr_value = 0.0
            loss_values.append(fivep_utr_value)

        w_3p_utr = self._loss_weights.get('threep_utr_hybridisation', 0.0)
        if w_3p_utr != 0.0:
            if 'threep_utr_hybridisation' in loss_components:
                threep_utr_value = loss_components['threep_utr_hybridisation']
            else:
                try:
                    threep_utr_value = self.threep_utr_hybridisation_penalty
                except Exception:
                    threep_utr_value = 0.0
            loss_values.append(threep_utr_value)

        w_hairpin = self._loss_weights.get('initial_hybridisation', 0.0)
        if w_hairpin != 0.0:
            if 'initial_hybridisation' in loss_components:
                hairpin_value = loss_components['initial_hybridisation']
            else:
                try:
                    hairpin_value = self.initial_hybridisation_penalty
                except (ValueError, AttributeError):
                    hairpin_value = 0.0
            loss_values.append(hairpin_value)
        
        w_restriction = self._loss_weights.get('restriction_sites', 0.0)
        if w_restriction != 0.0:
            if 'restriction_sites' in loss_components:
                restriction_value = loss_components['restriction_sites']
            else:
                try:
                    restriction_value = self.restriction_sites_count()
                except Exception:
                    restriction_value = 0.0
            loss_values.append(restriction_value)
        
        w_codon_divergence = self._loss_weights.get('codon_divergence', 0.0)
        if w_codon_divergence != 0.0:
            if 'codon_divergence' in loss_components:
                codon_divergence_value = loss_components['codon_divergence']
            else:
                try:
                    codon_divergence_value = self.codon_distribution_divergence
                except Exception:
                    codon_divergence_value = 0.0
            loss_values.append(codon_divergence_value)
        
        # Format with fixed-width columns (24 chars each) for alignment with header
        col_width = 24
        parts = [f"{step:<{col_width}}", f"{acceptance_rate:<{col_width}.4f}", f"{sequence_identity:<{col_width}.4f}"]
        for val in loss_values:
            parts.append(f"{val:<{col_width}.6e}")
        line = "".join(parts) + "\n"
        
        if isinstance(output_file, str):
            # If it's a file path, open in append mode
            with open(output_file, 'a+') as f:
                f.write(line)
        else:
            # If it's already a file object, write to it
            output_file.write(line)
            output_file.flush()
        
        # Save codon sequence to sequence_sampled.txt
        sequence_filename = f"cds_sequence_sampled.txt"
        try:
            with open(sequence_filename, 'a') as f:
                # Write the codon sequence as space-separated codons
                codon_sequence = ' '.join(self.codons)
                start_codon = self.fivep_utr_sequence[-3:]
                f.write( f"{step}: " + start_codon + " " + codon_sequence + '\n')
        except Exception as e:
            print(f"Warning: Could not save codon sequence to {sequence_filename}: {e}")

        # Save codon sequence to sequence_sampled.txt
        sequence_filename = f"fivep_utr_sequence_sampled.txt"
        try:
            with open(sequence_filename, 'a') as f:
                # Write the codon sequence as space-separated codons
                f.write( f"{step}: " + self.fivep_utr_sequence[:-3] + '\n')
        except Exception as e:
            print(f"Warning: Could not save codon sequence to {sequence_filename}: {e}")
        # Save codon sequence to sequence_sampled.txt
        sequence_filename = f"threep_utr_sequence_sampled.txt"
        try:
            with open(sequence_filename, 'a') as f:
                # Write the codon sequence as space-separated codons
                f.write( f"{step}: " + self.threep_utr_sequence + '\n')
        except Exception as e:
            print(f"Warning: Could not save codon sequence to {sequence_filename}: {e}")

    def save_prob_matrix(self, step):
        """
        Save the probability matrix to a file named prob_matrix.{step}.
        Format: i, j, prob_matrix(i,j) for each entry.
        
        Args:
            step: Current optimization step number
        """
        # Use prob_matrix property (will compute if needed)
        try:
            prob_matrix = self.prob_matrix
        except Exception as e:
            print(f"Warning: Could not compute prob_matrix at step {step}: {e}")
            return
        
        # Safety check
        if prob_matrix is None or not hasattr(prob_matrix, 'shape'):
            print(f"Warning: prob_matrix not available at step {step}, skipping save")
            return
        
        filename = f"prob_matrix.{step}"
        
        try:
            with open(filename, 'w') as f:
                # Write header
                f.write("# i, j, probability\n")
                
                # Write all entries (i, j, prob_matrix(i,j))
                n = prob_matrix.shape[0]
                m = prob_matrix.shape[1] if len(prob_matrix.shape) > 1 else n
                for i in range(n):
                    for j in range(m):
                        prob = prob_matrix[i, j]
                        if prob > 0:  # Only write non-zero probabilities
                            f.write(f"{i}, {j}, {prob:.6e}\n")
        except Exception as e:
            print(f"Error saving prob_matrix at step {step}: {e}")
            return

    def save_structure(self, step):
        """
        Save the RNA secondary structure in dot-bracket notation to a file named structure.{step}.
        
        Args:
            step: Current optimization step number
        """
        # Use structure property (will compute if needed)
        structure = self.structure
        
        filename = f"structure.{step}"
        
        with open(filename, 'w') as f:
            # Write the structure in dot-bracket notation
            f.write(f"{structure}\n")

    def save_structure_svg(self, step):
        """
        Save the RNA secondary structure in SVG format to a file named structure.{step}.svg.
        
        Note: RNA.svg_rna_plot() may cause segmentation faults in some ViennaRNA versions.
        If this occurs, consider using EPS format instead.
        
        Args:
            step: Current optimization step number
        """
        try:
            seq = self.sequence
            if not seq or len(seq) == 0:
                print(f"Warning: Empty sequence at step {step}, skipping structure plot save")
                return
            
            # Use structure property (will compute if needed)
            try:
                structure = self.structure
            except Exception as e:
                print(f"Warning: Could not compute structure for plot at step {step}: {e}")
                return
            
            # Validate structure exists and matches sequence length
            if structure is None:
                print(f"Warning: Structure not available at step {step}, skipping structure plot save")
                return
            
            if len(structure) != len(seq):
                print(f"Warning: Structure length {len(structure)} doesn't match sequence length {len(seq)} at step {step}, skipping structure plot save")
                return
            
            # Validate structure contains only valid characters
            valid_chars = set('().')
            if not all(c in valid_chars for c in structure):
                print(f"Warning: Invalid characters in structure at step {step}, skipping structure plot save")
                return
            
            filename = f"structure.{step}.svg"
            
            # Note: LinearPartition doesn't provide visualization functions
            # Save structure to text file for use with external visualization tools
            try:
                with open(filename + ".txt", 'w') as f:
                    f.write(f"{seq}\n{structure}\n")
            except Exception as e:
                print(f"Warning: Could not save structure plot at step {step}: {e}")
                return
        except Exception as e:
            print(f"Error saving structure plot at step {step}: {e}")
            import traceback
            traceback.print_exc()
            return

    def save_loss_components(self, step):
        """
        Save the different components of the loss function to a file named loss_components.{step}.
        Only saves components where the loss weight is not 0.0.
        Format: component_name, component_value
        
        Args:
            step: Current optimization step number
        """
        filename = f"loss_components.{step}"
        
        with open(filename, 'w') as f:
            # Write header
            f.write("# component_name, component_value\n")
            
            # Check each component and save if weight is not 0.0
            w_cai = self._loss_weights.get('cai', 0.0)
            if w_cai != 0.0:
                cai_value = np.exp(self.calculate_CAI_log())
                f.write(f"cai, {cai_value:.6e}\n")
            
            w_mfe = self._loss_weights.get('mfe', 0.0)
            if w_mfe != 0.0:
                try:
                    mfe_value = self.mfe
                    f.write(f"mfe, {mfe_value:.6e}\n")
                except Exception as e:
                    print(f"Warning: Could not compute mfe for loss components: {e}")
                    f.write(f"mfe, 0.000000e+00\n")
            
            w_fe = self._loss_weights.get('fe', 0.0)
            if w_fe != 0.0:
                try:
                    fe_value = self.free_energy
                    f.write(f"fe, {fe_value:.6e}\n")
                except Exception as e:
                    print(f"Warning: Could not compute free_energy for loss components: {e}")
                    f.write(f"fe, 0.000000e+00\n")
            
            w_cpg = self._loss_weights.get('cpg', 0.0)
            if w_cpg != 0.0:
                cpg_value = self.cpg_count
                f.write(f"cpg, {cpg_value:.6e}\n")
            
            w_stem = self._loss_weights.get('stem', 0.0)
            if w_stem != 0.0:
                try:
                    stem_value = self.stem_penalty
                    f.write(f"stem, {stem_value:.6e}\n")
                except Exception as e:
                    print(f"Warning: Could not compute stem_penalty for loss components: {e}")
                    f.write(f"stem, 0.000000e+00\n")
            
            w_fivep_utr = self._loss_weights.get('fivep_utr_hybridisation', 0.0)
            if w_fivep_utr != 0.0:
                try:
                    fivep_utr_value = self.fivep_utr_hybridisation_penalty
                    f.write(f"fivep_utr_hybridisation, {fivep_utr_value:.6e}\n")
                except Exception as e:
                    print(f"Warning: Could not compute fivep_utr_hybridisation_penalty for loss components: {e}")
                    f.write(f"fivep_utr_hybridisation, 0.000000e+00\n")
            
            w_hairpin = self._loss_weights.get('initial_hybridisation', 0.0)
            if w_hairpin != 0.0:
                try:
                    hairpin_value = self.initial_hybridisation_penalty
                    f.write(f"initial_hybridisation, {hairpin_value:.6e}\n")
                except (ValueError, AttributeError):
                    # Skip if initial_region_end_index is not set
                    pass
            
            w_restriction = self._loss_weights.get('restriction_sites', 0.0)
            if w_restriction != 0.0:
                # Note: restriction_sites component may not be implemented yet
                # If it exists as a property, add it here
                if hasattr(self, 'restriction_sites'):
                    restriction_value = self.restriction_sites
                    f.write(f"restriction_sites, {restriction_value:.6e}\n")
    
    def sample_structures_from_ensemble(self, num_samples=1, use_unique_ml=True):
        """
        Sample secondary structures from the Boltzmann ensemble.
        
        Note: LinearPartition doesn't directly support structure sampling.
        This method returns MEA structures as an approximation.
        
        Args:
            num_samples (int): Number of structures to sample. Default is 1.
            use_unique_ml (bool): Not used with LinearPartition, kept for compatibility.
        
        Returns:
            list or str: If num_samples == 1, returns a single structure string.
                        If num_samples > 1, returns a list of structure strings.
                        Returns None if sampling fails.
        
        Raises:
            ValueError: If sequence is empty or partition function cannot be computed.
        """
        try:
            # Validate sequence
            seq = self.sequence
            if not seq or len(seq) == 0:
                raise ValueError("Empty sequence for structure sampling")
            
            # LinearPartition doesn't support structure sampling directly
            # We return MEA structures as an approximation
            if num_samples == 1:
                structure, _ = self._linearpartition.calculate_mea_structure(seq)
                return structure
            else:
                # Return multiple MEA structures (they'll be the same, but this maintains API compatibility)
                structures = []
                for _ in range(num_samples):
                    structure, _ = self._linearpartition.calculate_mea_structure(seq)
                    structures.append(structure)
                return structures
                
        except Exception as e:
            print(f"Error sampling structures from ensemble: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    @staticmethod
    def _calculate_largest_stem_length(structure):
        """
        Helper method to calculate the length of the largest contiguous stem in a structure.
        
        Args:
            structure (str): RNA structure in dot-bracket notation
        
        Returns:
            int: Length of the largest contiguous stem (number of base pairs)
        """
        if not structure or len(structure) == 0:
            return 0
        
        n = len(structure)
        
        # Find all base pairs
        pairs = {}  # position -> paired position
        stack = []
        
        for i, char in enumerate(structure):
            if char == '(':
                stack.append(i)
            elif char == ')':
                if stack:
                    j = stack.pop()
                    pairs[j] = i
                    pairs[i] = j
        
        # Find contiguous stems
        # A stem is a contiguous region where every position is paired
        # and pairs are adjacent (i pairs with j, i+1 pairs with j-1, etc.)
        
        visited = set()
        best_length = 0
        
        for start in range(n):
            if start in visited or structure[start] != '(':
                continue
                
            # Try to extend a stem starting from this position
            i = start
            stem_pairs = []
            
            while i < n and i in pairs:
                j = pairs[i]
                if j < i:  # j should be to the right
                    break
                    
                # Check if this forms a contiguous stem
                # (i, j) should be adjacent to previous pair (i-1, j+1)
                if stem_pairs:
                    prev_i, prev_j = stem_pairs[-1]
                    if i != prev_i + 1 or j != prev_j - 1:
                        break
                
                stem_pairs.append((i, j))
                visited.add(i)
                visited.add(j)
                i += 1
            
            if len(stem_pairs) > best_length:
                best_length = len(stem_pairs)
        
        return best_length
    
    def optimize_codon_usage(
                            self, 
                            T_opt: float = 1.0, 
                            nsteps = 100, 
                            sample_frequency = 100, 
                            verbose = False,
                            output_filename: str = "opt_statistics.txt",
                            n_sample: int = None,
                            use_average_stem_length: bool = False,
                            average_stem_num_samples: int = 20
                            ):
        """
        Perform optimisation of the codon usage of the mRNA sequence.
        
        Args:
            T_opt: Initial temperature for simulated annealing
            nsteps: Number of optimization steps (multiplied by number of codons)
            sample_frequency: Frequency of printing progress
            verbose: Whether to print verbose output
            output_filename: Filename for statistics output (default: opt_statistics.txt)
            n_sample: Number of steps between statistics saves (default: same as sample_frequency)
            use_average_stem_length: If True, use average stem length from ensemble for stem penalty.
                                     If False, use longest stem from MFE structure (default, faster).
            average_stem_num_samples: Number of structures to sample for average stem calculation.
                                     Only used if use_average_stem_length=True. Default is 100.
        """
        # Perform optimisation of the codon usage of the mRNA sequence 
        nsteps = int(nsteps) * len(self.codons)
        T_schedule = T_opt * np.linspace(1, 0, nsteps)
        
        # Record start time for timing estimates
        start_time = time.time()
        first_n_sample_time = None
        
        # Set n_sample default to sample_frequency if not provided
        if n_sample is None:
            n_sample = sample_frequency
        
        # Set stem calculation method flags
        self._use_average_stem_length = use_average_stem_length
        self._average_stem_num_samples = average_stem_num_samples
        self._cached_average_stem_length = None  # Reset cache at start of optimization
        
        # Initialize statistics tracking
        # n_sample is the number of steps between statistics saves
        # Track (step_number, accepted) tuples for mutations in the last n_sample steps
        acceptance_history = []  # List of (step, accepted) tuples
        
        # Initialize statistics file with header
        # Build header with loss component columns for non-zero weights
        header_parts = ["# Step", "Acceptance_Rate", "Sequence_Identity_%"]
        
        # Add loss component names in the same order as save_statistics
        if self._loss_weights.get('cai', 0.0) != 0.0:
            header_parts.append("CAI")
        if self._loss_weights.get('mfe', 0.0) != 0.0:
            header_parts.append("MFE")
        if self._loss_weights.get('fe', 0.0) != 0.0:
            header_parts.append("Free_Energy")
        if self._loss_weights.get('cpg', 0.0) != 0.0:
            header_parts.append("CpG")
        if self._loss_weights.get('stem', 0.0) != 0.0:
            header_parts.append("Stem")
        if self._loss_weights.get('fivep_utr_hybridisation', 0.0) != 0.0:
            header_parts.append("Fivep_Utr_Hybridisation")
        if self._loss_weights.get('threep_utr_hybridisation', 0.0) != 0.0:
            header_parts.append("Threep_Utr_Hybridisation")
        if self._loss_weights.get('initial_hybridisation', 0.0) != 0.0:
            header_parts.append("Initial_Hybridisation")
        if self._loss_weights.get('restriction_sites', 0.0) != 0.0:
            header_parts.append("Restriction_Sites")
        if self._loss_weights.get('codon_divergence', 0.0) != 0.0:
            header_parts.append("Codon_Divergence")
        
        # Format header with fixed-width columns (24 chars each) to match data rows
        col_width = 24
        header_line = "".join([f"{h:<{col_width}}" for h in header_parts]) + "\n"
        with open(output_filename, 'w') as f:
            f.write(header_line)
        
        # Initialize loss and save initial state
        self._loss = None
        best_state = self.save_state()  # Lightweight state save instead of deepcopy
        best_loss = self.loss
        
        for i, T in enumerate( T_schedule ):
            # Save state before mutation
            prev_state = self.save_state()
            prev_loss = self.loss
            
            self.make_mutation()
            self.reset()
            accepted = False
            new_loss = self.loss
            delta = new_loss - prev_loss
            if delta < 0:
                # Always accept if loss decreases
                print(f"Change accepted, delta loss = {delta}")
                best_state = self.save_state()
                best_loss = new_loss
                accepted = True
            else:
                # Avoid division by zero or very small T
                if T > 1e-10:
                    accept_probability = np.exp(-delta/T)
                else:
                    # If T is very small, only accept if loss decreases
                    accept_probability = 0.0
                if np.random.rand() < accept_probability:
                    print(f"Change accepted, delta loss = {delta}")
                    # Accept with probability
                    best_state = self.save_state()
                    best_loss = new_loss
                    accepted = True
                else:
                    # Reject: restore previous state
                    print(f"Change rejected, delta loss = {delta}")
                    self.restore_state(prev_state)
                    if verbose:
                        print(f"Rejected mutation")
                
            # Track acceptance for statistics (with step number)
            acceptance_history.append((i + 1, accepted))
            
            # Save statistics every n_sample steps
            if i % n_sample == 0:
                # Calculate acceptance rate over last n_sample steps
                # Filter to only include mutations from the last n_sample steps
                current_step = i + 1
                recent_mutations = [(step, acc) for step, acc in acceptance_history 
                                   if step > current_step - n_sample]

                print(f"Step {i} of {nsteps} completed. Current T: {T}")
                energy = self.free_energy if self._pf_computed else self.mfe
                print(f"x nucleotide quantities:")
                n_nts = len( self.sequence )
                print(f"Current loss: {self.loss/n_nts}, energy: {energy/n_nts}, cai (x codon): {self.calculate_CAI(form = 'linear',normalise=True)}")
                
                if len(recent_mutations) > 0:
                    # Extract just the accepted values
                    accepted_values = [acc for _, acc in recent_mutations]
                    acceptance_rate = np.mean(accepted_values)
                else:
                    acceptance_rate = 0.0
                
                # Clean up old entries (keep only last n_sample * 2 to be safe)
                acceptance_history = [(step, acc) for step, acc in acceptance_history 
                                     if step > current_step - n_sample * 2]
                
                # Save statistics
                try:
                    self.save_statistics(i, acceptance_rate, output_filename)
                except Exception as e:
                    print(f"Error saving statistics at step {i}: {e}")
                
                # Calculate and display timing information
                current_time = time.time()
                elapsed_time = current_time - start_time
                elapsed_hours = elapsed_time / 3600.0
                
                # Record time after first n_sample steps to calculate rate
                if first_n_sample_time is None and current_step >= n_sample:
                    first_n_sample_time = current_time
                    time_per_n_sample = first_n_sample_time - start_time
                elif first_n_sample_time is not None:
                    # Use the rate from the first n_sample steps
                    time_per_n_sample = first_n_sample_time - start_time
                else:
                    # Not enough steps yet, estimate based on current progress
                    time_per_n_sample = elapsed_time / max(1, current_step)
                
                # Calculate estimated time to complete
                if current_step > 0:
                    steps_remaining = nsteps - current_step
                    # Estimate based on time per step from first n_sample steps
                    time_per_step = time_per_n_sample / n_sample
                    estimated_time_remaining = time_per_step * steps_remaining
                    estimated_time_remaining_hours = estimated_time_remaining / 3600.0
                    print(f"Time elapsed: {elapsed_hours:.4f} hours ({elapsed_time:.1f} seconds)")
                    if steps_remaining > 0:
                        print(f"Estimated time to complete: {estimated_time_remaining_hours:.4f} hours ({estimated_time_remaining:.1f} seconds)")
                else:
                    print(f"Time elapsed: {elapsed_hours:.4f} hours ({elapsed_time:.1f} seconds)")
                

            # Save probability matrix and structure less frequently (every sample_frequency * 10 steps)
            # This is done less often because these files can be large and expensive to compute
            if i % (sample_frequency * 10) == 0:
                try:
                    self.save_prob_matrix(i)
                except Exception as e:
                    print(f"Error saving prob_matrix at step {i}: {e}")

                try:
                    self.save_structure(i)
                except Exception as e:
                    print(f"Error saving structure at step {i}: {e}")
                
        # Save final statistics
        # Calculate acceptance rate over last n_sample steps
        recent_mutations = [(step, acc) for step, acc in acceptance_history 
                           if step > nsteps - n_sample]
        
        if len(recent_mutations) > 0:
            accepted_values = [acc for _, acc in recent_mutations]
            final_acceptance_rate = np.mean(accepted_values)
        else:
            final_acceptance_rate = 0.0
        self.save_statistics(nsteps, final_acceptance_rate, output_filename)
        # Save final probability matrix
        self.save_prob_matrix(nsteps)
        # Save final structure
        self.save_structure(nsteps)
        # Save final structure in SVG format
        self.save_structure_svg(nsteps)
        # Save final loss components
        self.save_loss_components(nsteps)
        
        print(f"Final loss: {self._loss}, mfe: {self.mfe}, cai: {np.exp(self.calculate_CAI_log())}")
        print(f"Statistics saved to {output_filename}")

    def start_from_optimal_cai(self):
        """
        For the initial sequence, start from the optimal CAI.
        Only change though the part that is mutable
        """
        # Replace each codon with the most frequent codon for the amino acid,
        # i.e. the one for which CAI is highest.
        print(f"Changing the sequence to start from optimal CAI")
        for i, cod in enumerate(self.codons):
            if self.codon_mutability[i]:
                amino_acid = codon_table[cod]  # Get full amino acid name (e.g., "Alanine", "Methionine")
                most_frequent_codon = max(self.aa_to_codon_cai[amino_acid], key=self.aa_to_codon_cai[amino_acid].get)
                self.codons[i] = most_frequent_codon
        return 

if __name__ == "__main__":
    sequence = "AUUAAAGGGUUUUAGCAGGGGCCCCUUAAGGCGGCGGAGGACGGG"
    sequence = sequence * 2 
    print(f"Sequence length: {len(sequence)}")
    assert len(sequence) % 3 == 0, "Sequence length must be divisible by 3"
    weights = {
                'mfe': 1.0, 
                'cai': 10.0, 
                #'fe': 1.0, 
                'cpg': 1.0, 
                #'stem': 10.0, 
                #'utr_hybridisation': 1.0, 
                #'initial_hybridisation': 1.0, 
                #'restriction_sites': 0.0 
                }
    verbose = False
    system = mRNA( 
                sequence = sequence, 
                species = 'human', 
                aa_to_codon_cai=human_aa_to_codon_cai,
                verbose = verbose, 
                loss_weights=weights, 
                modify_fivep_utr=False,
                modify_threep_utr=False,
                initial_region_end_index = 30,
                T_K = 310,
                average_stem_num_samples = 20
                )
    print(f"Number of nucleotides: {3*len(system.codons)}" )
    print(f"McCaskill free-energy (*nucleotide): {system.free_energy} / {3*len(system.codons)} kcal/mol per nucleotide")
    system.optimize_codon_usage( 
        T_opt=1.0, 
        nsteps=100, 
        sample_frequency = 100, 
        verbose = False,
        output_filename = "opt_statistics.txt"
        )
    print(f"Optimized mRNA: {system.sequence}")
    print(f"Optimized codon usage: {system._codons}")
    print(f"Initial amino acid sequence: {system.initial_aminoacid}")
    print(f"Optimized amino acid usage: {system.codons_to_amino_acids()}")
    print(f"Optimized loss: {system.loss}")
    print(f"Optimized mfe (*codon): {system.mfe} kcal/mol per codon")
    print(f"Optimized mfe (total): {system.mfe*len(system.codons)} kcal/mol")
    print(f"McCaskill free-energy (*codon): {system.free_energy} kcal/mol per codon")
    print(f"Optimized cai: {np.exp(system.calculate_CAI_log())}")
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
