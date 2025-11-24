# This code implements the class mRNA and its methods, used to optimize 
# the codon usage of an mRNA sequence. 

import numpy as np
from codons import *
from linearpartition_wrapper import LinearPartitionWrapper

class mRNA:
    def __init__(
        self, 
        sequence:str, 
        species:str, 
        aa_to_codon_cai:dict, 
        loss_weights:dict = {'mfe': 1.0, 'cai': 1.0},
        T_K = 310, 
        modify_utr = False,
        verbose = False,
        initial_region_end_index = None,
        beamsize = 100,
        bpp_cutoff = 0.0,
        ):
        # Retrieve codon adaptation index dictionary for the species
        assert species in ['human', 'mouse', 'rat', 'yeast', 'e_coli', 'other'], "Invalid species"

        assert len(sequence) % 3 == 0, "Sequence length must be divisible by 3"
        if 'T' in sequence:
            if 'U' in sequence:
                raise ValueError("Sequence contains both T and U, which is not allowed")
            sequence = sequence.replace('T', 'U')
            print(f"Sequence converted from T to U: {sequence}")        
        
        # Check that the sequence is divisible by 3
        assert len(sequence) % 3 == 0, "Sequence length must be divisible by 3"
        self.codons = sequence  # Sets ._codons and ._sequence
        self.start_index = None
        # Check that a start codon is present
        try:
            assert 'AUG' in self.codons, "Sequence must contain a start codon"
        except AssertionError:
            print( 'WARNING: Start codon added at beginning of sequence because it was not found' )
            sequence = 'AUG' + sequence
            self.codons = sequence
            self.start_index = 0

        # Check that a stop codon is present
        self.stop_index()

        if initial_region_end_index is not None:
            self.initial_region_end_index = initial_region_end_index
        # Check that the sequence contains only valid codons
        for codon in self.codons:
            assert codon in codon_table.keys(), f"Invalid codon: {codon}"

        self.initial_aminoacid = self.codons_to_amino_acids()
        self._initial_codons = self.codons.copy()  # Store initial codon sequence for statistics
        self._species = species
        self.aa_to_codon_cai = aa_to_codon_cai
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
        self._utr_index_range = None
        self._utr_sequence = None
        self.retrieve_utr() # This will initialize self.utr_index_range and self.utr_sequence
        self._largest_stem = None
        self.set_largest_stem() # This will initialize self.largest_stem and self.largest_stem_length
        self.modify_utr = modify_utr
        self._skip_svg = False  # SVG format is now used (may cause segfaults in some ViennaRNA versions)
        self._use_average_stem_length = False  # Flag to use average stem length from ensemble vs longest stem from MFE
        self._cached_average_stem_length = None  # Cache for average stem length
        self._average_stem_num_samples = 100  # Number of samples for average stem calculation

    @property
    def sequence(self):
        return ''.join(self.codons)
    
    def validate_weights(self):
        all_keys = set(self._loss_weights.keys())
        for k in all_keys:
            if k not in ['mfe', 'fe', 'cai','cpg', 'stem', 'utr_hybridisation', 'initial_hybridisation', 'restriction_sites']: 
                raise ValueError(f"Invalid loss weight key: {k}")
        
        if self._loss_weights.get('mfe',0.0) != 0.0 and self._loss_weights.get('fe',0.0) != 0.0:
            raise ValueError("Only one of mfe or fe can be non-zero")
        return


    def codons_to_amino_acids(self):
        codons = self.codons
        amino_acids = []
        for codon in codons:
            amino_acids.append(codon_to_amino_acid_1L(codon))
        return amino_acids

    def propose_codon_mutation(self):
        if self.modify_utr:
            index = np.random.randint(0, len(self.codons))
        else:
            index = np.random.randint(self.utr_index_range[1]+1, len(self.codons))
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
    def codons(self, sequence:str):
        codons = []
        for i in range(0, len(sequence), 3):
            codons.append(sequence[i:i+3])
        self._codons = np.array(codons, dtype=str)
        self._codons_changed = True  # Mark codons as changed
        self._cached_codons_string = None  # Invalidate cached string
        # Invalidate cached results when codons change
        if hasattr(self, '_pf_computed'):
            self._pf_computed = False
        return

    def calculate_CAI(self, form = 'log'):
        if form == 'log':
            return self.calculate_CAI_log
        elif form == 'linear':
            return np.exp(self.calculate_CAI_log)
        else:
            raise ValueError("Invalid form")

    @property
    def calculate_CAI_log(self):
        """
        Vectorized calculation of the average log CAI of the mRNA sequence.
        """
        if self._cai_log is None or self._codons_changed:
            caidict = codon_adaptation_index.get(self._species, None)
        
            # Get codons as numpy array of str
            codon_arr = self.codons
            # Get CAI values for each codon (use np.vectorize for clarity)
            cai_vec = np.vectorize(lambda codon: caidict.get(codon, 0.0))
            cai_values = cai_vec(codon_arr)

            valid = ( cai_values > 0 ).all()
            if not valid:
                raise ValueError(f"Some codons have a CAI of 0: {codon_arr[cai_values <= 0]}")
            self._cai_log = np.mean(np.log(cai_values))
            # Note: _codons_changed is reset in codons_string property

        return self._cai_log

    @property
    def mfe(self) -> float:
        # Check sequence-based cache
        seq = self.codons_string
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
                self._minimum_folding_energy = ensemble_energy / self.RT / len(self.codons)
                
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
        seq = self.codons_string
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
                    self._free_energy = ensemble_energy / self.RT / len(self.codons)
                    self._pf_computed = True
                else:
                    # If pf was already computed, we need to recompute it (shouldn't happen often)
                    ensemble_energy = self._linearpartition.calculate_partition_function(seq)
                    self._free_energy = ensemble_energy / self.RT / len(self.codons)
                
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
        seq = self.codons_string
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
        seq = self.codons_string
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
        if what is None or what == 'all' or (isinstance(what, list) and 'all' in what):
            # Reset everything (original behavior)
            self._structure = None
            self._minimum_folding_energy = None
            self._free_energy = None
            self._prob_matrix = None
            self._cai_log = None
            self._largest_stem = None
            self._cpg_count = None
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
            if self.modify_utr:
                self._utr_index_range = None
                self._utr_sequence = None
        else:
            # Selective reset based on what's needed
            if isinstance(what, str):
                what = [what]
            
            if 'energy' in what or 'fe' in what or 'mfe' in what:
                self._free_energy = None
                self._minimum_folding_energy = None
                if hasattr(self, '_cached_seq_hash'):
                    self._cached_seq_hash = None
                if hasattr(self, '_cached_mfe_seq_hash'):
                    self._cached_mfe_seq_hash = None
                # Note: Don't reset fold_compound, just reset pf_computed flag
                self._pf_computed = False
            
            if 'structure' in what:
                self._structure = None
                self._largest_stem = None
                if hasattr(self, '_cached_structure_seq_hash'):
                    self._cached_structure_seq_hash = None
            
            if 'bpp' in what or 'prob_matrix' in what:
                self._prob_matrix = None
                if hasattr(self, '_cached_bpp_seq_hash'):
                    self._cached_bpp_seq_hash = None
            
            if 'cai' in what:
                self._cai_log = None
            
            if 'cpg' in what:
                self._cpg_count = None
            
            if 'stem' in what:
                self._largest_stem = None
                self._cached_average_stem_length = None
        
        return

    @property
    def codons_string(self):
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

    @property
    def cpg_count(self) -> int:
        """
        Count the number of CpG dinucleotides in the codon sequence.
        CpG refers to a cytosine (C) followed by a guanine (G).
        
        Returns:
            int: Number of CpG dinucleotides in the sequence
        """
        if self._cpg_count is None or self._codons_changed:
            seq = self.codons_string
            count = 0
            for i in range(len(seq) - 1):
                if seq[i] == 'C' and seq[i+1] == 'G':
                    count += 1
            self._cpg_count = count
        return self._cpg_count

    @property
    def loss(self) -> float:
        loss = 0.0
        w_cai = self._loss_weights.get('cai', 0.0)
        if w_cai != 0.0:
            loss += -w_cai * self.calculate_CAI_log
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
        w_utr = self._loss_weights.get('utr_hybridisation', 0.0)
        if w_utr != 0.0:
            loss += w_utr * self.utr_hybridisation_penalty
        w_hairpin = self._loss_weights.get('initial_hybridisation', 0.0)
        if w_hairpin != 0.0:
            loss += w_hairpin * self.initial_hybridisation_penalty
        return loss

    def visualize_structure(self, filename="structure.svg", format="svg"):
        print(f"Visualizing structure with minimum energy in {filename}")
        try:
            seq = self.codons_string
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
    def utr_hybridisation_penalty(self):
        return self.hybridisation_penalty(self.utr_index_range[0], self.utr_index_range[1])

    @property
    def initial_hybridisation_penalty(self):
        if self.initial_region_end_index is None:
            print(f"Initial region end index not set")
            raise ValueError("Initial region end index not set but weight not equal to 0")
        return self.hybridisation_penalty(self.start_index, self.initial_region_end_index)

    def hybridisation_penalty(self, start_index, end_index ):
        if start_index == end_index:
            raise ValueError("Start index and end index are the same")
        
        # Use prob_matrix property (will compute if needed)
        try:
            prob_matrix = self.prob_matrix
        except Exception as e:
            print(f"Warning: Could not compute prob_matrix for hybridisation_penalty: {e}")
            return 0.0
        
        # Safety check
        if prob_matrix is None or not hasattr(prob_matrix, 'shape'):
            return 0.0
        
        # Get the hybridisation probability matrix for the initial part of the mRNA sequence
        # Remember P(i,j) is the probability of the base pair (i,j), not the codon pair (i,j)
        start_index = 3*start_index
        end_index = 3*end_index
        
        # Bounds checking
        seq_len = len(self.codons_string)
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

    def retrieve_utr(self):
        """Find the UTR of the mRNA sequence, corresponding to the 5'UTR until the start codon
        appears.
        """
        try:
            mask = self.codons == 'AUG'
            indices = np.arange( len( self.codons ) )
            #print( f"Indices: {indices}" )
            #print( f"Mask: {mask}" )
            start_index = indices[ mask ].item(0)
            #print( f"Start codon found at index {start_index}" )
        except IndexError:
            raise ValueError( "Start codon not found in the sequence" )

        self._start_index = start_index
        self._utr_sequence = self.codons[:start_index]
        self._utr_index_range = (0, start_index - 1)
        return

    @property
    def start_index(self):
        if self._start_index is None:
            self.retrieve_utr()
        return self._start_index

    @start_index.setter
    def start_index(self, value):
        self._start_index = value
        return

    def stop_index(self):
        mask = np.zeros(len(self.codons), dtype=bool)
        for codon in aminoacid_to_codon_table['Stop']:
            mask = mask | (self.codons == codon)
        indices = np.arange( len( self.codons ) )
        stop_indices = indices[ mask ]
        try:
            if len( stop_indices ) > 1:
                print( "WARNING: Multiple stop codons found in the sequence")
            _stop_index = stop_indices[0]
            assert _stop_index > self.start_index, "Stop codon found before start codon"
            print( f"Stop codon found at index {_stop_index}" )
        except IndexError:
            raise ValueError( f"No stop codon found in the sequence: mask is {mask}" )
        self._stop_index = _stop_index
        return _stop_index

    @property
    def utr_index_range(self):
        """Return the index range of the UTR of the mRNA sequence."""
        if self._utr_index_range is None:
            self.retrieve_utr()
        return self._utr_index_range

    @property
    def utr_sequence(self):
        """Return the sequence of the UTR of the mRNA sequence."""
        if self._utr_sequence is None:
            self.retrieve_utr()
        return self._utr_sequence

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
            seq = self.codons_string
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
            if 'cai' in loss_components:
                cai_value = loss_components['cai']
            else:
                cai_value = np.exp(self.calculate_CAI_log)
            loss_values.append(f"{cai_value:.6e}")
        
        w_mfe = self._loss_weights.get('mfe', 0.0)
        if w_mfe != 0.0:
            if 'mfe' in loss_components:
                mfe_value = loss_components['mfe']
            else:
                try:
                    mfe_value = self.mfe
                except Exception:
                    mfe_value = 0.0
            loss_values.append(f"{mfe_value:.6e}")
        
        w_fe = self._loss_weights.get('fe', 0.0)
        if w_fe != 0.0:
            if 'fe' in loss_components:
                fe_value = loss_components['fe']
            else:
                try:
                    fe_value = self.free_energy
                except Exception:
                    fe_value = 0.0
            loss_values.append(f"{fe_value:.6e}")
        
        w_cpg = self._loss_weights.get('cpg', 0.0)
        if w_cpg != 0.0:
            if 'cpg' in loss_components:
                cpg_value = loss_components['cpg']
            else:
                cpg_value = self.cpg_count
            loss_values.append(f"{cpg_value:.6e}")
        
        w_stem = self._loss_weights.get('stem', 0.0)
        if w_stem != 0.0:
            if 'stem' in loss_components:
                stem_value = loss_components['stem']
            else:
                try:
                    stem_value = self.stem_penalty
                except Exception:
                    stem_value = 0.0
            loss_values.append(f"{stem_value:.6e}")
        
        w_utr = self._loss_weights.get('utr_hybridisation', 0.0)
        if w_utr != 0.0:
            if 'utr_hybridisation' in loss_components:
                utr_value = loss_components['utr_hybridisation']
            else:
                try:
                    utr_value = self.utr_hybridisation_penalty
                except Exception:
                    utr_value = 0.0
            loss_values.append(f"{utr_value:.6e}")
        
        w_hairpin = self._loss_weights.get('initial_hybridisation', 0.0)
        if w_hairpin != 0.0:
            if 'initial_hybridisation' in loss_components:
                hairpin_value = loss_components['initial_hybridisation']
            else:
                try:
                    hairpin_value = self.initial_hybridisation_penalty
                except (ValueError, AttributeError):
                    hairpin_value = 0.0
            loss_values.append(f"{hairpin_value:.6e}")
        
        w_restriction = self._loss_weights.get('restriction_sites', 0.0)
        if w_restriction != 0.0:
            if 'restriction_sites' in loss_components:
                restriction_value = loss_components['restriction_sites']
            else:
                try:
                    restriction_value = self.restriction_sites_count
                except Exception:
                    restriction_value = 0.0
            loss_values.append(f"{restriction_value:.6e}")
        
        # Format: step, acceptance_rate, sequence_identity_%, [loss_component_values...]
        loss_str = "\t".join(loss_values) if loss_values else ""
        if loss_str:
            line = f"{step}\t{acceptance_rate:.4f}\t{sequence_identity:.4f}\t{loss_str}\n"
        else:
            line = f"{step}\t{acceptance_rate:.4f}\t{sequence_identity:.4f}\n"
        
        if isinstance(output_file, str):
            # If it's a file path, open in append mode
            with open(output_file, 'a') as f:
                f.write(line)
        else:
            # If it's already a file object, write to it
            output_file.write(line)
            output_file.flush()

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
            seq = self.codons_string
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
                cai_value = np.exp(self.calculate_CAI_log)
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
            
            w_utr = self._loss_weights.get('utr_hybridisation', 0.0)
            if w_utr != 0.0:
                try:
                    utr_value = self.utr_hybridisation_penalty
                    f.write(f"utr_hybridisation, {utr_value:.6e}\n")
                except Exception as e:
                    print(f"Warning: Could not compute utr_hybridisation_penalty for loss components: {e}")
                    f.write(f"utr_hybridisation, 0.000000e+00\n")
            
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
            seq = self.codons_string
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
        loss = self.loss
        
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
        if self._loss_weights.get('utr_hybridisation', 0.0) != 0.0:
            header_parts.append("UTR_Hybridisation")
        if self._loss_weights.get('initial_hybridisation', 0.0) != 0.0:
            header_parts.append("Initial_Hybridisation")
        if self._loss_weights.get('restriction_sites', 0.0) != 0.0:
            header_parts.append("Restriction_Sites")
        
        with open(output_filename, 'w') as f:
            f.write("\t".join(header_parts) + "\n")
        
        for i, T in enumerate( T_schedule ):
            current_codon, new_codon, index = self.propose_codon_mutation()
            accepted = False
            
            if new_codon != current_codon:
                # Determine what needs to be reset based on loss weights
                reset_what = []
                if self._loss_weights.get('fe', 0.0) != 0.0 or self._loss_weights.get('mfe', 0.0) != 0.0:
                    reset_what.append('energy')
                if self._loss_weights.get('stem', 0.0) != 0.0:
                    reset_what.append('structure')
                    reset_what.append('stem')
                if self._loss_weights.get('utr_hybridisation', 0.0) != 0.0 or \
                   self._loss_weights.get('initial_hybridisation', 0.0) != 0.0:
                    reset_what.append('bpp')
                if self._loss_weights.get('cai', 0.0) != 0.0:
                    reset_what.append('cai')
                if self._loss_weights.get('cpg', 0.0) != 0.0:
                    reset_what.append('cpg')
                
                # Selective reset instead of full reset
                if reset_what:
                    self.reset(what=reset_what)
                else:
                    # If no weights require energy/structure, minimal reset
                    self.reset(what=['cai'])  # CAI might change
                
                # Make the mutation
                self._codons[index] = new_codon
                self._codons_changed = True  # Mark codons as changed
                
                try:
                    new_loss = self.loss
                except Exception as e:
                    print(f"Error computing loss at step {i+1}: {e}")
                    import traceback
                    traceback.print_exc()
                    # Revert mutation on error
                    self._codons[index] = current_codon
                    self._codons_changed = True
                    self.reset(what=reset_what if reset_what else ['cai'])
                    continue
                if new_loss < loss:
                    # Always accept if loss decreases
                    loss = new_loss
                    accepted = True
                else:
                    delta = new_loss - loss
                    # Avoid division by zero or very small T
                    if T > 1e-10:
                        accept_probability = np.exp(-delta/T)
                    else:
                        # If T is very small, only accept if loss decreases
                        accept_probability = 0.0
                    if np.random.rand() < accept_probability:
                        # Accept with probability
                        loss = new_loss
                        accepted = True
                    else:
                        # Revert the mutation by restoring the original codon
                        self._codons[index] = current_codon
                        self._codons_changed = True
                        # Selective reset after reverting
                        if reset_what:
                            self.reset(what=reset_what)
                        accepted = False
                        if verbose:
                            print(f"Rejected mutation: {current_codon} to {new_codon} at index {index}")
                
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
                print(f"Current loss: {loss}, energy: {energy}, cai: {np.exp(self.calculate_CAI_log)}")
                
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
                
                # Save loss components
                try:
                    self.save_loss_components(i)
                except Exception as e:
                    print(f"Error saving loss components at step {i}: {e}")
                    #import traceback
                    #traceback.print_exc()

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
        
        print(f"Final loss: {loss}, mfe: {self.mfe}, cai: {np.exp(self.calculate_CAI_log)}")
        print(f"Statistics saved to {output_filename}")
            
        return 

if __name__ == "__main__":
    sequence = "AUUAAAGGGUUUUAGCAGGGGCCCCUUAAGGCGGCGGAGGACGGG"
    sequence = sequence * 2 
    print(f"Sequence length: {len(sequence)}")
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
                modify_utr=True,
                initial_region_end_index = 30,
                T_K = 310,
                average_stem_num_samples = 20
                )
    print(f"McCaskill free-energy (*codon): {system.free_energy} kcal/mol per codon")
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


def plot_optimization_statistics(statistics_file="opt_statistics.txt", output_file=None, figsize=(14, 10)):
    """
    Plot optimization statistics from opt_statistics.txt file.
    
    Creates subplots for:
    - Acceptance rate over steps
    - Sequence identity over steps
    - Loss components over steps (CAI, Free Energy, CpG, Stem, UTR Hybridisation, Initial Hybridisation)
    - Total loss (if computable from components)
    
    Args:
        statistics_file: Path to the statistics file (default: "opt_statistics.txt")
        output_file: Path to save the plot (default: None, displays interactively)
        figsize: Figure size tuple (default: (14, 10))
    
    Returns:
        matplotlib figure object
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError("matplotlib is required for plotting. Install with: pip install matplotlib")
    
    # Read the statistics file
    try:
        with open(statistics_file, 'r') as f:
            all_lines = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        raise FileNotFoundError(f"Statistics file not found: {statistics_file}")
    
    if not all_lines:
        raise ValueError(f"No data found in {statistics_file}")
    
    # Parse header (first line starting with '#')
    header_line = None
    data_lines = []
    for line in all_lines:
        if line.startswith('#'):
            if header_line is None:
                header_line = line.lstrip('#').strip()
        else:
            data_lines.append(line)
    
    if not header_line:
        raise ValueError("Could not find header in statistics file")
    
    if not data_lines:
        raise ValueError("No data lines found in statistics file")
    
    # Parse header
    header = [col.strip() for col in header_line.split('\t')]
    
    # Parse data
    data = {}
    for col in header:
        data[col] = []
    
    for line in data_lines:
        values = [v.strip() for v in line.split('\t')]
        for i, val in enumerate(values):
            if i < len(header):
                try:
                    data[header[i]].append(float(val))
                except ValueError:
                    data[header[i]].append(val)
    
    # Convert to numpy arrays for easier plotting
    steps = np.array(data.get('Step', []))
    if len(steps) == 0:
        raise ValueError("No data points found in statistics file")
    
    # Determine number of subplots needed
    # Always plot: Acceptance_Rate, Sequence_Identity_%
    # Then plot all loss components that are present
    loss_components = []
    for col in header:
        if col not in ['Step', 'Acceptance_Rate', 'Sequence_Identity_%']:
            loss_components.append(col)
    
    # Create subplots
    n_plots = 2 + len(loss_components)  # Acceptance rate + Sequence identity + loss components
    n_cols = 2
    n_rows = (n_plots + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_plots == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    plot_idx = 0
    
    # Plot 1: Acceptance Rate
    if 'Acceptance_Rate' in data and len(data['Acceptance_Rate']) > 0:
        ax = axes[plot_idx]
        ax.plot(steps, data['Acceptance_Rate'], 'b-', linewidth=2, marker='o', markersize=4)
        ax.set_xlabel('Step', fontsize=11)
        ax.set_ylabel('Acceptance Rate', fontsize=11)
        ax.set_title('Acceptance Rate Over Optimization', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_ylim([0, 1.1])
        plot_idx += 1
    
    # Plot 2: Sequence Identity
    if 'Sequence_Identity_%' in data and len(data['Sequence_Identity_%']) > 0:
        ax = axes[plot_idx]
        ax.plot(steps, data['Sequence_Identity_%'], 'g-', linewidth=2, marker='s', markersize=4)
        ax.set_xlabel('Step', fontsize=11)
        ax.set_ylabel('Sequence Identity (%)', fontsize=11)
        ax.set_title('Sequence Identity Over Optimization', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_ylim([0, 110])
        plot_idx += 1
    
    # Plot loss components
    colors = ['r', 'm', 'c', 'orange', 'purple', 'brown', 'pink', 'gray']
    for i, component in enumerate(loss_components):
        if component in data and len(data[component]) > 0:
            ax = axes[plot_idx]
            color = colors[i % len(colors)]
            values = np.array(data[component])
            
            # Handle different scales - some might be very small (scientific notation)
            ax.plot(steps, values, color=color, linewidth=2, marker='.', markersize=4, label=component)
            ax.set_xlabel('Step', fontsize=11)
            ax.set_ylabel(component, fontsize=11)
            ax.set_title(f'{component} Over Optimization', fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3)
            
            # Use scientific notation for very small or large values
            if np.any(np.abs(values) < 1e-3) or np.any(np.abs(values) > 1e3):
                ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
            
            plot_idx += 1
    
    # Hide unused subplots
    for idx in range(plot_idx, len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    
    # Save or show
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()
    
    return fig

