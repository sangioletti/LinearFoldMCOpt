"""
Python wrapper for LinearPartition to calculate partition function and base pair probabilities.
"""

import subprocess
import os
import tempfile
import re
import numpy as np
from pathlib import Path


class LinearPartitionWrapper:
    """
    Wrapper class to interface with LinearPartition executable.
    """
    
    def __init__(self, linearpartition_path=None, use_vienna=True, beamsize=100, verbose=False):
        """
        Initialize the LinearPartition wrapper.
        
        Args:
            linearpartition_path: Path to linearpartition executable. If None, uses './linearpartition'
            use_vienna: Whether to use Vienna parameters (default True)
            beamsize: Beam size for LinearPartition (default 100)
            verbose: Whether to print verbose output (default False)
        """
        if linearpartition_path is None:
            # Try to find linearpartition in the current directory
            script_dir = Path(__file__).parent.absolute()
            linearpartition_path = script_dir / "linearpartition"
            if not os.path.exists(linearpartition_path):
                # Try relative path
                linearpartition_path = "./linearpartition"
        
        # Convert Path to string if needed
        if isinstance(linearpartition_path, Path):
            self.linearpartition_path = str(linearpartition_path)
        else:
            self.linearpartition_path = linearpartition_path
        self.use_vienna = use_vienna
        self.beamsize = beamsize
        self.verbose = verbose
        
        # Check if executable exists
        if not os.path.exists(self.linearpartition_path):
            raise FileNotFoundError(f"LinearPartition executable not found at {self.linearpartition_path}")
    
    def calculate_partition_function(self, sequence):
        """
        Calculate the partition function (ensemble free energy) for a sequence.
        
        Args:
            sequence: RNA sequence string
            
        Returns:
            float: Ensemble free energy in kcal/mol
        """
        # Build command
        cmd = [self.linearpartition_path, "-b", str(self.beamsize), "-p"]  # -p for partition function only
        if self.use_vienna:
            cmd.append("-V")
        if self.verbose:
            cmd.append("--verbose")
        
        # Run LinearPartition
        try:
            result = subprocess.run(
                cmd,
                input=sequence,
                text=True,
                capture_output=True,
                check=True
            )
            
            # Parse stderr for free energy
            # Format: "Free Energy of Ensemble: -XX.XX kcal/mol"
            stderr_lines = result.stderr.split('\n')
            for line in stderr_lines:
                if "Free Energy of Ensemble:" in line:
                    match = re.search(r'Free Energy of Ensemble:\s+([-\d.]+)\s+kcal/mol', line)
                    if match:
                        return float(match.group(1))
            
            # If not found in stderr, try stdout
            stdout_lines = result.stdout.split('\n')
            for line in stdout_lines:
                if "Free Energy of Ensemble:" in line:
                    match = re.search(r'Free Energy of Ensemble:\s+([-\d.]+)\s+kcal/mol', line)
                    if match:
                        return float(match.group(1))
            
            raise ValueError("Could not parse free energy from LinearPartition output")
            
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"LinearPartition failed: {e.stderr}")
    
    def calculate_bpp_matrix(self, sequence, cutoff=0.0):
        """
        Calculate the base pair probability matrix for a sequence.
        
        Args:
            sequence: RNA sequence string
            cutoff: Minimum probability threshold (default 0.0)
            
        Returns:
            numpy.ndarray: Base pair probability matrix (n x n, 0-indexed)
        """
        # Create temporary file for output
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bpp') as tmp_file:
            tmp_filename = tmp_file.name
        
        try:
            # Build command
            cmd = [self.linearpartition_path, "-b", str(self.beamsize), "-r", tmp_filename]
            if self.use_vienna:
                cmd.append("-V")
            if cutoff > 0.0:
                cmd.extend(["-c", str(cutoff)])
            if self.verbose:
                cmd.append("--verbose")
            
            # Run LinearPartition
            try:
                result = subprocess.run(
                    cmd,
                    input=sequence,
                    text=True,
                    capture_output=True,
                    check=True
                )
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"LinearPartition failed: {e.stderr}")
            
            # Parse BPP file
            # Format: "i j prob" (1-indexed)
            seq_len = len(sequence)
            bpp_matrix = np.zeros((seq_len, seq_len), dtype=np.float64)
            
            if os.path.exists(tmp_filename):
                with open(tmp_filename, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith('#'):
                            continue
                        
                        parts = line.split()
                        if len(parts) >= 3:
                            try:
                                i = int(parts[0]) - 1  # Convert to 0-indexed
                                j = int(parts[1]) - 1  # Convert to 0-indexed
                                prob = float(parts[2])
                                
                                if 0 <= i < seq_len and 0 <= j < seq_len:
                                    bpp_matrix[i, j] = prob
                                    bpp_matrix[j, i] = prob  # Symmetric matrix
                            except (ValueError, IndexError):
                                continue
                
                # Clean up
                os.unlink(tmp_filename)
            
            return bpp_matrix
            
        except Exception as e:
            # Clean up on error
            if os.path.exists(tmp_filename):
                os.unlink(tmp_filename)
            raise
    
    def calculate_mea_structure(self, sequence, gamma=3.0):
        """
        Calculate the MEA (Maximum Expected Accuracy) structure for a sequence.
        
        Args:
            sequence: RNA sequence string
            gamma: MEA gamma parameter (default 3.0)
            
        Returns:
            tuple: (structure_string, ensemble_free_energy)
                structure_string: Dot-bracket notation structure
                ensemble_free_energy: Ensemble free energy in kcal/mol
        """
        # Build command
        cmd = [self.linearpartition_path, "-b", str(self.beamsize), "-M", "-g", str(gamma)]
        if self.use_vienna:
            cmd.append("-V")
        if self.verbose:
            cmd.append("--verbose")
        
        # Run LinearPartition
        try:
            result = subprocess.run(
                cmd,
                input=sequence,
                text=True,
                capture_output=True,
                check=True
            )
            
            # Parse output
            # Format:
            # Free Energy of Ensemble: -XX.XX kcal/mol
            # SEQUENCE
            # STRUCTURE
            
            lines = result.stdout.split('\n')
            ensemble_energy = None
            
            # Try to find free energy in stderr first
            stderr_lines = result.stderr.split('\n')
            for line in stderr_lines:
                if "Free Energy of Ensemble:" in line:
                    match = re.search(r'Free Energy of Ensemble:\s+([-\d.]+)\s+kcal/mol', line)
                    if match:
                        ensemble_energy = float(match.group(1))
                        break
            
            # If not found, try stdout
            if ensemble_energy is None:
                for line in lines:
                    if "Free Energy of Ensemble:" in line:
                        match = re.search(r'Free Energy of Ensemble:\s+([-\d.]+)\s+kcal/mol', line)
                        if match:
                            ensemble_energy = float(match.group(1))
                            break
            
            # Find structure (should be after sequence line)
            structure = None
            found_sequence = False
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                if line == sequence:
                    found_sequence = True
                    continue
                if found_sequence and all(c in '().' for c in line):
                    structure = line
                    break
            
            if structure is None:
                raise ValueError("Could not parse MEA structure from LinearPartition output")
            
            if ensemble_energy is None:
                raise ValueError("Could not parse ensemble free energy from LinearPartition output")
            
            return structure, ensemble_energy
            
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"LinearPartition failed: {e.stderr}")
    
    def calculate_all(self, sequence, cutoff=0.0, gamma=3.0):
        """
        Calculate partition function, BPP matrix, and MEA structure in one call.
        
        Args:
            sequence: RNA sequence string
            cutoff: Minimum probability threshold for BPP (default 0.0)
            gamma: MEA gamma parameter (default 3.0)
            
        Returns:
            dict: Dictionary with keys:
                'ensemble_energy': Ensemble free energy in kcal/mol
                'bpp_matrix': Base pair probability matrix (n x n)
                'mea_structure': MEA structure in dot-bracket notation
        """
        # Calculate partition function
        ensemble_energy = self.calculate_partition_function(sequence)
        
        # Calculate BPP matrix
        bpp_matrix = self.calculate_bpp_matrix(sequence, cutoff=cutoff)
        
        # Calculate MEA structure
        mea_structure, _ = self.calculate_mea_structure(sequence, gamma=gamma)
        
        return {
            'ensemble_energy': ensemble_energy,
            'bpp_matrix': bpp_matrix,
            'mea_structure': mea_structure
        }

