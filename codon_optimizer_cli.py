#!/usr/bin/env python3
"""
Codon Optimization CLI

A command-line interface for mRNA codon optimization using fire.
Supports YAML configuration files for batch processing.

Usage:
    python codon_optimizer_cli.py optimize --config config.yaml
    python codon_optimizer_cli.py fetch_pdb --pdb_id 1BA3
    python codon_optimizer_cli.py batch --csv_file sequences.csv --config config.yaml --output_dir jobs
"""

import os
import sys
import yaml
import numpy as np
import fire
import string
import random
from typing import Optional, Dict, Any, List, Union
from pathlib import Path

# Import from the local modules
from codons import (
    human_aa_to_codon_cai, 
    ecoli_aa_to_codon_cai,
    codon_table
)

from utils import (
    codon_to_amino_acid_1L,
    codon_to_amino_acid_3L,
    aa_1L_to_codons,
    aa_3L_to_codons,
    aa_to_codon_sequence
)

from mRNA import mRNA, aa_to_codon_sequence


def fetch_sequence_from_pdb(pdb_id: str, chain: Optional[str] = None, target_path: str = ".") -> str:
    """
    Fetch amino acid sequence from PDB database using biotite.
    
    Args:
        pdb_id: PDB identifier (e.g., "1BA3")
        chain: Optional chain identifier. If None, returns the first chain.
        target_path: Directory to download FASTA file to
        
    Returns:
        Amino acid sequence as a string (1-letter code)
    """
    try:
        from biotite.database import rcsb
        from biotite.sequence.io import fasta
    except ImportError:
        raise ImportError(
            "biotite is required for PDB fetching. "
            "Install with: pip install biotite"
        )
    
    # Download the FASTA for the PDB entry
    fasta_file = rcsb.fetch(pdb_id, "fasta", target_path=target_path)
    
    # Read the sequence from the FASTA file
    fasta_file_obj = fasta.FastaFile.read(fasta_file)
    
    if len(fasta_file_obj) == 0:
        raise ValueError(f"No sequences found in FASTA file for PDB ID: {pdb_id}")
    
    # Get first header and corresponding sequence, or specific chain
    if chain is not None:
        # Look for specific chain
        for header, seq in fasta_file_obj.items():
            if f"Chain {chain}" in header or f"chain {chain.lower()}" in header.lower():
                return str(seq)
        raise ValueError(f"Chain {chain} not found in PDB {pdb_id}")
    else:
        # Get first sequence
        first_header = list(fasta_file_obj.keys())[0]
        return str(fasta_file_obj[first_header])


def generate_unique_id(length: int = 10) -> str:
    """Generate a random unique identifier of specified length."""
    chars = string.ascii_uppercase + string.digits
    return ''.join(random.choice(chars) for _ in range(length))


def load_config(config_path: str) -> Dict[str, Any]:
    """Load and validate YAML configuration file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Set defaults for missing values
    defaults = {
        'species': 'human',
        'verbose': False,
        'T_K': 300,
        'modify_fivep_utr': False,
        'modify_threep_utr': False,
        'initial_region_end_index': 0,
        'beamsize': 100,
        'bpp_cutoff': 0.01,
        'start_from_optimal_cai': True,
        'loss_weights': {
            'cai': 10.0,
            'fe': 1.0,
            'cpg': 3.0,
            'stem': 10.0,
            'fivep_utr_hybridisation': 3.0,
            'threep_utr_hybridisation': 0.0,
            'initial_hybridisation': 0.0,
            'codon_divergence': 10.0,
        },
        'optimization': {
            'T_opt': 1.0,
            'nsteps': 100,
            'sample_frequency': 20,
            'output_filename': 'opt_statistics.txt',
            'n_sample': None,
            'use_average_stem_length': False,
            'average_stem_num_samples': 20,
        },
        'output': {
            'save_structure_svg': True,
            'save_prob_matrix': True,
            'output_dir': '.',
        },
        # Unique identifier for output files
        'identifier': None,  # If None, a random 10-char ID will be generated
        # Sequence source options
        'sequence': None,  # Direct nucleotide sequence
        'aa_sequence': None,  # Direct amino acid sequence
        'pdb_id': None,  # PDB ID to fetch
        'pdb_chain': None,  # Optional PDB chain
        'fasta_file': None,  # FASTA file path
        # UTR options
        'five_prime_utr': 'GCCACCAUG',  # Default 5' UTR
        'three_prime_utr': 'A'*201,  # Default 3' UTR (after stop codon)
        'additional_car_codons': '',
        'binder_linker': '',
    }
    
    # Merge defaults with config (config takes precedence)
    def deep_merge(base: dict, override: dict) -> dict:
        result = base.copy()
        for key, value in override.items():
            if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                result[key] = deep_merge(result[key], value)
            else:
                result[key] = value
        return result
    
    return deep_merge(defaults, config)


def get_sequence_from_config(config: Dict[str, Any]) -> str:
    """
    Get nucleotide sequence from config following priority:
    1. Direct nucleotide sequence
    2. Direct amino acid sequence (converted)
    3. PDB ID (fetched and converted)
    4. FASTA file (read and converted)
    
    Returns:
        Nucleotide sequence ready for optimization
    """
    sequence = None
    
    # Priority 1: Direct nucleotide sequence
    if config.get('cds_sequence'):
        sequence = config['cds_sequence']
        print(f"Using provided nucleotide sequence (length: {len(sequence)})")
    
    # Priority 2: Direct amino acid sequence
    elif config.get('aa_sequence'):
        aa_seq = config['aa_sequence']
        sequence = aa_to_codon_sequence(aa_seq)
        print(f"Converted amino acid sequence to nucleotides (length: {len(sequence)})")
    
    # Priority 3: PDB ID
    elif config.get('pdb_id'):
        pdb_id = config['pdb_id']
        chain = config.get('pdb_chain')
        print(f"Fetching sequence from PDB: {pdb_id}" + (f" chain {chain}" if chain else ""))
        aa_seq = fetch_sequence_from_pdb(pdb_id, chain)
        sequence = aa_to_codon_sequence(aa_seq)
        print(f"Fetched and converted PDB sequence (length: {len(sequence)})")
    
    # Priority 4: FASTA file
    elif config.get('fasta_file'):
        from biotite.sequence.io import fasta
        fasta_path = config['fasta_file']
        fasta_file_obj = fasta.FastaFile.read(fasta_path)
        if len(fasta_file_obj) == 0:
            raise ValueError(f"No sequences found in FASTA file: {fasta_path}")
        first_header = list(fasta_file_obj.keys())[0]
        aa_seq = str(fasta_file_obj[first_header])
        sequence = aa_to_codon_sequence(aa_seq)
        print(f"Read FASTA file and converted to nucleotides (length: {len(sequence)})")
    
    else:
        raise ValueError(
            "No sequence source provided. Specify one of: "
            "sequence, aa_sequence, pdb_id, or fasta_file"
        )

    # Add binder linker, if specified
    binder_linker = config.get('binder_linker', '')
    if binder_linker:
        sequence += binder_linker

    # Add additional CAR codons if specified
    extra_cds = config.get('extra_cds', '')
    if extra_cds:
        sequence += extra_cds

    # Add 5' UTR if specified
    five_prime_utr = config.get('five_prime_utr', '')
    
    # Add 3' UTR if specified (after stop codon)
    three_prime_utr = config.get('three_prime_utr', '')
    
    return five_prime_utr, sequence, three_prime_utr

def get_cai_dict(species: str) -> dict:
    """Get the appropriate CAI dictionary for the species."""
    if species == 'human':
        return human_aa_to_codon_cai
    elif species == 'e_coli':
        return ecoli_aa_to_codon_cai
    else:
        # Default to human for unsupported species
        print(f"Warning: CAI data for '{species}' not available. Using human CAI.")
        return human_aa_to_codon_cai


class CodonOptimizerCLI:
    """
    Codon Optimization Command Line Interface.
    
    Commands:
        optimize    Run codon optimization using a YAML config file
        fetch_pdb   Fetch amino acid sequence from PDB database
        batch       Generate PBS batch scripts for multiple sequences
        example     Generate an example YAML configuration file
    """
    
    def optimize(
        self, 
        config: str,
        cds_sequence: Optional[str] = None,
        fivep_utr: Optional[str] = None,
        threep_utr: Optional[str] = None,
        output_dir: Optional[str] = None,
        job_name: Optional[str] = None,
        identifier: Optional[str] = None,
    ):
        """
        Run codon optimization using configuration from a YAML file.
        
        Args:
            config: Path to YAML configuration file
            cds_sequence: Optional sequence override (nucleotide or amino acid)
            fivep_utr: Optional 5' UTR override
            threep_utr: Optional 3' UTR override
            output_dir: Optional output directory override
            job_name: Optional job name for output files
            identifier: Optional unique identifier for output files
        """
        # Load configuration
        cfg = load_config(config)
        
        # Override sequence if provided via command line
        if cds_sequence:
            # Detect if it's amino acid or nucleotide
            if all(c in 'ACDEFGHIKLMNPQRSTVWY*X' for c in cds_sequence.upper()):
                if all(c in 'AUGC' for c in cds_sequence.upper().replace( 'T', 'U')):
                    cfg['cds_sequence'] = cds_sequence
                else:
                    # If it's amino acid sequence, convert to nucleotide
                    cds_sequence = aa_to_codon_sequence(cds_sequence)
                    cfg['cds_sequence'] = cds_sequence
            else:
                cfg['cds_sequence'] = cds_sequence
        
        if fivep_utr:
            cfg['fivep_utr'] = fivep_utr
        if threep_utr:
            cfg['threep_utr'] = threep_utr
        
        # Determine unique identifier
        # Priority: CLI arg > config file > auto-generate
        if identifier:
            unique_id = identifier
        elif cfg.get('identifier'):
            unique_id = cfg['identifier']
        else:
            unique_id = generate_unique_id(10)
        
        print(f"Using identifier: {unique_id}")
        
        # Override output directory if provided
        if output_dir:
            base_out_dir = output_dir
        else:
            base_out_dir = cfg['output']['output_dir']
        
        # Create result folder with unique identifier
        result_dir = os.path.join(base_out_dir, f"result_{unique_id}")
        os.makedirs(result_dir, exist_ok=True)
        
        # Change to result directory
        original_dir = os.getcwd()
        os.chdir(result_dir)
        
        try:
            # Get sequence
            fivep_utr, cds_sequence, threep_utr = get_sequence_from_config(cfg)
            
            # Get CAI dictionary
            cai_dict = get_cai_dict(cfg['species'])

            all_sequence = fivep_utr + cds_sequence + threep_utr
            
            # Create mRNA system
            print(f"\nInitializing mRNA system...")
            print(f"  Sequence length: {len(all_sequence)} nucleotides")
            print(f"  Species: {cfg['species']}")
            print(f"  Loss weights: {cfg['loss_weights']}")
            
            system = mRNA(
                cds_sequence=cds_sequence,
                fivep_utr=fivep_utr,
                threep_utr=threep_utr,
                species=cfg['species'],
                aa_to_codon_cai=cai_dict,
                start_from_optimal_cai=cfg['start_from_optimal_cai'],
                verbose=cfg['verbose'],
                loss_weights=cfg['loss_weights'],
                modify_fivep_utr=cfg['modify_fivep_utr'],
                modify_threep_utr=cfg['modify_threep_utr'],
                mutable_range=cfg['mutable_range'],
                initial_region_end_index=cfg['initial_region_end_index'],
                T_K=cfg['T_K'],
                bpp_cutoff=cfg['bpp_cutoff'],
                beamsize=cfg['beamsize'],
            )
            
            print(f"  Number of codons: {system.n_cds}")
            print(f"  Initial free energy (x nucleotide): {system.free_energy/system.n_nts:.4f} kcal/mol")
            print(f"  Initial CAI (x codon): {system.calculate_CAI(form='linear',normalise=True):.4f}")
            
            # Run optimization
            opt_cfg = cfg['optimization']
            output_filename = job_name + "_statistics.txt" if job_name else opt_cfg['output_filename']
            
            print(f"\nStarting optimization...")
            print(f"  T_opt: {opt_cfg['T_opt']}")
            print(f"  nsteps: {opt_cfg['nsteps']}")
            
            system.optimize_codon_usage(
                T_opt=opt_cfg['T_opt'],
                nsteps=opt_cfg['nsteps'],
                sample_frequency=opt_cfg['sample_frequency'],
                verbose=cfg['verbose'],
                output_filename=output_filename,
                n_sample=opt_cfg.get('n_sample'),
                use_average_stem_length=opt_cfg['use_average_stem_length'],
                average_stem_num_samples=opt_cfg['average_stem_num_samples'],
            )
            
            # Print results
            print(f"\n{'='*60}")
            print("OPTIMIZATION COMPLETE")
            print(f"{'='*60}")
            print(f"Final loss: {system.loss:.6f}")
            print(f"Final MFE: {system.mfe:.4f} kcal/mol per codon")
            print(f"Final free energy: {system.free_energy:.4f} kcal/mol per codon")
            print(f"Final CAI: {np.exp(system.calculate_CAI_log()):.4f}")
            print(f"Structure: {system._structure[:50]}..." if len(system._structure) > 50 else f"Structure: {system._structure}")
            
            # Save optimized sequence
            seq_filename = job_name + "_optimized.txt" if job_name else "optimized_sequence.txt"
            with open(seq_filename, 'w') as f:
                f.write(f"# Optimized mRNA sequence\n")
                f.write(f"# Length: {len(system.codons_string)} nucleotides\n")
                f.write(f"# CAI: {np.exp(system.calculate_CAI_log()):.4f}\n")
                f.write(f"# MFE: {system.mfe:.4f} kcal/mol\n")
                f.write(f"\n{system.codons_string}\n")
            print(f"\nOptimized sequence saved to: {seq_filename}")
            
            # Save structure visualization if enabled
            if cfg['output'].get('save_structure_svg', True):
                svg_filename = job_name + "_structure.svg" if job_name else "structure.svg"
                try:
                    system.visualize_structure(filename=svg_filename, format="svg")
                    print(f"Structure visualization saved to: {svg_filename}")
                except Exception as e:
                    print(f"Warning: Could not save structure SVG: {e}")
            
            return {
                'final_loss': system.loss,
                'final_cai': np.exp(system.calculate_CAI_log()),
                'final_mfe': system.mfe,
                'sequence': system.codons_string,
                'structure': system._structure,
            }
            
        finally:
            os.chdir(original_dir)
    
    def fetch_pdb(
        self, 
        pdb_id: str, 
        chain: Optional[str] = None,
        output: Optional[str] = None,
        convert_to_nucleotide: bool = False,
    ):
        """
        Fetch amino acid sequence from PDB database.
        
        Args:
            pdb_id: PDB identifier (e.g., "1BA3")
            chain: Optional chain identifier
            output: Optional output file path
            convert_to_nucleotide: If True, convert to nucleotide sequence
        """
        aa_sequence = fetch_sequence_from_pdb(pdb_id, chain)
        
        print(f"\nPDB ID: {pdb_id}")
        if chain:
            print(f"Chain: {chain}")
        print(f"Amino acid sequence length: {len(aa_sequence)}")
        print(f"\nAmino acid sequence:")
        print(aa_sequence)
        
        if convert_to_nucleotide:
            nuc_sequence = aa_to_codon_sequence(aa_sequence)
            print(f"\nNucleotide sequence length: {len(nuc_sequence)}")
            print(f"\nNucleotide sequence:")
            print(nuc_sequence)
            
            if output:
                with open(output, 'w') as f:
                    f.write(f">{pdb_id}" + (f"_{chain}" if chain else "") + "\n")
                    f.write(nuc_sequence + "\n")
                print(f"\nSaved to: {output}")
        elif output:
            with open(output, 'w') as f:
                f.write(f">{pdb_id}" + (f"_{chain}" if chain else "") + "\n")
                f.write(aa_sequence + "\n")
            print(f"\nSaved to: {output}")
        
        return aa_sequence
    
    def batch(
        self,
        csv_file: str,
        config: str,
        output_dir: str = "batch_jobs",
        pbs_template: Optional[str] = None,
        queue: str = "normal",
        walltime: str = "24:00:00",
        ncpus: int = 1,
        mem: str = "8GB",
        python_path: Optional[str] = None,
        id_column: Optional[str] = None,
    ):
        """
        Generate PBS batch scripts for multiple sequences from a CSV file.
        
        The CSV file should have columns:
        - name: Job name / identifier (first column)
        - identifier: Unique identifier for output (optional, second column or specified by id_column)
        - sequence: Nucleotide sequence OR amino acid sequence
        - type: 'nucleotide' or 'aminoacid' (optional, auto-detected if not present)
        
        If the CSV has more than one column, the second column can be used as the identifier.
        
        Args:
            csv_file: Path to CSV file with sequences
            config: Path to base YAML configuration file
            output_dir: Directory to write PBS scripts and job configs
            pbs_template: Optional custom PBS template file
            queue: PBS queue name
            walltime: PBS walltime (HH:MM:SS)
            ncpus: Number of CPUs per job
            mem: Memory per job (e.g., "8GB")
            python_path: Path to Python interpreter (default: current python)
            id_column: Name of column to use as identifier (default: 'identifier' or second column)
        """
        import csv
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Load base config
        base_config = load_config(config)
        
        # Read CSV file
        sequences = []
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            columns = reader.fieldnames or []
            for row in reader:
                sequences.append(row)
        
        print(f"Found {len(sequences)} sequences in {csv_file}")
        print(f"CSV columns: {columns}")
        
        # Determine identifier column
        # Priority: id_column arg > 'identifier' column > second column > auto-generate
        if id_column and id_column in columns:
            id_col = id_column
        elif 'identifier' in columns:
            id_col = 'identifier'
        elif len(columns) >= 2:
            id_col = columns[1]  # Use second column as identifier
            print(f"Using second column '{id_col}' as identifier")
        else:
            id_col = None
        
        # Get Python path
        if python_path is None:
            python_path = sys.executable
        
        # Get script path
        script_path = os.path.abspath(__file__)
        
        # Generate jobs
        job_scripts = []
        for i, seq_info in enumerate(sequences):
            name = seq_info.get('name', f'job_{i:04d}')
            sequence = seq_info.get('sequence', '')
            seq_type = seq_info.get('type', 'auto')
            
            # Get identifier from CSV or generate one
            if id_col and seq_info.get(id_col):
                job_identifier = seq_info[id_col]
            else:
                job_identifier = generate_unique_id(10)
            
            # Auto-detect sequence type
            if seq_type == 'auto':
                seq_upper = sequence.upper() if sequence else ''
                if seq_upper and all(c in 'ACDEFGHIKLMNPQRSTVWY*X' for c in seq_upper):
                    seq_type = 'aminoacid'
                else:
                    seq_type = 'nucleotide'
            
            # Create job-specific config (deep copy to avoid issues)
            import copy
            job_config = copy.deepcopy(base_config)
            
            # Handle PDB ID from CSV
            if seq_info.get('pdb_id'):
                job_config['pdb_id'] = seq_info['pdb_id']
                job_config['pdb_chain'] = seq_info.get('chain', seq_info.get('pdb_chain'))
                job_config['sequence'] = None
                job_config['aa_sequence'] = None
            elif seq_type == 'aminoacid':
                job_config['aa_sequence'] = sequence
                job_config['sequence'] = None
                job_config['pdb_id'] = None
            else:
                job_config['sequence'] = sequence
                job_config['aa_sequence'] = None
                job_config['pdb_id'] = None
            
            job_config['fasta_file'] = None
            
            # Set identifier in config
            job_config['identifier'] = job_identifier
            
            # Set job-specific output directory
            job_dir = os.path.join(output_dir, name)
            os.makedirs(job_dir, exist_ok=True)
            job_config['output']['output_dir'] = job_dir
            
            # Write job config
            job_config_path = os.path.join(job_dir, 'config.yaml')
            with open(job_config_path, 'w') as f:
                yaml.dump(job_config, f, default_flow_style=False)
            
            # Create PBS script
            pbs_script = f"""#!/bin/bash
#PBS -N codon_opt_{name}
#PBS -q {queue}
#PBS -l walltime={walltime}
#PBS -l ncpus={ncpus}
#PBS -l mem={mem}
#PBS -o {job_dir}/job.out
#PBS -e {job_dir}/job.err
#PBS -j oe

# Change to working directory
cd $PBS_O_WORKDIR

# Activate environment if needed (uncomment and modify as needed)
# source /path/to/conda/etc/profile.d/conda.sh
# conda activate myenv

# Run codon optimization
{python_path} {script_path} optimize \\
    --config {os.path.abspath(job_config_path)} \\
    --job_name {name}

echo "Job {name} completed at $(date)"
"""
            
            pbs_path = os.path.join(job_dir, f'{name}.pbs')
            with open(pbs_path, 'w') as f:
                f.write(pbs_script)
            
            job_scripts.append(pbs_path)
            print(f"  Created job: {name}")
        
        # Create master submission script
        master_script = f"""#!/bin/bash
# Master script to submit all codon optimization jobs
# Generated by codon_optimizer_cli.py

cd {os.path.abspath(output_dir)}

echo "Submitting {len(job_scripts)} jobs..."

"""
        for pbs_path in job_scripts:
            master_script += f"qsub {os.path.abspath(pbs_path)}\n"
        
        master_script += f"""
echo "All jobs submitted!"
echo "Monitor with: qstat -u $USER"
"""
        
        master_path = os.path.join(output_dir, 'submit_all.sh')
        with open(master_path, 'w') as f:
            f.write(master_script)
        os.chmod(master_path, 0o755)
        
        print(f"\n{'='*60}")
        print(f"Batch job generation complete!")
        print(f"{'='*60}")
        print(f"Generated {len(job_scripts)} PBS scripts in: {output_dir}")
        print(f"\nTo submit all jobs, run:")
        print(f"  bash {master_path}")
        print(f"\nOr submit individual jobs with:")
        print(f"  qsub {job_scripts[0]}")
        
        return job_scripts
    
    def example(self, output: str = "example_config.yaml"):
        """
        Generate an example YAML configuration file.
        
        Args:
            output: Output path for the example config file
        """
        example_config = """# Codon Optimization Configuration File
# =====================================

# Unique Identifier
# -----------------
# Used for naming output folder (result_{identifier}) and files
# If not provided, a random 10-character ID will be generated
identifier: null  # e.g., "luciferase_v1" or leave null for auto-generation

# Sequence Source (use ONE of the following options)
# --------------------------------------------------

# Option 1: Direct nucleotide sequence
# sequence: "AUGUUUAAAGGG..."

# Option 2: Direct amino acid sequence (1-letter code)
# aa_sequence: "MFKG..."

# Option 3: PDB ID to fetch sequence from
pdb_id: "1BA3"
pdb_chain: null  # Optional: specify chain (e.g., "A")

# Option 4: FASTA file path
# fasta_file: "/path/to/sequence.fasta"

# Sequence Processing
# -------------------
five_prime_utr: "AUUUGGUGGAGG"  # 5' UTR to prepend (before start codon)
three_prime_utr: ""             # 3' UTR to append (after stop codon)
add_start_codon: true           # Add AUG if not present
add_stop_codon: true            # Add UAA if not present

# Species Configuration
# --------------------
species: "human"  # Options: human, e_coli, mouse, rat, yeast, other

# mRNA System Parameters
# ----------------------
T_K: 310                        # Temperature in Kelvin
modify_utr: false               # Allow mutations in UTR region
initial_region_end_index: 30    # End of initial region for hybridization penalty
beamsize: 100                   # LinearPartition beam size
bpp_cutoff: 0.01               # Base pair probability cutoff
start_from_optimal_cai: true    # Start from optimal CAI sequence
verbose: false                  # Verbose output

# Loss Function Weights
# --------------------
# Set weight to 0.0 to disable a component
loss_weights:
  cai: 10.0                     # Codon Adaptation Index (maximize)
  fe: 1.0                       # Free energy (minimize)
  cpg: 3.0                      # CpG content (minimize)
  stem: 3.0                     # Stem/hairpin penalty
  utr_hybridisation: 3.0        # UTR hybridization penalty
  initial_hybridisation: 3.0    # Initial region hybridization penalty
  # mfe: 0.0                    # Minimum free energy (alternative to fe)
  # restriction_sites: 0.0      # Restriction site penalty

# Optimization Parameters
# ----------------------
optimization:
  T_opt: 0.01                   # Initial temperature for simulated annealing
  nsteps: 100                   # Number of steps (multiplied by number of codons)
  sample_frequency: 100         # Frequency of progress output
  output_filename: "opt_statistics.txt"
  n_sample: null                # Steps between statistics saves (default: sample_frequency)
  use_average_stem_length: false
  average_stem_num_samples: 20

# Output Configuration
# -------------------
# Output will be saved to: {output_dir}/result_{identifier}/
output:
  output_dir: "."               # Base output directory
  save_structure_svg: true
  save_prob_matrix: true
"""
        
        with open(output, 'w') as f:
            f.write(example_config)
        
        print(f"Example configuration saved to: {output}")
        print("\nEdit this file to customize your optimization, then run:")
        print(f"  python codon_optimizer_cli.py optimize --config {output}")


def main():
    """Main entry point for the CLI."""
    fire.Fire(CodonOptimizerCLI)


if __name__ == "__main__":
    main()
