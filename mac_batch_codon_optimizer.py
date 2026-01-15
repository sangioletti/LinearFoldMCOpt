#!/usr/bin/env python3
"""
Mac/Local Batch Codon Optimization Script

This script reads a CSV file containing multiple sequences and generates
local bash scripts for parallel or sequential processing on a Mac or local workstation.
Derived from batch_codon_optimizer.py but simplified for local execution.

Usage:
    python mac_batch_codon_optimizer.py sequences.csv config.yaml --output_dir jobs
    
CSV Format:
    name,sequence,type
    protein1,MFKGV...,aminoacid
    protein2,AUGUUU...,nucleotide
    
The 'type' column is optional - if omitted, will auto-detect.
"""

import os
import sys
import csv
import yaml
import argparse
import shutil
from pathlib import Path
from typing import Optional, List, Dict, Any


def load_config(config_path: str) -> Dict[str, Any]:
    """Load YAML configuration file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def detect_sequence_type(sequence: str) -> str:
    """
    Auto-detect if a sequence is amino acid or nucleotide.
    
    Returns:
        'aminoacid' or 'nucleotide'
    """
    seq_upper = sequence.upper().replace(' ', '').replace('\n', '')
    
    # Amino acid characters (excluding nucleotides ACGU)
    aa_only_chars = set('DEFHIKLMNPQRSVWY*X')
    
    # If sequence contains any amino-acid-only characters, it's amino acid
    if any(c in aa_only_chars for c in seq_upper):
        return 'aminoacid'
    
    # If only contains ACGUT, it's nucleotide
    if all(c in 'ACGUT' for c in seq_upper):
        return 'nucleotide'
    
    # Default to amino acid for ambiguous cases
    return 'aminoacid'


def read_sequences_csv(csv_path: str) -> List[Dict[str, str]]:
    """
    Read sequences from CSV file.
    
    Expected columns:
    - name: Job identifier (required)
    - sequence: The sequence string (required)
    - type: 'aminoacid' or 'nucleotide' (optional, auto-detected)
    - pdb_id: PDB ID to fetch (optional, alternative to sequence)
    - chain: PDB chain (optional)
    
    Returns:
        List of dictionaries with sequence information
    """
    sequences = []
    
    if not os.path.exists(csv_path):
        print(f"Error: CSV file not found: {csv_path}")
        return []

    with open(csv_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        
        for i, row in enumerate(reader):
            # Normalize column names (lowercase, strip whitespace), skip None keys
            row = {k.lower().strip(): v.strip() if isinstance(v, str) else v 
                   for k, v in row.items() if k is not None}
            
            # Check for required fields
            if 'name' not in row:
                row['name'] = f'seq_{i:04d}'
            
            # Accept either 'sequence' or 'seq' column
            if 'sequence' not in row and 'seq' in row:
                row['sequence'] = row['seq']
            
            # Validate row has some sequence source
            if not row.get('sequence') and not row.get('pdb_id'):
                print(f"Warning: Row {i+1} has no sequence or pdb_id, skipping")
                continue
            
            # Auto-detect type if not specified
            if 'type' not in row or not row['type']:
                if row.get('sequence'):
                    row['type'] = detect_sequence_type(row['sequence'])
                else:
                    row['type'] = 'aminoacid'  # PDB sequences are amino acids
            
            sequences.append(row)
    
    return sequences


def generate_local_script(
    job_name: str,
    config_path: str,
    job_dir: str,
    script_path: str,
    python_path: str,
    conda_env: Optional[str] = None,
) -> str:
    """Generate a local bash job script."""
    
    script = f"""#!/bin/bash
# Local job script for {job_name}

# Print job information
echo "Job {job_name} started at $(date)"
echo "Running on host: $(hostname)"
echo "Working directory: {job_dir}"
echo ""

# Change to job directory
cd {job_dir}

"""
    
    # Add conda activation if specified
    if conda_env:
        script += f"""# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate {conda_env}
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

"""
    
    script += f"""# Run codon optimization
echo "Starting codon optimization for {job_name}..."
{python_path} {script_path} optimize \\
    --config {config_path} \\
    --job_name {job_name}

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "Job {job_name} completed successfully at $(date)"
else
    echo ""
    echo "Job {job_name} FAILED at $(date)"
    exit 1
fi
"""
    
    return script


def main():
    parser = argparse.ArgumentParser(
        description="Generate batch job scripts for local codon optimization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
CSV Format Example:
    name,sequence,type
    luciferase,MFKVLS...,aminoacid
    gfp,AUGAGUAAA...,nucleotide
    1ba3_protein,,pdb_id  # Will fetch from PDB

Or with PDB IDs:
    name,pdb_id,chain
    luciferase,1BA3,A
    gfp,1EMA,
"""
    )
    
    parser.add_argument("--csv_file", "-c", help="CSV file containing sequences")
    parser.add_argument("--config", "-g", default="config.yaml",
                        help="YAML configuration file (default: config.yaml)")
    parser.add_argument("--output_dir", "-o", default="local_batch_jobs",
                        help="Output directory for job files (default: local_batch_jobs)")
    parser.add_argument("--parallel", "-p", type=int, default=None,
                        help="Number of jobs to run in parallel in the master script (default: 1 or value in config)")
    parser.add_argument("--conda_env", default=None,
                        help="Conda environment to activate")
    parser.add_argument("--python_path", default=None,
                        help="Path to Python interpreter")
    parser.add_argument("--dry_run", action="store_true",
                        help="Print what would be done without creating files")
    
    args = parser.parse_args()
    
    # Validate config
    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}")
        sys.exit(1)
    
    # Load base config
    base_config = load_config(args.config)
    
    # Get CSV file path
    csv_file = args.csv_file or base_config.get('sequences_csv')
    if not csv_file:
        print("Error: No CSV file specified. Provide it via --csv_file or set 'sequences_csv' in config.yaml")
        sys.exit(1)
        
    if not os.path.exists(csv_file):
        print(f"Error: CSV file not found: {csv_file}")
        sys.exit(1)

    # Get parallel cores
    parallel_cores = args.parallel
    if parallel_cores is None:
        parallel_cores = base_config.get('parallel_cores', 1)

    # Get paths
    python_path = args.python_path or sys.executable
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cli_script = os.path.join(script_dir, "codon_optimizer_cli.py")
    
    if not os.path.exists(cli_script):
        print(f"Error: CLI script not found: {cli_script}")
        sys.exit(1)
    
    # Read sequences
    sequences = read_sequences_csv(csv_file)
    print(f"Found {len(sequences)} sequences in {csv_file}")
    
    if len(sequences) == 0:
        print("No valid sequences found. Exiting.")
        sys.exit(1)
    
    if args.dry_run:
        print("\n[DRY RUN] Would generate the following jobs:")
        for seq in sequences:
            print(f"  - {seq['name']}: {seq['type']} sequence")
        sys.exit(0)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Generate job scripts
    job_scripts = []
    
    for seq_info in sequences:
        name = seq_info['name']
        
        # Create job directory
        job_dir = os.path.join(args.output_dir, name)
        os.makedirs(job_dir, exist_ok=True)
        
        # Create job-specific config
        job_config = base_config.copy() if base_config else {}
        
        # Set sequence based on type
        if seq_info.get('pdb_id'):
            job_config['pdb_id'] = seq_info['pdb_id']
            job_config['pdb_chain'] = seq_info.get('chain')
            job_config['sequence'] = None
            job_config['aa_sequence'] = None
        elif seq_info['type'] == 'aminoacid':
            job_config['aa_sequence'] = seq_info.get('sequence', '')
            job_config['sequence'] = None
            job_config['pdb_id'] = None
        else:
            job_config['sequence'] = seq_info.get('sequence', '')
            job_config['aa_sequence'] = None
            job_config['pdb_id'] = None
        
        job_config['fasta_file'] = None
        
        # Set output directory
        if 'output' not in job_config:
            job_config['output'] = {}
        job_config['output']['output_dir'] = os.path.abspath(job_dir)
        
        # Merge all other CSV columns into job config
        for k, v in seq_info.items():
            if k not in ['name', 'sequence', 'seq', 'type', 'pdb_id', 'chain'] and v:
                job_config[k] = v

        # Copy plasmid file if specified in the config
        if job_config.get('plasmid_file'):
            plasmid_filename = job_config['plasmid_file']
            # Try to find the file relative to the base config file directory,
            # relative to the CSV file directory, relative to current directory, 
            # or as an absolute path
            config_dir = os.path.dirname(os.path.abspath(args.config))
            csv_dir = os.path.dirname(os.path.abspath(csv_file))
            
            potential_paths = [
                os.path.join(config_dir, plasmid_filename),
                os.path.join(csv_dir, plasmid_filename),
                os.path.abspath(plasmid_filename),
                plasmid_filename
            ]
            
            src_plasmid = None
            for p in potential_paths:
                if os.path.exists(p):
                    src_plasmid = p
                    break
            
            if src_plasmid:
                dst_plasmid = os.path.join(job_dir, os.path.basename(plasmid_filename))
                # Only copy if it doesn't already exist or if we want to overwrite
                shutil.copy2(src_plasmid, dst_plasmid)
                # Ensure the job-specific config uses the base filename
                job_config['plasmid_file'] = os.path.basename(plasmid_filename)
            else:
                print(f"Warning: Plasmid file '{plasmid_filename}' not found at any of: {potential_paths}")

        # Write job config
        config_path = os.path.join(job_dir, 'config.yaml')
        with open(config_path, 'w') as f:
            yaml.dump(job_config, f, default_flow_style=False)
        
        # Generate local script
        script = generate_local_script(
            job_name=name,
            config_path=os.path.abspath(config_path),
            job_dir=os.path.abspath(job_dir),
            script_path=cli_script,
            python_path=python_path,
            conda_env=args.conda_env,
        )
        
        script_path = os.path.join(job_dir, f'run_{name}.sh')
        with open(script_path, 'w') as f:
            f.write(script)
        os.chmod(script_path, 0o755)
        
        job_scripts.append(script_path)
        print(f"  Created: {name}")
    
    # Create master execution script
    master_script = f"""#!/bin/bash
# Master script to run all codon optimization jobs locally
# Generated by mac_batch_codon_optimizer.py

cd {os.path.abspath(args.output_dir)}

echo "Starting {len(job_scripts)} codon optimization jobs..."
echo "Output directory: $(pwd)"
echo ""

"""
    
    if parallel_cores > 1:
        master_script += f"""# Running in parallel (max {parallel_cores} jobs)
# Using a basic backgrounding and wait mechanism
MAX_JOBS={parallel_cores}
COUNT=0

for script in {" ".join([os.path.abspath(s) for s in job_scripts])}; do
    echo "Launching $(basename $script)..."
    bash $script > "${{script%.sh}}.log" 2>&1 &
    COUNT=$((COUNT+1))
    
    # Simple rate limiting for parallel jobs
    if (( COUNT % MAX_JOBS == 0 )); then
        echo "Waiting for batch of $MAX_JOBS jobs to complete..."
        wait
    fi
done

wait
echo ""
echo "All parallel jobs finished execution."
"""
    else:
        for i, script_path in enumerate(job_scripts):
            name = os.path.basename(os.path.dirname(script_path))
            master_script += f'echo "Running job {i+1}/{len(job_scripts)}: {name}"\n'
            master_script += f"bash {os.path.abspath(script_path)}\n"
            master_script += "echo \"\"\n"
        
        master_script += "echo \"All sequential jobs finished execution.\"\n"

    master_path = os.path.join(args.output_dir, f'run_all.sh')
    with open(master_path, 'w') as f:
        f.write(master_script)
    os.chmod(master_path, 0o755)
    
    # Print summary
    print(f"\n{'='*60}")
    print("LOCAL BATCH JOB GENERATION COMPLETE")
    print(f"{'='*60}")
    print(f"Generated {len(job_scripts)} job scripts in: {args.output_dir}")
    if parallel_cores > 1:
        print(f"Parallel mode enabled: {parallel_cores} simultaneous jobs")
    else:
        print(f"Sequential mode enabled")
    print(f"\nTo run all jobs, execute:")
    print(f"  bash {master_path}")
    print(f"\nOr run individual jobs:")
    print(f"  bash {job_scripts[0]}")


if __name__ == "__main__":
    main()

