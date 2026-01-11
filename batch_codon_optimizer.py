#!/usr/bin/env python3
"""
Batch Codon Optimization Script

This script reads a CSV file containing multiple sequences and generates
PBS job scripts for parallel processing on an HPC cluster.

Usage:
    python batch_codon_optimizer.py sequences.csv config.yaml --output_dir jobs
    
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
    
    with open(csv_path, 'r', newline='') as f:
        # Use standard comma-delimited CSV (don't rely on sniffer - it's unreliable)
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

def generate_pbs_script(
    job_name: str,
    config_path: str,
    job_dir: str,
    script_path: str,
    python_path: str,
    queue: str = "hx",
    walltime: str = "72:00:00",
    ncpus: int = 1,
    mem: str = "8GB",
    modules: Optional[List[str]] = None,
    conda_env: Optional[str] = None,
) -> str:
    """Generate a PBS job script."""
    
    script = f"""#!/bin/bash
#PBS -N codon_{job_name[:12]}
#PBS -q {queue}
#PBS -l walltime={walltime}
#PBS -l ncpus={ncpus}
#PBS -l mem={mem}
#PBS -o {job_dir}/job.out
#PBS -e {job_dir}/job.err
#PBS -j oe

# Print job information
echo "Job started at $(date)"
echo "Running on host: $(hostname)"
echo "Working directory: $PBS_O_WORKDIR"
echo ""

# Change to submission directory
cd $PBS_O_WORKDIR

"""
    
    # Add module loads if specified
    if modules:
        script += "# Load modules\n"
        for module in modules:
            script += f"module load {module}\n"
        script += "\n"
    
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


def generate_slurm_script(
    job_name: str,
    config_path: str,
    job_dir: str,
    script_path: str,
    python_path: str,
    partition: str = "hx",
    time: str = "72:00:00",
    cpus: int = 1,
    mem: str = "8G",
    modules: Optional[List[str]] = None,
    conda_env: Optional[str] = None,
) -> str:
    """Generate a SLURM job script (alternative to PBS)."""
    
    script = f"""#!/bin/bash
#SBATCH --job-name=codon_{job_name[:12]}
#SBATCH --partition={partition}
#SBATCH --time={time}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --output={job_dir}/job.out
#SBATCH --error={job_dir}/job.err

# Print job information
echo "Job started at $(date)"
echo "Running on host: $(hostname)"
echo "Working directory: $SLURM_SUBMIT_DIR"
echo ""

# Change to submission directory
cd $SLURM_SUBMIT_DIR

"""
    
    if modules:
        script += "# Load modules\n"
        for module in modules:
            script += f"module load {module}\n"
        script += "\n"
    
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
        description="Generate batch job scripts for codon optimization",
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
    
    parser.add_argument("csv_file", help="CSV file containing sequences")
    parser.add_argument("config", help="YAML configuration file")
    parser.add_argument("--output_dir", "-o", default="batch_jobs",
                        help="Output directory for job files (default: batch_jobs)")
    parser.add_argument("--scheduler", choices=["pbs", "slurm"], default="pbs",
                        help="Job scheduler type (default: pbs)")
    parser.add_argument("--queue", "-q", default="hx",
                        help="Queue/partition name (default: hx)")
    parser.add_argument("--walltime", "-t", default="72:00:00",
                        help="Wall time limit (default: 72:00:00)")
    parser.add_argument("--ncpus", "-n", type=int, default=1,
                        help="Number of CPUs per job (default: 1)")
    parser.add_argument("--mem", "-m", default="8GB",
                        help="Memory per job (default: 8GB)")
    parser.add_argument("--modules", nargs="+", default=None,
                        help="Modules to load (space-separated)")
    parser.add_argument("--conda_env", default=None,
                        help="Conda environment to activate")
    parser.add_argument("--python_path", default=None,
                        help="Path to Python interpreter")
    parser.add_argument("--dry_run", action="store_true",
                        help="Print what would be done without creating files")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.csv_file):
        print(f"Error: CSV file not found: {args.csv_file}")
        sys.exit(1)
    
    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}")
        sys.exit(1)
    
    # Get paths
    python_path = args.python_path or sys.executable
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cli_script = os.path.join(script_dir, "codon_optimizer_cli.py")
    
    if not os.path.exists(cli_script):
        print(f"Error: CLI script not found: {cli_script}")
        sys.exit(1)
    
    # Read sequences
    sequences = read_sequences_csv(args.csv_file)
    print(f"Found {len(sequences)} sequences in {args.csv_file}")
    
    if len(sequences) == 0:
        print("No valid sequences found. Exiting.")
        sys.exit(1)
    
    if args.dry_run:
        print("\n[DRY RUN] Would generate the following jobs:")
        for seq in sequences:
            print(f"  - {seq['name']}: {seq['type']} sequence")
        sys.exit(0)
    
    # Load base config
    base_config = load_config(args.config)
    
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
        
        # Write job config
        config_path = os.path.join(job_dir, 'config.yaml')
        with open(config_path, 'w') as f:
            yaml.dump(job_config, f, default_flow_style=False)
        
        # Generate scheduler script
        if args.scheduler == "pbs":
            script = generate_pbs_script(
                job_name=name,
                config_path=os.path.abspath(config_path),
                job_dir=os.path.abspath(job_dir),
                script_path=cli_script,
                python_path=python_path,
                queue=args.queue,
                walltime=args.walltime,
                ncpus=args.ncpus,
                mem=args.mem,
                modules=args.modules,
                conda_env=args.conda_env,
            )
            script_ext = ".pbs"
        else:
            script = generate_slurm_script(
                job_name=name,
                config_path=os.path.abspath(config_path),
                job_dir=os.path.abspath(job_dir),
                script_path=cli_script,
                python_path=python_path,
                partition=args.queue,
                time=args.walltime,
                cpus=args.ncpus,
                mem=args.mem.replace('GB', 'G'),  # SLURM uses 'G' not 'GB'
                modules=args.modules,
                conda_env=args.conda_env,
            )
            script_ext = ".slurm"
        
        script_path = os.path.join(job_dir, f'{name}{script_ext}')
        with open(script_path, 'w') as f:
            f.write(script)
        os.chmod(script_path, 0o755)
        
        job_scripts.append(script_path)
        print(f"  Created: {name}")
    
    # Create master submission script
    submit_cmd = "qsub" if args.scheduler == "pbs" else "sbatch"
    
    master_script = f"""#!/bin/bash
# Master script to submit all codon optimization jobs
# Generated by batch_codon_optimizer.py
# Scheduler: {args.scheduler.upper()}

cd {os.path.abspath(args.output_dir)}

echo "Submitting {len(job_scripts)} codon optimization jobs..."
echo ""

"""
    
    for i, script_path in enumerate(job_scripts):
        master_script += f'echo "Submitting job {i+1}/{len(job_scripts)}: {os.path.basename(os.path.dirname(script_path))}"\n'
        master_script += f"{submit_cmd} {os.path.abspath(script_path)}\n"
        master_script += "sleep 0.5  # Small delay between submissions\n\n"
    
    master_script += f"""
echo ""
echo "All {len(job_scripts)} jobs submitted!"
echo "Monitor with: {'qstat -u $USER' if args.scheduler == 'pbs' else 'squeue -u $USER'}"
"""
    
    master_path = os.path.join(args.output_dir, f'submit_all.sh')
    with open(master_path, 'w') as f:
        f.write(master_script)
    os.chmod(master_path, 0o755)
    
    # Print summary
    print(f"\n{'='*60}")
    print("BATCH JOB GENERATION COMPLETE")
    print(f"{'='*60}")
    print(f"Generated {len(job_scripts)} job scripts in: {args.output_dir}")
    print(f"Scheduler: {args.scheduler.upper()}")
    print(f"\nTo submit all jobs, run:")
    print(f"  bash {master_path}")
    print(f"\nOr submit individual jobs:")
    print(f"  {submit_cmd} {job_scripts[0]}")


if __name__ == "__main__":
    main()
