# LinearFold Monte Carlo Codon Optimizer

A Python package for optimizing mRNA codon usage using Monte Carlo simulated annealing, with RNA secondary structure prediction powered by LinearPartition.

## Overview

This tool optimizes mRNA sequences to balance multiple objectives:
- **Codon Adaptation Index (CAI)**: Maximize expression in the target organism
- **Minimum Free Energy (MFE)**: Control RNA secondary structure stability
- **CpG Content**: Minimize immunogenic CpG dinucleotides
- **Stem Length**: Avoid long hairpin structures that may impede translation
- **UTR Hybridization**: Prevent coding sequence from interfering with UTR elements

## Features

- 🧬 **Multi-objective optimization** with configurable loss weights
- ⚡ **Fast RNA folding** using LinearPartition (O(n²) complexity)
- 📊 **Detailed statistics** tracking during optimization
- 🔧 **YAML-based configuration** for reproducible experiments
- 🖥️ **HPC batch processing** with PBS/SLURM script generation
- 🔗 **PDB integration** to fetch sequences directly from protein structures

## Installation

### Prerequisites

- Python 3.8+
- Conda (recommended for environment management)

### Quick Install

```bash
# Clone the repository
git clone https://github.com/sangiole/protein_design.git
cd protein_design/LinearFoldMCOpt

# Create and activate conda environment
conda create -n codon-opt python=3.10 numpy matplotlib pyyaml -y
conda activate codon-opt

# Install additional dependencies
pip install fire biotite

# Optional: Install as editable package
pip install -e .
```

### Using Existing Environment

If you have the `bagel-recipes` conda environment:

```bash
conda activate bagel-recipes
```

## Quick Start

### 1. Generate Example Configuration

```bash
python codon_optimizer_cli.py example --output my_config.yaml
```

### 2. Run Optimization

```bash
# Using a PDB ID
python codon_optimizer_cli.py optimize --config my_config.yaml

# Using a custom sequence
python codon_optimizer_cli.py optimize --config my_config.yaml --sequence "MFKVL..."

# Specify output directory
python codon_optimizer_cli.py optimize --config my_config.yaml --output_dir results/
```

### 3. Fetch PDB Sequence

```bash
# Get amino acid sequence from PDB
python codon_optimizer_cli.py fetch_pdb --pdb_id 1BA3

# Convert to nucleotides
python codon_optimizer_cli.py fetch_pdb --pdb_id 1BA3 --convert_to_nucleotide

# Save to file
python codon_optimizer_cli.py fetch_pdb --pdb_id 1BA3 --output sequence.fasta
```

## Configuration File

The YAML configuration file controls all aspects of the optimization. Here's a complete example:

```yaml
# Sequence Source (use ONE option)
# --------------------------------
# Option 1: Direct nucleotide sequence
sequence: "AUGUUUAAAGGG..."

# Option 2: Direct amino acid sequence
aa_sequence: "MFKG..."

# Option 3: PDB ID
pdb_id: "1BA3"
pdb_chain: "A"  # Optional

# Option 4: FASTA file
fasta_file: "/path/to/sequence.fasta"

# Sequence Processing
# -------------------
five_prime_utr: "AUUUGGUGGAGG"
add_start_codon: true
add_stop_codon: true

# Species
# -------
species: "human"  # Options: human, e_coli

# mRNA Parameters
# ---------------
T_K: 310                      # Temperature (Kelvin)
modify_utr: false             # Allow UTR mutations
initial_region_end_index: 30  # For hybridization penalty
beamsize: 100                 # LinearPartition accuracy
bpp_cutoff: 0.01             # Base pair probability threshold
start_from_optimal_cai: true  # Initialize with optimal codons
verbose: false

# Loss Weights
# ------------
loss_weights:
  cai: 10.0                   # Codon Adaptation Index
  fe: 1.0                     # Free energy
  cpg: 3.0                    # CpG dinucleotides
  stem: 3.0                   # Stem/hairpin penalty
  utr_hybridisation: 3.0      # UTR interference
  initial_hybridisation: 3.0  # Initial region interference
  # mfe: 0.0                  # Alternative to fe
  # restriction_sites: 0.0    # Restriction enzyme sites

# Optimization Parameters
# -----------------------
optimization:
  T_opt: 0.01                 # Initial annealing temperature
  nsteps: 100                 # Steps per codon
  sample_frequency: 100       # Progress output frequency
  output_filename: "opt_statistics.txt"
  use_average_stem_length: false
  average_stem_num_samples: 20

# Output Configuration
# --------------------
output:
  output_dir: "./optimization_results"
  save_structure_svg: true
  save_prob_matrix: true
```

## Batch Processing

For processing multiple sequences on an HPC cluster:

### 1. Create a CSV file with sequences

```csv
name,sequence,type
luciferase,,pdb_id
gfp,MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLK...,aminoacid
my_gene,AUGUUUAAAGGGCCCCCC...,nucleotide
```

Or with PDB IDs:

```csv
name,pdb_id,chain
luciferase,1BA3,A
gfp,1EMA,
```

### 2. Generate PBS/SLURM scripts

```bash
# Generate PBS scripts
python batch_codon_optimizer.py sequences.csv config.yaml \
    --output_dir batch_jobs \
    --scheduler pbs \
    --queue normal \
    --walltime 24:00:00 \
    --ncpus 1 \
    --mem 8GB \
    --conda_env bagel-recipes

# Generate SLURM scripts
python batch_codon_optimizer.py sequences.csv config.yaml \
    --output_dir batch_jobs \
    --scheduler slurm \
    --queue normal \
    --walltime 24:00:00
```

### 3. Submit all jobs

```bash
bash batch_jobs/submit_all.sh
```

### Batch Script Options

| Option | Description | Default |
|--------|-------------|---------|
| `--output_dir` | Output directory for jobs | `batch_jobs` |
| `--scheduler` | Job scheduler (pbs/slurm) | `pbs` |
| `--queue` | Queue/partition name | `normal` |
| `--walltime` | Wall time limit | `24:00:00` |
| `--ncpus` | CPUs per job | `1` |
| `--mem` | Memory per job | `8GB` |
| `--modules` | Modules to load | None |
| `--conda_env` | Conda environment | None |
| `--dry_run` | Preview without creating | False |

## Output Files

After optimization, you'll find:

| File | Description |
|------|-------------|
| `opt_statistics.txt` | Step-by-step optimization statistics |
| `optimized_sequence.txt` | Final optimized mRNA sequence |
| `structure.svg` | RNA secondary structure visualization |
| `prob_matrix.{step}` | Base pair probability matrices |
| `structure.{step}` | Structure in dot-bracket notation |
| `loss_components.{step}` | Individual loss components |
| `sequence_{step}.txt` | Codon sequence at each checkpoint |

## CLI Reference

### Main Commands

```bash
# Show help
python codon_optimizer_cli.py --help

# Run optimization
python codon_optimizer_cli.py optimize --config CONFIG [--sequence SEQ] [--output_dir DIR] [--job_name NAME]

# Fetch PDB sequence
python codon_optimizer_cli.py fetch_pdb --pdb_id ID [--chain CHAIN] [--output FILE] [--convert_to_nucleotide]

# Generate batch scripts
python codon_optimizer_cli.py batch --csv_file CSV --config CONFIG [--output_dir DIR] [OPTIONS]

# Generate example config
python codon_optimizer_cli.py example [--output FILE]
```

## Loss Function

The optimization minimizes a weighted sum of objectives:

```
L = -w_cai × log(CAI) + w_fe × FE + w_cpg × N_cpg + w_stem × L_stem + w_utr × P_utr + w_init × P_init
```

Where:
- **CAI**: Codon Adaptation Index (logarithm, so we negate to maximize)
- **FE**: Free energy from LinearPartition
- **N_cpg**: Number of CpG dinucleotides
- **L_stem**: Length of longest stem in predicted structure
- **P_utr**: Sum of base pair probabilities between coding and UTR regions
- **P_init**: Sum of base pair probabilities in initial region

## Algorithm

The optimizer uses **Simulated Annealing** with a linear temperature schedule:

1. Start with optimal CAI sequence (fastest initial translation)
2. At each step:
   - Randomly select a mutable codon position
   - Propose a synonymous codon change
   - Accept if loss decreases, or with probability exp(-ΔL/T)
3. Temperature decreases linearly from T_opt to 0

## Performance Tips

- **Increase beamsize** (e.g., 200) for more accurate folding at the cost of speed
- **Reduce nsteps** for faster iterations during testing
- **Use `start_from_optimal_cai: true`** to begin with a good initial sequence
- **Set unused weights to 0** to skip unnecessary calculations
- **Increase bpp_cutoff** (e.g., 0.05) to reduce BPP storage overhead

## Troubleshooting

### ImportError: No module named 'fire'

```bash
pip install fire
```

### ImportError: No module named 'biotite'

```bash
pip install biotite
```

### LinearPartition not found

Ensure the `linearpartition_wrapper.py` is in the same directory and LinearPartition is properly compiled.

### PBS job fails immediately

Check that:
1. The conda environment exists on compute nodes
2. Paths in the PBS script are absolute
3. Required modules are loaded

## License

MIT License - See [LICENSE](LICENSE) for details.

## Citation

If you use this tool, please cite:

```bibtex
@software{linearfold_mcopt,
  author = {Sangiole},
  title = {LinearFold Monte Carlo Codon Optimizer},
  year = {2025},
  url = {https://github.com/sangiole/protein_design}
}
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.