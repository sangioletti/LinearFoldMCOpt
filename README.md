# LinearFold Monte Carlo Codon Optimizer

A Python package for optimizing mRNA codon usage using Monte Carlo simulated annealing, with RNA secondary structure prediction powered by LinearPartition.

## Overview

This tool optimizes mRNA sequences to balance multiple objectives:
- **Codon Adaptation Index (CAI)**: Maximize expression in the target organism
- **Minimum Free Energy (MFE)**: Control RNA secondary structure stability
- **CpG Content**: Minimize immunogenic CpG dinucleotides
- **Stem Length**: Avoid long hairpin structures that may impede translation
- **UTR Hybridization**: Prevent coding sequence from interfering with UTR elements

## Directory Structure

```
LinearFoldMCOpt/
├── codon_optimizer_cli.py      # Main CLI for single sequence optimization
├── batch_codon_optimizer.py    # Batch processing for HPC clusters
├── config.yaml                 # Default configuration
├── example_config.yaml         # Example configuration with comments
├── *.csv                       # Sequence data files
│
├── src/                        # Core Python modules
│   ├── __init__.py
│   ├── mRNA.py                 # Main mRNA class
│   ├── codons.py               # Codon tables and CAI values
│   ├── utils.py                # Utility functions
│   ├── linearpartition_wrapper.py
│   ├── linearpartition_native.py
│   └── gflags.py
│
├── bin/                        # Compiled binaries
│   ├── linearpartition         # Wrapper script
│   ├── linearpartition_c       # LinearPartition (CONTRAfold params)
│   └── linearpartition_v       # LinearPartition (Vienna params)
│
├── linearpartition/            # LinearPartition C++ source
│
├── scripts/                    # Utility scripts
│   ├── plot_statistics.py
│   └── clean_all_output.x
│
└── analysis_scripts/           # Analysis utilities
```

## Installation

### Prerequisites

- Python 3.8+
- C++ compiler (g++ >= 4.8)
- Conda (recommended for environment management)

### Quick Install

```bash
# Clone the repository
git clone https://github.com/sangiole/protein_design.git
cd protein_design/LinearFoldMCOpt

# Create and activate conda environment
conda create -n codon-opt python=3.10 -y
conda activate codon-opt

# Install Python dependencies
pip install -e .

# Build LinearPartition binaries
make
```

### Build LinearPartition

The optimization requires LinearPartition binaries. Build them with:

```bash
# Build both CONTRAfold and Vienna parameter versions
make

# Or build individually:
make linearpartition_c    # CONTRAfold parameters
make linearpartition_v    # Vienna parameters (recommended)

# Clean build artifacts
make clean
```

The compiled binaries will be placed in `bin/`.

### Verify Installation

```bash
# Check that LinearPartition works
./bin/linearpartition_v -h

# Run a test optimization
python codon_optimizer_cli.py optimize --config example_config.yaml
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

The YAML configuration file controls all aspects of optimization. See `example_config.yaml` for a complete example with comments.

Key options:

```yaml
# Sequence source (use one)
aa_sequence: "MFKG..."           # Amino acid sequence
pdb_id: "1BA3"                   # Or fetch from PDB

# UTR sequences
five_prime_utr: "GCCACCATG"      # 5' UTR (must end with Kozak + AUG)
three_prime_utr: "..."           # 3' UTR

# Optimization
species: "human"
loss_weights:
  cai: 10.0                      # Codon Adaptation Index
  fe: 1.0                        # Free energy
  cpg: 3.0                       # CpG content
  stem: 3.0                      # Stem penalty
  codon_divergence: 10.0         # Match expected codon usage

optimization:
  T_opt: 0.01                    # Annealing temperature
  nsteps: 100                    # Steps per codon
```

## Batch Processing

For processing multiple sequences on an HPC cluster:

```bash
# Generate PBS scripts
python batch_codon_optimizer.py sequences.csv config.yaml \
    --output_dir batch_jobs \
    --scheduler pbs \
    --queue hx \
    --walltime 72:00:00 \
    --conda_env codon-opt

# Submit all jobs
bash batch_jobs/submit_all.sh
```

## Dependencies

All dependencies are listed in `pyproject.toml`:

- **numpy** >= 1.20.0
- **fire** >= 0.4.0 (CLI framework)
- **pyyaml** >= 5.4.0 (configuration)
- **matplotlib** >= 3.3.0 (plotting)
- **requests** >= 2.25.0 (PDB fetching)
- **biotite** >= 0.30.0 (optional, for PDB parsing)

Install all with:
```bash
pip install -e ".[all]"
```

## License

MIT License - See [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.