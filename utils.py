import numpy as np

from codons import (
    amino_acid_1L_to_full_name,
    amino_acid_3L_to_full_name,
    codon_table,
    amino_acid_to_1_letter
)

def get_sequence_as_dna(sequence):
        """
        Convert the mRNA sequence to DNA sequence.
        
        Returns:
            str: DNA sequence
        """
        # Make hard copy to make sure you are not changing the original 
        sequence = sequence.copy()
        return sequence.replace('U', 'T')

def codon_to_amino_acid_3L(codon):
    return amino_acid_to_1_letter[codon_table[codon]]

def codon_to_amino_acid_1L(codon):
    return amino_acid_to_1_letter[codon_table[codon]]

def aa_3L_to_codons(amino_acid):
    try:
        assert len(amino_acid) == 3, "Amino acid in 3-letter code expected"
    except AssertionError:
        assert amino_acid == "Stop", "Amino acid must be a 3-letter code"
        raise ValueError
    # Transform 3-letter code to full name
    amino_acid = amino_acid_3L_to_full_name[amino_acid]
    return tuple(codon for codon, aa in codon_table.items() if aa == amino_acid)

def aa_1L_to_codons(amino_acid):
    try:
        assert len(amino_acid) == 1, "Amino acid in 1-letter code expected"
    except AssertionError:
        assert amino_acid == "Stop", "Amino acid must be a 1-letter code"
        raise ValueError
    # Transform 1-letter code to full name
    amino_acid = amino_acid_1L_to_full_name[amino_acid]
    return tuple(codon for codon, aa in codon_table.items() if aa == amino_acid)


def aa_to_codon_sequence(aa_sequence):
    # Convert amino acid 1-letter code to nucleotide sequence
    # Use aa_1L_to_codons to get possible codons for each amino acid
    # Randomly select one of the possible codons for each amino acid
    nucleotide_seq = ""
    for aa in aa_sequence:
        if aa == "*" or aa == "X":
            # Stop codon or unknown, use UAA
            codon = "UAA"
        else:
            # Get possible codons for this amino acid
            possible_codons = aa_1L_to_codons(aa)
            if not possible_codons:
                raise ValueError(f"No codon found for amino acid: {aa}")
            codon = np.random.choice(possible_codons)
        nucleotide_seq += codon
    return nucleotide_seq

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
        ax.plot(steps, data['Acceptance_Rate'], 'bo', linestyle='none', markersize=4)
        ax.set_xlabel('Step', fontsize=11)
        ax.set_ylabel('Acceptance Rate', fontsize=11)
        ax.set_title('Acceptance Rate Over Optimization', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_ylim([0, 1.1])
        plot_idx += 1
    
    # Plot 2: Sequence Identity
    if 'Sequence_Identity_%' in data and len(data['Sequence_Identity_%']) > 0:
        ax = axes[plot_idx]
        ax.plot(steps, data['Sequence_Identity_%'], 'gs', linestyle='none', markersize=4)
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
            ax.plot(steps, values, color=color, marker='.', linestyle='none', markersize=4, label=component)
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
