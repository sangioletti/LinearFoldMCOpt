#!/usr/bin/env python3
"""
Standalone script to plot optimization statistics from opt_statistics.txt

Usage:
    python plot_statistics.py [statistics_file] [output_file]
    
Examples:
    python plot_statistics.py
    python plot_statistics.py opt_statistics.txt
    python plot_statistics.py opt_statistics.txt optimization_plot.png
"""

import sys
from mRNA import plot_optimization_statistics

if __name__ == "__main__":
    # Parse command line arguments
    statistics_file = "opt_statistics.txt"
    output_file = None
    
    if len(sys.argv) > 1:
        statistics_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    
    try:
        plot_optimization_statistics(statistics_file, output_file)
        if output_file:
            print(f"\n✓ Plot successfully saved to {output_file}")
        else:
            print("\n✓ Plot displayed (close window to continue)")
    except Exception as e:
        print(f"Error plotting statistics: {e}", file=sys.stderr)
        sys.exit(1)

