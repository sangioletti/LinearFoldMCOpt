"""
Plasmid file parser module.

Provides functions to read and write plasmid files in GenBank format,
compatible with ApE (A plasmid Editor) and other sequence editors.

Note: ApE can save files in two formats:
1. Native ApE binary format (.dna) - NOT supported by this module
2. GenBank text format (.dna or .gb) - Supported

To convert ApE binary to GenBank text:
- Open the file in ApE
- File → Export → GenBank format
"""

from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, field
from typing import Optional, Union


try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


def read_genbank(file_path: str | Path) -> tuple[str, dict]:
    """
    Parse a GenBank format file and extract sequence and features.
    
    Args:
        file_path: Path to the GenBank file (.dna, .gb, .gbk)
    
    Returns:
        tuple: (sequence_string, features_dict)
            - sequence_string: The DNA sequence as a string
            - features_dict: Dictionary of features with structure:
                {feature_name: {
                    'type': feature_type,
                    'start': start_position (0-based),
                    'end': end_position,
                    'strand': 1 (forward) or -1 (reverse),
                    'sequence': feature_sequence,
                    'qualifiers': dict of additional qualifiers
                }}
    
    Raises:
        ImportError: If BioPython is not installed
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file is in binary format (not GenBank text)
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError(
            "BioPython is required to parse GenBank files. "
            "Install with: pip install biopython"
        )
    
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    # Check if file is binary (ApE native format)
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            first_line = f.readline()
            if not first_line.startswith('LOCUS'):
                raise ValueError(
                    f"File '{file_path}' appears to be in ApE binary format, not GenBank text. "
                    "Please export from ApE as GenBank format: File → Export → GenBank"
                )
    except UnicodeDecodeError:
        raise ValueError(
            f"File '{file_path}' appears to be in binary format, not GenBank text. "
            "Please export from ApE as GenBank format: File → Export → GenBank"
        )
    
    # Parse the GenBank file
    record = SeqIO.read(file_path, "genbank")
    
    # Extract sequence
    sequence = str(record.seq)
    
    # Extract features
    features_dict = {}
    for feature in record.features:
        # Get feature name from label qualifier, or use type as fallback
        label = feature.qualifiers.get('label', [feature.type])[0]
        
        # Handle duplicate feature names by appending a number
        original_label = label
        counter = 1
        while label in features_dict:
            label = f"{original_label}_{counter}"
            counter += 1
        
        features_dict[label] = {
            'type': feature.type,
            'start': int(feature.location.start),  # 0-based
            'end': int(feature.location.end),
            'strand': feature.location.strand,  # 1 for forward, -1 for reverse
            'sequence': str(feature.extract(record.seq)),
            'qualifiers': dict(feature.qualifiers)
        }
    
    return sequence, features_dict


def read_plasmid(file_path: str | Path) -> str:
    """
    Read just the sequence from a GenBank file.
    
    Args:
        file_path: Path to the GenBank file
    
    Returns:
        str: The DNA sequence
    """
    sequence, _ = read_genbank(file_path)
    return sequence


def get_features(file_path: str | Path) -> dict:
    """
    Read just the features from a GenBank file.
    
    Args:
        file_path: Path to the GenBank file
    
    Returns:
        dict: Dictionary of features
    """
    _, features = read_genbank(file_path)
    return features


@dataclass(frozen=True)
class Plasmid:
    """Immutable representation of a plasmid with sequence and features."""
    sequence: str
    locus: str = "plasmid"
    definition: str = ""
    organism: str = "synthetic"
    circular: bool = True
    features: Optional[dict] = field(default=None)
    
    def __len__(self) -> int:
        """Return the length of the plasmid sequence."""
        return len(self.sequence)
    
    @classmethod
    def from_genbank(cls, file_path: Union[str, Path]) -> "Plasmid":
        """
        Create a Plasmid from a GenBank file.
        
        Args:
            file_path: Path to the GenBank file
            
        Returns:
            Plasmid: A new Plasmid instance
        """
        sequence, features = read_genbank(file_path)
        return cls(
            sequence=sequence,
            features=features,
            locus=Path(file_path).stem
        )
    
    def to_genbank(self, output_file: Union[str, Path]) -> None:
        """
        Write the plasmid to a GenBank file.
        
        Args:
            output_file: Output file path
        """
        write_genbank(
            sequence=self.sequence,
            output_file=output_file,
            features=self.features,
            locus=self.locus,
            definition=self.definition,
            organism=self.organism,
            circular=self.circular
        )
    
    def get_feature_sequence(self, feature_name: str) -> Optional[str]:
        """Get the sequence of a specific feature."""
        if self.features and feature_name in self.features:
            return self.features[feature_name].get('sequence')
        return None
    
    def print_info(self) -> None:
        """Print a summary of the plasmid."""
        print(f"Plasmid: {self.locus}")
        print(f"Total Length: {len(self)} bp")
        print(f"Circular: {self.circular}")
        if self.features:
            print(f"Number of Features: {len(self.features)}")

    def merge_plasmid(self, insert_plasmid: "Plasmid", start_cut: int, end_cut: int ) -> "Plasmid": 
        """Merge a plasmid with an insert at a specific location."""
        if start_cut < 0 or end_cut > len(self.sequence):
            raise ValueError(f"Start and end cut positions must be within the plasmid sequence length. Start: {start_cut}, End: {end_cut}, Length: {len(self.sequence)}")
        if start_cut > end_cut:
            raise ValueError(f"Start cut position must be less than end cut position. Start: {start_cut}, End: {end_cut}")
        if start_cut == end_cut:
            raise ValueError(f"Start and end cut positions must be different. Start: {start_cut}, End: {end_cut}")
        
        # Check for features that partially overlap the cut region
        if self.features:
            for feature_name, feature_dict in self.features.items():
                feat_start = feature_dict['start']
                feat_end = feature_dict['end']
                # Check if feature is completely within cut region (OK to remove)
                completely_within = feat_start >= start_cut and feat_end <= end_cut
                # Check if feature partially overlaps (ERROR)
                partially_overlaps = (
                    (feat_start < start_cut < feat_end) or  # Starts before, ends within
                    (feat_start < end_cut < feat_end) or    # Starts within, ends after
                    (feat_start < start_cut and feat_end > end_cut)  # Spans entire cut region
                )
                if partially_overlaps:
                    raise ValueError(
                        f"Feature '{feature_name}' (position {feat_start}..{feat_end}) partially overlaps "
                        f"with the cut region ({start_cut}..{end_cut}). Features must either be completely "
                        f"within the cut region (and will be removed) or completely outside it."
                    )
        
        new_sequence = self.sequence[:start_cut] + insert_plasmid.sequence + self.sequence[end_cut:]
        new_features = {}
        
        # Calculate the shift amount for features after the cut region
        cut_length = end_cut - start_cut
        insert_length = len(insert_plasmid.sequence)
        shift_amount = insert_length - cut_length
        
        # Process features from the base plasmid
        if self.features:
            for feature_name, feature_dict in self.features.items():
                feat_start = feature_dict['start']
                feat_end = feature_dict['end']
                
                # Remove features completely within the cut region
                if feat_start >= start_cut and feat_end <= end_cut:
                    continue  # Skip this feature (it's being removed)
                
                # Adjust positions for features after the cut region
                if feat_start >= end_cut:
                    # Feature is completely after the cut region, shift it
                    adjusted_feature = feature_dict.copy()
                    adjusted_feature['start'] = feat_start + shift_amount
                    adjusted_feature['end'] = feat_end + shift_amount
                    new_features[feature_name] = adjusted_feature
                else:
                    # Feature is before the cut region, keep it as is
                    new_features[feature_name] = feature_dict.copy()

        # Add features from the insert plasmid
        if insert_plasmid.features:
            for feature_name, feature_dict in insert_plasmid.features.items():
                new_feature_name = feature_name
                if new_feature_name in new_features:
                    counter = 1
                    while f"{feature_name}_{counter}" in new_features:
                        counter += 1
                    new_feature_name = f"{feature_name}_{counter}"
                new_features[new_feature_name] = feature_dict.copy()
        
        return Plasmid(
            sequence=new_sequence,
            features=new_features if new_features else None,
            locus=self.locus
        )

def write_genbank(
    sequence: str,
    output_file: str | Path,
    features: Optional[dict] = None,
    locus: str = "plasmid",
    definition: str = "Designed by Stefano Angioletti-Uberti",
    organism: str = "synthetic",
    circular: bool = True
) -> None:
    """
    Write a sequence and features to a GenBank format file.
    
    Args:
        sequence: DNA sequence string
        output_file: Output file path
        features: Optional dictionary of features (same format as read_genbank output)
        locus: Locus name (appears in header)
        definition: Description of the plasmid
        organism: Source organism
        circular: Whether the plasmid is circular (default True)
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError(
            "BioPython is required to write GenBank files. "
            "Install with: pip install biopython"
        )
    
    # Clean sequence (convert U to T for DNA)
    sequence = sequence.upper().replace('U', 'T')
    
    # Create SeqRecord
    record = SeqRecord(
        Seq(sequence),
        id=locus,
        name=locus,
        description=definition,
        annotations={
            "molecule_type": "DNA",
            "topology": "circular" if circular else "linear",
            "organism": organism,
            "date": datetime.now().strftime("%d-%b-%Y").upper()
        }
    )
    
    # Add features
    if features:
        for name, info in features.items():
            feature = SeqFeature(
                location=FeatureLocation(
                    start=info.get('start', 0),
                    end=info.get('end', 0),
                    strand=info.get('strand', 1)
                ),
                type=info.get('type', 'misc_feature'),
                qualifiers={'label': [name]}
            )
            record.features.append(feature)
    
    # Write to file
    output_file = Path(output_file)
    with open(output_file, 'w') as f:
        SeqIO.write(record, f, "genbank")


def print_plasmid_info(file_path: str | Path) -> None:
    """
    Print a summary of a plasmid file.
    
    Args:
        file_path: Path to the GenBank file
    """
    sequence, features = read_genbank(file_path)
    
    print(f"Plasmid: {file_path}")
    print(f"Total Length: {len(sequence)} bp")
    print(f"Number of Features: {len(features)}")
    print("\nFeatures:")
    print("-" * 60)
    
    for name, info in features.items():
        strand = "+" if info['strand'] == 1 else "-"
        print(f"  [{name}]")
        print(f"    Type: {info['type']}")
        print(f"    Location: {info['start']}..{info['end']} ({strand})")
        print(f"    Length: {info['end'] - info['start']} bp")

def define_inserted_car_features( fivep_utr_seq, threep_utr_seq, car_cds_seq):
    car_features = {}
    car_features['car_fivep_utr'] = {
        'type': 'misc_feature',
        'start': 0,
        'end': len(fivep_utr_seq),
        'strand': 1,
        'sequence': fivep_utr_seq,
        'color': '#FFFF00',  # Yellow color
        'qualifiers': {'label': '5p UTR'}
    }
    car_features['kozak_car'] = {
        'type': 'misc_feature',
        'start': len(fivep_utr_seq)-9,
        'end': len(fivep_utr_seq),
        'strand': 1,
        'sequence': fivep_utr_seq[-9:],
        'color': '#FFFFFF',  
        'qualifiers': {'label': 'kozak_car'}
    }
    car_features['car_start_codon'] = {
        'type': 'misc_feature',
        'start': len(fivep_utr_seq)-3,
        'end': len(fivep_utr_seq),
        'strand': 1,
        'sequence': fivep_utr_seq[-3:],
        'color': '#FF0000',  # Red color
        'qualifiers': {'label': 'car start codon'}
    }
    car_features['car_CDS'] = {
        'type': 'misc_feature',
        'start': len(fivep_utr_seq)-3, # We put the start codon at the end of the 5p UTR
        'end': len(fivep_utr_seq) + len(car_cds_seq),
        'strand': 1,
        'sequence': car_cds_seq,
        'color': '#0000FF',  # Red color
        'qualifiers': {'label': 'car_CDS'}
    }
    car_features['car_stop_codon'] = {
        'type': 'misc_feature',
        'start': len(fivep_utr_seq) + len(car_cds_seq) -3,
        'end': len(fivep_utr_seq) + len(car_cds_seq),
        'strand': 1,
        'sequence': car_cds_seq[-3:],
        'color': '#FF0000',  # Red color
        'qualifiers': {'label': 'kozak_car'}
    }
    car_features['threep_utr'] = {
        'type': 'misc_feature',
        'start': len(fivep_utr_seq) + len(car_cds_seq),
        'end': len(fivep_utr_seq) + len(car_cds_seq) + len(threep_utr_seq),
        'strand': 1,
        'sequence': threep_utr_seq,
        'color': '#FFDD00',  # Yellow color
        'qualifiers': {'label': '3p UTR'}
    }
    return car_features

def move_features_to_insertion_point( start_nt, all_features ):
    """Redefine the position of all features by moving them to the new position
    after the insertion point defined by start_nt"""
    features_dict = all_features.copy()
    for feature in features_dict:
        features_dict[feature]['start'] += start_nt
        features_dict[feature]['end'] += start_nt
    return features_dict


# Main block for testing
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python plasmid.py <genbank_file>")
        print("\nExample:")
        print("  python plasmid.py my_plasmid.gb")
        sys.exit(1)
    
    file_path = sys.argv[1]
    
    try:
        print_plasmid_info(file_path)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
