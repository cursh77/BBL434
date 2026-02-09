import sys
import os
import re

def read_fasta(file_path):
    """Reads a FASTA file and returns the sequence string."""
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    seq = []
    for line in lines:
        if not line.startswith('>'):
            seq.append(line.strip())
    return "".join(seq).upper()

def parse_markers(marker_file):
    """
    Parses the messy markers.tab.txt file to map names to sequences.
    
    Since the provided marker file only lists RECOGNITION sites for enzymes 
    but only NAMES for genes (like AmpR), we will:
    1. Parse Restriction Enzymes sites (e.g. BamHI -> GGATCC)
    2. Use a fallback dictionary for common genes (AmpR, LacZ) since 
       their sequences are not explicitly in the text description of the table.
    """
    markers_db = {}
    
    # 1. Fallback/Built-in sequences for Genes (since tab file lacks full gene sequences)
    # These are standard sequences found in pUC vectors
    markers_db['AmpR_gene'] = (
        "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAAC"
        "GCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGA"
        "TCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCC"
        "CGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCAC"
        "AGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCA"
        "ACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTT"
        "GATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAAC"
        "GTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAG"
        "TTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCT"
        "CGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAAC"
        "TATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAA"
    )
    
    markers_db['lacZ_alpha'] = (
        "ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCACTGGC"
        "CGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCA"
        "GCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGC"
    )

    # 2. Parse the text file for Restriction Enzymes
    # Pattern looks for: Name (e.g. BamHI) ... Recognizes SEQUENCE
    if os.path.exists(marker_file):
        with open(marker_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
            
            # Regex to find enzyme names and their recognition sites
            # Example hit: BamHI ... Recognizes GGATCC
            pattern = re.compile(r"\|\s*([A-Za-z0-9]+)\s*\|.*?Recognizes\s*([ATCG]+)", re.DOTALL | re.IGNORECASE)
            matches = pattern.findall(content)
            
            for name, seq in matches:
                # Add to DB, map e.g., "BamHI_site" -> sequence
                # The design file uses "BamHI_site", the table has "BamHI"
                markers_db[f"{name}_site"] = seq.upper()
                markers_db[name] = seq.upper() # Store plain name too just in case

    return markers_db

def calculate_gc_skew(sequence, window_size=100):
    """
    Calculates the Cumulative GC Skew to find the ORI.
    
    Theory:
    In bacteria/plasmids, the Origin of Replication (ORI) is typically located 
    at the global minimum of the cumulative GC skew.
    Skew = (G - C) / (G + C)
    """
    skew = [0]
    curr_skew = 0
    
    # We will traverse the genome and calculate cumulative skew
    # To save time on large files, we step through found chars
    for i, base in enumerate(sequence):
        if base == 'G':
            curr_skew += 1
        elif base == 'C':
            curr_skew -= 1
        skew.append(curr_skew)
    
    return skew

def find_ori(sequence):
    """
    Finds the likely ORI sequence by locating the GC skew minimum.
    Returns the ORI sequence (approx 500bp region).
    """
    print("Analyzing Genome for Origin of Replication (GC Skew Method)...")
    
    # Calculate Skew
    skew_profile = calculate_gc_skew(sequence)
    
    # Find the index of the minimum skew (The theoretical ORI)
    # The minimum value in the cumulative skew diagram points to the Start of Replication
    min_skew_index = skew_profile.index(min(skew_profile))
    
    print(f"ORI location detected around base pair: {min_skew_index}")
    
    # Extract a standard ORI size (e.g., ~600bp centered on the minimum)
    # Handling circular genome boundaries
    genome_len = len(sequence)
    start = min_skew_index - 300
    end = min_skew_index + 300
    
    if start < 0:
        ori_seq = sequence[start:] + sequence[:end]
    elif end > genome_len:
        ori_seq = sequence[start:] + sequence[:end-genome_len]
    else:
        ori_seq = sequence[start:end]
        
    return ori_seq

def main():
    # 1. Handle Arguments
    if len(sys.argv) < 3:
        print("Usage: python plasmid_designer.py <dna_genome.fa> <design_file.txt>")
        sys.exit(1)
        
    genome_file = sys.argv[1]
    design_file = sys.argv[2]
    markers_file = "markers.tab.txt" # Assumed to be in local dir

    # 2. Read Input Data
    print(f"Reading genome from {genome_file}...")
    genome_seq = read_fasta(genome_file)
    
    print(f"Parsing marker database from {markers_file}...")
    markers_db = parse_markers(markers_file)
    
    # 3. Find ORI
    # This fulfills the assignment requirement to "make an ori finder"
    ori_sequence = find_ori(genome_seq)
    
    # 4. Construct Plasmid from Design
    print(f"Reading design from {design_file}...")
    final_plasmid_seq = ""
    
    # We assume the user wants the ORI first, then the designed parts
    # (Or we check if ORI is listed in the design. The prompt's Design file lists 'ori_pMB1')
    
    with open(design_file, 'r') as df:
        design_lines = df.readlines()
    
    components_added = []
    
    for line in design_lines:
        line = line.strip()
        if not line or ',' not in line:
            continue
            
        feature_type, feature_name = [x.strip() for x in line.split(',', 1)]
        
        # Determine sequence to add
        seq_to_add = ""
        
        # Check if it is an ORI request
        if "ori" in feature_type.lower() or "replication" in feature_name.lower():
            seq_to_add = ori_sequence
            components_added.append(f"ORI (Found at {feature_type})")
        
        # Check if it is in our parsed markers DB
        elif feature_type in markers_db:
            seq_to_add = markers_db[feature_type]
            components_added.append(feature_type)
            
        elif feature_name in markers_db:
            seq_to_add = markers_db[feature_name]
            components_added.append(feature_name)
            
        else:
            print(f"Warning: Marker '{feature_name}' (Type: {feature_type}) not found in database. Skipping.")
            continue
            
        final_plasmid_seq += seq_to_add

    # 5. Output
    output_filename = "Output.Fa"
    with open(output_filename, 'w') as out:
        out.write(">Constructed_Plasmid [Design based on " + design_file + "]\n")
        # Format lines of 80 chars for standard FASTA
        for i in range(0, len(final_plasmid_seq), 80):
            out.write(final_plasmid_seq[i:i+80] + "\n")
            
    print("-" * 30)
    print("Plasmid construction complete.")
    print(f"Features added: {', '.join(components_added)}")
    print(f"Total length: {len(final_plasmid_seq)} bp")
    print(f"Output saved to: {output_filename}")

if __name__ == "__main__":
    main()