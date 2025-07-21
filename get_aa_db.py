import pandas as pd

def convert_to_fasta(df):
    """
    Convert DataFrame with three-letter amino acid codes to one-letter FASTA format
    
    Args:
        df: DataFrame where first column is sequence ID and remaining columns are amino acids
    
    Returns:
        String in FASTA format
    """
    
    # Mapping from three-letter to one-letter amino acid codes
    aa_mapping = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    }
    
    fasta_sequences = []
    
    for index, row in df.iterrows():
        # First column is sequence ID
        seq_id = row.iloc[0]
        
        # Remaining columns are amino acids
        amino_acids = row.iloc[1:]
        
        # Convert to one-letter codes
        one_letter_sequence = ""
        for aa in amino_acids:
            if pd.isna(aa) or aa == '' or aa not in aa_mapping:
                one_letter_sequence += "-"
            else:
                one_letter_sequence += aa_mapping[aa]
        
        # Create FASTA entry
        fasta_sequences.append(f">{seq_id}")
        fasta_sequences.append(one_letter_sequence)
    
    return "\n".join(fasta_sequences)



df = pd.read_csv('Substrate_search.txt', sep='\t', header=None, encoding='latin-1',
                               on_bad_lines='skip', engine='python')

df_subset = df.iloc[:, [1] + list(range(4, 12))]

df_subset = df_subset.replace("'", "", regex=True)

    
a = convert_to_fasta(df_subset)

outfile = open("cleavage_db.fasta", "w")
outfile.write(a)