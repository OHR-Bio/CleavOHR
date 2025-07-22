# Protein Cleavage Analysis Pipeline

A bioinformatics pipeline for analyzing protein cleavage sites using BLAST similarity searches and custom cleavage prediction tools based on MEROPS database (https://www.ebi.ac.uk/merops/index.shtml).

## Prerequisites

- Python 3.x
- Docker (for BLAST analysis)
- CleaveOHR tool (ensure it's in your PATH or current directory)

## Directory Structure

```
project/
├── CleaveOHR.py                    # Main cleavage analysis script
├── db/
│   ├── Substrate_search.txt        # Substrate database for cleavage analysis
│   └── cleavage_db_clean.fasta    # Clean cleavage database for BLAST
├── data/
│   └── WT_with_tag.fasta          # Alternative format (check filename)
└── results/
    ├── results_WT_with_taag_all_species.csv
    └── blast_results_with_headers_WT_with_tag.txt
```

## Usage

### Step 1: Cleavage Site Prediction

Run the CleaveOHR analysis to predict cleavage sites:

```bash
python3 CleaveOHR.py db/Substrate_search.txt --fasta data/WT._with_tagfasta -o results_WT_with_taag_all_species.csv
```

**Parameters:**
- `db/Substrate_search.txt`: Substrate database file
- `--fasta data/WT._with_tagfasta`: Input FASTA file with target sequences
- `-o results_WT_with_taag_all_species.csv`: Output file for results

### Step 2: BLAST Similarity Search

#### Create Header for BLAST Results

```bash
echo -e "Query_ID\tSubject_ID\tPercent_Identity\tAlignment_Length\tMismatches\tGap_Opens\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_value\tBit_Score\tQuery_Sequence\tSubject_Sequence" > blast_results_with_headers_WT_with_tag.txt
```

#### Run BLAST Analysis

```bash
docker run --rm -v $(pwd):/blast/blastdb ncbi/blast:latest blastp \
  -query /blast/blastdb/db/cleavage_db_clean.fasta \
  -subject /blast/blastdb/data/WT_with_tag.fasta \
  -evalue 1000 \
  -word_size 2 \
  -task blastp-short \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
  -max_target_seqs 100000 \
  -qcov_hsp_perc 87.5 >> blast_results_with_headers_WT_with_tag.txt
```

**BLAST Parameters:**
- `-query`: Database of known cleavage sequences
- `-subject`: Target protein sequences
- `-evalue 1000`: Relaxed E-value threshold for comprehensive search
- `-word_size 2`: Small word size for short sequence matching
- `-task blastp-short`: Optimized for short sequences
- `-qcov_hsp_perc 87.5`: Minimum query coverage percentage
- `-max_target_seqs 100000`: Maximum number of target sequences

## Output Files

### CleaveOHR Results
- **File:** `results_WT_with_taag_all_species.csv`
- **Content:** Predicted cleavage sites with scores and positions

### BLAST Results
- **File:** `blast_results_with_headers_WT_with_tag.txt`
- **Content:** Detailed alignment results with the following columns:
  - Query_ID, Subject_ID
  - Percent_Identity, Alignment_Length
  - Mismatches, Gap_Opens
  - Query_Start, Query_End, Subject_Start, Subject_End
  - E_value, Bit_Score
  - Query_Sequence, Subject_Sequence

## File Format Requirements

- **FASTA files:** Standard protein FASTA format
- **Database files:** Ensure proper formatting for respective tools
- **Output:** Tab-separated values with headers

## Notes

- The BLAST search uses relaxed parameters (`-evalue 1000`, `-word_size 2`) to capture short, potentially significant matches
- High query coverage requirement (87.5%) ensures meaningful alignments
- Docker containerization ensures reproducible BLAST results across environments

## Troubleshooting

1. **File path issues:** Ensure all input files exist in specified locations
2. **Docker permissions:** Make sure Docker has access to your working directory
3. **Memory requirements:** Large databases may require additional system resources
4. **File naming:** Check for typos in filenames (note: `WT._with_tagfasta` vs `WT_with_tag.fasta`)

## Dependencies

- Python 3.x
- Docker
- NCBI BLAST+ (via Docker container)
- CleaveOHR tool and its dependencies
