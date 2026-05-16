# Phylogenetic Analysis Scripts

This directory contains scripts for phylogenetic analysis of maize AGC gene family.

## Script Overview

### Main Scripts

| Script | Language | Purpose |
|--------|----------|---------|
| `prepare_phylogeny_data.py` | Python | Extract exemplar sequences from pan-genome data |
| `run_mafft.ps1` | PowerShell | Run MAFFT multiple sequence alignment |
| `run_iqtree.py` | Python | Run trimAI and IQ-TREE tree construction |

### Usage

```bash
# Step 1: Prepare sequence data
python prepare_phylogeny_data.py

# Step 2: MAFFT alignment
powershell run_mafft.ps1

# Step 3: IQ-TREE phylogenetic tree
python run_iqtree.py
```

### Output Files

- `phylogeny/maize.pep_fulllength.fasta` - Maize AGC protein sequences
- `phylogeny/all_sequences_for_tree.fasta` - Combined maize + Arabidopsis sequences
- `phylogeny/all_sequences_aligned.fasta` - MAFFT alignment output
- `phylogeny/all_sequences_aligned_trimmed.fasta` - Trimmed alignment
- `phylogeny/full_length.treefile` - Maximum likelihood tree
- `phylogeny/full_length.iqtree` - IQ-TREE report

### Requirements

- Python 3.x with pandas
- MAFFT (for alignment)
- IQ-TREE (for tree construction)
- trimAI (for alignment trimming)

### Notes

- Scripts assume relative paths from the project root
- IQ-TREE may take several hours to run with MFP model selection
- Intermediate files are saved in the `phylogeny/` directory