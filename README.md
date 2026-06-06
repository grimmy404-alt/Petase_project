# ML-Guided Discovery of Novel PETase-Like Enzymes from Marine Metagenomes

A computational pipeline integrating BLAST homology screening, Random Forest machine learning, and structural bioinformatics to identify putative PET-degrading enzymes from marine metagenomic sequences.

> **Status: COMPLETE** — Pipeline finished

---

## Overview

Polyethylene terephthalate (PET) plastic enters the ocean at ~8 million tonnes/year. Enzymatic degradation offers a path to circular recycling, but known PETase diversity is narrow. This project screens 21,128 marine metagenomic protein sequences to find novel PETase-like candidates using a multi-stage computational approach.

**Key results:**
- Random Forest classifier: **96.08% accuracy** (precision/recall >0.95)
- **580 / 672** marine candidates predicted as PETase-like
- **8 high-confidence candidates** selected for structural validation
- **7 / 8** candidates confirmed Ser–Asp–His catalytic triad by MAFFT alignment
- Top 2 candidates (A0A0F9UIZ8, A0A0F9X315) achieve **RMSD 0.714–0.715 Å** vs IsPETase — near-identical folds
- 3 candidates reclassified from "lipase / alpha-beta hydrolase" → high-confidence PETase-like

---

## Repository Structure

```
petase_project/
├── data/
│   ├── known_petases.fasta          # 179 confirmed PETases (PANDA dataset, positives)
│   ├── non_petases.fasta            # 74 non-PETase cutinases (UniProt, negatives)
│   ├── dataset.pkl                  # Combined labeled dataset (253 sequences)
│   ├── features.npy                 # Feature matrix (253 × 24)
│   ├── labels.npy                   # Labels (1=PETase, 0=non-PETase)
│   ├── marine_candidates.fasta      # 780 marine candidates post-BLAST
│   ├── aligned_candidates.fasta     # MAFFT alignment (8 candidates + IsPETase)
│   ├── alignment_input.fasta        # MAFFT input
│   ├── top_candidates_final2.fasta  # 8 final candidates (UniProt sequences)
│   └── external/
│       ├── marine_alpha_beta_hydrolases.fasta  # 7,123 seqs
│       ├── marine_hydrolases.fasta             # 13,871 seqs
│       ├── marine_pet_hydrolases.fasta         # 134 seqs
│       └── all_marine.fasta                    # Combined (21,128 seqs)
├── models/
│   └── petase_classifier.pkl        # Trained Random Forest model
├── scripts/
│   ├── fasta_parser.py              # FASTA parsing utilities (shared)
│   ├── prepare_dataset.py           # Phase 1: combine FASTAs, label, pickle
│   ├── extract_features.py          # Phase 2: extract 24 features per sequence
│   ├── train_model.py               # Phase 2: train Random Forest, save model
│   ├── predict_candidates.py        # Phase 3: predict on marine candidates
│   ├── analyze_structures.py        # Phase 4: parse pLDDT from PDB/CIF files
│   └── catalytic_triad.py           # Phase 4: check Ser/Asp/His conservation
├── results/
│   ├── blast_*.tsv                  # Raw BLAST output (3 datasets)
│   ├── high_confidence_*.tsv        # Filtered BLAST hits (identity≥35%, e≤1e-10)
│   ├── all_candidate_ids.txt        # 672 unique candidate IDs
│   ├── predictions.txt              # ML scores for all 672 candidates (ranked)
│   ├── top_candidates.txt           # 8 final candidate IDs
│   ├── catalytic_triad_analysis.txt # Catalytic triad conservation results
│   ├── structure_analysis.txt       # pLDDT scores per candidate
│   ├── structures/                  # Predicted PDB/CIF files (ESMFold, AlphaFold2)
│   └── figures/                     # PyMOL figures (superpositions, active sites)
├── docs/
│   ├── pipeline_documentation.md    # Full technical documentation
│   └── final_result_documentation.md
└── archive/
    └── alphafold_raw/               # Raw ColabFold/AlphaFold Server outputs
```

---

## Pipeline

### Phase 1 — Training Data
- **Positives:** 179 PETases from [PANDA dataset](https://doi.org/10.1093/protein/gzae026) (Ahituv et al., 2025)
- **Negatives:** 74 cutinases from UniProt Swiss-Prot (`cutinase NOT PET NOT polyester`, reviewed only)
- Hard negatives chosen deliberately: cutinases share the alpha/beta hydrolase fold with PETases

### Phase 2 — Machine Learning
- **Features:** 24 per sequence (20 amino acid frequencies, length, aromatic/hydrophobic/charged content)
- **Algorithm:** Random Forest, 100 trees, 80/20 stratified train/test split, `random_state=42`
- **Performance:**

| Class | Precision | Recall | F1-Score |
|---|---|---|---|
| Non-PETase | 0.95 | 0.95 | 0.95 |
| PETase | 0.97 | 0.97 | 0.97 |
| **Overall** | **0.9608** | | |

### Phase 3 — Marine Metagenome Screening
```bash
# Build BLAST database
makeblastdb -in data/known_petases.fasta -dbtype prot -out data/petase_db

# Search each marine dataset
blastp -query data/external/marine_alpha_beta_hydrolases.fasta \
       -db data/petase_db -out results/blast_alpha_beta.tsv \
       -outfmt 6 -evalue 0.001 -num_threads 4

# Filter: identity >= 35%, E-value <= 1e-10
awk '$3 >= 35 && $11 <= 1e-10' results/blast_alpha_beta.tsv \
    > results/high_confidence_alpha_beta.tsv

# Deduplicate and extract
cut -f1 results/high_confidence_*.tsv | sort -u > results/all_candidate_ids.txt
seqkit grep -f results/all_candidate_ids.txt data/external/all_marine.fasta \
    > data/marine_candidates.fasta
```
→ 21,128 sequences → 672 unique candidates → 580 ML-predicted PETases → **8 top candidates selected**

### Phase 4 — Structural Validation
- Structure prediction: ESMFold (candidates ≤400 aa), AlphaFold2 Server (candidate 07, 588 aa)
- pLDDT assessment: all candidates mean pLDDT 88.6–94.2 (high confidence)
- Catalytic triad: MAFFT alignment → mapped IsPETase Ser160/Asp206/His237 → 7/8 conserved
- Structural superposition: PyMOL `align` command

---

## Top 8 Candidates

| Candidate | Source | Annotation | ML Prob | Triad | RMSD vs IsPETase |
|---|---|---|---|---|---|
| A0A365H682 | *Actinomadura craniellae* | PET hydrolase | 1.00 | ✓ | — |
| A0A365H6Y8 | *Actinomadura craniellae* | PET hydrolase | 1.00 | ✓ | — |
| A0A0F9UIZ8 | Marine sediment metagenome | PET hydrolase/cutinase-like | 1.00 | ✓ | **0.714 Å** |
| A0A0F9X315 | Marine sediment metagenome | PET hydrolase/cutinase-like | 1.00 | ✓ | **0.715 Å** |
| A0A1E7LM55 | *Streptomyces nanshensis* | Lipase ⚠️ | 0.99 | ✓ | — |
| A0A386WI52 | *Micromonospora tulbaghiae* | Alpha/beta hydrolase ⚠️ | 0.99 | ✓ | — |
| A0A553PS56 | *Tigriopus californicus* | Carboxylic ester hydrolase | 0.99 | ✗ | 7.68 Å (negative control) |
| A0A679PDB4 | *Streptomyces* sp. SM14 | Alpha/beta hydrolase ⚠️ | 0.99 | ✓ | — |

⚠️ = annotation potentially incomplete; pipeline identifies as PETase-like

---

## Requirements

```bash
conda create -n petase python=3.11
conda activate petase
pip install scikit-learn numpy biopython
# External tools: BLAST+, MAFFT, SeqKit, PyMOL
```

---

## Limitations

- All structures computationally predicted (ESMFold / AlphaFold2), not experimentally determined
- RMSD measured against ColabFold-predicted IsPETase, not the experimental crystal structure (PDB: 6EQE)
- ML trained on 253 sequences — small dataset, predominantly terrestrial PETases
- No wet-lab validation — all candidates are *putative* pending biochemical assay
- BLAST 35% identity threshold may miss highly divergent marine PETases (HMMER noted as future work)
