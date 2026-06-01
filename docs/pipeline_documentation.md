# PETase Discovery Pipeline — Complete Documentation

## Project Overview

**Goal:** ML-guided computational discovery of novel PET-degrading enzymes from marine metagenomes.

**Approach:** Two-stage filtering pipeline — BLAST homology search to reduce search space, followed by machine learning classification — to identify marine PETase candidates from ocean metagenomic datasets. Candidates are then validated via catalytic triad analysis and structural alignment against the reference enzyme *I. sakaiensis* IsPETase.

**Rationale for marine focus:** Ocean environments cover ~70% of Earth's surface and are heavily polluted with PET plastic, yet remain largely unexplored for plastic-degrading microorganisms. *I. sakaiensis* itself cannot survive in 3% NaCl (marine conditions), meaning marine PETases, if they exist, would represent genuinely novel enzymes adapted to a distinct environment.

---

## Environment

| Item | Detail |
|------|--------|
| OS | Windows 10 + WSL2 (Ubuntu) |
| IDE | VS Code with Remote-WSL extension |
| Python | 3.11 (conda environment: `petase`) |

### Installed Packages

| Tool | Purpose | Install |
|------|---------|---------|
| Biopython | Sequence handling | `pip install biopython` |
| NumPy | Numerical arrays | `pip install numpy` |
| Scikit-learn | ML library | `pip install scikit-learn` |
| BLAST+ | Sequence similarity search | `conda install -c bioconda blast` |
| SeqKit | FASTA manipulation | `conda install -c bioconda seqkit` |
| MAFFT | Multiple sequence alignment | `conda install -c bioconda mafft` |
| PyMOL | Structure visualization and alignment | Separate install |

---

## Project Folder Structure

```
petase_project/
├── README.md
├── archive/
│   └── alphafold_raw/                        # Archived raw AlphaFold/ColabFold outputs
│       ├── IsPETase_reference_58458/         # ColabFold output for IsPETase reference
│       │   ├── *.pdb                         # 5 predicted structures (ranks 1–5)
│       │   ├── *.json                        # Scores and PAE data
│       │   ├── *.png                         # Coverage, PAE, pLDDT plots
│       │   └── IsPETase_reference_58458.a3m  # MSA file
│       └── candidate_07_a0a553ps56/          # AlphaFold Server output for Candidate 07
│           ├── *.cif                         # 5 predicted structures
│           ├── *.json                        # Full data and confidence summaries
│           └── msas/ templates/              # MSA and template files
├── data/
│   ├── known_petases.fasta                   # 179 confirmed PETase sequences (training positives)
│   ├── non_petases.fasta                     # 74 non-PETase cutinases/esterases (training negatives)
│   ├── marine_candidates.fasta               # 780 marine candidate sequences (post-BLAST)
│   ├── dataset.pkl                           # Combined labeled dataset (253 sequences)
│   ├── features.npy                          # Feature matrix (253 × 24)
│   ├── labels.npy                            # Labels (1=PETase, 0=non-PETase)
│   ├── alignment_input.fasta                 # Input for MAFFT alignment
│   ├── aligned_candidates.fasta              # MAFFT-aligned output
│   ├── top_candidates.fasta                  # Initial top candidate sequences
│   ├── top_candidates_unique.fasta           # Deduplicated top candidates
│   ├── top_candidates_final.fasta            # Final filtered candidate set
│   ├── top_candidates_final2.fasta           # Final candidate set (alternate iteration)
│   ├── petase_db.phr / .pin / .psq           # BLAST database files (built from known_petases.fasta)
│   └── external/
│       ├── marine_alpha_beta_hydrolases.fasta  # 7,123 sequences
│       ├── marine_hydrolases.fasta             # 13,871 sequences
│       ├── marine_pet_hydrolases.fasta         # 134 sequences
│       └── all_marine.fasta                    # Combined marine sequences (21,128 total)
├── models/
│   ├── petase_classifier.pkl                 # Trained Random Forest classifier
├── results/
│   ├── blast_alpha_beta.tsv                  # Raw BLAST results — alpha/beta hydrolases
│   ├── blast_hydrolases.tsv                  # Raw BLAST results — hydrolases
│   ├── blast_pet_hydrolases.tsv              # Raw BLAST results — PET hydrolases
│   ├── blast_hits.tsv                        # Combined BLAST hits
│   ├── high_confidence_alpha_beta.tsv        # Filtered hits (≥35% identity, e-value ≤1e-10)
│   ├── high_confidence_hydrolases.tsv        # Filtered hits
│   ├── high_confidence_pet_hydrolases.tsv    # Filtered hits
│   ├── all_candidate_ids.txt                 # 672 unique candidate IDs (deduplicated)
│   ├── predictions.txt                       # ML classifier output on 780 candidates
│   ├── top_candidates.txt                    # Final ranked candidate list
│   ├── catalytic_triad_analysis.txt          # Catalytic residue conservation results
│   ├── structure_analysis.txt                # ESMFold/AlphaFold pLDDT statistics
│   ├── petase_overlay.png                    # PyMOL whole-structure overlay figure
│   ├── candidate03_overlay.png               # PyMOL overlay — Candidate 03
│   ├── catalytic_triad.png                   # PyMOL catalytic triad visualization
│   └── structures/
│       ├── IsPETase_reference.pdb            # Reference structure
│       ├── candidate_01_A0A365H682.pdb       # ESMFold predicted structure
│       ├── candidate_02_A0A365H6Y8.pdb
│       ├── candidate_03_A0A0F9UIZ8.pdb
│       ├── candidate_04_A0A0F9X315.pdb
│       ├── candidate_05_A0A1E7LM55.pdb
│       ├── candidate_06_A0A386WI52.pdb
│       ├── candidate_07_A0A553PS56.cif       # AlphaFold Server output (CIF format)
│       └── candidate_08_A0A679PDB4.pdb
├── scripts/
│   ├── fasta_parser.py                       # FASTA parsing utility
│   ├── prepare_dataset.py                    # Dataset preparation (Phase 1)
│   ├── extract_features.py                   # Feature extraction (Phase 2)
│   ├── train_model.py                        # ML model training (Phase 2)
│   ├── predict_candidates.py                 # Predict marine candidates (Phase 3)
│   ├── catalytic_triad.py                    # Catalytic triad analysis (Phase 4)
│   └── analyze_structures.py                 # Structure confidence analysis (Phase 4)
├── docs/
│   └── pipeline_documentation.md            # This file
├── petase_analysis.pse                       # PyMOL session — full candidate analysis
├── petase_comparison.pse                     # PyMOL session — IsPETase comparison
├── logs/
└── notebooks/
```

---

## Phase 1: Training Data Collection

**Goal:** Collect known PETase sequences (positive examples) and non-PETase esterases (negative examples) to train the ML classifier.

### 1A — Positive Examples (Known PETases)

**Source:** PANDA dataset (Ahituv et al., 2025, *Protein Science*)

> Ahituv et al. "The diversity of PET degrading enzymes: A systematic review of sequence, structure, and function." *Protein Science* (2025). PMC12432417.

- Downloaded Supplementary Table 1 (Excel file) from the paper
- Extracted FASTA sequences manually
- Reference sequence: UniProt **A0A0K8P6T7** (*I. sakaiensis* IsPETase, 290 aa)

**Result:** 179 PETase sequences → `data/known_petases.fasta`

Verified with:
```python
# scripts/fasta_parser.py
python scripts/fasta_parser.py
# Output: Total: 179
```

### 1B — Negative Examples (Non-PETases)

**Source:** UniProt (SwissProt — reviewed entries only)

**Search query:** `cutinase NOT PET NOT polyester`

**Why cutinases as negatives?**
Cutinases share the same alpha/beta hydrolase fold as PETases. Using them as negatives forces the model to learn *functional* differences rather than just fold-level structural patterns — this produces a much more discriminating classifier.

**Result:** 74 sequences → `data/non_petases.fasta`

### FASTA Parser

```python
# scripts/fasta_parser.py
def parse_fasta(file):
    with open(file, "r") as f:
        data = f.read().splitlines()

    sequences = []
    current_sequence = ""

    for line in data:
        if line.startswith(">"):
            if current_sequence:
                sequences.append(current_sequence)
            current_sequence = ""
        elif line.strip():
            current_sequence += line.strip()

    if current_sequence:
        sequences.append(current_sequence)
    return sequences
```

---

## Phase 2: ML Model Training

**Goal:** Train a binary classifier to distinguish PETases from non-PETases based on sequence-derived numerical features.

### 2A — Dataset Preparation

```python
# scripts/prepare_dataset.py
import pickle

petases = parse_fasta("data/known_petases.fasta")      # 179 sequences
non_petases = parse_fasta("data/non_petases.fasta")    # 74 sequences

labels = [1] * len(petases) + [0] * len(non_petases)
all_seqs = petases + non_petases                        # 253 total

with open("data/dataset.pkl", "wb") as f:
    pickle.dump({"sequences": all_seqs, "labels": labels}, f)
```

Output: `data/dataset.pkl` — combined labeled dataset (253 sequences)

> **Why pickle?** Saves Python objects directly to disk. Faster than re-parsing FASTA files on every run.

### 2B — Feature Extraction

ML models require numerical input. Each sequence is converted into 24 numerical features:

| Feature | Count | Biological Rationale |
|---------|-------|----------------------|
| Amino acid frequencies (20 aa) | 20 | Sequence composition fingerprint |
| Sequence length | 1 | PETases are typically 250–310 aa |
| Aromatic content (F, Y, W) | 1 | Direct PET ring binding |
| Hydrophobic content (A, I, L, M, F, V, P, W, G) | 1 | PET surface binding |
| Charged content (D, E, K, R) | 1 | pH stability |
| **Total** | **24** | |

```python
# scripts/extract_features.py
import pickle
import numpy as np

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

def extract_features(seq):
    features = []
    for aa in AMINO_ACIDS:
        features.append(seq.count(aa) / len(seq))               # 20 compositional features
    features.append(len(seq))                                    # length
    aromatic = sum(seq.count(aa) for aa in 'FYW')
    features.append(aromatic / len(seq))                         # aromatic content
    hydrophobic = sum(seq.count(aa) for aa in 'AILMFVPWG')
    features.append(hydrophobic / len(seq))                      # hydrophobic content
    charged = sum(seq.count(aa) for aa in 'DEKR')
    features.append(charged / len(seq))                          # charged content
    return features

with open("data/dataset.pkl", "rb") as f:
    data = pickle.load(f)

X = np.array([extract_features(seq) for seq in data["sequences"]])
y = np.array(data["labels"])

np.save("data/features.npy", X)   # shape: (253, 24)
np.save("data/labels.npy", y)     # shape: (253,)
```

### 2C — Model Training

**Algorithm:** Random Forest Classifier

**Why Random Forest?**
- Performs well on small datasets (~253 sequences)
- Does not require large training sets unlike deep learning
- Interpretable via feature importance scores
- Robust to overfitting
- Standard approach for tabular biological classification tasks

**Train/test split:** 80/20 (202 training, 51 testing)

```python
# scripts/train_model.py
import numpy as np
import pickle
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report

X = np.load("data/features.npy")
y = np.load("data/labels.npy")

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

y_pred = model.predict(X_test)
print(f"Accuracy: {accuracy_score(y_test, y_pred):.2%}")
print(classification_report(y_test, y_pred, target_names=["Non-PETase", "PETase"]))

with open("models/petase_classifier.pkl", "wb") as f:
    pickle.dump(model, f)
```

**Results:**

| Metric | Non-PETase | PETase |
|--------|-----------|--------|
| Precision | 0.95 | 0.97 |
| Recall | 0.95 | 0.97 |
| F1-score | 0.95 | 0.97 |
| **Overall Accuracy** | **96.08%** | |

Output: `models/petase_classifier.pkl`

An HMMER profile was also built for the PETase family: `models/petase.hmm`

---

## Phase 3: Marine Metagenome Search and Candidate Prediction

**Goal:** Search marine metagenomic sequences for PETase-like candidates using BLAST homology, then classify with the trained ML model.

### 3A — Marine Dataset Download

**Source:** UniProt (metagenome-derived sequences, filtered for marine origin)

| File | Sequences | Search Query |
|------|-----------|-------------|
| `marine_alpha_beta_hydrolases.fasta` | 7,123 | alpha/beta hydrolase marine metagenome |
| `marine_hydrolases.fasta` | 13,871 | hydrolase marine metagenome |
| `marine_pet_hydrolases.fasta` | 134 | PET hydrolase marine metagenome |
| **`all_marine.fasta`** (combined) | **21,128** | |

### 3B — BLAST Database Build

```bash
makeblastdb -in data/known_petases.fasta \
            -dbtype prot \
            -out data/petase_db
```

Generates: `petase_db.phr`, `petase_db.pin`, `petase_db.psq`

### 3C — BLAST Search

```bash
# Alpha-beta hydrolases
blastp -query data/external/marine_alpha_beta_hydrolases.fasta \
       -db data/petase_db \
       -out results/blast_alpha_beta.tsv \
       -outfmt 6 -evalue 0.001 -num_threads 4

# General hydrolases
blastp -query data/external/marine_hydrolases.fasta \
       -db data/petase_db \
       -out results/blast_hydrolases.tsv \
       -outfmt 6 -evalue 0.001 -num_threads 4

# PET hydrolases
blastp -query data/external/marine_pet_hydrolases.fasta \
       -db data/petase_db \
       -out results/blast_pet_hydrolases.tsv \
       -outfmt 6 -evalue 0.001 -num_threads 4
```

**BLAST output format 6 columns:**

| Column | Field | Description |
|--------|-------|-------------|
| 1 | qseqid | Query sequence ID (marine sequence) |
| 2 | sseqid | Subject ID (known PETase matched) |
| 3 | pident | % sequence identity |
| 4 | length | Alignment length |
| 5 | mismatch | Number of mismatches |
| 6 | gapopen | Number of gap openings |
| 7–10 | — | Alignment positions |
| 11 | evalue | Statistical significance (lower = better) |
| 12 | bitscore | Alignment score |

### 3D — Filter High-Confidence Hits

**Thresholds:** Identity ≥ 35%, E-value ≤ 1e-10

> **Why 35% identity?**
> - \>70% = rediscovering known enzymes, not novel sequences
> - <20% = too distant, likely not functionally related
> - 35% = "twilight zone" threshold — captures genuinely novel but functionally related enzymes

```bash
awk '$3 >= 35 && $11 <= 1e-10' results/blast_alpha_beta.tsv \
    > results/high_confidence_alpha_beta.tsv

awk '$3 >= 35 && $11 <= 1e-10' results/blast_hydrolases.tsv \
    > results/high_confidence_hydrolases.tsv

awk '$3 >= 35 && $11 <= 1e-10' results/blast_pet_hydrolases.tsv \
    > results/high_confidence_pet_hydrolases.tsv
```

| Dataset | High-Confidence Hits |
|---------|---------------------|
| Alpha-beta hydrolases | 2,036 |
| Hydrolases | 5,789 |
| PET hydrolases | 4,700 |
| **Total** | **12,525** |

### 3E — Deduplication and Sequence Extraction

```bash
# Extract unique candidate IDs
cut -f1 results/high_confidence_alpha_beta.tsv \
    results/high_confidence_hydrolases.tsv \
    results/high_confidence_pet_hydrolases.tsv \
    | sort -u > results/all_candidate_ids.txt

wc -l results/all_candidate_ids.txt
# Result: 672 unique candidate IDs

# Combine marine FASTA files
cat data/external/marine_alpha_beta_hydrolases.fasta \
    data/external/marine_hydrolases.fasta \
    data/external/marine_pet_hydrolases.fasta \
    > data/external/all_marine.fasta

# Extract candidate sequences by ID
seqkit grep -f results/all_candidate_ids.txt \
    data/external/all_marine.fasta \
    > data/marine_candidates.fasta

grep ">" data/marine_candidates.fasta | wc -l
# Result: 780 sequences
```

> **Why 780 sequences from 672 IDs?** Some UniProt IDs map to multiple isoforms or variants in the combined database. This is expected behavior.

### 3F — ML Prediction on Candidates

```python
# scripts/predict_candidates.py
```

The trained classifier extracted the same 24 features from all 780 candidates and predicted PETase probability (0–1 score). Sequences were ranked by confidence score and filtered to yield the final high-confidence candidates.

**Output:** `results/predictions.txt`, `data/top_candidates.fasta`

**Result:** 8 high-confidence marine PETase candidates selected for structural validation.

---

## Phase 4: Structural Validation

**Goal:** Validate the 8 ML-selected candidates through catalytic residue conservation analysis and 3D structure prediction and alignment.

### 4A — Catalytic Triad Analysis

PETase activity depends on a conserved catalytic triad: **Ser–Asp–His**.

```python
# scripts/catalytic_triad.py
```

**Results:** (`results/catalytic_triad_analysis.txt`)

| Candidate | UniProt ID | Triad Residues | Retained? |
|-----------|-----------|----------------|-----------|
| Candidate 01 | A0A365H682 | S–D–H | ✅ |
| Candidate 02 | A0A365H6Y8 | S–D–H | ✅ |
| Candidate 03 | A0A0F9UIZ8 | S–D–H | ✅ |
| Candidate 04 | A0A0F9X315 | S–D–H | ✅ |
| Candidate 05 | A0A1E7LM55 | S–D–H | ✅ |
| Candidate 06 | A0A386WI52 | S–D–H | ✅ |
| Candidate 07 | A0A553PS56 | S–N–E | ❌ |
| Candidate 08 | A0A679PDB4 | S–D–H | ✅ |

**7 of 8 candidates retained the canonical Ser–Asp–His triad.** Candidate 07 (A0A553PS56) showed a non-canonical S–N–E substitution, suggesting it may not be a functional PETase despite ML classification.

### 4B — Structure Prediction

Structures were predicted using two complementary tools:

**ESMFold** (via API/web): Candidates 01–06 and 08 → `.pdb` files
> ESMFold uses a language-model-based folding approach; fast and suitable for single-sequence prediction without MSA.

**AlphaFold Server** (web): Candidate 07 (A0A553PS56) → `.cif` files (archived in `archive/alphafold_raw/`)
> AlphaFold was used for Candidate 07 due to its lower confidence under ESMFold; full MSA-based folding was applied for a more reliable prediction.

**Reference:** IsPETase was folded using **ColabFold** (archived in `archive/alphafold_raw/IsPETase_reference_58458/`) to provide a consistent computational reference structure.

**Confidence scores** (mean pLDDT, from `results/structure_analysis.txt`):

| Candidate | UniProt ID | Mean pLDDT |
|-----------|-----------|-----------|
| Candidate 01 | A0A365H682 | 93.92 |
| Candidate 02 | A0A365H6Y8 | 93.58 |
| Candidate 03 | A0A0F9UIZ8 | 92.24 |
| Candidate 04 | A0A0F9X315 | 92.91 |
| Candidate 05 | A0A1E7LM55 | 94.21 |
| Candidate 06 | A0A386WI52 | 93.84 |
| Candidate 07 | A0A553PS56 | 88.60 |
| Candidate 08 | A0A679PDB4 | 91.92 |

> pLDDT >90 = high confidence. All candidates except Candidate 07 exceeded this threshold.

```python
# scripts/analyze_structures.py
```

### 4C — Structural Alignment (PyMOL)

All candidate structures were aligned against the IsPETase reference structure using PyMOL's `align` command.

**Sessions saved as:** `petase_analysis.pse` and `petase_comparison.pse`

**Key results:**

| Candidate | UniProt ID | RMSD vs IsPETase (Å) |
|-----------|-----------|----------------------|
| Candidate 03 | A0A0F9UIZ8 | **0.714** |
| Candidate 04 | A0A0F9X315 | **0.715** |

> RMSD < 1 Å indicates near-identical global fold. Values of 0.714–0.715 Å represent exceptional structural conservation relative to IsPETase.

**Figures generated** (`results/`):
- `petase_overlay.png` — whole-structure overlay of all candidates
- `candidate03_overlay.png` — detailed overlay of Candidate 03 vs IsPETase
- `catalytic_triad.png` — catalytic triad residue visualization

---

## Final Results Summary

**Pipeline:** 21,128 marine sequences → BLAST filter → 780 candidates → ML classification → 8 high-confidence candidates → catalytic triad analysis → structural prediction → structural alignment

**Key numbers:**

| Metric | Value |
|--------|-------|
| Marine sequences screened | 21,128 |
| BLAST high-confidence hits | 12,525 |
| Unique candidate IDs | 672 |
| Sequences extracted | 780 |
| ML-selected candidates | 8 |
| Catalytic triad retained (S–D–H) | 7 / 8 |
| Mean pLDDT > 90 | 7 / 8 |
| Top candidates by RMSD | A0A0F9UIZ8 (0.714 Å), A0A0F9X315 (0.715 Å) |

**Main finding:**

> Among eight machine-learning-selected marine candidates, **A0A0F9UIZ8** and **A0A0F9X315** emerged as the most promising PETase-like enzymes. Both retain the conserved Ser–Asp–His catalytic triad and exhibit near-identical three-dimensional structures relative to IsPETase (RMSD ≈ 0.71 Å), indicating strong conservation of the structural framework required for PET hydrolysis.

---

## Limitations

1. No experimental wet-lab validation performed — enzymatic activity on PET substrate remains to be confirmed.
2. Structure predictions are computational; crystallographic or cryo-EM structures would provide higher confidence.
3. ML predictions are dependent on training data quality and size (253 sequences total).
4. Structural alignment was performed using predicted (not experimental) structures for both reference and candidates.
5. Sequence identity between top candidates and IsPETase remains in the twilight zone (~35–40%); distant homology does not guarantee equivalent function.

---

## Future Work

1. Molecular docking with PET oligomers (MHET, BHET) to assess active site compatibility
2. Molecular dynamics simulations to evaluate thermostability and flexibility
3. Recombinant expression in *E. coli* and purification
4. Experimental PET degradation assays (turbidity, HPLC product detection)
5. Protein engineering of top candidates to improve thermostability (mimicking beta 6-7 loop modifications in engineered IsPETase variants)

---

## Key Concepts Glossary

| Term | Definition |
|------|-----------|
| FASTA | Standard text format for biological sequences |
| Metagenome | DNA extracted and sequenced from an entire environmental community without culturing |
| BLAST | Basic Local Alignment Search Tool — finds sequences similar to a query in a database |
| E-value | Probability that a BLAST match is due to random chance; lower values = more significant |
| pLDDT | Per-residue confidence score from AlphaFold/ESMFold (0–100; >90 = high confidence) |
| RMSD | Root-Mean-Square Deviation — measures structural difference between two aligned proteins (Å) |
| Random Forest | Ensemble ML algorithm that aggregates predictions from many decision trees |
| Catalytic triad | Three residues (Ser–Asp–His) required for serine hydrolase enzymatic activity |
| Alpha/beta hydrolase | Protein fold family that PETases, cutinases, and related esterases belong to |
| Precision | Of all sequences predicted as PETase, what fraction actually are |
| Recall | Of all true PETases in the test set, what fraction were correctly predicted |
| ColabFold | Open-source implementation of AlphaFold2 optimized for speed using MMseqs2 MSA |
| ESMFold | Meta's language-model-based protein structure prediction tool (no MSA required) |

---
