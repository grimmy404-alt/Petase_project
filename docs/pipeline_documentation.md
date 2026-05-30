# PETase Discovery Pipeline - Complete Documentation

## Project Overview

**Goal:** ML-guided discovery of novel PET-degrading enzymes from marine metagenomes

Integrating sequence-based ML + BLAST homology search to identify marine PETase candidates from ocean metagenomes — environments underexplored for plastic-degrading enzymes.

---

## Environment Setup

### System
- OS: Windows 10 + WSL2 (Ubuntu)
- IDE: VS Code with Remote-WSL extension

### Tools & Packages

| Tool | Purpose | Install |
|------|---------|---------|
| Conda | Environment manager | Miniconda |
| Python 3.11 | Core language | `conda create -n petase python=3.11` |
| Biopython | Sequence handling | `pip install biopython` |
| NumPy | Numerical arrays | `pip install numpy` |
| Scikit-learn | ML library | `pip install scikit-learn` |
| BLAST+ | Sequence similarity search | `conda install -c bioconda blast` |
| HMMER | Hidden Markov Models | `conda install -c bioconda hmmer` |
| SeqKit | FASTA manipulation | `conda install -c bioconda seqkit` |
| MAFFT | Multiple sequence alignment | `conda install -c bioconda mafft` |

### Project Folder Structure

```
petase_project/
├── data/
│   ├── known_petases.fasta       # 179 confirmed PETase sequences
│   ├── non_petases.fasta         # 74 non-PETase cutinases/esterases
│   ├── marine_candidates.fasta   # 780 marine candidate sequences
│   ├── dataset.pkl               # Combined labeled dataset
│   ├── features.npy              # Feature matrix (253 × 24)
│   ├── labels.npy                # Labels (1=PETase, 0=non-PETase)
│   └── external/
│       ├── marine_alpha_beta_hydrolases.fasta  # 7,123 sequences
│       ├── marine_hydrolases.fasta             # 13,871 sequences
│       ├── marine_pet_hydrolases.fasta         # 134 sequences
│       └── all_marine.fasta                    # Combined marine sequences
├── scripts/
│   ├── fasta_parser.py           # FASTA parsing utility
│   ├── prepare_dataset.py        # Dataset preparation
│   ├── extract_features.py       # Feature extraction
│   ├── train_model.py            # ML model training
│   └── predict_candidates.py     # Predict marine candidates
├── results/
│   ├── blast_alpha_beta.tsv
│   ├── blast_hydrolases.tsv
│   ├── blast_pet_hydrolases.tsv
│   ├── high_confidence_alpha_beta.tsv   # 2,036 hits
│   ├── high_confidence_hydrolases.tsv  # 5,789 hits
│   ├── high_confidence_pet_hydrolases.tsv # 4,700 hits
│   ├── all_candidate_ids.txt     # 672 unique candidate IDs
│   └── predictions.txt           # ML predictions (Phase 3 output)
├── models/
│   └── petase_classifier.pkl     # Trained Random Forest model
├── notebooks/
├── logs/
├── docs/
└── README.md
```

---

## Phase 1: Training Data Collection

### Goal
Collect known PETase sequences (positive examples) and non-PETase esterases (negative examples) to train ML classifier.

### 1A: Positive Examples (Known PETases)

**Source:** PANDA dataset (Ahituv et al., 2025, Protein Science)
- Paper: "The diversity of PET degrading enzymes: A systematic review of sequence, structure, and function"
- URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC12432417/
- Downloaded Supplementary Table 1 (Excel file)
- Extracted FASTA sequences manually
- Reference sequence: UniProt A0A0K8P6T7 (I. sakaiensis IsPETase)

**Result:** 179 PETase sequences saved as `data/known_petases.fasta`

**Count verification:**
```bash
python scripts/fasta_parser.py
# Output: Total: 179
```

### 1B: Negative Examples (Non-PETases)

**Source:** UniProt
- Search query: `cutinase NOT PET NOT polyester`
- Filter: Reviewed (SwissProt only)
- Download: FASTA format

**Why these negatives?**
- Cutinases are structurally similar to PETases (same alpha/beta hydrolase fold)
- Hard negatives = better ML model
- Model must learn functional differences, not just structural

**Result:** 74 sequences saved as `data/non_petases.fasta`

### FASTA Parser Script

```python
# scripts/fasta_parser.py
def parse_fasta(file):
    with open(file,"r") as f:
        data = f.read().splitlines()

    DNA_array = []
    current_sequence = ""

    for line in data:
        if line.startswith(">"):
            if current_sequence:
                DNA_array.append(current_sequence)
            current_sequence = ""
        elif line.strip():
            current_sequence += line.strip()

    if current_sequence:
        DNA_array.append(current_sequence)
    return DNA_array

DNA_array = parse_fasta("data/known_petases.fasta")
print(f"Total: {len(DNA_array)}")
```

---

## Phase 2: ML Model Training

### Goal
Train a classifier to distinguish PETases from non-PETases based on sequence features.

### Step 2A: Dataset Preparation

```python
# scripts/prepare_dataset.py
import pickle

def parse_fasta(file):
    # [same parser as above]

petases = parse_fasta("data/known_petases.fasta")      # 179
non_petases = parse_fasta("data/non_petases.fasta")    # 74

labels = [1]*len(petases) + [0]*len(non_petases)
all_seqs = petases + non_petases                        # 253 total

with open("data/dataset.pkl", "wb") as f:
    pickle.dump({"sequences": all_seqs, "labels": labels}, f)
```

**Output:** `data/dataset.pkl` — combined labeled dataset

**Why pickle?**
Saves Python objects directly to disk. Faster than re-reading FASTA files every time.

### Step 2B: Feature Extraction

**Problem:** ML models need numbers, not amino acid letters.

**Solution:** Convert each sequence into 24 numerical features.

| Feature | Count | Biological Rationale |
|---------|-------|---------------------|
| Amino acid frequencies | 20 | Protein fingerprint |
| Sequence length | 1 | PETases ~250-310 aa |
| Aromatic content (F,Y,W) | 1 | Direct PET binding |
| Hydrophobic content (A,I,L,M,F,V,P,W,G) | 1 | PET surface binding |
| Charged content (D,E,K,R) | 1 | pH stability |

```python
# scripts/extract_features.py
import pickle
import numpy as np

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

def extract_features(seq):
    features = []
    for aa in AMINO_ACIDS:
        features.append(seq.count(aa) / len(seq))     # 20 features
    features.append(len(seq))                          # length
    aromatic = sum(seq.count(aa) for aa in 'FYW')
    features.append(aromatic / len(seq))               # aromatic
    hydrophobic = sum(seq.count(aa) for aa in 'AILMFVPWG')
    features.append(hydrophobic / len(seq))            # hydrophobic
    charged = sum(seq.count(aa) for aa in 'DEKR')
    features.append(charged / len(seq))                # charged
    return features

with open("data/dataset.pkl", "rb") as f:
    data = pickle.load(f)

X = np.array([extract_features(seq) for seq in data["sequences"]])
y = np.array(data["labels"])

np.save("data/features.npy", X)
np.save("data/labels.npy", y)
```

**Output:**
- `data/features.npy` — matrix shape (253, 24)
- `data/labels.npy` — array shape (253,)

### Step 2C: Model Training

**Algorithm:** Random Forest Classifier

**Why Random Forest?**
- Works well with small datasets (~253 sequences)
- No need for large data like deep learning
- Interpretable (feature importance available)
- Robust to overfitting
- Industry standard for tabular biological data

**Train/test split:** 80/20
- Training: 202 sequences (model learns)
- Testing: 51 sequences (unseen, model evaluated)

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
print(classification_report(y_test, y_pred,
      target_names=["Non-PETase", "PETase"]))

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

---

## Phase 3: Marine Metagenome Search

### Goal
Search marine microbial sequences for novel PETase candidates using BLAST homology search, then score with trained ML model.

### Why Marine Metagenomes?

- Ocean = 70% of Earth's surface, mostly unexplored microbially
- I. sakaiensis cannot survive in 3% NaCl (marine conditions)
- Marine microbes exposed to ocean plastic = evolutionary pressure
- Marine PETases = genuinely novel sequences
- Direct relevance to ocean plastic pollution crisis

### Step 3A: Download Marine Sequences

**Source:** UniProt (metagenome sequences)

Downloaded three datasets filtered for marine origin:

| File | Sequences | Search Query |
|------|-----------|-------------|
| marine_alpha_beta_hydrolases.fasta | 7,123 | alpha/beta hydrolase marine metagenome |
| marine_hydrolases.fasta | 13,871 | hydrolase marine metagenome |
| marine_pet_hydrolases.fasta | 134 | PET hydrolase marine metagenome |

### Step 3B: BLAST Search

**What is BLAST?**
BLAST (Basic Local Alignment Search Tool) finds sequences similar to known PETases in marine databases.

**Why BLAST before ML?**
Two-stage filtering:
1. BLAST = fast similarity filter (millions → thousands)
2. ML = functional prediction on remaining candidates

**Build reference database:**
```bash
makeblastdb -in data/known_petases.fasta \
            -dbtype prot \
            -out data/petase_db
```

**Run BLAST on each marine dataset:**
```bash
# Alpha-beta hydrolases
blastp -query data/external/marine_alpha_beta_hydrolases.fasta \
       -db data/petase_db \
       -out results/blast_alpha_beta.tsv \
       -outfmt 6 -evalue 0.001 -num_threads 4

# Hydrolases
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
```
1. qseqid    = query sequence ID (marine sequence)
2. sseqid    = subject ID (known PETase matched)
3. pident    = % sequence identity
4. length    = alignment length
5. mismatch  = number of mismatches
6. gapopen   = number of gap openings
7-10.        = alignment positions
11. evalue   = statistical significance
12. bitscore = alignment score
```

### Step 3C: Filter High Confidence Hits

**Thresholds:**
- Identity ≥ 35% (distant but related — novel enzyme range)
- E-value ≤ 1e-10 (statistically very significant)

**Why 35% identity?**
- >70% = find same enzymes again, not novel
- <20% = too distant, probably not functional
- 35% = sweet spot for finding novel but related enzymes

```bash
awk '$3 >= 35 && $11 <= 1e-10' results/blast_alpha_beta.tsv \
    > results/high_confidence_alpha_beta.tsv

awk '$3 >= 35 && $11 <= 1e-10' results/blast_hydrolases.tsv \
    > results/high_confidence_hydrolases.tsv

awk '$3 >= 35 && $11 <= 1e-10' results/blast_pet_hydrolases.tsv \
    > results/high_confidence_pet_hydrolases.tsv
```

**Results:**

| Dataset | High Confidence Hits |
|---------|---------------------|
| Alpha-beta hydrolases | 2,036 |
| Hydrolases | 5,789 |
| PET hydrolases | 4,700 |
| **Total** | **12,525** |

### Step 3D: Deduplicate Candidates

```bash
cut -f1 results/high_confidence_alpha_beta.tsv \
    results/high_confidence_hydrolases.tsv \
    results/high_confidence_pet_hydrolases.tsv \
    | sort -u > results/all_candidate_ids.txt

wc -l results/all_candidate_ids.txt
# Result: 672 unique candidates
```

### Step 3E: Extract Candidate Sequences

```bash
# Combine all marine FASTA files
cat data/external/marine_alpha_beta_hydrolases.fasta \
    data/external/marine_hydrolases.fasta \
    data/external/marine_pet_hydrolases.fasta \
    > data/external/all_marine.fasta

# Extract candidate sequences by ID
seqkit grep -f results/all_candidate_ids.txt \
    data/external/all_marine.fasta \
    > data/marine_candidates.fasta

# Verify
grep ">" data/marine_candidates.fasta | wc -l
# Result: 780 sequences
```

**Why 780 sequences from 672 IDs?**
Some IDs match multiple isoforms/variants in the database. Normal behavior.

### Step 3F: ML Prediction (NEXT STEP)

Script: `scripts/predict_candidates.py`

**What it will do:**
1. Load trained classifier
2. Extract same 24 features from 780 candidates
3. Predict probability of being PETase (0-1 score)
4. Rank candidates by confidence
5. Output top candidates for structure prediction

---

## Key Concepts Glossary

| Term | Definition |
|------|-----------|
| FASTA | Text format for biological sequences |
| Metagenome | DNA from entire environmental community |
| BLAST | Tool to find similar sequences in databases |
| E-value | Probability match is due to chance (lower = better) |
| Random Forest | ML algorithm using many decision trees voting |
| Feature extraction | Converting sequences to numbers for ML |
| Precision | Of predicted PETases, how many actually are? |
| Recall | Of all real PETases, how many did model find? |
| Alpha/beta hydrolase | Protein fold family PETases belong to |

---

## Current Status

| Phase | Status |
|-------|--------|
| Phase 1: Data collection | ✅ Complete |
| Phase 2: ML training (96% accuracy) | ✅ Complete |
| Phase 3: Marine search + BLAST | ✅ Complete |
| Phase 3: ML prediction on candidates | ⏳ Next step |
| Phase 4: Structure prediction (AlphaFold) | ⏳ Pending |
| Phase 5: Writing + publication | ⏳ Pending |

