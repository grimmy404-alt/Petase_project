# PETase Discovery Pipeline — Complete Documentation

## Project Overview

**Goal:** ML-guided computational discovery of novel PET-degrading enzymes from marine metagenomes.

**Approach:** Two-stage filtering pipeline — BLAST homology search to reduce search space, followed by machine learning classification — to identify marine PETase candidates from ocean metagenomic datasets. Candidates are then validated via catalytic triad analysis and structural alignment against the reference enzyme *I. sakaiensis* IsPETase.

**Rationale for marine focus:** Ocean environments cover ~70% of Earth's surface and are heavily polluted with PET plastic, yet remain largely unexplored for plastic-degrading microorganisms. *I. sakaiensis* itself cannot survive in 3% NaCl (marine conditions), meaning marine PETases, if they exist, would represent genuinely novel enzymes adapted to a distinct environment.

---

## Environment

| Item | Detail |
|---|---|
| OS | Windows 10 + WSL2 (Ubuntu) |
| IDE | VS Code with Remote-WSL extension |
| Python | 3.11 (conda environment: `petase`) |

### Installed Packages

| Tool | Purpose | Install |
|---|---|---|
| NumPy | Numerical arrays | `pip install numpy` |
| scikit-learn | ML library | `pip install scikit-learn` |
| Biopython | Sequence handling | `pip install biopython` |
| BLAST+ | Sequence similarity search | `conda install -c bioconda blast` |
| SeqKit | FASTA manipulation | `conda install -c bioconda seqkit` |
| MAFFT | Multiple sequence alignment | `conda install -c bioconda mafft` |
| PyMOL | Structure visualization and alignment | Separate install |

---

## Phase 1: Training Data Collection

**Goal:** Collect known PETase sequences (positive examples) and non-PETase esterases (negative examples) to train the ML classifier.

### 1A — Positive Examples (Known PETases)

**Source:** PANDA dataset (Ahituv et al., 2025, *Protein Science*)

> Ahituv et al. "The diversity of PET degrading enzymes: A systematic review of sequence, structure, and function." *Protein Science* (2025). PMC12432417.

- Downloaded Supplementary Table 1 (Excel file) from the paper
- Extracted FASTA sequences manually
- Reference sequence: UniProt A0A0K8P6T7 (*I. sakaiensis* IsPETase, 290 aa)

**Result:** 179 PETase sequences → `data/known_petases.fasta`

### 1B — Negative Examples (Non-PETases)

**Source:** UniProt Swiss-Prot (reviewed entries only)

Search query: `cutinase NOT PET NOT polyester`

**Why cutinases as negatives?** Cutinases share the same alpha/beta hydrolase fold as PETases. Using them as negatives forces the model to learn functional differences rather than just fold-level structural patterns — this produces a more discriminating classifier.

**Result:** 74 sequences → `data/non_petases.fasta`

### 1C — FASTA Parsing

All FASTA files throughout the pipeline are parsed with a consistent approach: sequential line-by-line reading, header detection via `>` prefix, sequence accumulation across continuation lines.

```python
# Pattern used across scripts
def parse_fasta(file):
    sequences = {}
    current_header = ""
    current_seq = ""
    with open(file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_seq:
                    sequences[current_header] = current_seq
                current_header = line[1:]
                current_seq = ""
            elif line:
                current_seq += line
    if current_seq:
        sequences[current_header] = current_seq
    return sequences
```

---

## Phase 2: ML Model Training

**Goal:** Train a binary classifier to distinguish PETases from non-PETases based on sequence-derived numerical features.

### 2A — Dataset Preparation (`prepare_dataset.py`)

```python
import pickle

petases = parse_fasta("data/known_petases.fasta")      # 179 sequences
non_petases = parse_fasta("data/non_petases.fasta")    # 74 sequences

labels = [1] * len(petases) + [0] * len(non_petases)
all_seqs = petases + non_petases                        # 253 total

with open("data/dataset.pkl", "wb") as f:
    pickle.dump({"sequences": all_seqs, "labels": labels}, f)
```

**Output:** `data/dataset.pkl` — combined labeled dataset (253 sequences)

**Why pickle?** Saves Python objects directly to disk. Faster than re-parsing FASTA files on every run.

### 2B — Feature Extraction (`extract_features.py`)

ML models require numerical input. Each sequence is converted into 24 numerical features:

| Feature | Count | Biological Rationale |
|---|---|---|
| Amino acid frequencies (20 aa) | 20 | Sequence composition fingerprint |
| Sequence length | 1 | PETases are typically 250–310 aa |
| Aromatic content (F, Y, W) | 1 | Direct PET aromatic ring binding |
| Hydrophobic content (A, I, L, M, F, V, P, W, G) | 1 | PET surface binding |
| Charged content (D, E, K, R) | 1 | pH stability and electrostatics |
| **Total** | **24** | |

```python
import numpy as np

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

def extract_features(seq):
    features = []
    for aa in AMINO_ACIDS:
        features.append(seq.count(aa) / len(seq))
    features.append(len(seq))
    aromatic = sum(seq.count(aa) for aa in 'FYW')
    features.append(aromatic / len(seq))
    hydrophobic = sum(seq.count(aa) for aa in 'AILMFVPWG')
    features.append(hydrophobic / len(seq))
    charged = sum(seq.count(aa) for aa in 'DEKR')
    features.append(charged / len(seq))
    return features

X = np.array([extract_features(seq) for seq in data["sequences"]])
y = np.array(data["labels"])

np.save("data/features.npy", X)   # shape: (253, 24)
np.save("data/labels.npy", y)     # shape: (253,)
```

### 2C — Model Training (`train_model.py`)

**Algorithm:** Random Forest Classifier

**Why Random Forest?**
- Performs well on small datasets (~253 sequences)
- Does not require large training sets unlike deep learning
- Robust to overfitting via ensemble averaging
- Standard approach for tabular biological classification

**Train/test split:** 80/20 stratified (202 training, 51 testing)

```python
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
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
|---|---|---|
| Precision | 0.95 | 0.97 |
| Recall | 0.95 | 0.97 |
| F1-score | 0.95 | 0.97 |
| **Overall Accuracy** | **96.08%** | |

**Output:** `models/petase_classifier.pkl`

---

## Phase 3: Marine Metagenome Search and Candidate Prediction

**Goal:** Search marine metagenomic sequences for PETase-like candidates using BLAST homology, then classify with the trained ML model.

### 3A — Marine Dataset Download

**Source:** UniProt (metagenome-derived sequences, filtered for marine origin)

| File | Sequences | Search Query |
|---|---|---|
| `marine_alpha_beta_hydrolases.fasta` | 7,123 | alpha/beta hydrolase marine metagenome |
| `marine_hydrolases.fasta` | 13,871 | hydrolase marine metagenome |
| `marine_pet_hydrolases.fasta` | 134 | PET hydrolase marine metagenome |
| `all_marine.fasta` (combined) | 21,128 | |

### 3B — BLAST Database Build

```bash
makeblastdb -in data/known_petases.fasta \
            -dbtype prot \
            -out data/petase_db
```

Generates: `petase_db.phr`, `petase_db.pin`, `petase_db.psq`

### 3C — BLAST Search

```bash
blastp -query data/external/marine_alpha_beta_hydrolases.fasta \
       -db data/petase_db \
       -out results/blast_alpha_beta.tsv \
       -outfmt 6 -evalue 0.001 -num_threads 4

blastp -query data/external/marine_hydrolases.fasta \
       -db data/petase_db \
       -out results/blast_hydrolases.tsv \
       -outfmt 6 -evalue 0.001 -num_threads 4

blastp -query data/external/marine_pet_hydrolases.fasta \
       -db data/petase_db \
       -out results/blast_pet_hydrolases.tsv \
       -outfmt 6 -evalue 0.001 -num_threads 4
```

**BLAST outfmt 6 columns:**

| Column | Field | Description |
|---|---|---|
| 1 | qseqid | Query sequence ID (marine sequence) — **always use col 1, not col 2** |
| 2 | sseqid | Subject ID (known PETase matched) |
| 3 | pident | % sequence identity |
| 11 | evalue | Statistical significance |
| 12 | bitscore | Alignment score |

> **Important:** Column 1 (`qseqid`) contains the marine sequence accession. Column 2 is the reference PETase it matched against. Always extract col 1 for candidate IDs.

### 3D — Filter High-Confidence Hits

**Thresholds:** Identity ≥ 35%, E-value ≤ 1e-10

**Why 35% identity?** Above the sequence twilight zone (20–35%). Captures genuinely novel but functionally related enzymes without rediscovering already-known ones.

```bash
awk '$3 >= 35 && $11 <= 1e-10' results/blast_alpha_beta.tsv \
    > results/high_confidence_alpha_beta.tsv

awk '$3 >= 35 && $11 <= 1e-10' results/blast_hydrolases.tsv \
    > results/high_confidence_hydrolases.tsv

awk '$3 >= 35 && $11 <= 1e-10' results/blast_pet_hydrolases.tsv \
    > results/high_confidence_pet_hydrolases.tsv
```

| Dataset | High-Confidence Hits |
|---|---|
| Alpha-beta hydrolases | 2,036 |
| Hydrolases | 5,789 |
| PET hydrolases | 4,700 |
| **Total** | **12,525** |

### 3E — Deduplication and Sequence Extraction

```bash
cut -f1 results/high_confidence_alpha_beta.tsv \
    results/high_confidence_hydrolases.tsv \
    results/high_confidence_pet_hydrolases.tsv \
    | sort -u > results/all_candidate_ids.txt

# Result: 672 unique candidate IDs

seqkit grep -r -p -f results/all_candidate_ids.txt \
    data/external/all_marine.fasta \
    > data/marine_candidates.fasta

# Result: 780 sequences
```

> **Note:** `-r -p` regex flags are required for SeqKit to correctly match against full UniProt FASTA headers. Without them, pattern matching fails silently.

**Why 780 from 672 IDs?** Some UniProt IDs map to multiple isoforms or fragment records in the combined database. Expected behavior.

### 3F — ML Prediction on Candidates (`predict_candidates.py`)

The trained classifier extracts identical 24 features from all 780 candidates and returns `predict_proba` scores (0–1). Candidates are ranked by PETase probability. The top 8 (probability 1.00 or 0.99) were selected for structural validation.

**Output:** `results/predictions.txt`, `data/top_candidates_final2.fasta`

---

## Phase 4: Structural Validation

**Goal:** Validate the 8 ML-selected candidates through catalytic residue conservation analysis and 3D structure prediction and alignment.

### 4A — Catalytic Triad Analysis (`catalytic_triad.py`)

PETase activity depends on a conserved Ser–Asp–His catalytic triad. This script identifies whether each candidate retains these residues at positions equivalent to IsPETase Ser160, Asp206, His237.

**Input:** `data/aligned_candidates.fasta` — all 8 candidates + IsPETase reference, aligned by MAFFT

```bash
mafft --auto data/alignment_input.fasta > data/aligned_candidates.fasta
```

**How it works:**

1. Load the MAFFT alignment (all candidates + IsPETase reference in one file)
2. Walk through the aligned IsPETase sequence, counting only non-gap characters until reaching residue positions 160, 206, and 237 — recording the alignment column index for each
3. For every other sequence in the alignment, extract the residue at those same column indices
4. A candidate is classified as triad-positive if it has S, D, H at those three positions

```python
def find_triad_positions(ref_seq_aligned):
    """Map IsPETase residue numbers to alignment column indices."""
    positions = {}
    residue_count = 0
    for align_pos, aa in enumerate(ref_seq_aligned):
        if aa != "-":
            residue_count += 1
            if residue_count == 160:
                positions["Ser160"] = align_pos
            elif residue_count == 206:
                positions["Asp206"] = align_pos
            elif residue_count == 237:
                positions["His237"] = align_pos
    return positions
```

**Why this approach?** Residue numbers differ across candidates due to insertions and deletions. Mapping through the alignment columns gives positionally equivalent residues regardless of individual sequence numbering.

**Output:** `results/catalytic_triad_analysis.txt`

**Results:**

| Candidate | UniProt ID | Ser | Asp | His | Triad? |
|---|---|---|---|---|---|
| IsPETase (ref) | A0A0K8P6T7 | S160 | D206 | H237 | ✓ |
| Candidate 01 | A0A365H682 | S158 | D204 | H236 | ✓ |
| Candidate 02 | A0A365H6Y8 | S159 | D205 | H237 | ✓ |
| Candidate 03 | A0A0F9UIZ8 | S174 | D218 | H250 | ✓ |
| Candidate 04 | A0A0F9X315 | S170 | D214 | H246 | ✓ |
| Candidate 05 | A0A1E7LM55 | S179 | D225 | H257 | ✓ |
| Candidate 06 | A0A386WI52 | S166 | D212 | H244 | ✓ |
| Candidate 07 | A0A553PS56 | S206 | **N244** | **E288** | ✗ |
| Candidate 08 | A0A679PDB4 | S156 | D202 | H234 | ✓ |

Candidate 07 (A0A553PS56, *Tigriopus californicus*) shows Asp→Asn and His→Glu substitutions — non-canonical, inconsistent with classical serine hydrolase mechanism. Treated as negative control.

### 4B — Structure Prediction

Structures were predicted using two complementary tools:

**ESMFold** (Candidates 01–06, 08 — all ≤400 aa) → `.pdb` files
- Single-sequence prediction, no MSA required
- Fast; suitable for short globular proteins
- pLDDT output in B-factor field, **scale 0–1** (normalized to 0–100 by `analyze_structures.py`)

**AlphaFold Server** (Candidate 07, 588 aa; IsPETase reference) → `.cif` files
- Full MSA-based prediction
- Used for Candidate 07 due to length exceeding ESMFold limit
- pLDDT output in mmCIF column 15, **scale 0–100**
- Raw outputs archived in `archive/alphafold_raw/`

> **Key difference:** ESMFold outputs pLDDT on a 0–1 scale; AlphaFold outputs 0–100. `analyze_structures.py` auto-detects this and normalizes before reporting.

### 4C — Structure Quality Assessment (`analyze_structures.py`)

This script automatically processes all `.pdb` and `.cif` files in `results/structures/`, extracts per-residue pLDDT scores, detects the scale, normalizes if needed, and reports summary statistics.

**How it works:**

**For PDB files (ESMFold):** pLDDT is stored in the B-factor column (characters 60–66) of `ATOM` records. Only Cα atoms are used (one score per residue).

```python
def parse_pdb_plddt(filepath):
    plddt_scores = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                if atom_name == "CA":
                    plddt = float(line[60:66].strip())
                    plddt_scores.append(plddt)
    return plddt_scores
```

**For CIF files (AlphaFold):** pLDDT is in column index 14 (0-indexed, i.e. column 15) of `ATOM` records. Column position was verified manually against known coordinate values before use.

```python
def parse_cif_plddt(filepath):
    plddt_scores = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("ATOM"):
                parts = line.split()
                if parts[3] == "CA":
                    plddt = float(parts[14])
                    if 0 <= plddt <= 100:
                        plddt_scores.append(plddt)
    return plddt_scores
```

**Scale auto-detection:** If mean pLDDT < 2.0, the script assumes ESMFold 0–1 scale and multiplies all scores by 100.

```python
if mean < 2.0:
    scores = [x * 100 for x in scores]
    mean = mean * 100
```

**File filtering:** Only files containing `candidate_0` or `IsPETase_reference` in their filename are processed — raw AlphaFold output files in `archive/` are excluded automatically.

**Output:** `results/structure_analysis.txt`

**Results:**

| Candidate | Length (aa) | Tool | Mean pLDDT | Residues >90 | Residues <50 |
|---|---|---|---|---|---|
| IsPETase (reference) | 290 | AlphaFold2 | 92.93 | 260 | 19 |
| A0A365H682 | 288 | ESMFold | 93.92 | 251 | 0 |
| A0A365H6Y8 | 290 | ESMFold | 93.58 | 242 | 3 |
| A0A0F9UIZ8 | 305 | ESMFold | 92.24 | 223 | 0 |
| A0A0F9X315 | 300 | ESMFold | 92.91 | 224 | 0 |
| A0A1E7LM55 | 310 | ESMFold | 94.21 | 251 | 0 |
| A0A386WI52 | 295 | ESMFold | 93.84 | 229 | 1 |
| A0A553PS56 | 588 | AlphaFold2 | 88.60 | 404 | 38 |
| A0A679PDB4 | 284 | ESMFold | 91.92 | 204 | 3 |

All candidates exceed mean pLDDT 88 — sufficient confidence for structural comparison.

### 4D — Structural Alignment (PyMOL)

All candidate structures aligned against IsPETase reference using PyMOL `align` command (sequence-guided superposition + iterative outlier rejection).

Sessions saved: `petase_analysis.pse`, `petase_comparison.pse`

**Key RMSD results:**

| Candidate | RMSD vs IsPETase | Atoms Aligned | Interpretation |
|---|---|---|---|
| A0A365H682 | 0.657 Å | — | Near-identical fold |
| A0A0F9UIZ8 | 0.714 Å | 1,412 | Near-identical fold |
| A0A0F9X315 | 0.715 Å | 1,430 | Near-identical fold |
| A0A553PS56 | 7.68 Å | — | Divergent fold (negative control) |

RMSD < 1 Å = near-identical global fold. Values of 0.657–0.715 Å represent exceptional structural conservation.

**Figures generated** (`results/figures/`):
- `Figure1A_UIZ8_Superposition.png` — A0A0F9UIZ8 vs IsPETase global overlay
- `Figure1B_X315_Superposition.png` — A0A0F9X315 vs IsPETase global overlay
- `Figure1C_PS56_NegativeControl.png` — A0A553PS56 vs IsPETase (negative control)
- `Figure2A_UIZ8_Triad.png` — A0A0F9UIZ8 active site, catalytic triad sticks
- `Figure2B_X315_Triad.png` — A0A0F9X315 active site, catalytic triad sticks

---

## Final Results Summary

**Pipeline:** 21,128 marine sequences → BLAST filter → 780 candidates → ML classification → 8 high-confidence candidates → catalytic triad analysis → structural prediction → structural alignment

| Metric | Value |
|---|---|
| Marine sequences screened | 21,128 |
| BLAST high-confidence hits | 12,525 |
| Unique candidate IDs | 672 |
| Sequences extracted | 780 |
| ML-selected candidates | 8 |
| Catalytic triad retained (S–D–H) | 7 / 8 |
| Mean pLDDT > 90 | 7 / 8 |
| Top candidates by RMSD | A0A365H682 (0.657 Å), A0A0F9UIZ8 (0.714 Å), A0A0F9X315 (0.715 Å) |

**Main finding:** Among eight ML-selected marine candidates, A0A0F9UIZ8 and A0A0F9X315 (marine sediment metagenomes) and A0A365H682 (*Actinomadura craniellae*) emerged as the strongest putative PETases. All retain the conserved Ser–Asp–His catalytic triad and exhibit near-identical three-dimensional structures relative to IsPETase (RMSD < 0.72 Å). Three candidates (A0A1E7LM55, A0A386WI52, A0A679PDB4) annotated only as lipase or alpha/beta hydrolase were reclassified as high-confidence PETase-like by the pipeline.

---

## Limitations

- No experimental wet-lab validation — enzymatic activity on PET substrate unconfirmed
- Structure predictions are computational; crystallographic structures would provide higher confidence
- ML trained on 253 sequences — small dataset, predominantly terrestrial PETases
- Structural alignment performed using predicted (not experimental) structures for both reference and candidates
- BLAST 35% identity threshold may miss highly divergent marine PETases (HMMER flagged as future work)

---

## Future Work

- Molecular docking with PET oligomers (MHET, BHET) to assess active site compatibility
- Molecular dynamics simulations for thermostability and flexibility assessment
- Recombinant expression in *E. coli* and purification
- Experimental PET degradation assays (turbidity, HPLC product detection)
- HMMER profile-based search to recover divergent candidates below BLAST twilight zone
- Expansion to TARA Oceans dataset for broader marine coverage

---

## Key Concepts Glossary

| Term | Definition |
|---|---|
| FASTA | Standard text format for biological sequences |
| Metagenome | DNA sequenced from an environmental community without culturing |
| BLAST | Basic Local Alignment Search Tool — finds similar sequences in a database |
| E-value | Probability a BLAST match is due to chance; lower = more significant |
| pLDDT | Per-residue confidence score from AlphaFold/ESMFold (0–100; >90 = high confidence) |
| RMSD | Root-Mean-Square Deviation — structural difference between two aligned proteins (Å) |
| Random Forest | Ensemble ML algorithm aggregating many decision trees |
| Catalytic triad | Ser–Asp–His residues required for serine hydrolase activity |
| Alpha/beta hydrolase | Protein fold family containing PETases, cutinases, and related esterases |
| Precision | Of all sequences predicted as PETase, fraction that actually are |
| Recall | Of all true PETases in test set, fraction correctly predicted |
| ColabFold | Open-source AlphaFold2 implementation using MMseqs2 MSA |
| ESMFold | Meta's language-model-based structure prediction (no MSA required) |
| Twilight zone | 20–35% sequence identity range where function is uncertain from sequence alone |