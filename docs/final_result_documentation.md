# Results

## Candidate Discovery

The machine learning pipeline was applied to marine protein sequences to identify putative PETase-like enzymes. Following feature extraction, classification, confidence filtering, and sequence-based validation, eight high-confidence candidate proteins were selected for downstream analysis.

| Candidate ID | Source Organism | Annotation |
|--------------|----------------|------------|
| A0A365H682 | *Actinomadura craniellae* | Poly(ethylene terephthalate) hydrolase |
| A0A365H6Y8 | *Actinomadura craniellae* | Poly(ethylene terephthalate) hydrolase |
| A0A0F9UIZ8 | Marine sediment metagenome | PET hydrolase/cutinase-like domain-containing protein |
| A0A0F9X315 | Marine sediment metagenome | PET hydrolase/cutinase-like domain-containing protein |
| A0A1E7LM55 | *Streptomyces nanshensis* | Lipase |
| A0A386WI52 | *Micromonospora tulbaghiae* | Alpha/beta hydrolase |
| A0A553PS56 | *Tigriopus californicus* | Carboxylic ester hydrolase |
| A0A679PDB4 | *Streptomyces* sp. SM14 | Alpha/beta hydrolase |

Reference PETase: IsPETase from *Ideonella sakaiensis* 201-F6

These candidates originated from marine and marine-associated organisms and represented the highest-confidence predictions produced by the pipeline.

---

## Catalytic Triad Conservation

PETases belong to the serine hydrolase family and rely on a conserved catalytic triad composed of:

* Serine (Ser)
* Aspartate (Asp)
* Histidine (His)

Multiple sequence alignment was performed against the reference PETase from *Ideonella sakaiensis* (IsPETase). Conserved catalytic positions corresponding to Ser160, Asp206, and His237 were examined in each candidate.

| Candidate | Catalytic Ser | Catalytic Asp | Catalytic His | Conserved Triad |
|-----------|--------------|--------------|--------------|-----------------|
| IsPETase (ref) | S160 | D206 | H237 | ✓ |
| A0A365H682 | S158 | D204 | H236 | ✓ |
| A0A365H6Y8 | S159 | D205 | H237 | ✓ |
| A0A0F9UIZ8 | S174 | D218 | H250 | ✓ |
| A0A0F9X315 | S170 | D214 | H246 | ✓ |
| A0A1E7LM55 | S179 | D225 | H257 | ✓ |
| A0A386WI52 | S166 | D212 | H244 | ✓ |
| A0A679PDB4 | S156 | D202 | H234 | ✓ |
| A0A553PS56 | S206 | N244 | E288 | ✗ |

Seven of the eight candidates retained the complete Ser–Asp–His catalytic triad, suggesting preservation of the catalytic machinery characteristic of PET-degrading enzymes. Candidate A0A553PS56 showed non-canonical substitutions at the Asp and His positions and was retained as a negative structural control.

---

## Structure Prediction

Three-dimensional structures were predicted using deep-learning-based protein structure prediction methods. Prediction quality was assessed using the predicted Local Distance Difference Test (pLDDT).

| Protein | Length (aa) | Mean pLDDT |
|---------|-------------|------------|
| IsPETase Reference | 290 | 92.93 |
| A0A365H682 | 288 | 93.92 |
| A0A365H6Y8 | 290 | 93.58 |
| A0A0F9UIZ8 | 305 | 92.24 |
| A0A0F9X315 | 300 | 92.91 |
| A0A1E7LM55 | 310 | 94.21 |
| A0A386WI52 | 295 | 93.84 |
| A0A553PS56 | 588 | 88.60 |
| A0A679PDB4 | 284 | 91.92 |

Most candidates achieved mean pLDDT scores above 90, indicating highly reliable structural predictions.

---

## Structural Comparison with IsPETase

Structural superposition was performed in PyMOL using the predicted IsPETase structure as the reference model. All seven catalytic triad-positive candidates were superposed against IsPETase. Candidate A0A553PS56, which lacks the canonical triad, served as a negative structural control.

### Global Structural Similarity

All seven triad-positive candidates aligned with exceptional fidelity to the IsPETase reference structure.

| Candidate | Source | RMSD (Å) | Interpretation |
|-----------|--------|-----------|----------------|
| A0A365H6Y8 | *Actinomadura craniellae* | 0.585 | Near-identical fold |
| A0A386WI52 | *Micromonospora tulbaghiae* | 0.585 | Near-identical fold |
| A0A1E7LM55 | *Streptomyces nanshensis* | 0.645 | Near-identical fold |
| A0A365H682 | *Actinomadura craniellae* | 0.657 | Near-identical fold |
| A0A0F9UIZ8 | Marine sediment metagenome | 0.714 | Near-identical fold |
| A0A0F9X315 | Marine sediment metagenome | 0.715 | Near-identical fold |
| A0A679PDB4 | *Streptomyces* sp. SM14 | 0.753 | Near-identical fold |
| A0A553PS56 | *Tigriopus californicus* | 7.68 | Divergent fold (negative control) |

RMSD values below 1 Å indicate exceptional three-dimensional conservation. All seven triad-positive candidates fell within this threshold (range: 0.585–0.753 Å). Structural overlays showed that the overall α/β-hydrolase architecture, central β-sheet, surrounding α-helices, and catalytic loop regions were highly conserved relative to IsPETase.

Notably, two candidates annotated only as "lipase" or "alpha/beta hydrolase" in public databases — A0A365H6Y8 and A0A386WI52 — achieved the lowest RMSD values in the entire dataset (0.585 Å each), matching or exceeding the structural similarity of candidates with explicit PETase annotations.

### Figure 1. Structural Superposition of Candidates with IsPETase

**Figure 1A.** Structural superposition of IsPETase (gray) and A0A1E7LM55 (magenta). RMSD = 0.645 Å.

**Figure 1B.** Structural superposition of IsPETase (gray) and A0A386WI52 (magenta). RMSD = 0.585 Å.

**Figure 1C.** Structural superposition of IsPETase (gray) and A0A679PDB4 (magenta). RMSD = 0.753 Å.

**Figure 1D.** Structural superposition of IsPETase (gray) and A0A0F9UIZ8 (magenta). RMSD = 0.714 Å.

**Figure 1E.** Structural superposition of IsPETase (gray) and A0A0F9X315 (magenta). RMSD = 0.715 Å.

The nearly complete structural overlap visually confirms the low RMSD values and conservation of the α/β-hydrolase fold across all candidates.

---

### Catalytic Triad Spatial Conservation

Although residue numbering differs across candidates due to sequence insertions and deletions, the catalytic residues occupy equivalent structural positions within the active site in all seven triad-positive candidates.

| Protein | Catalytic Ser | Catalytic Asp | Catalytic His |
|---------|--------------|--------------|--------------|
| IsPETase | S160 | D206 | H237 |
| A0A365H682 | S158 | D204 | H236 |
| A0A365H6Y8 | S159 | D205 | H237 |
| A0A0F9UIZ8 | S174 | D218 | H250 |
| A0A0F9X315 | S170 | D214 | H246 |
| A0A1E7LM55 | S179 | D225 | H257 |
| A0A386WI52 | S166 | D212 | H244 |
| A0A679PDB4 | S156 | D202 | H234 |

### Figure 2. Catalytic Triad Active-Site Overlays

Active-site overlays for triad-positive candidates demonstrate convergence of the Ser–Asp–His geometry at nearly identical spatial positions relative to IsPETase, supporting preservation of the classical serine hydrolase catalytic mechanism.

---

### Negative Control: A0A553PS56

Candidate A0A553PS56 was included as a negative structural control because catalytic triad analysis revealed non-canonical substitutions:

| Position | IsPETase | A0A553PS56 |
|----------|----------|------------|
| Catalytic Ser | S160 | S206 |
| Catalytic Asp | D206 | N244 |
| Catalytic His | H237 | E288 |

Structural alignment produced an RMSD of **7.68 Å**, indicating substantial divergence from the PETase fold. Combined with the catalytic residue substitutions, these observations confirm that A0A553PS56 is unlikely to function as a classical PETase despite its high ML confidence score, validating the importance of structural follow-up after sequence-based classification.

### Figure 3. Negative Control — A0A553PS56 vs IsPETase

Structural overlay of IsPETase and A0A553PS56 demonstrates major deviations in overall fold architecture and active-site organization, in contrast to the near-identical overlaps observed for the triad-positive candidates.

---

## Top Candidates

By combining machine learning predictions, catalytic triad conservation, structural confidence, and structural similarity, seven proteins were confirmed as high-confidence PETase-like candidates. Three are of particular interest because their public database annotations (lipase or alpha/beta hydrolase) do not reflect PETase function — representing potential novel functional discoveries by this pipeline.

### Candidates with Explicit PETase Annotation

#### A0A365H682 — *Actinomadura craniellae*
- Conserved Ser158–Asp204–His236 catalytic triad
- Mean pLDDT: 93.92
- RMSD to IsPETase: 0.657 Å
- Annotated as poly(ethylene terephthalate) hydrolase

#### A0A365H6Y8 — *Actinomadura craniellae*
- Conserved Ser159–Asp205–His237 catalytic triad
- Mean pLDDT: 93.58
- RMSD to IsPETase: 0.585 Å (joint lowest in dataset)
- Annotated as poly(ethylene terephthalate) hydrolase

#### A0A0F9UIZ8 — Marine sediment metagenome
- Conserved Ser174–Asp218–His250 catalytic triad
- Mean pLDDT: 92.24
- RMSD to IsPETase: 0.714 Å
- Annotated as PET hydrolase/cutinase-like domain-containing protein

#### A0A0F9X315 — Marine sediment metagenome
- Conserved Ser170–Asp214–His246 catalytic triad
- Mean pLDDT: 92.91
- RMSD to IsPETase: 0.715 Å
- Annotated as PET hydrolase/cutinase-like domain-containing protein

### Candidates with Potentially Incomplete Annotation ⚠️

#### A0A1E7LM55 — *Streptomyces nanshensis*
- Conserved Ser179–Asp225–His257 catalytic triad
- Mean pLDDT: 94.21 (highest in dataset)
- RMSD to IsPETase: 0.645 Å
- Annotated only as **lipase** in UniProt — PETase activity not recorded

#### A0A386WI52 — *Micromonospora tulbaghiae*
- Conserved Ser166–Asp212–His244 catalytic triad
- Mean pLDDT: 93.84
- RMSD to IsPETase: 0.585 Å (joint lowest in dataset)
- Annotated only as **alpha/beta hydrolase** in UniProt — PETase activity not recorded

#### A0A679PDB4 — *Streptomyces* sp. SM14
- Conserved Ser156–Asp202–His234 catalytic triad
- Mean pLDDT: 91.92
- RMSD to IsPETase: 0.753 Å
- Annotated only as **alpha/beta hydrolase** in UniProt — PETase activity not recorded

---

## Key Findings

* 8 high-confidence PETase-like candidates were identified from 21,128 marine sequences.
* 7 of 8 candidates retained the canonical Ser–Asp–His catalytic triad.
* All predicted structures showed high confidence (mean pLDDT 88.6–94.2).
* All 7 triad-positive candidates exhibited near-identical structural folds relative to IsPETase, with RMSD values ranging from **0.585 to 0.753 Å**.
* The lowest RMSD values (0.585 Å) were observed for A0A365H6Y8 and A0A386WI52 — the latter currently annotated only as a lipase.
* 3 candidates (A0A1E7LM55, A0A386WI52, A0A679PDB4) carry no PETase annotation in public databases despite satisfying all computational criteria, representing the primary novel discovery contribution of this pipeline.
* The negative control (A0A553PS56, RMSD 7.68 Å) confirms that ML confidence alone is insufficient — structural validation is essential.

These results demonstrate that the pipeline successfully identified previously uncharacterized proteins with strong sequence and structural characteristics of known PET-degrading enzymes, including candidates whose functional potential is not captured by current database annotations.