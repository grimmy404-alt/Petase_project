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

| Candidate  | Ser160 | Asp206 | His237 | Conserved Triad |
| ---------- | ------ | ------ | ------ | --------------- |
| A0A365H682 | S      | D      | H      | ✓               |
| A0A365H6Y8 | S      | D      | H      | ✓               |
| A0A0F9UIZ8 | S      | D      | H      | ✓               |
| A0A0F9X315 | S      | D      | H      | ✓               |
| A0A1E7LM55 | S      | D      | H      | ✓               |
| A0A386WI52 | S      | D      | H      | ✓               |
| A0A679PDB4 | S      | D      | H      | ✓               |
| A0A553PS56 | S      | N      | E      | ✗               |

Seven of the eight candidates retained the complete catalytic triad, suggesting preservation of the catalytic machinery characteristic of PET-degrading enzymes.

---

## Structure Prediction

Three-dimensional structures were predicted using deep-learning-based protein structure prediction methods. Prediction quality was assessed using the predicted Local Distance Difference Test (pLDDT).

| Protein            | Length (aa) | Mean pLDDT |
| ------------------ | ----------- | ---------- |
| IsPETase Reference | 290         | 92.93      |
| A0A365H682         | 288         | 93.92      |
| A0A365H6Y8         | 290         | 93.58      |
| A0A0F9UIZ8         | 305         | 92.24      |
| A0A0F9X315         | 300         | 92.91      |
| A0A1E7LM55         | 310         | 94.21      |
| A0A386WI52         | 295         | 93.84      |
| A0A553PS56         | 588         | 88.60      |
| A0A679PDB4         | 284         | 91.92      |

Most candidates achieved mean pLDDT scores above 90, indicating highly reliable structural predictions.

---

## Structural Comparison with IsPETase

Structural superposition was performed in PyMOL using experimentally characterized IsPETase as the reference structure.

The two strongest candidates, A0A0F9UIZ8 and A0A0F9X315, were aligned against IsPETase and evaluated using Root Mean Square Deviation (RMSD).

| Candidate  | RMSD (Å) |
| ---------- | -------- |
| A0A0F9UIZ8 | 0.714    |
| A0A0F9X315 | 0.715    |

Both candidates exhibited RMSD values below 1 Å, demonstrating exceptionally strong structural similarity to the reference PETase.

Structural overlays revealed conservation of the overall α/β-hydrolase fold and catalytic region architecture.

---

## Top Candidates

By combining machine learning predictions, catalytic triad conservation, structural confidence, and structural similarity analyses, two proteins emerged as the strongest PETase-like candidates.

### A0A0F9UIZ8
**Source:** Marine sediment metagenome

- Conserved Ser-Asp-His catalytic triad
- Mean pLDDT: 92.24
- RMSD to IsPETase: 0.714 Å
- Annotated as a PET hydrolase/cutinase-like domain-containing protein

### A0A0F9X315
**Source:** Marine sediment metagenome

- Conserved Ser-Asp-His catalytic triad
- Mean pLDDT: 92.91
- RMSD to IsPETase: 0.715 Å
- Annotated as a PET hydrolase/cutinase-like domain-containing protein

---

## Key Findings

* 8 high-confidence PETase-like candidates were identified.
* 7 of 8 candidates retained the canonical Ser-Asp-His catalytic triad.
* Predicted structures showed high confidence (mean pLDDT generally >90).
* Two marine metagenomic proteins (A0A0F9UIZ8 and A0A0F9X315) exhibited near-identical structural folds to IsPETase.
* Structural alignment produced RMSD values of approximately 0.71 Å, supporting strong evolutionary and functional similarity.

These results demonstrate that the pipeline successfully identified previously uncharacterized marine proteins with strong sequence and structural characteristics of known PET-degrading enzymes.
