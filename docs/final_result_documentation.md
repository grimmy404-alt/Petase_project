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

Structural superposition was performed in PyMOL using the predicted IsPETase structure as the reference model. The two highest-confidence candidates, A0A0F9UIZ8 and A0A0F9X315, were selected for detailed structural comparison because they combined:

- Maximum machine learning confidence (1.00)
- Complete catalytic triad conservation
- High-confidence structural predictions
- Marine metagenomic origin

### Global Structural Similarity

Both candidates aligned extremely well with the IsPETase reference structure.

| Candidate | RMSD (Å) | Interpretation |
|------------|-----------|---------------|
| A0A0F9UIZ8 | 0.714 | Near-identical fold |
| A0A0F9X315 | 0.715 | Near-identical fold |

RMSD values below 1 Å indicate exceptional three-dimensional conservation. Structural overlays showed that the overall α/β-hydrolase architecture, central β-sheet, surrounding α-helices, and catalytic loop regions were highly conserved relative to IsPETase.

### Figure 1. Structural Superposition of Top Candidates with IsPETase

**Figure 1A.** Structural superposition of IsPETase (cyan) and candidate A0A0F9UIZ8 (magenta).

**Figure 1B.** Structural superposition of IsPETase (cyan) and candidate A0A0F9X315 (magenta).

The nearly complete overlap between structures visually confirms the low RMSD values obtained during alignment.

---

### Catalytic Triad Spatial Conservation

Catalytic triad analysis identified the following active-site residues:

| Protein | Catalytic Ser | Catalytic Asp | Catalytic His |
|----------|--------------|--------------|--------------|
| IsPETase | S160 | D206 | H237 |
| A0A0F9UIZ8 | S174 | D218 | H250 |
| A0A0F9X315 | S170 | D214 | H246 |

Although the residue numbering differs due to sequence insertions and deletions, the catalytic residues occupy equivalent structural positions within the active site.

### Figure 2. Catalytic Triad Superposition

**Figure 2A.** Active-site overlay of IsPETase and A0A0F9UIZ8 showing conservation of catalytic Ser-Asp-His geometry.

**Figure 2B.** Active-site overlay of IsPETase and A0A0F9X315 showing conservation of catalytic Ser-Asp-His geometry.

The catalytic residues converge in nearly identical spatial arrangements, supporting preservation of the classical serine hydrolase catalytic mechanism.

---

### Negative Control: A0A553PS56

Candidate A0A553PS56 was included as a negative structural control because catalytic triad analysis revealed non-canonical substitutions:

| Position | IsPETase | A0A553PS56 |
|----------|----------|------------|
| Catalytic Ser | S160 | S206 |
| Catalytic Asp | D206 | N244 |
| Catalytic His | H237 | E288 |

Structural alignment produced an RMSD of approximately **7.68 Å**, indicating substantial divergence from the PETase fold.

### Figure 3. Structural Superposition of A0A553PS56 with IsPETase

Structural overlay of IsPETase and A0A553PS56 demonstrates major deviations in overall fold architecture and active-site organization.

Combined with the catalytic residue substitutions, these observations suggest that A0A553PS56 is unlikely to function as a classical PETase despite receiving a high machine learning confidence score.

---

### Summary

The combined sequence and structural analyses identify A0A0F9UIZ8 and A0A0F9X315 as the strongest PETase-like candidates recovered from marine metagenomic datasets. Both proteins exhibit:

- Complete catalytic triad conservation
- High-confidence predicted structures
- Near-identical global folds relative to IsPETase
- Sub-angstrom RMSD values
- Conserved active-site geometry

These results provide strong computational evidence that both candidates belong to the PETase family and warrant experimental validation.

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
* Structural alignment produced RMSD values of approximately 0.71 Å, supporting strong evolutionary and probable functional similarity.

These results demonstrate that the pipeline successfully identified previously uncharacterized marine proteins with strong sequence and structural characteristics of known PET-degrading enzymes.
