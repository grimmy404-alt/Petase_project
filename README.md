# Machine Learning-Guided Discovery of Novel PETases from Marine Metagenomes

## Overview

Plastic pollution is a major environmental challenge affecting marine ecosystems worldwide. Certain microbial enzymes known as PETases can degrade PET polymers, offering potential biotechnological solutions for plastic waste remediation.

The goal of this project is to explore whether machine learning and sequence-based bioinformatics approaches can help identify potential PETase-like enzymes from marine microbial datasets.

The workflow combines:
- sequence dataset construction,
- homology-based screening,
- feature extraction,
- machine learning classification,
- and structural prediction tools.

The project is being developed as a reproducible Linux-based bioinformatics pipeline using WSL, Python, Conda, and Git/GitHub.

---

## Objectives

- Build a curated dataset of known PETases and non-PETase homologs
- Develop preprocessing and sequence validation workflows
- Extract sequence-derived biochemical and structural features
- Train machine learning models to distinguish PETase-like enzymes
- Screen marine metagenomic datasets for novel candidate enzymes
- Prioritize candidates for downstream structural and functional analysis

---

## Biological Background

Polyethylene terephthalate (PET) is a synthetic polyester widely used in plastic bottles, textiles, and packaging materials. Due to its high stability and large-scale production, PET accumulates extensively in terrestrial and marine environments.

Recent discoveries of PET-degrading enzymes such as PETase from *Ideonella sakaiensis* have demonstrated the possibility of enzymatic PET biodegradation. However, many currently known PETases suffer from limitations including:
- low thermal stability,
- reduced efficiency against highly crystalline PET,
- and poor tolerance to saline marine conditions.

Marine microbial ecosystems remain relatively underexplored as a potential source of novel polyester-degrading enzymes.

---

## Current Project Status

### Completed
- Linux/WSL bioinformatics environment setup
- Conda-based project environment
- FASTA dataset parsing workflow
- Collection of known PETase sequences
- Collection of non-PETase homolog sequences
- Git/GitHub version control integration

### In Progress
- Metadata extraction pipeline
- Sequence validation
- Feature engineering

### Planned
- Multiple sequence alignment
- HMMER profile construction
- Machine learning model training
- Marine metagenome screening
- Structural prediction and candidate prioritization

---

## Project Structure

```text
petase_project/
├── data/
├── docs/
├── logs/
├── models/
├── notebooks/
├── results/
└── scripts/

## Technologies and Tools

- Python
- Biopython
- Conda
- Git/GitHub
- WSL (Windows Subsystem for Linux)
- SeqKit
- BLAST
- HMMER
- MAFFT
- Scikit-learn

## Future Directions

Future work may include:
- advanced feature engineering,
- integration of protein language models,
- AlphaFold structural analysis,
- molecular docking,
- and experimental candidate prioritization.

## Disclaimer

This project is an independent educational and research-oriented bioinformatics study currently under active development.This project is being developed as a learning-oriented bioinformatics pipeline combining sequence analysis, Linux workflows, and machine learning methods.

