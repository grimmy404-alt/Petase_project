import os
import glob

def parse_pdb_plddt(filepath):
    """Extract pLDDT from PDB file"""
    plddt_scores = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                if atom_name == "CA":
                    try:
                        plddt = float(line[60:66].strip())
                        plddt_scores.append(plddt)
                    except:
                        pass
    return plddt_scores

def parse_cif_plddt(filepath):
    """Extract pLDDT from CIF file - column 15"""
    plddt_scores = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("ATOM"):
                parts = line.split()
                try:
                    if parts[3] == "CA":
                        plddt = float(parts[14])  # Column 15 (0-indexed = 14)
                        if 0 <= plddt <= 100:
                            plddt_scores.append(plddt)
                except:
                    pass
    return plddt_scores

def analyze_structure(filepath, name):
    """Analyze structure and return metrics"""
    if filepath.endswith(".cif"):
        scores = parse_cif_plddt(filepath)
    else:
        scores = parse_pdb_plddt(filepath)

    if not scores:
        return None

    mean = sum(scores)/len(scores)

    # Detect scale (ESMFold uses 0-1, AlphaFold uses 0-100)
    if mean < 2.0:
        scores = [x * 100 for x in scores]
        mean = mean * 100

    high = sum(1 for x in scores if x > 90)
    low = sum(1 for x in scores if x < 50)

    return {
        "name": name,
        "length": len(scores),
        "mean_plddt": round(mean, 2),
        "high_conf": high,
        "low_conf": low
    }

def main():
    structures_dir = "results/structures"
    results = []

    # Get all PDB and CIF files
    files = glob.glob(f"{structures_dir}/*.pdb") + \
            glob.glob(f"{structures_dir}/*.cif")

    # Only keep candidate and reference files
    files = [f for f in files if 
         ("candidate_0" in os.path.basename(f) or "IsPETase_reference" in f)
         and "fold_candidate" not in os.path.basename(f)]

    for filepath in sorted(files):
        name = os.path.basename(filepath).replace(".pdb","").replace(".cif","")
        result = analyze_structure(filepath, name)
        if result:
            results.append(result)

    # Print table
    print(f"\n{'Name':<35} {'Length':>8} {'Mean pLDDT':>12} {'>90%':>8} {'<50%':>8}")
    print("-"*75)
    for r in results:
        print(f"{r['name']:<35} {r['length']:>8} {r['mean_plddt']:>12} {r['high_conf']:>8} {r['low_conf']:>8}")

    # Save results
    with open("results/structure_analysis.txt", "w") as f:
        f.write(f"{'Name':<35} {'Length':>8} {'Mean pLDDT':>12} {'>90%':>8} {'<50%':>8}\n")
        f.write("-"*75 + "\n")
        for r in results:
            f.write(f"{r['name']:<35} {r['length']:>8} {r['mean_plddt']:>12} {r['high_conf']:>8} {r['low_conf']:>8}\n")

    print("\nSaved to results/structure_analysis.txt")

if __name__ == "__main__":
    main()