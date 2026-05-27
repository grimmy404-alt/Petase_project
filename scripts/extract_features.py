import pickle
import numpy as np

# Amino acid groups

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

def extract_features(seq):
    """Extract simple features from sequence"""
    features = []
    
    # 1. Amino acid composition (20 features)
    for aa in AMINO_ACIDS:
        count = seq.count(aa)
        freq = count / len(seq)
        features.append(freq)
    
    # 2. Length
    features.append(len(seq))
    
    # 3. Aromatic residue content [Aromatic residues directly interact with PET polymer]
    aromatic = sum(seq.count(aa) for aa in 'FYW')
    features.append(aromatic / len(seq))
    
    # 4. Hydrophobic content [Enzymes that bind PET need hydrophobic surface]
    hydrophobic = sum(seq.count(aa) for aa in 'AILMFVPWG')
    features.append(hydrophobic / len(seq))
    
    # 5. Charged residue content [Charged residues affect pH stability and enzyme behavior]
    charged = sum(seq.count(aa) for aa in 'DEKR')
    features.append(charged / len(seq))
    
    return features

# Load dataset

with open("data/dataset.pkl", "rb") as f:
    data = pickle.load(f)

sequences = data["sequences"]
labels = data["labels"]

# Extract features for all sequences

print("Extracting features...")
X = []
for i, seq in enumerate(sequences):
    if i % 50 == 0:
        print(f"  {i}/{len(sequences)}")
    features = extract_features(seq)
    X.append(features)

X = np.array(X)
y = np.array(labels)

print(f"\nFeature matrix shape: {X.shape}")
print(f"Labels shape: {y.shape}")

# Save

np.save("data/features.npy", X)
np.save("data/labels.npy", y)
print("Features saved!")