import pickle
import numpy as np

# Parse data
def parse_fasta(file):
    with open(file,"r") as f:
        data = f.read().splitlines()
    sequences = []
    headers = []
    current_seq = ""
    current_header = ""
    for line in data:
        if line.startswith(">"):
            if current_seq:
                sequences.append(current_seq)
                headers.append(current_header)
            current_seq = ""
            current_header = line[1:]
        elif line.strip():
            current_seq += line.strip()
    if current_seq:
        sequences.append(current_seq)
        headers.append(current_header)
    return headers, sequences

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

# Load model
with open("models/petase_classifier.pkl", "rb") as f:
    model = pickle.load(f)

# Load candidates
headers, sequences = parse_fasta("data/marine_candidates.fasta")
print(f"Loaded {len(sequences)} candidates")

# Deduplicate by header ID
seen = set()
unique_headers = []
unique_sequences = []

for h, s in zip(headers, sequences):
    uid = h.split('|')[1] if '|' in h else h.split()[0]
    if uid not in seen:
        seen.add(uid)
        unique_headers.append(h)
        unique_sequences.append(s)

headers, sequences = unique_headers, unique_sequences
print(f"After deduplication: {len(sequences)} unique candidates")

# Extract features
X = np.array([extract_features(seq) for seq in sequences])

# Predict
predictions = model.predict(X)
probabilities = model.predict_proba(X)[:,1]  # Probability of being PETase

# Combine results
results = list(zip(headers, predictions, probabilities))

# Sort by probability (highest first)
results.sort(key=lambda x: x[2], reverse=True)

# Show top 20
print(f"\nTop 20 PETase candidates:")
print(f"{'Rank':<5} {'Probability':<12} {'Header'}")
print("-"*80)
for i, (header, pred, prob) in enumerate(results[:20]):
    print(f"{i+1:<5} {prob:.4f}       {header[:60]}")

# Save all results
with open("results/predictions.txt", "w") as f:
    f.write("Rank\tProbability\tPrediction\tHeader\n")
    for i, (header, pred, prob) in enumerate(results):
        f.write(f"{i+1}\t{prob:.4f}\t{pred}\t{header}\n")

print(f"\nAll results saved to results/predictions.txt")
print(f"Total predicted PETases: {sum(predictions)}")