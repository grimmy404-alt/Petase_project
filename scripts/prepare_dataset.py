#fasta parser

def parse_fasta(file):
    with open(file,"r") as f:
        data= f.read().splitlines()

    sequences= []
    current_sequence=""

    for line in data:
        if line.startswith(">"):
            if current_sequence:
                sequences.append(current_sequence)
            current_sequence=""
        elif line.strip():
            current_sequence += line.strip()

    if current_sequence:
        sequences.append(current_sequence)
    return sequences

#loading data

petases = parse_fasta("data/known_petases.fasta")
non_petases = parse_fasta("data/non_petases.fasta")

#creating labels

labels= [1]*len(petases) + [0]*len(non_petases)
all_seqs= petases + non_petases

print(f"PETases: {len(petases)}")
print(f"Non-PETases: {len(non_petases)}")
print(f"Total: {len(all_seqs)}")

# Saving

import pickle
with open("data/dataset.pkl", "wb") as f:
    pickle.dump({"sequences": all_seqs, "labels": labels}, f)

print("Dataset saved!")