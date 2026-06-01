def parse_fasta(file):
    sequences = {}
    current_header = ""
    current_seq = ""
    with open(file,"r") as f:
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

def find_triad_positions(ref_seq_aligned):
    """Find alignment positions of catalytic triad in IsPETase"""
    # Count to residue 160, 206, 237 ignoring gaps
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

# Load aligned sequences
seqs = parse_fasta("data/aligned_candidates.fasta")

# Find IsPETase reference
ref_seq = None
ref_key = None
for header, seq in seqs.items():
    if "IsPETase_reference" in header:
        ref_seq = seq
        ref_key = header
        break

if not ref_seq:
    print("ERROR: IsPETase reference not found in alignment")
    exit()

# Find triad positions in alignment
positions = find_triad_positions(ref_seq)
ser_pos = positions["Ser160"]
asp_pos = positions["Asp206"]
his_pos = positions["His237"]

print(f"\nCatalytic triad alignment positions: Ser={ser_pos}, Asp={asp_pos}, His={his_pos}")
print(f"\n{'Candidate':<25} {'Ser(160)':>10} {'Asp(206)':>10} {'His(237)':>10} {'Triad?':>8}")
print("-"*67)

for header, seq in seqs.items():
    if "IsPETase_reference" in header:
        uid = "IsPETase_ref"
    else:
        uid = header.split("|")[1] if "|" in header else header.split()[0]

    s = seq[ser_pos]
    d = seq[asp_pos]
    h = seq[his_pos]

    triad = "YES" if s=="S" and d=="D" and h=="H" else "NO"
    print(f"{uid:<25} {s:>10} {d:>10} {h:>10} {triad:>8}")

    # Save results
with open("results/catalytic_triad_analysis.txt", "w") as f:
    f.write(f"Catalytic triad alignment positions: Ser={ser_pos}, Asp={asp_pos}, His={his_pos}\n\n")
    f.write(f"{'Candidate':<25} {'Ser(160)':>10} {'Asp(206)':>10} {'His(237)':>10} {'Triad?':>8}\n")
    f.write("-"*67 + "\n")
    for header, seq in seqs.items():
        if "IsPETase_reference" in header:
            uid = "IsPETase_ref"
        else:
            uid = header.split("|")[1] if "|" in header else header.split()[0]
        s = seq[ser_pos]
        d = seq[asp_pos]
        h = seq[his_pos]
        triad = "YES" if s=="S" and d=="D" and h=="H" else "NO"
        f.write(f"{uid:<25} {s:>10} {d:>10} {h:>10} {triad:>8}\n")

print("\nSaved to results/catalytic_triad_analysis.txt")