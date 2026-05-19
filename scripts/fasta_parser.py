def parse_fasta(file):

    with open(file,"r") as f:
        data= f.read().splitlines()

    DNA_array= []
    current_sequence=""

    for line in data:
        if line.startswith(">"):
            if current_sequence:
                DNA_array.append(current_sequence)
            current_sequence=""
        elif line.strip():
            current_sequence += line.strip()

    if current_sequence:
        DNA_array.append(current_sequence)
    return DNA_array




DNA_array = parse_fasta("data/known_petases.fasta")

print (DNA_array)