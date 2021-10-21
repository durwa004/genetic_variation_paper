with open("GCA_000002305.1_EquCab2.0_genomic.fna", "r") as input_file, open("GCA_000002305.1_EquCab2.0_genomic_renamed.fna", "w") as output_file:
    for line in input_file:
        line = line.rstrip("\n").split()
        if ">" in line[0]:
            print(line[0], file = output_file)
        else:
            print("".join(line), file = output_file)
