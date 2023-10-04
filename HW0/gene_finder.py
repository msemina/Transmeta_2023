 
import re

def find_genes(dna_sequence, start_codon, stop_codons):
    gene_coordinates = []
    start_positions = [match.start() for match in re.finditer(start_codon, dna_sequence)]

    for start_pos in start_positions:
        for stop_codon in stop_codons:
            stop_match = re.search(stop_codon, dna_sequence[start_pos:])
            if stop_match:
                stop_pos = start_pos + stop_match.start() + 3 # triplet code
                gene_coordinates.append((start_pos, stop_pos))
                break  # if stop codon found

    return gene_coordinates


# Merge overlapping gene coordinates
def merge_overlapping(coordinates):
    merged = []
    coordinates.sort(key=lambda x: x[0])
    for start, end in coordinates:
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def write_bed(filename, coordinates, test_seq):
    with open(filename, "w") as bed_file:
        for i, (start, end) in enumerate(coordinates, start=1):
            name = test_seq[start:end]
            bed_file.write(f"chr1\t{start}\t{end}\tGene{i}_{name}\n")


# CONST
start_cdn = 'ATG'
stop_cdn = ['TAG', 'TGA', 'TAA']


# MAIN
text_file = open("/home/msem/Documents/Study/ITMO/transcript/0/data/gcvP.fna", "r")

test_seq = text_file.read()
predicted_genes = find_genes(test_seq, start_cdn, stop_cdn)
merged_genes = merge_overlapping(predicted_genes)

write_bed('predicted_genes.bed', merged_genes, test_seq)

for i, (start, end) in enumerate(merged_genes, start=1):
    name = test_seq[start:end]
    print(f'Predicted Gene {i}: Start: {start}, End: {end}, Name: {name}')

text_file.close()