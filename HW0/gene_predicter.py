import subprocess
import pandas as pd

def read_fasta(file_path: str) -> Dict[str, str]:
    """
    Read a FASTA file and return a dictionary of sequences.
    """
    sequences = {}
    current_sequence_id = None
    with open(file_path, "r") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"): #read only seq lines
                current_sequence_id = line[1:]
                sequences[current_sequence_id] = ""
            else:
                sequences[current_sequence_id] += line
    return sequences


def find_orfs(sequence: str, min_len: int = 10) -> List[tuple]:
    """
    Find open reading frames (ORFs) in a DNA sequence.
    """
    orfs = []
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    
    for i in range(len(sequence) - 2):
        if sequence[i:i + 3] == start_codon:
            for j in range(i + 3, len(sequence) - 2, 3):
                codon = sequence[j:j + 3]
                if codon in stop_codons and j - i >= min_len:
                    orfs.append((i, j + 2))
                    break
    
    return orfs

def run_blastx(sequence: str, db_path: str) -> pd.DataFrame:
    """
    Run BLASTX on a DNA sequence against a specified database.
    """
    cmd = ["blastx", "-query", sequence, "-db", db_path, "-outfmt", "6 qseqid sseqid stitle"]
    blastx_output = subprocess.check_output(cmd, universal_newlines=True)
    blastx_df = pd.read_csv(pd.compat.StringIO(blastx_output), sep='\t', header=0)
    return blastx_df

def write_bed(out_bed: str, sequences: Dict[str, str]) -> None:
    with open(out_bed, "w") as bed_file:
        for seq_id, sequence in sequences.items():
            orfs = find_orfs(sequence)
            for i, j in orfs:
                blastx_df = run_blastx(sequence, "nr-prot") # Run BLASTX
                description = blastx_df.iloc[0, 2] # Get protein id
                description = description.replace(' ', '_')
                bed_file.write(f"{seq_id}\t{i}\t{j}\t{description}\n")


if __name__ == "__main__":
    in_fasta = "gcvP.fna"
    out_bed = "predicted_genes.bed"
    sequences = read_fasta(in_fasta)
    print(sequences)
    write_bed(out_bed, sequences)   
    