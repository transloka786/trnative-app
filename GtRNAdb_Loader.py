from Bio import SeqIO

def load_trna_from_gtrnadb(target_aa, fasta_file="Hsapi19-tRNAs.fa"):
    """
    Load real human tRNAs for a given amino acid from GtRNAdb FASTA.
    `target_aa`: e.g. "Gln", "Arg", "Tyr"
    """
    aa_map = {
        "Gln": "Q", "Arg": "R", "Lys": "K", "Asn": "N", "Glu": "E", "Asp": "D", "Tyr": "Y", "Cys": "C", "Trp": "W",
        "Ser": "S", "Leu": "L", "Gly": "G", "Ala": "A", "Val": "V", "Ile": "I", "Thr": "T", "Phe": "F", "Pro": "P", "Met": "M", "His": "H"
    }
    aa_code = aa_map.get(target_aa)
    trnas = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        if aa_code and f"AminoAcid={aa_code}" in record.description:
            trnas.append({
                "name": record.id,
                "sequence": str(record.seq)
            })
    return trnas
