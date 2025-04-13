from Bio import SeqIO

def load_trna_from_gtrnadb(target_aa, fasta_file="Hsapi19-tRNAs.fa"):
    aa_map = {
        "Gln": "Q", "Arg": "R", "Lys": "K", "Asn": "N", "Glu": "E", "Asp": "D", "Tyr": "Y", "Cys": "C", "Trp": "W",
        "Ser": "S", "Leu": "L", "Gly": "G", "Ala": "A", "Val": "V", "Ile": "I", "Thr": "T", "Phe": "F", "Pro": "P", "Met": "M", "His": "H"
    }
    aa_code = aa_map.get(target_aa)
    trnas = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        desc = record.description
        if aa_code in desc:
            trnas.append({
                "name": record.id,
                "sequence": str(record.seq)
            })
    return trnas
