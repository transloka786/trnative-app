from Bio import SeqIO

def load_trna_from_gtrnadb(target_aa, fasta_file="Hsapi19-tRNAs.fa"):
    aa_map = {
        "Gln": "Q", "Arg": "R", "Lys": "K", "Asn": "N", "Glu": "E", "Asp": "D", "Tyr": "Y", "Cys": "C", "Trp": "W",
        "Ser": "S", "Leu": "L", "Gly": "G", "Ala": "A", "Val": "V", "Ile": "I", "Thr": "T", "Phe": "F", "Pro": "P", "Met": "M", "His": "H"
    }
    aa_code = aa_map.get(target_aa)
    trnas = []

    print(f"🧪 Looking for tRNAs matching amino acid: {target_aa} ({aa_code})")

    for record in SeqIO.parse(fasta_file, "fasta"):
        # ✅ Print first few headers for inspection
        if len(trnas) < 5:
            print("FASTA HEADER:", record.description)

        if aa_code in record.description:
            trnas.append({
                "name": record.id,
                "sequence": str(record.seq)
            })

    print(f"✅ Found {len(trnas)} tRNAs for {target_aa}")
    return trnas
