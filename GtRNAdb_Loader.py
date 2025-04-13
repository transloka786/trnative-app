from Bio import SeqIO

def load_trna_from_gtrnadb(target_aa, fasta_file="Hsapi19-tRNAs.fa"):
    trnas = []
    print(f"🧪 Looking for tRNAs matching: {target_aa}")

    for record in SeqIO.parse(fasta_file, "fasta"):
        # Show some headers to debug
        if len(trnas) < 5:
            print("FASTA HEADER:", record.description)

        # ✅ Match the 3-letter amino acid name (e.g., "Trp", "Val")
        if target_aa in record.description:
            trnas.append({
                "name": record.id,
                "sequence": str(record.seq)
            })

    print(f"✅ Found {len(trnas)} tRNAs for {target_aa}")
    return trnas
