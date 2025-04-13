import random

def codon_to_anticodon(stop_codon):
    """Convert a stop codon (RNA) to its suppressor anticodon (DNA)"""
    complement = {'A': 'T', 'U': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[b] for b in stop_codon[::-1])

def modify_anticodon(seq, new_anticodon, anticodon_pos=(34, 37)):
    """
    Replace the anticodon in the tRNA sequence.
    Assumes anticodon is at position 34-36 (0-based indexing).
    """
    return seq[:anticodon_pos[0]] + new_anticodon + seq[anticodon_pos[1]:]

def extract_features_for_variants(seq, stop_codon="UGA"):
    """
    Build a suppressor tRNA by modifying the anticodon,
    then simulate single-nucleotide variants and extract features.
    """
    suppressor_seq = modify_anticodon(seq, codon_to_anticodon(stop_codon))
    variants = []

    for i in range(len(suppressor_seq)):
        for nt in "ACGT":
            if nt != suppressor_seq[i]:
                new_seq = suppressor_seq[:i] + nt + suppressor_seq[i+1:]
                gc = (new_seq.count("G") + new_seq.count("C")) / len(new_seq)
                dG = round(random.uniform(-30, -10), 2)
                variants.append({
                    "position": i,
                    "ref": suppressor_seq[i],
                    "alt": nt,
                    "sequence": new_seq,
                    "gc_content": round(gc, 3),
                    "dummy_deltaG": dG,
                    "is_valid_structure": 1
                })

    return variants
