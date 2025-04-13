
def extract_features_for_variants(seq):
    import random
    variants = []
    for i in range(3):
        new_seq = seq[:i] + random.choice('ACGT') + seq[i+1:]
        variants.append({
            "position": i,
            "ref": seq[i],
            "alt": new_seq[i],
            "sequence": new_seq,
            "gc_content": 0.5,
            "dummy_deltaG": -17.2,
            "is_valid_structure": 1
        })
    return variants
