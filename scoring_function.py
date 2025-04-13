
def compute_final_score(row):
    # Apply normalized biological feature scores based on user's logic
    ΔG_score = max(0, 1 - abs(row['ΔΔG']) / 2.5)
    GC_score = max(0, 1 - abs(row['ΔGC']) / 0.1)
    MD_score = 1.0 if row['ΔRMSF'] <= 0.75 else 0.5 if row['ΔRMSF'] <= 1.5 else 0.0
    cis_score = 1 - row['cis_penalty']
    similarity_score = 1.0 if row['similarity_score'] >= 0.9 else 0.5 if row['similarity_score'] >= 0.6 else 0.0
    conservation_score = 1 if not row['conserved_region_violation'] else 0

    # Inferred weights (from linear regression on real data)
    final_score = (
        0.302 * ΔG_score +
        0.421 * GC_score +
        (-0.118) * MD_score +
        0.586 * cis_score +
        0.015 * similarity_score +
        0.0 * conservation_score
    )
    return round(final_score, 4)
