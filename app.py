
import streamlit as st
import pandas as pd
import joblib
import requests
import re
import matplotlib.pyplot as plt
from tRNA_Feature_Extractor import extract_features_for_variants
from GtRNAdb_Loader import load_trna_from_gtrnadb
from vienna_folding import get_vienna_dG

st.set_page_config(page_title="tRNative - Final Pre-GNN Version", layout="wide")
st.title("🔬 tRNative – Suppressor tRNA Predictor (Final Random Forest Version)")

# Inputs
gene = st.text_input("🧬 Gene Symbol (e.g., CFTR)")
mutation = st.text_input("🧬 Mutation (e.g., W1282X)")
stop_codon = st.selectbox("🛑 Stop Codon Introduced", ["UAA", "UAG", "UGA"])

aa_1_to_3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
    "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"
}

def infer_target_aa_from_mutation(mut):
    match = re.match(r"([A-Z])\d+[X*]", mut.upper())
    if match:
        return aa_1_to_3.get(match.group(1))
    return None

def get_rnacomposer_link(seq):
    return f"https://rnacomposer.cs.put.poznan.pl/#sequence={seq}"

def mock_3d_match_score(seq):
    return round(sum([1 for m in ["TTCGA", "GGG", "CCG"] if m in seq]) / 3, 2)

def mock_md_rmsf(sequence):
    return {
        "rmsf_acceptor": 1.3,
        "rmsf_dloop": 0.9,
        "rmsf_anticodon": 2.5,
        "rmsf_tloop": 1.6
    }

try:
    model = joblib.load("tRNA_readthrough_model.pkl")
    st.success("✅ ML model loaded.")
except Exception as e:
    st.error(f"❌ Failed to load ML model: {e}")
    st.stop()

target_aa = infer_target_aa_from_mutation(mutation)
mutation_pos = int(re.findall(r"\d+", mutation)[0]) if mutation else None

if target_aa:
    st.write(f"🧠 Inferred amino acid to restore: **{target_aa}**")

if gene and mutation and stop_codon and target_aa and mutation_pos:
    try:
        tRNAs = load_trna_from_gtrnadb(target_aa)
    except Exception as e:
        st.error(f"Failed to load tRNAs: {e}")
        tRNAs = []

    tRNAs = tRNAs[:2]
    all_variants = []

    for entry in tRNAs:
        features = extract_features_for_variants(entry["sequence"], stop_codon=stop_codon)
        for row in features:
            row["parent_tRNA"] = entry["name"]
            row["deltaG"] = get_vienna_dG(row["sequence"]) or -20.0
            row["similarity_to_known_suppressors"] = 0.85
            row["cis_penalty"] = 0.15
            rmsf = mock_md_rmsf(row["sequence"])
            row.update(rmsf)
            row["domain_fluctuation_score"] = round(sum(rmsf.values()) / 4, 2)
            row["RNAComposer_Link"] = get_rnacomposer_link(row["sequence"])
            row["fold_3D_score"] = mock_3d_match_score(row["sequence"])
            all_variants.append(row)

    df = pd.DataFrame(all_variants)
    df.rename(columns={"gc_content": "GC_content"}, inplace=True)
    df["predicted_readthrough"] = model.predict(df[[
        "GC_content", "deltaG", "cis_penalty", "domain_fluctuation_score", "similarity_to_known_suppressors"
    ]])

    top2 = df.sort_values("predicted_readthrough", ascending=False).head(2)

    for _, row in top2.iterrows():
        st.markdown("----")
        st.write(f"**tRNA Parent:** `{row['parent_tRNA']}`")
        st.write(f"**Sequence:** `{row['sequence']}`")
        st.write(f"🎯 **Predicted Readthrough:** `{row['predicted_readthrough']:.3f}`")
        with st.expander("🔬 Score Breakdown"):
            st.write(f"GC Content: {row['GC_content']}")
            st.write(f"ΔG: {row['deltaG']} kcal/mol")
            st.write(f"Cis Penalty: {row['cis_penalty']}")
            st.write(f"Similarity to Known Suppressors: {row['similarity_to_known_suppressors']}")
            st.write(f"Domain Fluctuation Score: {row['domain_fluctuation_score']}")
            st.markdown(f"[📎 RNAComposer Structure]({row['RNAComposer_Link']})")

        st.markdown("📊 MD Fluctuation Plot")
        fig, ax = plt.subplots()
        ax.bar(["Acceptor", "D-loop", "Anticodon", "T-loop"], [
            row["rmsf_acceptor"], row["rmsf_dloop"], row["rmsf_anticodon"], row["rmsf_tloop"]
        ])
        ax.set_ylabel("RMSF (Å)")
        ax.set_title("MD Domain Fluctuation (Mocked)")
        st.pyplot(fig)

    csv = top2.to_csv(index=False).encode("utf-8")
    st.download_button("📥 Download Top 2 Suppressors", data=csv,
                       file_name=f"{gene}_{mutation}_top2_suppressors.csv", mime="text/csv")
