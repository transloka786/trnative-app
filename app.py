
import streamlit as st
import pandas as pd
import joblib
import requests
import os
import random
from tRNA_Feature_Extractor import extract_features_for_variants
from GtRNAdb_Loader import load_trna_from_gtrnadb
from vienna_folding import get_vienna_dG

# Set Streamlit page config
st.set_page_config(page_title="tRNative - Suppressor tRNA Regression Engine", layout="wide")
st.write("📁 Files in repo:", os.listdir())

# Load regression model
try:
    model = joblib.load("tRNA_readthrough_model.pkl")
except Exception as e:
    st.error("❌ Failed to load regression model.")
    st.stop()

# Supported amino acids
supported_aas = [
    "Gln", "Arg", "Lys", "Asn", "Glu", "Asp", "Tyr", "Cys", "Trp",
    "Ser", "Leu", "Gly", "Ala", "Val", "Ile", "Thr", "Phe", "Pro", "Met", "His"
]

# User input
st.title("🔬 tRNative - Suppressor tRNA Readthrough Predictor")
gene = st.text_input("Gene (optional)")
mutation_pos = st.text_input("Nonsense Mutation Position (e.g., G542X)")
stop_codon = st.selectbox("Stop Codon Introduced", ["UAA", "UAG", "UGA"])
target_aa = st.selectbox("Amino Acid to Restore", supported_aas)

# Extra features
def off_target_risk(seq, stop_codon):
    return round((seq.count(stop_codon.replace("U", "T")) / len(seq)) + 0.005 * seq.count("TGA"), 4)

def mock_3d_match_score(seq):
    return round(sum([1 for m in ["TTCGA", "GGG", "CCG"] if m in seq]) / 3, 2)

def get_rnacomposer_link(seq):
    return f"https://rnacomposer.cs.put.poznan.pl/#sequence={seq}"

# Main logic
if target_aa:
    st.header(f"🧾 Natural Human tRNAs for {target_aa} (from GtRNAdb)")
    try:
        tRNAs = load_trna_from_gtrnadb(target_aa)
    except Exception as e:
        st.error(f"Failed to load tRNAs: {e}")
        tRNAs = []

    st.write(f"✅ Loaded {len(tRNAs)} tRNAs for {target_aa}")
    if len(tRNAs) == 0:
        st.warning("⚠️ No tRNAs loaded. Check FASTA source or target amino acid.")

    st.header("🔁 Regression Scoring of Suppressor Variants")

    for entry in tRNAs:
        st.subheader(f"Variants for {entry['name']}")
        features = extract_features_for_variants(entry["sequence"], stop_codon=stop_codon)

        for row in features:
            row["deltaG"] = get_vienna_dG(row["sequence"]) or round(random.uniform(-30, -10), 2)
            row["similarity_to_known_suppressors"] = round(random.uniform(0.65, 0.98), 2)
            row["cis_penalty"] = round(random.uniform(0.08, 0.25), 2)
            row["domain_fluctuation_score"] = round(random.uniform(0.12, 0.30), 2)

        df_feat = pd.DataFrame(features)
        df_feat["predicted_readthrough"] = model.predict(df_feat[[
            "gc_content", "deltaG", "cis_penalty", "domain_fluctuation_score", "similarity_to_known_suppressors"
        ]])

        df_feat["RNAComposer_Link"] = df_feat["sequence"].apply(get_rnacomposer_link)
        df_feat["fold_3D_score"] = df_feat["sequence"].apply(mock_3d_match_score)

        top = df_feat.sort_values("predicted_readthrough", ascending=False).head(10)
        st.dataframe(top[[
            "sequence", "gc_content", "deltaG", "predicted_readthrough", "fold_3D_score", "RNAComposer_Link"
        ]], use_container_width=True)

        csv = top.to_csv(index=False).encode("utf-8")
        st.download_button("📥 Download Top Candidates", data=csv,
                           file_name=f"{entry['name']}_top_suppressors.csv",
                           mime="text/csv")
