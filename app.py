import streamlit as st
import pandas as pd
import joblib
import requests
import os
import random
from tRNA_Feature_Extractor import extract_features_for_variants
from GtRNAdb_Loader import load_trna_from_gtrnadb
from vienna_folding import get_vienna_dG

# ✅ Streamlit settings
st.set_page_config(page_title="tRNative tRNA Engine", layout="wide")
st.write("📁 Files in repo:", os.listdir())

# ✅ Load regression model
try:
    model = joblib.load("tRNA_readthrough_model.pkl")
except Exception as e:
    st.error("❌ Failed to load ML model. Make sure 'tRNA_readthrough_model.pkl' is in your repo.")
    st.stop()

# Supported AAs
supported_aas = [
    "Gln", "Arg", "Lys", "Asn", "Glu", "Asp", "Tyr", "Cys", "Trp",
    "Ser", "Leu", "Gly", "Ala", "Val", "Ile", "Thr", "Phe", "Pro", "Met", "His"
]

# Inputs
st.title("🔬 tRNative - Suppressor tRNA Prediction Engine")
st.header("🧬 Mutation Details")
gene = st.text_input("Gene (optional)")
mutation_pos = st.text_input("Nonsense Mutation Position (e.g., G542X)")
stop_codon = st.selectbox("Stop Codon Introduced", ["UAA", "UAG", "UGA"])
target_aa = st.selectbox("Amino Acid to Restore", supported_aas)

# Utility functions
def off_target_risk(seq, stop_codon):
    count = seq.count(stop_codon.replace("U", "T"))
    penalty = 0.005 * seq.count("TGA")
    return round((count / len(seq)) + penalty, 4)

def mock_3d_match_score(seq):
    core_motifs = ["TTCGA", "GGG", "CCG"]
    score = sum([1 for motif in core_motifs if motif in seq]) / len(core_motifs)
    return round(score, 2)

def get_rnacomposer_link(seq):
    return f"https://rnacomposer.cs.put.poznan.pl/#sequence={seq}"

# Main workflow
if target_aa:
    st.markdown("---")
    st.header(f"🧾 Natural Human tRNAs for {target_aa} (from GtRNAdb)")

    try:
        tRNAs = load_trna_from_gtrnadb(target_aa)
    except Exception as e:
        st.error(f"❌ Failed to load tRNAs: {e}")
        tRNAs = []

    st.write(f"✅ Loaded {len(tRNAs)} tRNAs for {target_aa}")
    if len(tRNAs) > 0:
        st.write(tRNAs[:2])
    else:
        st.warning("⚠️ No tRNAs loaded. Check FASTA file or matching logic.")

    st.markdown("---")
    st.header("🔁 Suppressor Candidate Evaluation (Regression)")

    for entry in tRNAs:
        st.subheader(f"Variants for {entry['name']}")
        features = extract_features_for_variants(entry["sequence"], stop_codon=stop_codon)

        for row in features:
            row["dummy_deltaG"] = get_vienna_dG(row["sequence"]) or round(random.uniform(-30, -10), 2)
            row["cis_penalty"] = round(random.uniform(0.08, 0.25), 3)  # placeholder until integrated
            row["domain_fluctuation_score"] = round(random.uniform(0.12, 0.3), 3)  # placeholder RMSF
            row["similarity_to_known_suppressors"] = round(random.uniform(0.7, 0.98), 3)  # proxy similarity

        df_feat = pd.DataFrame(features)

        df_feat["predicted_readthrough"] = model.predict(df_feat[[
            "gc_content", "dummy_deltaG", "cis_penalty",
            "domain_fluctuation_score", "similarity_to_known_suppressors"
        ]])

        df_feat = df_feat.rename(columns={"dummy_deltaG": "vienna_dG"})

        df_feat["off_target_risk"] = df_feat["sequence"].apply(lambda x: off_target_risk(x, stop_codon))
        df_feat["fold_3D_score"] = df_feat["sequence"].apply(mock_3d_match_score)
        df_feat["RNAComposer_Link"] = df_feat["sequence"].apply(get_rnacomposer_link)

        top_candidates = df_feat.sort_values("predicted_readthrough", ascending=False).head(10)[[
            "sequence", "gc_content", "vienna_dG", "predicted_readthrough",
            "off_target_risk", "fold_3D_score", "RNAComposer_Link"
        ]]

        st.dataframe(top_candidates, use_container_width=True)

        csv = top_candidates.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download Top Candidates as CSV",
            data=csv,
            file_name=f"{entry['name']}_top_tRNA_candidates.csv",
            mime="text/csv"
        )
