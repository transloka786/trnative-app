import streamlit as st
import pandas as pd
import joblib
import requests
import os
from tRNA_Feature_Extractor import extract_features_for_variants
from GtRNAdb_Loader import load_trna_from_gtrnadb


# Load model
try:
    model = joblib.load("tRNA_structure_model.pkl")
except Exception as e:
    st.error("❌ Failed to load ML model. Make sure 'tRNA_structure_model.pkl' is in your repo.")
    st.stop()

supported_aas = [
    "Gln", "Arg", "Lys", "Asn", "Glu", "Asp", "Tyr", "Cys", "Trp",
    "Ser", "Leu", "Gly", "Ala", "Val", "Ile", "Thr", "Phe", "Pro", "Met", "His"
]

st.set_page_config(page_title="tRNative tRNA Engine", layout="wide")
st.title("🔬 tRNative - Suppressor tRNA Prediction Engine")

st.header("🧬 Mutation Details")
gene = st.text_input("Gene (optional)")
mutation_pos = st.text_input("Nonsense Mutation Position (e.g., G542X)")
stop_codon = st.selectbox("Stop Codon Introduced", ["UAA", "UAG", "UGA"])
target_aa = st.selectbox("Amino Acid to Restore", supported_aas)

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

if target_aa:
    st.markdown("---")
    st.header(f"🧾 Natural Human tRNAs for {target_aa} (from GtRNAdb)")
    try:
        tRNAs = load_trna_from_gtrnadb(target_aa)
        df = pd.DataFrame(tRNAs)
        st.dataframe(df, use_container_width=True)

        st.markdown("---")
        st.header("🔁 Mutated Candidate Scoring")

        for entry in tRNAs[:2]:
            st.subheader(f"Variants for {entry['name']}")
            features = extract_features_for_variants(entry["sequence"])
            df_feat = pd.DataFrame(features)
            df_feat["prediction"] = model.predict(df_feat[["position", "gc_content", "dummy_deltaG"]])
            df_filtered = df_feat[df_feat["prediction"] == 1]

            df_filtered["off_target_risk"] = df_filtered["sequence"].apply(lambda x: off_target_risk(x, stop_codon))
            df_filtered["fold_3D_score"] = df_filtered["sequence"].apply(mock_3d_match_score)
            df_filtered["RNAComposer_Link"] = df_filtered["sequence"].apply(get_rnacomposer_link)

            top_candidates = df_filtered.head(10)[[
                "sequence", "gc_content", "dummy_deltaG",
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
    except Exception as e:
        st.error(f"❌ Failed to load tRNAs: {e}")
