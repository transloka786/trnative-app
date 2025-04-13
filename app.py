import streamlit as st
import pandas as pd
import joblib
import requests
import os
from tRNA_Feature_Extractor import extract_features_for_variants
from GtRNAdb_Loader import load_trna_from_gtrnadb
from vienna_folding import get_vienna_dG

# ✅ Must be first Streamlit command
st.set_page_config(page_title="tRNative tRNA Engine", layout="wide")
st.write("📁 Files in repo:", os.listdir())

# ✅ Load ML model
try:
    model = joblib.load("tRNA_structure_model.pkl")
except Exception as e:
    st.error("❌ Failed to load ML model. Make sure 'tRNA_structure_model.pkl' is in your repo.")
    st.stop()

# ✅ Supported amino acids
supported_aas = [
    "Gln", "Arg", "Lys", "Asn", "Glu", "Asp", "Tyr", "Cys", "Trp",
    "Ser", "Leu", "Gly", "Ala", "Val", "Ile", "Thr", "Phe", "Pro", "Met", "His"
]

# ✅ User Inputs
st.title("🔬 tRNative - Suppressor tRNA Prediction Engine")
st.header("🧬 Mutation Details")
gene = st.text_input("Gene (optional)")
mutation_pos = st.text_input("Nonsense Mutation Position (e.g., G542X)")
stop_codon = st.selectbox("Stop Codon Introduced", ["UAA", "UAG", "UGA"])
target_aa = st.selectbox("Amino Acid to Restore", supported_aas)

# ✅ Utility functions
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

# ✅ Load and process tRNAs
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
    st.header("🔁 Mutated Suppressor Candidate Scoring")

    for entry in tRNAs:
        st.subheader(f"Variants for {entry['name']}")
        features = extract_features_for_variants(entry["sequence"], stop_codon=stop_codon)

        # Apply ViennaRNA folding ΔG
        for row in features:
            row["dummy_deltaG"] = get_vienna_dG(row["sequence"]) or round(random.uniform(-30, -10), 2)

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
