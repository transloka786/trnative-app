
import streamlit as st
import re
import joblib

st.set_page_config(page_title="tRNative Suppressor tRNA Engine", layout="wide")
st.title("🔬 tRNative – Final Integrated App")

# Inputs
gene = st.text_input("🧬 Gene Symbol (e.g., CFTR)")
mutation = st.text_input("🧬 Mutation (e.g., W1282X)")
stop_codon = st.selectbox("🛑 Stop Codon Introduced", ["UAA", "UAG", "UGA"])

# Mutation parser
def infer_target_aa_from_mutation(mut):
    match = re.match(r"([A-Z])\d+[X*]", mut.upper())
    if match:
        aa_1_to_3 = {
            "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
            "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
            "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
            "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"
        }
        return aa_1_to_3.get(match.group(1))
    return None

# Load model
try:
    model = joblib.load("tRNA_readthrough_model.pkl")
    st.success("✅ ML model loaded.")
except Exception as e:
    st.error(f"❌ Failed to load ML model: {e}")
    st.stop()

# Logic
if gene and mutation and stop_codon:
    aa = infer_target_aa_from_mutation(mutation)
    if aa:
        st.info(f"🧠 Inferred amino acid to restore: **{aa}**")
    else:
        st.warning("⚠️ Unable to parse amino acid from mutation.")
else:
    st.info("Please enter gene, mutation, and stop codon to begin.")
