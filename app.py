
# tRNative app final integration placeholder with correct mutation parser
import re

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

# The rest of the app logic would go here...
