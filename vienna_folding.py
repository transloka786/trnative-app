import requests
import re

def get_vienna_dG(seq):
    try:
        url = "http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi"
        headers = {"Content-Type": "application/x-www-form-urlencoded"}
        payload = {"SCREEN": "1", "SEQUENCE": seq}

        response = requests.post(url, data=payload, headers=headers, timeout=15)

        if response.ok:
            # Try to extract ΔG value using regex
            match = re.search(r"dG =\s*(-?\d+\.\d+)", response.text)
            if match:
                return float(match.group(1))
    except Exception as e:
        print(f"ViennaRNA call failed: {e}")
    
    return None  # fallback
