#!/usr/bin/env python3
import json
import re
import hashlib
from pathlib import Path

# --- DIRECTORY SETUP ---
# Script is OUTSIDE, looking into the 'drugs' folder
BASE_DIR = Path(__file__).resolve().parent
INPUT_PATH = BASE_DIR / "drugs" / "blood_cancer_medicines_detailed.txt"
OUTPUT_DIR = BASE_DIR / "cleaned_output"

def clean_medical_text(text):
    """Simple regex to strip junk while keeping medical symbols."""
    text = re.sub(r'\s+', ' ', text)
    return re.sub(r'[^a-zA-Z0-9\s%/+\-\(\)\.,±μ≥≤²³]', '', text).strip()

def process_blood_cancer():
    # 1. Create output folder if it doesn't exist
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    if not INPUT_PATH.exists():
        print(f"❌ Error: Cannot find {INPUT_PATH}")
        return

    seen_hashes = set()
    cleaned_data = []

    print(f"Reading: {INPUT_PATH}")

    # 2. Process the text file line by line
    with open(INPUT_PATH, "r", encoding="utf-8") as f:
        for line in f:
            para = line.strip()
            if len(para.split()) < 15: continue # Skip fragments

            cleaned = clean_medical_text(para)
            
            # Deduplicate using MD5 hash
            para_hash = hashlib.md5(cleaned.lower().encode()).hexdigest()
            if para_hash in seen_hashes: continue
            seen_hashes.add(para_hash)

            # 3. Structure for LLM (Instruction, Input, Output)
            # Tagging clinical markers like 'phase', 'trial', or 'survival'
            has_clinical = any(m in cleaned.lower() for m in ['trial', 'phase', 'survival', 'pfs'])
            
            cleaned_data.append({
                "instruction": "Clean and categorize blood cancer medical text.",
                "input": cleaned,
                "output": f"Contains clinical data: {has_clinical}"
            })

    # 4. Save as JSONL (Training ready)
    output_file = OUTPUT_DIR / "blood_cancer_train.jsonl"
    with open(output_file, "w", encoding="utf-8") as f:
        for entry in cleaned_data:
            f.write(json.dumps(entry) + "\n")

    print(f"✅ Success: {len(cleaned_data)} records saved to {output_file}")

if __name__ == "__main__":
    process_blood_cancer()