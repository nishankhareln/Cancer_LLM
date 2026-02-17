#!/usr/bin/env python3
import json
import re
import hashlib
from pathlib import Path
from typing import List, Dict, Set
from datetime import datetime

BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR / "prevention"
OUTPUT_DIR = BASE_DIR / "prevention_cleaned_output"

CLINICAL_MARKERS = {
    "trial_phase": r'\b(phase\s*[1-4]|clinical\s*trial|randomized|RCT)\b',
    "survival_metrics": r'\b(overall\s*survival|OS|progression-free\s*survival|PFS|DFS|median|HR)\b',
    "drug_types": r'\b(chemotherapy|targeted|immunotherapy|kinase\s*inhibitor|monoclonal|antibody)\b',
    "genetics": r'\b(mutation|fusion|translocation|BCR-ABL|FLT3|TP53|WGS|sequencing)\b',
    "cancer_types": r'\b(cancer|carcinoma|lymphoma|myeloma|leukemia|sarcoma|melanoma)\b',
    "side_effects": r'\b(toxicity|adverse|safety|thrombocytopenia|infection|cytopenias)\b'
}

class BloodCancerProcessor:
    def __init__(self):
        self.seen_hashes: Set[str] = set()
        self.cleaned_data: List[Dict] = []
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    def clean_medical_text(self, text: str) -> str:
        text = re.sub(r'\s+', ' ', text.strip())
        text = re.sub(r'[^a-zA-Z0-9\s%/+\-\(\)\.,±μ≥≤²³\'\"]', '', text)
        return text.strip()

    def extract_clinical_features(self, text: str) -> Dict[str, bool]:
        text_lower = text.lower()
        return {
            name: bool(re.search(pattern, text_lower, re.IGNORECASE))
            for name, pattern in CLINICAL_MARKERS.items()
        }

    def get_hash(self, text: str) -> str:
        normalized = re.sub(r'\s+', ' ', text.lower().strip())
        return hashlib.md5(normalized.encode()).hexdigest()

    def process_text_file(self, filepath: Path):
        if not filepath.exists():
            print("❌ File not found:", filepath)
            return

        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                text = line.strip()
                if text and len(text.split()) >= 15:
                    text_hash = self.get_hash(text)
                    if text_hash in self.seen_hashes:
                        continue
                    self.seen_hashes.add(text_hash)

                    cleaned = self.clean_medical_text(text)
                    features = self.extract_clinical_features(cleaned)

                    self.cleaned_data.append({
                        "text": cleaned,
                        "clinical_features": features,
                        "has_clinical_data": any(features.values())
                    })

    def save_outputs(self):
        # JSONL
        jsonl_path = OUTPUT_DIR / "prostate_cancer_train.jsonl"
        with open(jsonl_path, "w", encoding="utf-8") as f:
            for entry in self.cleaned_data:
                llm_entry = {
                    "instruction": "Analyze blood cancer prevention medical text.",
                    "input": entry["text"],
                    "output": "Clinical features detected"
                }
                f.write(json.dumps(llm_entry) + "\n")

        # TXT
        txt_path = OUTPUT_DIR / "prostate_cancer_train.txt"
        with open(txt_path, "w", encoding="utf-8") as f:
            for entry in self.cleaned_data:
                f.write(entry["text"] + "\n\n")

        print("✅ Processed:", len(self.cleaned_data))
        print("📁 Saved in:", OUTPUT_DIR)

    def run(self):
        print("\n🔬 prostate CANCER PREVENTION PROCESSOR\n")

        file_path = INPUT_DIR / "prostate_cancer_prevention_detailed.txt"
        self.process_text_file(file_path)

        if self.cleaned_data:
            self.save_outputs()
        else:
            print("⚠️ No valid data found.")


if __name__ == "__main__":
    processor = BloodCancerProcessor()
    processor.run()
