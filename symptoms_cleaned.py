#!/usr/bin/env python3
import json
import re
import hashlib
from pathlib import Path
from typing import List, Dict, Set
from datetime import datetime

BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR / "symptoms"
OUTPUT_DIR = BASE_DIR / "symptoms_cleaned_output"


class SymptomsProcessor:
    def __init__(self, cancer_type: str):
        self.cancer_type = cancer_type
        self.seen_hashes: Set[str] = set()
        self.cleaned_data: List[Dict] = []
        self.metadata = {
            "cancer_type": cancer_type,
            "processed_at": datetime.now().isoformat(),
            "deduplicated": 0,
            "total_records": 0
        }

        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        INPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ------------------------
    # Cleaning
    # ------------------------
    def clean_text(self, text: str) -> str:
        text = re.sub(r'\s+', ' ', text.strip())
        text = re.sub(r'[^a-zA-Z0-9\s%/+\-\(\)\.,±μ≥≤²³\'\"]', '', text)
        return text.strip()

    def get_hash(self, text: str) -> str:
        normalized = re.sub(r'\s+', ' ', text.lower().strip())
        return hashlib.md5(normalized.encode()).hexdigest()

    def is_duplicate(self, text: str) -> bool:
        text_hash = self.get_hash(text)
        if text_hash in self.seen_hashes:
            self.metadata["deduplicated"] += 1
            return True
        self.seen_hashes.add(text_hash)
        return False

    # ------------------------
    # File Processing
    # ------------------------
    def process_json_file(self, filepath: Path) -> int:
        count = 0
        if not filepath.exists():
            return 0

        try:
            with open(filepath, "r", encoding="utf-8") as f:
                data = json.load(f)

            records = data if isinstance(data, list) else [data]

            for record in records:
                text = None
                if isinstance(record, dict):
                    text = (
                        record.get("paragraph")
                        or record.get("text")
                        or record.get("content")
                        or record.get("description")
                    )
                else:
                    text = str(record)

                if text and len(text.split()) >= 10:
                    if not self.is_duplicate(text):
                        cleaned = self.clean_text(text)
                        self.cleaned_data.append({
                            "text": cleaned,
                            "source": "json"
                        })
                        count += 1

        except json.JSONDecodeError:
            pass

        return count

    def process_text_file(self, filepath: Path) -> int:
        count = 0
        if not filepath.exists():
            return 0

        try:
            with open(filepath, "r", encoding="utf-8") as f:
                for line in f:
                    text = line.strip()
                    if text and len(text.split()) >= 10:
                        if not self.is_duplicate(text):
                            cleaned = self.clean_text(text)
                            self.cleaned_data.append({
                                "text": cleaned,
                                "source": "text"
                            })
                            count += 1
        except Exception:
            pass

        return count

    # ------------------------
    # Save Outputs
    # ------------------------
    def save_outputs(self):
        jsonl_path = OUTPUT_DIR / f"{self.cancer_type}_symptoms_train.jsonl"
        txt_path = OUTPUT_DIR / f"{self.cancer_type}_symptoms_train.txt"

        # JSONL (LLM ready)
        with open(jsonl_path, "w", encoding="utf-8") as f:
            for entry in self.cleaned_data:
                llm_entry = {
                    "instruction": f"Analyze {self.cancer_type} cancer symptoms.",
                    "input": entry["text"],
                    "output": f"{self.cancer_type.capitalize()} cancer symptom information."
                }
                f.write(json.dumps(llm_entry) + "\n")

        # TXT
        with open(txt_path, "w", encoding="utf-8") as f:
            for entry in self.cleaned_data:
                f.write(entry["text"] + "\n\n")

        print(f"✅ JSONL saved: {jsonl_path}")
        print(f"✅ TXT saved: {txt_path}")
        print(f"📊 Total Records: {len(self.cleaned_data)}")

    def process(self):
        print(f"\n{'='*60}")
        print(f"{self.cancer_type.upper()} CANCER SYMPTOMS PROCESSOR")
        print(f"{'='*60}\n")

        json_file = INPUT_DIR / f"{self.cancer_type}_cancer_symptoms_ONLY.json"
        text_file = INPUT_DIR / f"{self.cancer_type}_cancer_symptoms_ONLY.txt"

        json_count = self.process_json_file(json_file)
        text_count = self.process_text_file(text_file)

        self.metadata["total_records"] = len(self.cleaned_data)

        print(f"JSON Records: {json_count}")
        print(f"TEXT Records: {text_count}")
        print(f"Removed Duplicates: {self.metadata['deduplicated']}")
        print(f"Total Cleaned: {len(self.cleaned_data)}\n")

        if self.cleaned_data:
            self.save_outputs()
        else:
            print("⚠️ No valid data found.")


# ------------------------
# Run
# ------------------------
if __name__ == "__main__":
    processor = SymptomsProcessor("mouth")   # change to blood, lung, oral, etc.
    processor.process()
