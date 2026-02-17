#!/usr/bin/env python3
import json
import re
import unicodedata
from pathlib import Path

# --- DIRECTORY SETUP ---
BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR / "drugs"
OUTPUT_DIR = BASE_DIR / "cleaned_output"

def clean_text_for_llm(text: str) -> str:
    """Standardized cleaning for medical training data."""
    text = unicodedata.normalize("NFKC", text)
    # Strip emojis and artifacts while keeping medical symbols (%, +, -, etc)
    text = re.sub(r"[^\w\s.,;:%()\-+/±μ≥≤²³]", " ", text)
    text = re.sub(r"\[\d+\]", "", text)
    return " ".join(text.split()).strip()

def process_all_cancers():
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    # Target every detailed.txt file in the folder
    target_files = list(INPUT_DIR.glob("*_detailed.txt"))
    
    for file_path in target_files:
        cancer_name = file_path.stem.replace("_medicines_detailed", "").replace("_", " ").title()
        print(f"📖 Reading: {file_path.name}")

        file_data = []
        current_drug = "Unknown Medicine"
        capturing_para = False
        para_buffer = []

        # Read line-by-line to avoid memory/splitting errors
        with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                # 1. Detect Drug Name (Anchor: 💊)
                if "💊" in line:
                    drug_match = re.search(r"(?:💊|DRUG:)\s*([A-Z\s\-()]{3,})", line, re.I)
                    if drug_match:
                        current_drug = drug_match.group(1).split('(')[0].strip().title()
                    continue

                # 2. Detect Paragraph Start (Anchor: 📌)
                if "📌" in line or "PARAGRAPH:" in line.upper():
                    capturing_para = True
                    para_buffer = [] # Reset buffer for new snippet
                    # Add content if any exists on the same line as the emoji
                    line_content = line.split("📌")[-1].replace("PARAGRAPH:", "").strip()
                    if line_content:
                        para_buffer.append(line_content)
                    continue

                # 3. Detect Section End (Anchor: === or ---)
                if capturing_para and (re.match(r"^[=\-]{5,}", line.strip())):
                    if para_buffer:
                        full_text = " ".join(para_buffer)
                        cleaned = clean_text_for_llm(full_text)
                        
                        if len(cleaned.split()) > 10:
                            file_data.append({
                                "instruction": f"Extract clinical outcomes for {current_drug} in {cancer_name}.",
                                "input": cleaned,
                                "output": f"Clinical findings for {current_drug} treatment in {cancer_name}."
                            })
                    capturing_para = False
                    para_buffer = []
                    continue

                # 4. Collect paragraph lines
                if capturing_para:
                    para_buffer.append(line.strip())

        # Save Results
        if file_data:
            out_stem = file_path.stem
            with open(OUTPUT_DIR / f"{out_stem}.jsonl", "w", encoding="utf-8") as f_out:
                for entry in file_data:
                    f_out.write(json.dumps(entry, ensure_ascii=False) + "\n")
            
            with open(OUTPUT_DIR / f"{out_stem}_pure.txt", "w", encoding="utf-8") as f_txt:
                f_txt.write("\n\n".join([e["input"] for e in file_data]))
                
            print(f"   ✅ SUCCESS: {len(file_data)} snippets extracted for {cancer_name}.")
        else:
            print(f"   ⚠️ FAILED: No data found in {file_path.name}. Check if 📌 is followed by text.")

if __name__ == "__main__":
    process_all_cancers()