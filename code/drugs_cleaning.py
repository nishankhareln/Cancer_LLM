#!/usr/bin/env python3
"""
════════════════════════════════════════════════════════════════
  Drug Data Cleaner — TXT + JSONL
  
  Cleans ALL drug files in your drugs/ folder:
    ├── liver_cancer_medicines_detailed.txt    → cleaned .txt
    ├── liver_cancer_medicines_detailed.json   → cleaned .jsonl
    ├── prostate_cancer_medicines_detailed.txt → cleaned .txt
    ├── prostate_cancer_medicines_detailed.json→ cleaned .jsonl
    └── ... (any cancer)
  
  Output goes to:  cleaned_output/
    ├── liver_cancer_medicines_detailed_pure.txt   ← clean text
    ├── liver_cancer_medicines_detailed.jsonl      ← clean jsonl
    ├── prostate_cancer_medicines_detailed_pure.txt
    ├── prostate_cancer_medicines_detailed.jsonl
    └── ...
  
  Run: python clean_drug_data.py
════════════════════════════════════════════════════════════════
"""

import json
import re
import hashlib
from pathlib import Path

# ─────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────

BASE_DIR   = Path(__file__).resolve().parent
INPUT_DIR  = BASE_DIR / "drugs"           # ← your drugs folder
OUTPUT_DIR = BASE_DIR / "cleaned_drugs"  # ← output folder


# ═══════════════════════════════════════════════════════════════
# TEXT CLEANING
# ═══════════════════════════════════════════════════════════════

def clean_medical_text(text: str) -> str:
    """
    Clean medical text while keeping all important medical symbols.
    Removes: special unicode junk, excess whitespace, empty lines
    Keeps  : letters, numbers, medical symbols, punctuation
    """
    if not text:
        return ""

    # Normalize whitespace
    text = re.sub(r'\s+', ' ', text)

    # Keep medical-relevant characters:
    # letters, digits, spaces, medical symbols, punctuation
    text = re.sub(
        r'[^a-zA-Z0-9\s%/+\-\(\)\.,;:\'\"\!?±μαβγ≥≤²³°×→↑↓#&@]',
        '', text
    )

    # Remove lines that are just separators (===, ---, ███)
    text = re.sub(r'^[=\-_\s\|#█─]{3,}$', '', text, flags=re.MULTILINE)

    # Remove extra blank lines
    text = re.sub(r'\n{3,}', '\n\n', text)

    return text.strip()


def is_good_paragraph(text: str, min_words: int = 15) -> bool:
    """Check if a paragraph has enough content to be useful."""
    if not text or len(text.split()) < min_words:
        return False
    # Skip separator lines
    if re.match(r'^[=\-_\s\|#█─]{5,}$', text):
        return False
    # Skip lines that are just labels like "MECHANISM:" or "DOSING:"
    if re.match(r'^[A-Z\s]{3,20}:$', text.strip()):
        return False
    return True


def get_hash(text: str) -> str:
    return hashlib.md5(text.lower().strip().encode()).hexdigest()


def detect_cancer_type(filename: str) -> str:
    """Extract cancer type from filename."""
    name = filename.lower()
    for cancer in ["liver", "blood", "brain", "lung", "mouth",
                   "ovarian", "prostate", "skin"]:
        if cancer in name:
            return cancer.title() + " Cancer"
    return "Cancer"


def has_clinical_data(text: str) -> bool:
    """Check if text contains clinical trial or outcome data."""
    markers = [
        'trial', 'phase', 'survival', 'pfs', 'os', 'orr',
        'response rate', 'median', 'hazard ratio', 'hr ',
        'confidence interval', 'p value', 'p=', 'randomized',
        'efficacy', 'safety', 'adverse', 'toxicity', 'dose'
    ]
    text_lower = text.lower()
    return any(m in text_lower for m in markers)


# ═══════════════════════════════════════════════════════════════
# PROCESS TXT FILES
# ═══════════════════════════════════════════════════════════════

def process_txt_file(input_path: Path, output_dir: Path):
    """
    Reads a drug .txt file → cleans → saves:
      1. _pure.txt  : clean plain text (one paragraph per line)
      2. .jsonl     : structured {instruction, input, output}
    """
    cancer_type  = detect_cancer_type(input_path.name)
    stem         = input_path.stem
    out_txt      = output_dir / f"{stem}_pure.txt"
    out_jsonl    = output_dir / f"{stem}.jsonl"

    print(f"\n  📄 {input_path.name}")

    seen_hashes  = set()
    clean_paras  = []
    jsonl_records= []

    raw_text = input_path.read_text(encoding="utf-8", errors="ignore")

    # Split on double newlines or section separators
    paragraphs = re.split(r'\n{2,}|={5,}|─{5,}', raw_text)

    for para in paragraphs:
        para = para.strip()
        if not para:
            continue

        cleaned = clean_medical_text(para)

        if not is_good_paragraph(cleaned):
            continue

        # Deduplicate
        h = get_hash(cleaned)
        if h in seen_hashes:
            continue
        seen_hashes.add(h)

        clean_paras.append(cleaned)

        # Build JSONL record
        has_clinical = has_clinical_data(cleaned)
        instruction  = (
            f"Extract and analyze {cancer_type} drug treatment information."
            if not has_clinical else
            f"Extract clinical outcomes for drug treatment in {cancer_type}."
        )
        jsonl_records.append({
            "instruction": instruction,
            "input":       cleaned,
            "output":      f"{'Clinical trial data' if has_clinical else 'Drug information'} for {cancer_type} treatment."
        })

    # Save pure text
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(f"# {cancer_type} — Cleaned Drug Data\n")
        f.write(f"# Source: {input_path.name}\n")
        f.write(f"# Records: {len(clean_paras)}\n\n")
        f.write("\n\n".join(clean_paras))

    # Save JSONL
    with open(out_jsonl, "w", encoding="utf-8") as f:
        for record in jsonl_records:
            f.write(json.dumps(record, ensure_ascii=False) + "\n")

    print(f"     Paragraphs found   : {len(paragraphs)}")
    print(f"     After cleaning     : {len(clean_paras)}")
    print(f"     Saved pure .txt    : {out_txt.name}")
    print(f"     Saved .jsonl       : {out_jsonl.name}")

    return len(clean_paras)


# ═══════════════════════════════════════════════════════════════
# PROCESS JSON / JSONL FILES
# ═══════════════════════════════════════════════════════════════

def process_json_file(input_path: Path, output_dir: Path):
    """
    Reads a drug .json/.jsonl file → cleans → saves:
      1. _pure.txt  : clean plain text extracted from json
      2. .jsonl     : cleaned structured records
    """
    cancer_type  = detect_cancer_type(input_path.name)
    stem         = input_path.stem
    out_txt      = output_dir / f"{stem}_pure.txt"
    out_jsonl    = output_dir / f"{stem}_clean.jsonl"

    print(f"\n  📋 {input_path.name}")

    seen_hashes  = set()
    clean_paras  = []
    jsonl_records= []
    raw_records  = []

    # Load all records
    with open(input_path, encoding="utf-8", errors="ignore") as f:
        content = f.read().strip()

    # Handle both JSON array and JSONL (one per line)
    if content.startswith('['):
        try:
            raw_records = json.loads(content)
        except json.JSONDecodeError:
            raw_records = []
    else:
        for line in content.splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                raw_records.append(json.loads(line))
            except json.JSONDecodeError:
                continue

    print(f"     Raw records loaded : {len(raw_records)}")

    for rec in raw_records:
        # Extract text from common field names
        text = (
            rec.get("input") or
            rec.get("paragraph") or
            rec.get("text") or
            rec.get("content") or
            rec.get("description") or
            rec.get("body") or ""
        )

        # Also grab drug details if present
        drug_details = rec.get("drug_details", [])
        extra_text   = ""
        for d in drug_details:
            if isinstance(d, dict):
                parts = []
                for field in ["intro", "mechanism", "special_notes"]:
                    val = d.get(field, "")
                    if val:
                        parts.append(val)
                extra_text += " ".join(parts) + " "

        # Combine paragraph text with drug detail text
        full_text = (text + " " + extra_text).strip()
        cleaned   = clean_medical_text(full_text)

        if not is_good_paragraph(cleaned):
            continue

        # Deduplicate
        h = get_hash(cleaned)
        if h in seen_hashes:
            continue
        seen_hashes.add(h)

        clean_paras.append(cleaned)

        # Preserve original instruction/output if present,
        # otherwise build a new one
        instruction = rec.get("instruction") or (
            f"Extract clinical outcomes for drug treatment in {cancer_type}."
            if has_clinical_data(cleaned) else
            f"Extract and analyze {cancer_type} drug treatment information."
        )
        output_val = rec.get("output") or (
            f"{'Clinical trial data' if has_clinical_data(cleaned) else 'Drug information'} "
            f"for {cancer_type} treatment."
        )

        # Add drug names if available
        drugs_found = rec.get("drugs_found", [])
        if drugs_found:
            instruction = f"Extract clinical outcomes for {', '.join(drugs_found[:3])} in {cancer_type}."
            output_val  = f"Clinical findings for {', '.join(drugs_found[:3])} treatment in {cancer_type}."

        jsonl_records.append({
            "instruction": instruction,
            "input":       cleaned,
            "output":      output_val,
        })

    # Save pure text
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(f"# {cancer_type} — Cleaned Drug Data (from JSON)\n")
        f.write(f"# Source: {input_path.name}\n")
        f.write(f"# Records: {len(clean_paras)}\n\n")
        f.write("\n\n".join(clean_paras))

    # Save cleaned JSONL
    with open(out_jsonl, "w", encoding="utf-8") as f:
        for record in jsonl_records:
            f.write(json.dumps(record, ensure_ascii=False) + "\n")

    print(f"     After cleaning     : {len(clean_paras)}")
    print(f"     Saved pure .txt    : {out_txt.name}")
    print(f"     Saved clean .jsonl : {out_jsonl.name}")

    return len(clean_paras)


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    print("═" * 62)
    print("  Drug Data Cleaner — TXT + JSON/JSONL")
    print("  Input : drugs/")
    print("  Output: cleaned_output/")
    print("═" * 62)

    if not INPUT_DIR.exists():
        print(f"\n❌ drugs/ folder not found at: {INPUT_DIR}")
        print("   Make sure your drug extractor outputs are in drugs/")
        return

    OUTPUT_DIR.mkdir(exist_ok=True)

    # Find only LIVER CANCER files
    txt_files  = list(INPUT_DIR.glob("liver_cancer*.txt"))
    json_files = (list(INPUT_DIR.glob("liver_cancer*.json")) +
                  list(INPUT_DIR.glob("liver_cancer*.jsonl")))

    print(f"\n  Found: {len(txt_files)} .txt files")
    print(f"  Found: {len(json_files)} .json/.jsonl files")

    if not txt_files and not json_files:
        print("\n  ❌ No files found in drugs/ folder")
        return

    total_records = 0

    # Process TXT files
    if txt_files:
        print(f"\n  ── Processing TXT files ──")
        for f in sorted(txt_files):
            # Skip summary files
            if "SUMMARY" in f.name or "statistics" in f.name.lower():
                print(f"  ⏭️  Skipping: {f.name}")
                continue
            total_records += process_txt_file(f, OUTPUT_DIR)

    # Process JSON files
    if json_files:
        print(f"\n  ── Processing JSON/JSONL files ──")
        for f in sorted(json_files):
            if "statistics" in f.name.lower():
                print(f"  ⏭️  Skipping: {f.name}")
                continue
            total_records += process_json_file(f, OUTPUT_DIR)

    # Final summary
    out_files = list(OUTPUT_DIR.glob("*"))
    print(f"\n{'═'*62}")
    print(f"  ✅ CLEANING COMPLETE!")
    print(f"  Total clean records : {total_records:,}")
    print(f"  Output files        : {len(out_files)}")
    print(f"  Output folder       : {OUTPUT_DIR}")
    print(f"\n  Files created:")
    for f in sorted(out_files):
        size_kb = f.stat().st_size / 1024
        print(f"    {f.name:<50} {size_kb:>8.1f} KB")
    print(f"{'═'*62}\n")


if __name__ == "__main__":
    main()