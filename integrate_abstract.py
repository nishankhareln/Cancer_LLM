"""
════════════════════════════════════════════════════════════════
  Abstract Integrator — Updated for Your Exact Folder Structure
  
  SOURCE (abstracts):
  LLM/datasets/
    ├── Bloodcancer/abstract/     → 1.txt, 2.txt, 3.txt ...
    ├── brain cancer/abstract/
    ├── LiverCancer/abstract/
    ├── Lungscancer/abstract/
    ├── Mouth cancer/abstract/
    ├── Ovarian cancer/abstract/
    ├── prostate cancer/abstract/
    └── skin cancer/abstract/
  
  DESTINATION (main training data):
  Cleaned_Cancer_LLM_Data/
    ├── Blood_Cancer/
    ├── Brain_Cancer/
    ├── Liver_Cancer/
    ├── Lung_Cancer/
    ├── Mouth_Cancer/
    ├── Ovarian_Cancer/
    ├── Protaste_Cancer/
    └── Skin_Cancer/
  
  Run: python integrate_abstracts_v2.py
════════════════════════════════════════════════════════════════
"""

import shutil
from pathlib import Path

# ─────────────────────────────────────────────────────────────
# ✏️  EDIT THESE TWO PATHS ONLY
# ─────────────────────────────────────────────────────────────

LLM_DATASETS_ROOT  = Path("D:/bilogy/LLM/datasets")           # ← abstract source
MAIN_DATA_ROOT     = Path("D:/bilogy/Cleaned_Cancer_LLM_Data") # ← training destination

# ─────────────────────────────────────────────────────────────
# EXACT MAPPING:  LLM/datasets folder name  →  Cleaned_Cancer_LLM_Data folder name
# (handles all the name mismatches you have)
# ─────────────────────────────────────────────────────────────

FOLDER_MAP = {
    "Bloodcancer":    "Blood_Cancer",
    "blood cancer":   "Blood_Cancer",
    "brain cancer":   "Brain_Cancer",
    "Brain cancer":   "Brain_Cancer",
    "LiverCancer":    "Liver_Cancer",
    "Liver cancer":   "Liver_Cancer",
    "Lungscancer":    "Lung_Cancer",
    "Lungs cancer":   "Lung_Cancer",
    "Lung cancer":    "Lung_Cancer",
    "Mouth cancer":   "Mouth_Cancer",
    "Ovarian cancer": "Ovarian_Cancer",
    "Ovarian Cancer": "Ovarian_Cancer",
    "prostate cancer":"Protaste_Cancer",
    "Prostate cancer":"Protaste_Cancer",
    "skin cancer":    "Skin_Cancer",
    "Skin cancer":    "Skin_Cancer",
}

# ─────────────────────────────────────────────────────────────
# UPDATED CONFIG — includes Ovarian_Cancer and Skin_Cancer
# Copy this into biobert_finetune_pipeline.py CONFIG
# ─────────────────────────────────────────────────────────────

UPDATED_CANCER_TYPES = [
    "Blood_Cancer",
    "Brain_Cancer",
    "Liver_Cancer",
    "Lung_Cancer",
    "Mouth_Cancer",
    "Ovarian_Cancer",   # ← NEW
    "Protaste_Cancer",
    "Skin_Cancer",      # ← NEW
]


def integrate(dry_run=False):
    print("════════════════════════════════════════════════════════")
    print("  Abstract Integrator v2 — Exact Folder Mapping")
    print("════════════════════════════════════════════════════════")
    if dry_run:
        print("  🔍 DRY RUN — no files will be copied\n")
    else:
        print()

    if not LLM_DATASETS_ROOT.exists():
        print(f"❌ Not found: {LLM_DATASETS_ROOT}")
        return
    if not MAIN_DATA_ROOT.exists():
        print(f"❌ Not found: {MAIN_DATA_ROOT}")
        return

    total_copied = 0
    total_skipped = 0
    summary = []

    # Scan all folders inside LLM/datasets/
    for cancer_folder in sorted(LLM_DATASETS_ROOT.iterdir()):
        if not cancer_folder.is_dir():
            continue

        folder_name = cancer_folder.name

        # Match to Cleaned_Cancer_LLM_Data name
        main_cancer_name = FOLDER_MAP.get(folder_name)
        if not main_cancer_name:
            print(f"⚠️  No mapping for '{folder_name}' — add it to FOLDER_MAP")
            continue

        # Find the abstract subfolder
        abstract_dir = cancer_folder / "abstract"
        if not abstract_dir.exists():
            # Try alternate names
            for alt in ["abstracts", "Abstract", "Abstracts"]:
                alt_path = cancer_folder / alt
                if alt_path.exists():
                    abstract_dir = alt_path
                    break
            else:
                print(f"⚠️  No abstract/ subfolder in: {cancer_folder.name}")
                continue

        # Destination: Cleaned_Cancer_LLM_Data/{Cancer}/Abstracts/
        dest_dir = MAIN_DATA_ROOT / main_cancer_name / "Abstracts"
        if not dry_run:
            dest_dir.mkdir(parents=True, exist_ok=True)

        txt_files = sorted(abstract_dir.glob("*.txt"),
                           key=lambda x: int(x.stem) if x.stem.isdigit() else 0)

        if not txt_files:
            print(f"⚠️  No .txt files in {abstract_dir}")
            continue

        print(f"📁 {folder_name}  →  {main_cancer_name}/Abstracts/")

        copied = 0
        skipped = 0

        for txt_file in txt_files:
            text = txt_file.read_text(encoding="utf-8", errors="ignore").strip()
            words = len(text.split())

            if words < 30:
                print(f"   ⏭️  {txt_file.name}  ({words} words) — too short, skipped")
                skipped += 1
                total_skipped += 1
                continue

            # Rename to avoid conflicts: 1.txt → blood_cancer_abstract_1.txt
            safe_name = f"{main_cancer_name.lower()}_abstract_{txt_file.name}"
            dest_file = dest_dir / safe_name

            if dry_run:
                print(f"   📄 {txt_file.name} → {safe_name}  ({words} words)")
            else:
                shutil.copy2(txt_file, dest_file)
                print(f"   ✅ {txt_file.name} → {safe_name}  ({words} words)")

            copied += 1
            total_copied += 1

        summary.append((main_cancer_name, len(txt_files), copied, skipped))
        print()

    # ── Summary table ──
    print("════════════════════════════════════════════════════════")
    print("  SUMMARY")
    print("════════════════════════════════════════════════════════")
    print(f"\n  {'Cancer Type':<25} {'Total':>7} {'Copied':>8} {'Skipped':>9}")
    print(f"  {'─'*25} {'─'*7} {'─'*8} {'─'*9}")
    for name, total, copied, skipped in summary:
        print(f"  {name:<25} {total:>7} {copied:>8} {skipped:>9}")
    print(f"  {'─'*25} {'─'*7} {'─'*8} {'─'*9}")
    print(f"  {'TOTAL':<25} {'':>7} {total_copied:>8} {total_skipped:>9}")

    if not dry_run and total_copied > 0:
        print(f"""
  ✅  {total_copied} abstracts copied!

  ══ NEXT: Update biobert_finetune_pipeline.py ═════════════

  1. Update CONFIG cancer_types (add Ovarian + Skin):

     "cancer_types": [
         "Blood_Cancer",
         "Brain_Cancer",
         "Liver_Cancer",
         "Lung_Cancer",
         "Mouth_Cancer",
         "Ovarian_Cancer",    ← ADD
         "Protaste_Cancer",
         "Skin_Cancer",       ← ADD
     ],

  2. Add "Abstracts" to CONFIG categories:

     "categories": [
         "Biology", "Diagnostics",
         "Prognosis", "Treatment",
         "Abstracts",             ← ADD
     ],

  3. Run:  python biobert_finetune_pipeline.py
  ══════════════════════════════════════════════════════════
""")


def verify():
    """Check final state after copying."""
    print("\n════════════════════════════════════════════════════════")
    print("  VERIFICATION — Abstract counts in training folder")
    print("════════════════════════════════════════════════════════\n")

    for cancer in sorted(UPDATED_CANCER_TYPES):
        cancer_path = MAIN_DATA_ROOT / cancer
        abstracts_path = cancer_path / "Abstracts"
        if abstracts_path.exists():
            files = list(abstracts_path.glob("*.txt"))
            words = sum(
                len(f.read_text(encoding="utf-8", errors="ignore").split())
                for f in files
            )
            print(f"  ✅  {cancer:<25} → {len(files):>3} abstracts  |  ~{words:>6,} words")
        else:
            print(f"  ❌  {cancer:<25} → Abstracts/ folder not found")
    print()


if __name__ == "__main__":

    # ── Step 1: Preview ──
    integrate(dry_run=True)

    # ── Step 2: Confirm ──
    print("─" * 55)
    confirm = input("\n  Type 'yes' to copy all abstracts: ").strip().lower()

    if confirm == "yes":
        print()
        integrate(dry_run=False)
        verify()
    else:
        print("\n  Cancelled. No files copied.")