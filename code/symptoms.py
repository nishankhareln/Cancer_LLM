#!/usr/bin/env python3
"""
Enhanced PubMed Central Symptom Extractor - PROSTATE CANCER
Extracts ONLY symptom-related text with comprehensive keyword list
"""
from Bio import Entrez
import xml.etree.ElementTree as ET
import re
import time
import json

# ================= CONFIG =================
Entrez.email = "nkharel57@gmail.com"  # REQUIRED

SEARCH_QUERY = (
    "prostate cancer symptoms OR prostate cancer clinical presentation OR "
    "prostate cancer signs OR PSA elevated symptoms OR "
    "CRPC symptoms OR metastatic prostate cancer symptoms OR "
    "prostate cancer urinary symptoms OR prostate cancer bone pain OR "
    "advanced prostate cancer symptoms OR prostate cancer diagnosis signs"
)
MAX_ARTICLES = 1200
OUTPUT_TEXT = "prostate_cancer_symptoms_ONLY.txt"
OUTPUT_JSON = "prostate_cancer_symptoms_ONLY.json"

# ================= PROSTATE CANCER SYMPTOM KEYWORDS =================
SYMPTOM_KEYWORDS = [
    # Urinary symptoms (most common)
    "urinary frequency", "frequent urination", "nocturia", "urinary urgency",
    "urgency to urinate", "urinary hesitancy", "weak urine stream",
    "poor urine stream", "urinary retention", "incomplete bladder emptying",
    "straining to urinate", "difficulty urinating", "urinary obstruction",
    "bladder obstruction", "urinary incontinence", "leaking urine",
    "blood in urine", "hematuria", "urinary tract infection",
    "dysuria", "painful urination", "burning urination",
    "urinary dribbling", "interrupted urine stream",

    # Erectile and sexual symptoms
    "erectile dysfunction", "impotence", "sexual dysfunction",
    "decreased libido", "loss of libido", "ejaculatory dysfunction",
    "painful ejaculation", "blood in semen", "hematospermia",
    "decreased ejaculate", "ejaculation problems",

    # Pelvic and perineal symptoms
    "pelvic pain", "pelvic discomfort", "perineal pain",
    "perineal discomfort", "rectal pressure", "rectal pain",
    "lower abdominal pain", "suprapubic pain", "suprapubic discomfort",
    "lower pelvic pressure",

    # Bone pain (metastatic)
    "bone pain", "back pain", "lower back pain", "hip pain",
    "spinal pain", "bone metastasis", "bone metastases",
    "skeletal pain", "rib pain", "pelvic bone pain",
    "femur pain", "vertebral pain", "pathological fracture",
    "spinal cord compression", "bone fracture",

    # Leg symptoms (metastatic/lymph node)
    "leg swelling", "leg edema", "lower limb edema",
    "leg weakness", "leg pain", "limb weakness",
    "lower extremity weakness", "leg numbness",
    "deep vein thrombosis", "dvt",

    # Lymph node symptoms
    "lymph node enlargement", "lymphadenopathy",
    "inguinal swelling", "groin swelling", "swollen lymph nodes",
    "pelvic lymph node", "retroperitoneal lymph",

    # Constitutional / systemic symptoms
    "weight loss", "unexplained weight loss", "unintentional weight loss",
    "loss of appetite", "anorexia", "cachexia", "wasting",
    "fatigue", "tiredness", "weakness", "malaise",
    "exhaustion", "lethargy", "lack of energy",
    "fever", "night sweats", "chills",

    # PSA related
    "elevated psa", "rising psa", "psa elevation", "psa level",
    "psa progression", "biochemical recurrence", "psa recurrence",
    "psa doubling time",

    # Neurological (spinal cord compression)
    "numbness", "tingling", "paresthesia", "weakness in legs",
    "difficulty walking", "gait disturbance", "urinary incontinence",
    "bowel incontinence", "fecal incontinence", "saddle anesthesia",
    "paraplegia", "neurological deficit",

    # Rectal symptoms
    "rectal bleeding", "blood in stool", "rectal obstruction",
    "bowel habit change", "constipation", "rectal fullness",

    # Prostate-specific
    "enlarged prostate", "prostate enlargement", "hard prostate",
    "irregular prostate", "nodule on prostate", "prostate nodule",
    "abnormal digital rectal exam", "digital rectal examination",
    "prostatitis", "prostate inflammation",

    # Advanced disease
    "anemia", "shortness of breath", "dyspnea",
    "pleural effusion", "ascites", "jaundice",
    "renal failure", "ureteral obstruction", "hydronephrosis",
    "lymphedema",
]

SECTION_ALLOW = [
    "symptom", "clinical", "presentation", "sign",
    "manifestation", "features", "findings",
    "complaint", "patient", "diagnosis", "physical examination",
    "clinical feature", "case report", "case presentation",
    "urinary", "pelvic", "bone", "metastatic", "advanced"
]

# ================= FUNCTIONS =================

def search_pmc():
    try:
        handle = Entrez.esearch(db="pmc", term=SEARCH_QUERY, retmax=MAX_ARTICLES)
        record = Entrez.read(handle)
        return record["IdList"]
    except Exception as e:
        print(f"❌ Search error: {e}")
        return []

def fetch_xml(pmc_id):
    try:
        handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="full", retmode="xml")
        return handle.read()
    except Exception as e:
        print(f"⚠️ Fetch error for {pmc_id}: {e}")
        return None

def is_symptom_section(title):
    if not title:
        return False
    return any(key in title.lower() for key in SECTION_ALLOW)

def extract_symptom_text(xml_data):
    if not xml_data:
        return []
    try:
        root = ET.fromstring(xml_data)
    except:
        return []

    symptom_paragraphs = []
    current_section_is_symptom = False

    for elem in root.iter():
        tag = elem.tag.lower()

        if tag.endswith("title") and elem.text:
            current_section_is_symptom = is_symptom_section(elem.text)

        if tag.endswith("p") and elem.text:
            text = elem.text.strip().lower()
            if current_section_is_symptom:
                if any(keyword in text for keyword in SYMPTOM_KEYWORDS):
                    clean = re.sub(r'\s+', ' ', elem.text.strip())
                    if len(clean) > 150:
                        symptom_paragraphs.append(clean)

    return symptom_paragraphs

def count_symptoms_in_text(text):
    text_lower = text.lower()
    return [kw for kw in SYMPTOM_KEYWORDS if kw in text_lower]

# ================= MAIN =================

def main():
    print("\n" + "="*60)
    print("PROSTATE CANCER SYMPTOM EXTRACTOR")
    print("="*60 + "\n")

    print("🔍 Searching PubMed Central for PROSTATE CANCER symptoms...")
    pmc_ids = search_pmc()

    if not pmc_ids:
        print("❌ No articles found")
        return

    print(f"✅ Found {len(pmc_ids)} articles\n")

    all_symptoms = []

    for i, pmc_id in enumerate(pmc_ids, 1):
        print(f"📄 Processing {pmc_id} ({i}/{len(pmc_ids)})")
        try:
            xml = fetch_xml(pmc_id)
            if xml:
                symptoms = extract_symptom_text(xml)
                for s in symptoms:
                    found = count_symptoms_in_text(s)
                    all_symptoms.append({
                        "pmc_id":              pmc_id,
                        "text":                s,
                        "symptom_count":       len(found),
                        "symptoms_mentioned":  found,
                    })
            time.sleep(0.5)
        except Exception as e:
            print(f"⚠️ Error {pmc_id}: {e}")

    # Save text
    with open(OUTPUT_TEXT, "w", encoding="utf-8") as f:
        for i, item in enumerate(all_symptoms, 1):
            f.write(f"Context {i} (PMC:{item['pmc_id']}) - {item['symptom_count']} symptoms:\n")
            top = item['symptoms_mentioned'][:5]
            more = "..." if len(item['symptoms_mentioned']) > 5 else ""
            f.write(f"Symptoms: {', '.join(top)}{more}\n")
            f.write(item["text"] + "\n\n")

    # Save JSON
    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(all_symptoms, f, indent=2, ensure_ascii=False)

    total = len(all_symptoms)
    avg   = sum(i['symptom_count'] for i in all_symptoms) / total if total > 0 else 0

    print("\n" + "="*60)
    print("✅ EXTRACTION COMPLETE — PROSTATE CANCER SYMPTOMS")
    print("="*60)
    print(f"📝 Text : {OUTPUT_TEXT}")
    print(f"💾 JSON : {OUTPUT_JSON}")
    print(f"📊 Total symptom contexts : {total}")
    print(f"📈 Avg symptoms per context: {avg:.1f}")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()