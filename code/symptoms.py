#!/usr/bin/env python3
"""
Enhanced PubMed Central Symptom Extractor - MOUTH CANCER
Extracts ONLY symptom-related text with comprehensive keyword list
"""
from Bio import Entrez
import xml.etree.ElementTree as ET
import os
import re
import time
import json

# ================= CONFIG =================
Entrez.email = "nkharel57@gmail.com"  # REQUIRED

SEARCH_QUERY = (
    "mouth cancer symptoms OR oral cancer symptoms OR "
    "oral squamous cell carcinoma symptoms OR OSCC symptoms OR "
    "oropharyngeal cancer symptoms OR tongue cancer symptoms OR "
    "lip cancer symptoms OR clinical presentation OR signs"
)
MAX_ARTICLES = 1200
OUTPUT_TEXT = "mouth_cancer_symptoms_ONLY.txt"
OUTPUT_JSON = "mouth_cancer_symptoms_ONLY.json"

# ================= COMPREHENSIVE MOUTH CANCER SYMPTOM KEYWORDS =================
SYMPTOM_KEYWORDS = [
    # Oral ulcers/sores
    "oral ulcer", "mouth ulcer", "mouth sore", "oral sore",
    "ulceration", "ulcerations", "mucosal ulcer", "canker sore",
    "persistent sore", "non-healing ulcer", "wound",
    
    # Oral lesions/masses
    "oral lesion", "mouth lesion", "oral mass", "oral growth",
    "tongue lesion", "lip lesion", "palate lesion",
    "intraoral mass", "lump in mouth", "swelling in mouth",
    "nodule", "bump in mouth",
    
    # Pain in mouth
    "mouth pain", "oral pain", "tongue pain", "lip pain",
    "palatal pain", "oropharyngeal pain", "throat pain",
    "pain in mouth", "pain on tongue", "mouth discomfort",
    "pharyngeal pain", "gum pain", "dental pain",
    
    # Difficulty swallowing
    "difficulty swallowing", "dysphagia", "swallowing difficulty",
    "painful swallowing", "odynophagia", "trouble swallowing",
    "hard to swallow",
    
    # Voice changes
    "hoarseness", "hoarse voice", "voice changes",
    "speech changes", "slurred speech", "weak voice",
    "vocal changes", "dysphonia",
    
    # Difficulty chewing
    "difficulty chewing", "chewing pain", "trouble chewing",
    "impaired mastication", "chewing discomfort",
    
    # Mouth bleeding
    "bleeding from mouth", "mouth bleeding", "gum bleeding",
    "bleeding gums", "oral bleeding", "blood in mouth",
    "tongue bleeding", "lip bleeding",
    
    # Saliva changes
    "dry mouth", "xerostomia", "excessive saliva",
    "salivary changes", "difficulty salivating", "thick saliva",
    
    # Tongue symptoms
    "tongue swelling", "tongue enlargement", "glossitis",
    "tongue lesion", "tongue ulcer", "tongue mass",
    "tongue pain", "tongue numbness", "tongue discoloration",
    
    # Lip symptoms
    "lip swelling", "lip lesion", "lip ulcer", "lip sore",
    "lip thickening", "lip mass", "lip pain",
    "lip discoloration", "lip crust",
    
    # Jaw symptoms
    "jaw pain", "jaw stiffness", "restricted mouth opening",
    "trismus", "limited jaw mobility", "jaw swelling",
    "jaw lumps", "jaw numbness",
    
    # Throat symptoms
    "throat pain", "sore throat", "pharyngitis",
    "throat swelling", "throat mass", "throat lesion",
    "throat discomfort",
    
    # Gum symptoms
    "gum swelling", "gum enlargement", "gingival swelling",
    "gum ulcer", "gum sore", "gum pain", "gum bleeding",
    "gum lesion", "gum discoloration",
    
    # Tooth symptoms
    "loose tooth", "tooth pain", "tooth loss", "dental pain",
    "tooth mobility", "teeth loosening",
    
    # Swollen lymph nodes
    "swollen lymph nodes", "lymph node enlargement",
    "lymphadenopathy", "neck swelling", "cervical lymph",
    "submandibular swelling", "neck mass", "enlarged nodes",
    
    # Ear symptoms
    "ear pain", "otalgia", "ear discomfort", "referred pain",
    
    # Breathing difficulty
    "difficulty breathing", "dyspnea", "shortness of breath",
    "breathing difficulty", "airway obstruction",
    
    # Weight loss
    "weight loss", "unexplained weight loss", "unintentional weight loss",
    "loss of appetite", "anorexia", "cachexia", "wasting",
    
    # Fatigue
    "fatigue", "tiredness", "weakness", "malaise",
    "exhaustion", "lethargy", "lack of energy",
    
    # Fever
    "fever", "low grade fever", "febrile",
    "elevated temperature", "persistent fever",
    
    # Night symptoms
    "night sweats", "sweating at night",
    
    # Neck symptoms
    "neck pain", "neck stiffness", "neck mass",
    "neck enlargement", "cervical swelling",
    
    # Oral bleeding/bleeding disorders
    "hemoptysis", "blood in sputum", "bleeding",
    "hemorrhage", "hematuria",
    
    # Discoloration/Pigmentation
    "oral discoloration", "mouth discoloration",
    "white patch", "white lesion", "red patch", "red lesion",
    "leukoplakia", "erythroplakia", "pigmentation",
    "dark patch", "color change", "mucosal discoloration",
    
    # Odor symptoms
    "bad breath", "halitosis", "mouth odor", "foul odor",
    
    # Numbness/Tingling
    "numbness in mouth", "oral numbness", "paresthesia",
    "tingling in mouth", "burning sensation", "burning mouth",
    
    # Dental/Periodontal
    "periodontal disease", "periodontitis", "dental disease",
    "bad oral hygiene", "caries",
    
    # Nerve symptoms
    "facial numbness", "facial paralysis", "nerve paralysis",
    "nerve involvement", "cranial nerve",
    
    # Swallowing related
    "aspiration", "choking", "coughing while swallowing",
    
    # Systemic symptoms
    "general malaise", "feeling unwell", "systemic symptoms",
    
    # Smoking/Tobacco related
    "smokeless tobacco", "betel quid", "tobacco use",
    "smoking history", "tobacco chewing"
]

# Section headers that indicate symptom content
SECTION_ALLOW = [
    "symptom", "clinical", "presentation", "sign",
    "manifestation", "features", "findings",
    "complaint", "patient", "diagnosis", "physical examination",
    "clinical feature", "case report", "case presentation"
]

# ================= FUNCTIONS =================

def search_pmc():
    """Search PubMed Central for relevant articles"""
    try:
        handle = Entrez.esearch(
            db="pmc",
            term=SEARCH_QUERY,
            retmax=MAX_ARTICLES
        )
        record = Entrez.read(handle)
        return record["IdList"]
    except Exception as e:
        print(f"❌ Search error: {e}")
        return []


def fetch_xml(pmc_id):
    """Fetch full XML for a PMC article"""
    try:
        handle = Entrez.efetch(
            db="pmc",
            id=pmc_id,
            rettype="full",
            retmode="xml"
        )
        return handle.read()
    except Exception as e:
        print(f"⚠️ Fetch error for {pmc_id}: {e}")
        return None


def is_symptom_section(title):
    """Check if section title indicates symptom content"""
    if not title:
        return False
    title = title.lower()
    return any(key in title for key in SECTION_ALLOW)


def extract_symptom_text(xml_data):
    """Extract paragraphs containing symptom descriptions"""
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
        
        # Detect section titles
        if tag.endswith("title") and elem.text:
            current_section_is_symptom = is_symptom_section(elem.text)
        
        # Extract paragraphs only if in symptom section
        if tag.endswith("p") and elem.text:
            text = elem.text.strip().lower()
            
            if current_section_is_symptom:
                # Check if paragraph contains any symptom keywords
                if any(keyword in text for keyword in SYMPTOM_KEYWORDS):
                    # Clean up whitespace
                    clean = re.sub(r'\s+', ' ', elem.text.strip())
                    
                    # Only keep substantial paragraphs
                    if len(clean) > 150:
                        symptom_paragraphs.append(clean)
    
    return symptom_paragraphs


def count_symptoms_in_text(text):
    """Count how many different symptoms appear in a text"""
    text_lower = text.lower()
    found_symptoms = []
    
    for keyword in SYMPTOM_KEYWORDS:
        if keyword in text_lower:
            found_symptoms.append(keyword)
    
    return found_symptoms


# ================= MAIN =================

def main():
    print("\n" + "="*60)
    print("MOUTH CANCER SYMPTOM EXTRACTOR")
    print("="*60 + "\n")
    
    print("🔍 Searching PubMed Central for MOUTH CANCER symptoms...")
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
                    # Count which symptoms appear in this paragraph
                    found = count_symptoms_in_text(s)
                    
                    all_symptoms.append({
                        "pmc_id": pmc_id,
                        "text": s,
                        "symptom_count": len(found),
                        "symptoms_mentioned": found
                    })
            
            time.sleep(0.5)  # Rate limiting
            
        except Exception as e:
            print(f"⚠️ Error {pmc_id}: {e}")
    
    # Save text file
    with open(OUTPUT_TEXT, "w", encoding="utf-8") as f:
        for i, item in enumerate(all_symptoms, 1):
            f.write(f"Context {i} (PMC:{item['pmc_id']}) - {item['symptom_count']} symptoms:\n")
            f.write(f"Symptoms: {', '.join(item['symptoms_mentioned'][:5])}{'...' if len(item['symptoms_mentioned']) > 5 else ''}\n")
            f.write(item["text"] + "\n\n")
    
    # Save JSON
    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(all_symptoms, f, indent=2, ensure_ascii=False)
    
    # Statistics
    total_contexts = len(all_symptoms)
    avg_symptoms = sum(item['symptom_count'] for item in all_symptoms) / total_contexts if total_contexts > 0 else 0
    
    print("\n" + "="*60)
    print("✅ EXTRACTION COMPLETE")
    print("="*60)
    print(f"📝 Text file: {OUTPUT_TEXT}")
    print(f"💾 JSON file: {OUTPUT_JSON}")
    print(f"📊 Total symptom contexts: {total_contexts}")
    print(f"📈 Average symptoms per context: {avg_symptoms:.1f}")
    print("="*60 + "\n")


if __name__ == "__main__":
    main()