"""
Enhanced PubMed Central Medicine Extractor - PROSTATE CANCER
Extracts ALL medicines/drugs from 1000+ prostate cancer articles
Generates detailed TEXT, JSON, and STATISTICS for LLM training
Covers: Localized, Advanced, Metastatic, Hormone-Refractory, Castration-Resistant
"""
from Bio import Entrez
import xml.etree.ElementTree as ET
import os
import re
import time
import json
from collections import Counter

# ================= CONFIG =================
Entrez.email = "nkharel57@gmail.com"  # REQUIRED

SEARCH_QUERY = (
    "prostate cancer treatment drugs OR metastatic prostate cancer therapy OR "
    "castration-resistant prostate cancer drugs OR CRPC treatment OR "
    "androgen deprivation therapy prostate OR abiraterone enzalutamide OR "
    "docetaxel prostate cancer OR prostate cancer immunotherapy OR "
    "hormone-sensitive metastatic prostate cancer treatment OR "
    "PSA progression prostate cancer OR prostate cancer chemotherapy"
)
MAX_ARTICLES = 1000
OUTPUT_TEXT = "prostate_cancer_medicines_detailed.txt"
OUTPUT_JSON = "prostate_cancer_medicines_detailed.json"
OUTPUT_STATS = "prostate_cancer_medicines_statistics.json"
OUTPUT_SUMMARY = "prostate_cancer_medicines_SUMMARY.txt"

# ================= COMPREHENSIVE PROSTATE CANCER DRUG DATABASE =================

PROSTATE_CANCER_DRUG_DATABASE = {
    # ===== HORMONE THERAPY - GnRH AGONISTS =====
    "leuprolide": {
        "names": ["leuprolide", "lupron"],
        "brand": "Lupron, Lupron Depot",
        "generic": "Leuprolide Acetate",
        "category": "Hormone Therapy - GnRH Agonist",
        "cancer_types": ["Locally advanced prostate cancer", "Metastatic prostate cancer", "Recurrent prostate cancer"],
        "intro": "Leuprolide is a GnRH agonist that suppresses testosterone production, fundamental therapy for advanced prostate cancer.",
        "mechanism": "GnRH agonist; suppresses LH and FSH; reduces testosterone to castrate levels",
        "dose": "7.5 mg IM monthly or 22.5 mg every 3 months or 45 mg every 6 months",
        "frequency": "Monthly, quarterly, or 6-monthly depending on formulation",
        "duration": "Continuous indefinitely or until progression",
        "route": "Intramuscular injection",
        "onset": "Testosterone suppression occurs within 2-4 weeks",
        "usage_contexts": [
            "Androgen deprivation therapy (ADT) for locally advanced/metastatic prostate cancer",
            "First-line hormonal therapy for advanced disease",
            "Often combined with antiandrogen (bicalutamide) for maximal androgen blockade (MAB)",
            "Primary therapy in hormone-sensitive metastatic disease",
            "Used with radiation for intermediate/high-risk localized disease",
            "Neoadjuvant therapy before radiation"
        ],
        "indications": ["Advanced prostate cancer", "Metastatic hormone-sensitive prostate cancer", "Locally advanced prostate cancer"],
        "common_side_effects": ["Hot flashes (75%)", "Erectile dysfunction", "Decreased libido", "Fatigue", "Weight gain", "Gynecomastia"],
        "serious_side_effects": ["Testosterone flare (first 1-2 weeks)", "Spinal cord compression risk with flare", "Cardiovascular events", "Osteoporosis with long-term use", "Metabolic syndrome"],
        "monitoring": ["PSA levels every 3-6 months", "Testosterone levels at baseline and 2-4 weeks", "Bone density (DEXA) if long-term ADT", "Cardiovascular risk assessment"],
        "antiandrogen_use": "Often combined with bicalutamide 50 mg daily for first 2 weeks to block testosterone flare",
        "special_notes": "Testosterone flare phenomenon requires antiandrogen cover; castrate testosterone levels essential for efficacy"
    },
    
    "goserelin": {
        "names": ["goserelin", "zoladex"],
        "brand": "Zoladex",
        "generic": "Goserelin Acetate",
        "category": "Hormone Therapy - GnRH Agonist",
        "cancer_types": ["Prostate cancer", "Breast cancer"],
        "intro": "Goserelin is a GnRH agonist similar to leuprolide, providing testosterone suppression for prostate cancer.",
        "mechanism": "GnRH agonist; luteinizing hormone suppression",
        "dose": "3.6 mg subcutaneous every 28 days or 10.8 mg every 12 weeks",
        "frequency": "Monthly or 3-monthly",
        "duration": "Continuous",
        "route": "Subcutaneous implant",
        "usage_contexts": [
            "ADT for advanced prostate cancer",
            "Alternative to leuprolide",
            "Metastatic hormone-sensitive prostate cancer",
            "Locally advanced disease with radiation"
        ],
        "indications": ["Advanced prostate cancer", "Metastatic prostate cancer"]
    },
    
    "triptorelin": {
        "names": ["triptorelin", "trelstar"],
        "brand": "Trelstar",
        "generic": "Triptorelin Pamoate",
        "category": "Hormone Therapy - GnRH Agonist",
        "cancer_types": ["Prostate cancer"],
        "intro": "Triptorelin is a GnRH agonist with efficacy similar to leuprolide and goserelin.",
        "mechanism": "GnRH agonist",
        "dose": "3.75 mg IM every month or 11.25 mg every 3 months or 22.5 mg every 6 months",
        "frequency": "Monthly, quarterly, or 6-monthly",
        "duration": "Continuous",
        "route": "Intramuscular",
        "usage_contexts": [
            "ADT for advanced prostate cancer",
            "Metastatic hormone-sensitive prostate cancer",
            "Alternative GnRH agonist"
        ],
        "indications": ["Advanced prostate cancer"]
    },
    
    # ===== HORMONE THERAPY - ANTIANDROGENS =====
    "bicalutamide": {
        "names": ["bicalutamide", "casodex"],
        "brand": "Casodex",
        "generic": "Bicalutamide",
        "category": "Hormone Therapy - Nonsteroidal Antiandrogen",
        "cancer_types": ["Prostate cancer"],
        "intro": "Bicalutamide is a nonsteroidal antiandrogen blocking androgen receptor, used with GnRH agonists for maximal androgen blockade.",
        "mechanism": "Androgen receptor antagonist; blocks DHT and testosterone binding",
        "dose": "50 mg oral daily",
        "frequency": "Once daily",
        "duration": "Continuous with GnRH agonist or monotherapy",
        "route": "Oral tablet",
        "usage_contexts": [
            "Maximal androgen blockade (MAB): With GnRH agonist (7-14 days before or concurrent)",
            "Monotherapy in early prostate cancer (less common)",
            "Blocks testosterone flare when started before GnRH agonist",
            "Used with metastatic hormone-sensitive prostate cancer",
            "Combined with leuprolide, goserelin, or triptorelin"
        ],
        "indications": ["Advanced prostate cancer with GnRH agonist", "MAB therapy"],
        "common_side_effects": ["Gynecomastia (breast tissue growth)", "Breast pain", "Hot flashes (with GnRH agonist)", "Diarrhea", "Fatigue"],
        "serious_side_effects": ["Hepatotoxicity (rare)", "Interstitial lung disease (rare)"],
        "monitoring": ["Liver function tests baseline and periodic", "PSA levels", "Symptoms of toxicity"]
    },
    
    "flutamide": {
        "names": ["flutamide", "eulexin"],
        "brand": "Eulexin",
        "generic": "Flutamide",
        "category": "Hormone Therapy - Nonsteroidal Antiandrogen",
        "cancer_types": ["Prostate cancer"],
        "intro": "Flutamide is an older antiandrogen less commonly used now, replaced by bicalutamide and enzalutamide.",
        "mechanism": "Androgen receptor antagonist",
        "dose": "250 mg oral three times daily",
        "frequency": "Three times daily",
        "duration": "Continuous with GnRH agonist",
        "route": "Oral capsule",
        "usage_contexts": [
            "MAB: With GnRH agonist (less preferred than bicalutamide)",
            "Historical use; largely replaced by newer agents",
            "Liver toxicity concerns limit use"
        ],
        "indications": ["Advanced prostate cancer (less commonly used now)"],
        "serious_side_effects": ["Hepatotoxicity (more common than bicalutamide)"]
    },
    
    # ===== SECOND-GENERATION ANTIANDROGENS - AR ANTAGONISTS =====
    "enzalutamide": {
        "names": ["enzalutamide", "xtandi"],
        "brand": "Xtandi",
        "generic": "Enzalutamide",
        "category": "Targeted Therapy - Second-Generation AR Antagonist",
        "cancer_types": ["Castration-resistant prostate cancer (CRPC)", "Hormone-sensitive metastatic prostate cancer"],
        "intro": "Enzalutamide is a second-generation androgen receptor antagonist with enhanced binding and improved efficacy in CRPC.",
        "mechanism": "Androgen receptor antagonist; blocks AR binding and nuclear translocation; enhanced binding compared to first-generation",
        "dose": "160 mg oral daily (four 40 mg capsules)",
        "frequency": "Once daily",
        "duration": "Continuous until progression",
        "route": "Oral capsule",
        "usage_contexts": [
            "Castration-resistant prostate cancer (CRPC) - improves OS (AFFIRM trial)",
            "Hormone-sensitive metastatic prostate cancer (HSPC) - superior to ADT alone (ENZAMET)",
            "Metastatic CRPC first-line therapy",
            "Nonmetastatic CRPC to delay progression",
            "Combined with ADT (GnRH agonist/antagonist)"
        ],
        "indications": ["CRPC metastatic", "HSPC metastatic", "Non-metastatic CRPC"],
        "common_side_effects": ["Fatigue", "Diarrhea", "Hot flashes", "Headache", "Hypertension"],
        "serious_side_effects": ["Seizures (0.6%)", "Cardiac arrhythmias", "Ischemic heart disease", "Cerebrovascular events"],
        "monitoring": ["PSA levels every 3-6 months", "Neurologic symptoms (seizure risk)", "Blood pressure", "Cardiac assessment"],
        "special_notes": "Superior outcomes vs bicalutamide in HSPC; seizure risk higher than abiraterone; used in combination with ADT"
    },
    
    "abiraterone": {
        "names": ["abiraterone", "zytiga"],
        "brand": "Zytiga",
        "generic": "Abiraterone Acetate",
        "category": "Targeted Therapy - CYP17 Inhibitor",
        "cancer_types": ["Castration-resistant prostate cancer (CRPC)", "Hormone-sensitive metastatic prostate cancer (HSPC)"],
        "intro": "Abiraterone inhibits CYP17, blocking androgen synthesis at multiple levels, providing additional testosterone suppression in CRPC.",
        "mechanism": "CYP17 inhibitor; blocks androgen synthesis in testes, adrenal glands, and tumor; blocks 17α-hydroxylase and 17,20-lyase",
        "dose": "1000 mg oral daily (four 250 mg tablets)",
        "frequency": "Once daily on empty stomach (fasted 1 hour before, 2 hours after meals)",
        "duration": "Continuous until progression",
        "route": "Oral tablet",
        "requires_prednisone": "Requires prednisone 5 mg daily or dexamethasone to prevent mineralocorticoid excess",
        "usage_contexts": [
            "Metastatic CRPC - improves OS (COU-AA-302 trial)",
            "Hormone-sensitive metastatic prostate cancer (HSPC) - improves OS when added to ADT (LATITUDE/STAMPEDE)",
            "Nonmetastatic CRPC - delays progression to metastatic disease",
            "First-line or second-line depending on prior therapy",
            "Combined with ADT and prednisone"
        ],
        "indications": ["mCRPC", "mHSPC", "nmCRPC", "Biochemical-only recurrence"],
        "common_side_effects": ["Hypertension (31%)", "Hypokalemia", "Fluid retention", "Fatigue", "Diarrhea"],
        "serious_side_effects": ["Severe hypertension", "Hypokalemia with cardiac arrhythmia risk", "Hepatotoxicity", "Cardiac dysfunction"],
        "monitoring": ["Blood pressure weekly × 4 weeks, then monthly", "Serum potassium weekly × 4 weeks, then monthly", "Liver function tests monthly", "PSA levels"],
        "prednisone_essential": "Prednisone mandatory to prevent mineralocorticoid syndrome",
        "special_notes": "Superior to placebo in both CRPC and HSPC; requires prednisone co-administration; strict fasting requirements; careful electrolyte/BP monitoring"
    },
    
    "darolutamide": {
        "names": ["darolutamide", "nubeqa"],
        "brand": "Nubeqa",
        "generic": "Darolutamide",
        "category": "Targeted Therapy - AR Antagonist (Next-Generation)",
        "cancer_types": ["Non-metastatic CRPC", "Metastatic CRPC"],
        "intro": "Darolutamide is a next-generation AR antagonist with improved blood-brain barrier penetration and clinical efficacy in non-metastatic CRPC.",
        "mechanism": "Next-generation AR antagonist; high selectivity for AR; improved CNS penetration",
        "dose": "600 mg oral twice daily (300 mg BID) with food",
        "frequency": "Twice daily",
        "duration": "Continuous until progression",
        "route": "Oral tablet",
        "usage_contexts": [
            "Non-metastatic CRPC (nmCRPC) - improves metastasis-free survival (ARAMIS trial)",
            "Delays progression to metastatic disease",
            "Better CNS penetration may prevent brain metastases",
            "Excellent tolerability profile",
            "With ADT (GnRH agonist/antagonist)"
        ],
        "indications": ["Non-metastatic CRPC", "Metastatic CRPC"],
        "common_side_effects": ["Fatigue", "Rash", "Diarrhea", "Hot flashes", "Nausea"],
        "serious_side_effects": ["Seizures (rare)", "Cardiac toxicity (rare)"],
        "monitoring": ["PSA levels", "Metastasis-free survival assessment", "General safety monitoring"],
        "special_notes": "Excellent tolerability; lower seizure risk than enzalutamide; ARAMIS demonstrated benefit in nmCRPC"
    },
    
    # ===== CHEMOTHERAPY =====
    "docetaxel": {
        "names": ["docetaxel", "taxotere"],
        "brand": "Taxotere",
        "generic": "Docetaxel",
        "category": "Chemotherapy - Taxane",
        "cancer_types": ["Metastatic castration-resistant prostate cancer (mCRPC)", "Hormone-sensitive metastatic prostate cancer (mHSPC)"],
        "intro": "Docetaxel is a taxane chemotherapy showing OS benefit in both HSPC and CRPC prostate cancer.",
        "mechanism": "Taxane; stabilizes microtubules; prevents cell division",
        "dose": "75 mg/m² IV",
        "frequency": "Every 3 weeks",
        "duration": "6-10 cycles (18-30 weeks)",
        "route": "IV infusion",
        "usage_contexts": [
            "mHSPC: Added to ADT (STAMPEDE, CHAARTED trials) improves OS",
            "mCRPC: Improves OS vs mitoxantrone (TAXAN-301 trial)",
            "First-line chemotherapy for symptomatic CRPC",
            "Combined with prednisone/prednisolone",
            "High-volume metastatic disease benefits most"
        ],
        "indications": ["mHSPC", "mCRPC", "Symptomatic CRPC"],
        "common_side_effects": ["Myelosuppression", "Neutropenia", "Alopecia", "Nausea/vomiting", "Fatigue", "Peripheral neuropathy"],
        "serious_side_effects": ["Severe neutropenia with infection risk", "Sepsis", "Severe peripheral neuropathy", "Cardiac toxicity", "Fluid retention"],
        "monitoring": ["Complete blood count before each cycle", "Neurologic exam for neuropathy", "Renal function", "Liver function"],
        "prednisone_use": "Given with prednisone 5-10 mg daily to reduce fluid retention and hypersensitivity",
        "special_notes": "Improves OS in both HSPC and CRPC; benefit greater in high-volume disease; neuropathy cumulative with dose"
    },
    
    "cabazitaxel": {
        "names": ["cabazitaxel", "jevtana"],
        "brand": "Jevtana",
        "generic": "Cabazitaxel",
        "category": "Chemotherapy - Taxane",
        "cancer_types": ["mCRPC"],
        "intro": "Cabazitaxel is a second-generation taxane with activity in docetaxel-resistant mCRPC.",
        "mechanism": "Taxane; active against docetaxel-resistant tumors; less substrate for P-gp efflux pump",
        "dose": "25 mg/m² IV",
        "frequency": "Every 3 weeks",
        "duration": "Multiple cycles until progression",
        "route": "IV infusion",
        "usage_contexts": [
            "mCRPC after docetaxel progression (TROPIC trial) - improves OS",
            "Second-line chemotherapy",
            "Docetaxel-resistant disease",
            "Given with prednisone"
        ],
        "indications": ["mCRPC post-docetaxel", "Docetaxel-resistant disease"],
        "common_side_effects": ["Neutropenia", "Nausea/vomiting", "Diarrhea", "Fatigue"],
        "serious_side_effects": ["Severe neutropenia with sepsis", "Severe diarrhea", "Hepatotoxicity"]
    },
    
    # ===== IMMUNOTHERAPY =====
    "sipuleucel": {
        "names": ["sipuleucel-t", "provenge"],
        "brand": "Provenge",
        "generic": "Sipuleucel-T",
        "category": "Immunotherapy - Therapeutic Cancer Vaccine (Autologous)",
        "cancer_types": ["Asymptomatic or minimally symptomatic mCRPC"],
        "intro": "Sipuleucel-T is an autologous cellular immunotherapy using patient's own dendritic cells activated against PAP antigen.",
        "mechanism": "Autologous dendritic cells stimulated with PAP-GM-CSF fusion protein; induces CD8+ T-cell response against prostate cancer",
        "dose": "One dose of sipuleucel-T every 2 weeks × 3 doses",
        "frequency": "Three infusions over 6 weeks",
        "duration": "Single course",
        "route": "IV infusion",
        "usage_contexts": [
            "Asymptomatic or minimally symptomatic mCRPC",
            "Improves OS (IMPACT trial) by ~4.1 months",
            "No response in PSA or imaging required to benefit",
            "Immune response, not cytotoxic response",
            "Requires leukapheresis for cell collection"
        ],
        "indications": ["Asymptomatic/minimally symptomatic mCRPC"],
        "common_side_effects": ["Fever", "Chills", "Fatigue", "Nausea", "Headache"],
        "serious_side_effects": ["Leukapheresis complications", "Cerebrovascular events (rare)", "Myocardial infarction (rare)"],
        "special_notes": "OS benefit without tumor regression; works in asymptomatic disease; unique cellular therapy approach"
    },
    
    # ===== OTHER AGENTS =====
    "mitoxantrone": {
        "names": ["mitoxantrone", "novantrone"],
        "brand": "Novantrone",
        "generic": "Mitoxantrone",
        "category": "Chemotherapy - Anthracenedione",
        "cancer_types": ["mCRPC"],
        "intro": "Mitoxantrone is older chemotherapy with palliative benefits in symptomatic mCRPC, now largely replaced by docetaxel.",
        "mechanism": "Anthracenedione; topoisomerase II inhibitor",
        "dose": "12-14 mg/m² IV",
        "frequency": "Every 3 weeks",
        "duration": "Multiple cycles",
        "route": "IV infusion",
        "usage_contexts": [
            "Palliative therapy for symptomatic mCRPC (pain relief)",
            "Improves pain vs placebo but not OS",
            "Largely replaced by docetaxel/cabazitaxel",
            "Historical relevance"
        ],
        "indications": ["Symptomatic mCRPC (palliative)"],
        "special_notes": "No OS benefit; primarily for symptom palliation; replaced by superior taxanes"
    },
    
    # ===== RADIOTHERAPY AGENTS =====
    "radium_223": {
        "names": ["radium-223", "ra-223", "alpharadin"],
        "brand": "Xofigo",
        "generic": "Radium Ra 223 Dichloride",
        "category": "Radiotherapy - Alpha-Emitting Radiopharmaceutical",
        "cancer_types": ["mCRPC with bone metastases"],
        "intro": "Radium-223 is an alpha-emitting radiopharmaceutical that targets bone metastases, improving OS in mCRPC.",
        "mechanism": "Alpha-emitter; accumulates in bone metastases; delivers high local dose with minimal systemic toxicity",
        "dose": "55 kBq/kg (1.5 μCi/kg)",
        "frequency": "IV injection once every 4 weeks × 6 injections (24 weeks)",
        "duration": "Six doses over 6 months",
        "route": "Intravenous",
        "usage_contexts": [
            "mCRPC with symptomatic bone metastases (ALSYMPCA trial)",
            "Improves OS by ~3.6 months",
            "Delays skeletal-related events",
            "Can be combined with abiraterone or enzalutamide",
            "Excellent bone pain palliation"
        ],
        "indications": ["mCRPC with bone metastases", "Symptomatic bone disease"],
        "common_side_effects": ["Nausea", "Diarrhea", "Vomiting", "Bone pain flare"],
        "serious_side_effects": ["Myelosuppression", "Severe diarrhea"],
        "monitoring": ["Complete blood count (baseline and before each injection)", "Renal function", "Calcium levels"],
        "special_notes": "Alpha particle causes minimal external radiation exposure; delivers dose to bone with limited soft tissue damage"
    },
    
    # ===== SUPPORTIVE CARE =====
    "bisphosphonate": {
        "names": ["zoledronic acid", "alendronate"],
        "brand": "Zometa (zoledronic acid), Fosamax (alendronate)",
        "generic": "Bisphosphonates",
        "category": "Supportive Care - Bone Protective Agent",
        "cancer_types": ["mCRPC with bone metastases", "Patients on ADT"],
        "intro": "Bisphosphonates reduce skeletal-related events and bone loss in prostate cancer patients.",
        "mechanism": "Inhibit osteoclast function; reduce bone resorption",
        "dose": "Zoledronic acid 4 mg IV every 3-4 weeks (mCRPC) or every 12 weeks (ADT); Alendronate 70 mg weekly (ADT)",
        "frequency": "IV every 3-4 weeks or weekly oral",
        "duration": "Continuous",
        "route": "IV or oral",
        "usage_contexts": [
            "Prevention of skeletal-related events in mCRPC with bone metastases",
            "Bone loss prevention in men receiving ADT",
            "Reduces pathologic fractures",
            "Can be combined with radium-223"
        ],
        "indications": ["mCRPC with bone mets", "ADT-induced bone loss", "High fracture risk"],
        "common_side_effects": ["Acute phase reaction (fever, chills)", "Nausea", "Fatigue"],
        "serious_side_effects": ["Osteonecrosis of jaw", "Atypical fractures"],
        "monitoring": ["Renal function", "Calcium levels", "Dental exams"]
    },
    
    "denosumab": {
        "names": ["denosumab", "xgeva"],
        "brand": "Xgeva",
        "generic": "Denosumab",
        "category": "Supportive Care - RANKL Inhibitor",
        "cancer_types": ["mCRPC with bone metastases"],
        "intro": "Denosumab is a RANKL inhibitor that reduces skeletal-related events in mCRPC with bone metastases.",
        "mechanism": "RANKL inhibitor; prevents osteoclast activation",
        "dose": "120 mg subcutaneous",
        "frequency": "Every 4 weeks",
        "duration": "Continuous",
        "route": "Subcutaneous injection",
        "usage_contexts": [
            "Prevention of skeletal-related events in mCRPC (EMPOWER trial)",
            "Superior to zoledronic acid in some studies",
            "Reduces pathologic fractures and bone pain",
            "No renal dose adjustment needed"
        ],
        "indications": ["mCRPC with bone metastases"],
        "common_side_effects": ["Fatigue", "Nausea"],
        "serious_side_effects": ["Osteonecrosis of jaw", "Hypocalcemia"]
    }
}

# Flatten all drug keywords
ALL_DRUG_KEYWORDS = []
for drug, details in PROSTATE_CANCER_DRUG_DATABASE.items():
    ALL_DRUG_KEYWORDS.extend(details.get("names", []))
    if details.get("brand"):
        brand_words = details.get("brand").lower().replace("+", " ").split()
        ALL_DRUG_KEYWORDS.extend(brand_words)

SECTION_ALLOW = [
    "treatment", "therapy", "drug", "medication", "chemotherapy",
    "results", "management", "protocol", "efficacy", "safety",
    "adverse", "prostate", "psa", "crpc", "castration"
]

# ================= FUNCTIONS =================

def search_pmc():
    """Search PubMed Central"""
    try:
        print("🔍 Connecting to PubMed Central...")
        handle = Entrez.esearch(
            db="pmc",
            term=SEARCH_QUERY,
            retmax=MAX_ARTICLES,
            sort="relevance"
        )
        record = Entrez.read(handle)
        return record["IdList"]
    except Exception as e:
        print(f"❌ Search error: {e}")
        return []

def fetch_xml(pmc_id):
    """Fetch article XML"""
    handle = Entrez.efetch(
        db="pmc",
        id=pmc_id,
        rettype="full",
        retmode="xml"
    )
    return handle.read()

def extract_drug_paragraphs(xml_data):
    """Extract drug-related paragraphs"""
    root = ET.fromstring(xml_data)
    drug_paragraphs = []
    
    try:
        for elem in root.iter():
            tag = elem.tag.lower()
            if tag.endswith("p") and elem.text:
                text = elem.text.strip().lower()
                if any(keyword in text for keyword in ALL_DRUG_KEYWORDS):
                    clean = re.sub(r'\s+', ' ', elem.text.strip())
                    if len(clean) > 100:
                        drug_paragraphs.append(clean)
    except:
        pass
    
    return drug_paragraphs

def identify_drugs(text):
    """Identify drugs in text"""
    text_lower = text.lower()
    found_drugs = []
    
    for drug, details in PROSTATE_CANCER_DRUG_DATABASE.items():
        for name in details.get("names", []):
            if name in text_lower:
                found_drugs.append(drug)
                break
    
    return list(set(found_drugs))

def get_drug_details(drug_name):
    """Get drug details"""
    if drug_name in PROSTATE_CANCER_DRUG_DATABASE:
        return PROSTATE_CANCER_DRUG_DATABASE[drug_name]
    return None

# ================= MAIN EXECUTION =================

def main():
    print("="*150)
    print("PROSTATE CANCER MEDICINES EXTRACTOR - PUBMED CENTRAL")
    print("Processing up to 1000 articles for all prostate cancer drugs")
    print("="*150 + "\n")
    
    pmc_ids = search_pmc()
    print(f"✅ Found {len(pmc_ids)} relevant articles\n")
    
    if not pmc_ids:
        print("❌ No articles found.")
        return
    
    all_contexts = []
    processed = 0
    errors = 0
    
    # Process articles
    for i, pmc_id in enumerate(pmc_ids[:MAX_ARTICLES], 1):
        if i % 50 == 0:
            print(f"📊 Progress: {i}/{min(len(pmc_ids), MAX_ARTICLES)} articles...")
        
        try:
            xml = fetch_xml(pmc_id)
            paragraphs = extract_drug_paragraphs(xml)
            
            for para in paragraphs:
                drugs = identify_drugs(para)
                
                all_contexts.append({
                    "pmc_id": pmc_id,
                    "article_url": f"https://pubmed.ncbi.nlm.nih.gov/pmc/{pmc_id}",
                    "paragraph": para,
                    "drugs_found": drugs,
                    "drug_count": len(drugs),
                    "drug_details": [get_drug_details(d) for d in drugs if get_drug_details(d)]
                })
            
            processed += 1
            time.sleep(0.2)
        except:
            errors += 1
            continue
    
    print(f"\n✅ Processed: {processed} articles")
    print(f"❌ Errors: {errors}")
    print(f"📝 Drug contexts: {len(all_contexts)}\n")
    
    # Save TEXT
    print(f"💾 Saving TEXT file...")
    with open(OUTPUT_TEXT, "w", encoding="utf-8") as f:
        f.write("=" * 150 + "\n")
        f.write("PROSTATE CANCER MEDICINES - COMPREHENSIVE EXTRACTION\n")
        f.write(f"Contexts: {len(all_contexts)} | Articles: {processed}\n")
        f.write("=" * 150 + "\n\n")
        
        for i, ctx in enumerate(all_contexts, 1):
            f.write(f"\n{'='*150}\nCONTEXT #{i}\n{'='*150}\n")
            f.write(f"PMC: {ctx['pmc_id']}\nDrugs: {len(ctx['drugs_found'])}\n\n")
            
            for drug in ctx['drug_details']:
                if drug:
                    f.write(f"\n{'─'*150}\n💊 {drug['names'][0].upper()} ({drug['brand']})\n{'─'*150}\n")
                    f.write(f"Category: {drug['category']}\n")
                    f.write(f"Generic: {drug['generic']}\n")
                    f.write(f"Cancers: {', '.join(drug.get('cancer_types', ['N/A']))}\n\n")
                    f.write(f"📝 INTRO:\n{drug['intro']}\n\n")
                    f.write(f"🔬 MECHANISM:\n{drug['mechanism']}\n\n")
                    f.write(f"💉 DOSING:\n  Dose: {drug.get('dose', 'N/A')}\n")
                    f.write(f"  Frequency: {drug.get('frequency', 'N/A')}\n")
                    f.write(f"  Duration: {drug.get('duration', 'N/A')}\n")
                    f.write(f"  Route: {drug.get('route', 'N/A')}\n\n")
                    
                    if drug.get('usage_contexts'):
                        f.write(f"🎯 USAGE:\n")
                        for usage in drug['usage_contexts'][:4]:
                            f.write(f"  • {usage}\n")
                        f.write("\n")
                    
                    if drug.get('indications'):
                        f.write(f"✅ INDICATIONS:\n")
                        for ind in drug['indications']:
                            f.write(f"  • {ind}\n")
                        f.write("\n")
                    
                    if drug.get('common_side_effects'):
                        f.write(f"⚠️  SIDE EFFECTS:\n")
                        for se in drug['common_side_effects'][:5]:
                            f.write(f"  • {se}\n")
                        f.write("\n")
            
            f.write(f"\n📌 PARAGRAPH:\n{ctx['paragraph']}\n")
    
    # Save JSON
    print(f"💾 Saving JSON file...")
    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(all_contexts, f, indent=2, ensure_ascii=False)
    
    # Generate stats
    all_drugs = []
    for ctx in all_contexts:
        all_drugs.extend(ctx['drugs_found'])
    
    drug_freq = Counter(all_drugs)
    
    stats = {
        "date": str(time.ctime()),
        "articles_processed": processed,
        "errors": errors,
        "contexts": len(all_contexts),
        "unique_drugs": len(drug_freq),
        "top_drugs": dict(drug_freq.most_common(20)),
        "database_drugs": len(PROSTATE_CANCER_DRUG_DATABASE)
    }
    
    with open(OUTPUT_STATS, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
    
    # Summary report
    print(f"💾 Saving SUMMARY...")
    with open(OUTPUT_SUMMARY, "w", encoding="utf-8") as f:
        f.write("=" * 150 + "\n")
        f.write("PROSTATE CANCER MEDICINES - SUMMARY REPORT\n")
        f.write("=" * 150 + "\n\n")
        f.write(f"Articles: {processed} | Contexts: {len(all_contexts)} | Unique Drugs: {len(drug_freq)}\n\n")
        
        f.write(f"🏆 TOP 20 DRUGS:\n")
        for rank, (drug, count) in enumerate(drug_freq.most_common(20), 1):
            details = get_drug_details(drug)
            cat = details.get('category', 'Unknown') if details else 'Unknown'
            f.write(f"{rank:2d}. {drug:25s} - {count:4d} mentions - {cat}\n")
        f.write("\n")
        
        f.write(f"💊 CATEGORIES:\n")
        categories = {}
        for drug, details in PROSTATE_CANCER_DRUG_DATABASE.items():
            cat = details.get('category', 'Unknown')
            if cat not in categories:
                categories[cat] = 0
            categories[cat] += 1
        for cat in sorted(categories.keys()):
            f.write(f"  • {cat}: {categories[cat]} drugs\n")
        
        f.write(f"\n🔬 KEY TRIALS:\n")
        f.write(f"  • AFFIRM: Enzalutamide in mCRPC improves OS\n")
        f.write(f"  • STAMPEDE: Docetaxel + ADT in mHSPC improves OS\n")
        f.write(f"  • CHAARTED: Docetaxel + ADT in mHSPC (high volume benefit)\n")
        f.write(f"  • LATITUDE: Abiraterone + ADT in mHSPC improves OS\n")
        f.write(f"  • COU-AA-302: Abiraterone in mCRPC improves OS\n")
        f.write(f"  • ARAMIS: Darolutamide in nmCRPC improves MFS\n")
        f.write(f"  • ALSYMPCA: Radium-223 in mCRPC with bone mets improves OS\n")
        f.write(f"  • IMPACT: Sipuleucel-T in mCRPC improves OS\n")
    
    print("\n" + "🎉 " + "="*148)
    print("EXTRACTION COMPLETE - PROSTATE CANCER MEDICINES")
    print("=" * 150)
    print(f"\n📁 FILES: {OUTPUT_TEXT}, {OUTPUT_JSON}, {OUTPUT_STATS}, {OUTPUT_SUMMARY}")
    print(f"📊 Contexts: {len(all_contexts)} | Unique Drugs: {len(drug_freq)}")
    print(f"\n🏆 TOP 10:")
    for rank, (drug, count) in enumerate(drug_freq.most_common(10), 1):
        print(f"   {rank}. {drug}: {count}")
    print("\n" + "=" * 150)

if __name__ == "__main__":
    main()