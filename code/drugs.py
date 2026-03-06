"""
Enhanced PubMed Central Medicine Extractor - LIVER CANCER
Extracts ALL medicines/drugs from 1000+ liver cancer articles
Generates detailed TEXT, JSON, and STATISTICS for LLM training
Covers: HCC, Cholangiocarcinoma, Hepatoblastoma, Metastatic Liver Cancer,
        TACE, Ablation, Systemic Therapy, Immunotherapy
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
    "hepatocellular carcinoma treatment drugs OR liver cancer therapy OR "
    "HCC systemic therapy OR sorafenib liver cancer OR lenvatinib HCC OR "
    "regorafenib hepatocellular carcinoma OR atezolizumab bevacizumab HCC OR "
    "TACE transarterial chemoembolization OR liver cancer immunotherapy OR "
    "cholangiocarcinoma treatment OR hepatoblastoma chemotherapy OR "
    "nivolumab liver cancer OR cabozantinib HCC OR ramucirumab HCC"
)
MAX_ARTICLES = 1000
OUTPUT_TEXT    = "liver_cancer_medicines_detailed.txt"
OUTPUT_JSON    = "liver_cancer_medicines_detailed.json"
OUTPUT_STATS   = "liver_cancer_medicines_statistics.json"
OUTPUT_SUMMARY = "liver_cancer_medicines_SUMMARY.txt"

# ================= LIVER CANCER DRUG DATABASE =================

LIVER_CANCER_DRUG_DATABASE = {

    # ===== FIRST-LINE SYSTEMIC THERAPY =====
    "sorafenib": {
        "names": ["sorafenib", "nexavar"],
        "brand": "Nexavar",
        "generic": "Sorafenib",
        "category": "Targeted Therapy - Multi-Kinase Inhibitor (First-Line)",
        "cancer_types": ["Advanced HCC", "Unresectable HCC"],
        "intro": "Sorafenib was the first systemic therapy to improve OS in advanced HCC, remaining a standard first-line option for over a decade.",
        "mechanism": "Multi-kinase inhibitor; targets RAF/MEK/ERK pathway, VEGFR-2/3, PDGFR-β; anti-angiogenic and anti-proliferative",
        "dose": "400 mg oral twice daily",
        "frequency": "Twice daily (continuous)",
        "duration": "Until progression or intolerance",
        "route": "Oral tablet",
        "usage_contexts": [
            "First-line therapy for advanced/unresectable HCC (SHARP trial) — improved OS",
            "Child-Pugh A liver function preferred",
            "Previously the global standard of care for >10 years",
            "Now often replaced by atezolizumab + bevacizumab in eligible patients",
            "Alternative first-line when immunotherapy contraindicated",
            "Used in Barcelona Clinic Liver Cancer (BCLC) stage C disease"
        ],
        "indications": ["Advanced HCC", "Unresectable HCC", "BCLC stage C"],
        "common_side_effects": ["Hand-foot skin reaction (HFSR/HFSRD)", "Diarrhea", "Fatigue", "Hypertension", "Alopecia", "Anorexia"],
        "serious_side_effects": ["Severe HFSR", "Hepatotoxicity", "Bleeding/hemorrhage", "Cardiac ischemia", "QT prolongation"],
        "monitoring": ["Blood pressure weekly × 6 weeks", "Liver function tests monthly", "TSH every 3 months", "ECG if cardiac risk"],
        "special_notes": "SHARP trial showed 2.8-month OS improvement vs placebo; HFSR is dose-limiting toxicity; dose reduction to 400mg once daily for toxicity"
    },

    "lenvatinib": {
        "names": ["lenvatinib", "lenvima"],
        "brand": "Lenvima",
        "generic": "Lenvatinib",
        "category": "Targeted Therapy - Multi-Kinase Inhibitor (First-Line)",
        "cancer_types": ["Advanced HCC", "Unresectable HCC"],
        "intro": "Lenvatinib is non-inferior to sorafenib in first-line advanced HCC with potentially higher response rates.",
        "mechanism": "VEGFR1-3, FGFR1-4, PDGFRα, RET, KIT inhibitor; potent anti-angiogenic",
        "dose": "12 mg oral daily (≥60 kg) or 8 mg daily (<60 kg)",
        "frequency": "Once daily",
        "duration": "Until progression or intolerance",
        "route": "Oral capsule",
        "usage_contexts": [
            "First-line HCC — non-inferior to sorafenib (REFLECT trial)",
            "Higher objective response rate vs sorafenib (~24% vs 9%)",
            "Weight-based dosing important",
            "Excluded: main portal vein invasion, >50% liver involvement, bile duct invasion",
            "Alternative to sorafenib or atezolizumab + bevacizumab",
            "Child-Pugh A liver function"
        ],
        "indications": ["Unresectable HCC", "First-line advanced HCC"],
        "common_side_effects": ["Hypertension (42%)", "Fatigue", "Diarrhea", "Decreased appetite", "Weight loss", "Arthralgia"],
        "serious_side_effects": ["Severe hypertension", "Hepatic failure", "Arterial thromboembolic events", "Hemorrhage", "Proteinuria"],
        "monitoring": ["Blood pressure at least weekly", "Liver function tests", "Urinalysis for proteinuria", "Thyroid function"],
        "special_notes": "REFLECT trial: non-inferior OS to sorafenib; better ORR; weight-based dosing critical; excluded portal vein involvement"
    },

    "atezolizumab_bevacizumab": {
        "names": ["atezolizumab", "bevacizumab", "tecentriq", "avastin"],
        "brand": "Tecentriq (atezolizumab) + Avastin (bevacizumab)",
        "generic": "Atezolizumab + Bevacizumab",
        "category": "Immunotherapy + Anti-Angiogenic (First-Line Preferred)",
        "cancer_types": ["Advanced HCC", "Unresectable HCC"],
        "intro": "Atezolizumab + bevacizumab is now the preferred first-line regimen for advanced HCC, showing superior OS and PFS vs sorafenib.",
        "mechanism": "Atezolizumab: PD-L1 inhibitor (checkpoint inhibitor); Bevacizumab: VEGF-A inhibitor (anti-angiogenic); synergistic immunomodulatory effect",
        "dose": "Atezolizumab 1200 mg IV + Bevacizumab 15 mg/kg IV",
        "frequency": "Every 3 weeks",
        "duration": "Until progression or unacceptable toxicity",
        "route": "IV infusion",
        "usage_contexts": [
            "Preferred first-line therapy for advanced HCC (IMbrave150 trial)",
            "Superior OS and PFS vs sorafenib",
            "Requires upper GI endoscopy before starting (bleeding risk with bevacizumab)",
            "Child-Pugh A liver function",
            "ECOG PS 0-1",
            "Contraindicated with autoimmune disease or recent varices bleeding"
        ],
        "indications": ["Advanced/unresectable HCC", "First-line preferred therapy"],
        "common_side_effects": ["Fatigue", "Proteinuria", "Hypertension", "Rash", "Immune-related adverse events"],
        "serious_side_effects": ["GI bleeding (esophageal varices)", "Immune-related pneumonitis", "Hepatitis", "Arterial thrombosis", "GI perforation"],
        "monitoring": ["Endoscopy before initiation", "Blood pressure", "Urinalysis", "Liver function", "Thyroid function", "Immune-related toxicity screening"],
        "special_notes": "IMbrave150: OS 19.2 vs 13.4 months vs sorafenib; mandatory GI endoscopy before treatment; bevacizumab contraindicated with untreated varices"
    },

    "durvalumab_tremelimumab": {
        "names": ["durvalumab", "tremelimumab", "imfinzi"],
        "brand": "Imfinzi (durvalumab) + Tremelimumab",
        "generic": "Durvalumab + Tremelimumab",
        "category": "Immunotherapy - PD-L1 + CTLA-4 Inhibitor (First-Line)",
        "cancer_types": ["Advanced HCC"],
        "intro": "Durvalumab + tremelimumab (STRIDE regimen) is a dual checkpoint inhibitor combination approved for first-line advanced HCC.",
        "mechanism": "Durvalumab: PD-L1 inhibitor; Tremelimumab: CTLA-4 inhibitor; dual immune checkpoint blockade",
        "dose": "Tremelimumab 300 mg single dose (Cycle 1) + Durvalumab 1500 mg IV every 4 weeks",
        "frequency": "Tremelimumab single priming dose; Durvalumab every 4 weeks",
        "duration": "Until progression",
        "route": "IV infusion",
        "usage_contexts": [
            "First-line advanced HCC (HIMALAYA trial)",
            "Non-inferior OS vs sorafenib",
            "Alternative to atezo+bev when bevacizumab contraindicated",
            "Single priming dose of tremelimumab (STRIDE regimen)",
            "Child-Pugh A liver function"
        ],
        "indications": ["Advanced HCC first-line"],
        "common_side_effects": ["Fatigue", "Rash", "Diarrhea", "Immune-related AEs"],
        "serious_side_effects": ["Immune-mediated hepatitis", "Colitis", "Pneumonitis", "Endocrinopathies"],
        "monitoring": ["Liver function", "Thyroid", "Immune toxicity monitoring"],
        "special_notes": "HIMALAYA trial; STRIDE regimen — single tremelimumab priming dose; good alternative when bevacizumab contraindicated"
    },

    # ===== SECOND-LINE SYSTEMIC THERAPY =====
    "regorafenib": {
        "names": ["regorafenib", "stivarga"],
        "brand": "Stivarga",
        "generic": "Regorafenib",
        "category": "Targeted Therapy - Multi-Kinase Inhibitor (Second-Line)",
        "cancer_types": ["Advanced HCC post-sorafenib"],
        "intro": "Regorafenib is approved second-line HCC therapy after sorafenib progression in patients who tolerated sorafenib.",
        "mechanism": "Multi-kinase inhibitor; targets VEGFR1-3, TIE2, FGFR, PDGFR, RET, RAF, KIT; broader than sorafenib",
        "dose": "160 mg oral daily",
        "frequency": "Once daily for 21 days on / 7 days off (28-day cycle)",
        "duration": "Until progression or intolerance",
        "route": "Oral tablet (take with low-fat meal)",
        "usage_contexts": [
            "Second-line HCC after sorafenib progression (RESORCE trial)",
            "Must have tolerated sorafenib ≥400 mg/day for ≥20 of last 28 days",
            "Child-Pugh A liver function",
            "Improves OS vs placebo (10.6 vs 7.8 months)",
            "Shown effective even with sorafenib resistance"
        ],
        "indications": ["Advanced HCC post-sorafenib"],
        "common_side_effects": ["HFSR", "Fatigue", "Hypertension", "Diarrhea", "Rash"],
        "serious_side_effects": ["Severe HFSR", "Hepatotoxicity", "Hemorrhage", "Cardiac ischemia"],
        "monitoring": ["Liver function", "Blood pressure", "HFSR assessment"],
        "special_notes": "RESORCE trial; only for sorafenib-tolerant patients; 21-days on/7-days off schedule; take with low-fat meal"
    },

    "cabozantinib": {
        "names": ["cabozantinib", "cabometyx"],
        "brand": "Cabometyx",
        "generic": "Cabozantinib",
        "category": "Targeted Therapy - Multi-Kinase Inhibitor (Second-Line)",
        "cancer_types": ["Advanced HCC post-sorafenib"],
        "intro": "Cabozantinib is approved for second/third-line HCC, targeting VEGFR, MET, and AXL pathways.",
        "mechanism": "VEGFR1-3, MET, AXL, RET inhibitor; MET inhibition particularly relevant in HCC resistance",
        "dose": "60 mg oral daily",
        "frequency": "Once daily (fasted or low-fat meal)",
        "duration": "Until progression",
        "route": "Oral tablet",
        "usage_contexts": [
            "Second or third-line HCC after sorafenib (CELESTIAL trial)",
            "Improves OS vs placebo",
            "Can use after sorafenib AND one other systemic therapy",
            "MET overexpression in HCC may predict benefit",
            "Child-Pugh A"
        ],
        "indications": ["Advanced HCC after sorafenib", "2nd/3rd line HCC"],
        "common_side_effects": ["HFSR", "Hypertension", "Fatigue", "Diarrhea", "Decreased appetite"],
        "serious_side_effects": ["GI perforation/fistula", "Hemorrhage", "Thromboembolic events", "Hepatotoxicity"],
        "monitoring": ["Blood pressure", "Liver function", "HFSR"],
        "special_notes": "CELESTIAL trial; broadest kinase coverage including MET and AXL; fasting or low-fat meal requirement"
    },

    "ramucirumab": {
        "names": ["ramucirumab", "cyramza"],
        "brand": "Cyramza",
        "generic": "Ramucirumab",
        "category": "Targeted Therapy - VEGFR-2 Antibody (Second-Line, AFP-selected)",
        "cancer_types": ["Advanced HCC post-sorafenib with AFP ≥400 ng/mL"],
        "intro": "Ramucirumab is approved specifically for HCC patients with AFP ≥400 ng/mL after sorafenib — the first biomarker-selected HCC therapy.",
        "mechanism": "VEGFR-2 monoclonal antibody; blocks VEGF ligand binding; anti-angiogenic",
        "dose": "8 mg/kg IV",
        "frequency": "Every 2 weeks",
        "duration": "Until progression",
        "route": "IV infusion",
        "usage_contexts": [
            "Second-line HCC with AFP ≥400 ng/mL after sorafenib (REACH-2 trial)",
            "First biomarker-selected therapy in HCC",
            "AFP ≥400 ng/mL is mandatory selection criterion",
            "Child-Pugh A liver function",
            "Improves OS in AFP-high population"
        ],
        "indications": ["Advanced HCC post-sorafenib with AFP ≥400 ng/mL"],
        "common_side_effects": ["Hypertension", "Fatigue", "Peripheral edema", "Decreased appetite"],
        "serious_side_effects": ["Hemorrhage", "Arterial thrombosis", "GI perforation", "Severe hypertension"],
        "monitoring": ["AFP levels", "Blood pressure", "Liver function"],
        "special_notes": "REACH-2 trial; AFP ≥400 ng/mL mandatory; first biomarker-driven HCC therapy; IV every 2 weeks"
    },

    "nivolumab": {
        "names": ["nivolumab", "opdivo"],
        "brand": "Opdivo",
        "generic": "Nivolumab",
        "category": "Immunotherapy - PD-1 Inhibitor (Second-Line)",
        "cancer_types": ["Advanced HCC post-sorafenib"],
        "intro": "Nivolumab is a PD-1 checkpoint inhibitor with durable responses in HCC post-sorafenib.",
        "mechanism": "PD-1 monoclonal antibody; restores T-cell anti-tumor immunity",
        "dose": "240 mg IV every 2 weeks or 480 mg IV every 4 weeks",
        "frequency": "Every 2 or 4 weeks",
        "duration": "Until progression or unacceptable toxicity",
        "route": "IV infusion",
        "usage_contexts": [
            "Second-line HCC after sorafenib (CheckMate 040 trial)",
            "Accelerated approval based on response rate",
            "Durable responses in ~15-20% patients",
            "Combined with ipilimumab for enhanced response",
            "Child-Pugh A or B7"
        ],
        "indications": ["Advanced HCC post-sorafenib"],
        "common_side_effects": ["Fatigue", "Rash", "Pruritus", "Diarrhea", "Immune-related AEs"],
        "serious_side_effects": ["Immune-mediated hepatitis", "Pneumonitis", "Colitis", "Endocrinopathies", "Nephritis"],
        "monitoring": ["Liver function", "Thyroid function", "Immune toxicity panel"],
        "special_notes": "CheckMate 040; accelerated FDA approval; durable responses; combination with ipilimumab studied"
    },

    "nivolumab_ipilimumab": {
        "names": ["ipilimumab", "yervoy"],
        "brand": "Opdivo (nivolumab) + Yervoy (ipilimumab)",
        "generic": "Nivolumab + Ipilimumab",
        "category": "Immunotherapy - PD-1 + CTLA-4 Inhibitor Combination",
        "cancer_types": ["Advanced HCC post-sorafenib"],
        "intro": "Nivolumab + ipilimumab dual checkpoint blockade shows higher response rates and durable benefit in HCC.",
        "mechanism": "Nivolumab: PD-1 inhibitor; Ipilimumab: CTLA-4 inhibitor; dual checkpoint blockade enhances T-cell activation",
        "dose": "Nivolumab 1 mg/kg + Ipilimumab 3 mg/kg IV every 3 weeks × 4 doses → Nivolumab 240 mg every 2 weeks",
        "frequency": "Combination phase q3w × 4 then nivolumab maintenance",
        "duration": "Combination 4 cycles then maintenance",
        "route": "IV infusion",
        "usage_contexts": [
            "Second-line HCC after sorafenib (CheckMate 040 cohort 4)",
            "ORR ~32% — higher than nivolumab alone",
            "Durable responses in subset",
            "Child-Pugh A or B7",
            "Higher immune toxicity than monotherapy"
        ],
        "indications": ["Advanced HCC post-sorafenib"],
        "serious_side_effects": ["Higher immune toxicity than monotherapy", "Hepatitis", "Colitis", "Endocrinopathies"],
        "special_notes": "Higher ORR than nivolumab alone; more immune toxicity; CheckMate 040 cohort 4"
    },

    "pembrolizumab": {
        "names": ["pembrolizumab", "keytruda"],
        "brand": "Keytruda",
        "generic": "Pembrolizumab",
        "category": "Immunotherapy - PD-1 Inhibitor (Second-Line)",
        "cancer_types": ["Advanced HCC post-sorafenib"],
        "intro": "Pembrolizumab is a PD-1 inhibitor with established activity in HCC after sorafenib.",
        "mechanism": "PD-1 monoclonal antibody; restores T-cell mediated anti-tumor immunity",
        "dose": "200 mg IV every 3 weeks or 400 mg every 6 weeks",
        "frequency": "Every 3 or 6 weeks",
        "duration": "Until progression (max 2 years)",
        "route": "IV infusion",
        "usage_contexts": [
            "Second-line HCC after sorafenib (KEYNOTE-224, KEYNOTE-394)",
            "Accelerated approval in US",
            "ORR ~17% with durable responses",
            "Child-Pugh A",
            "Being studied in combination regimens"
        ],
        "indications": ["Advanced HCC post-sorafenib"],
        "common_side_effects": ["Fatigue", "Rash", "Pruritus", "Immune-related AEs"],
        "serious_side_effects": ["Immune-mediated hepatitis", "Pneumonitis", "Colitis", "Endocrinopathies"],
        "monitoring": ["Liver function", "Thyroid function", "Immune toxicity"],
        "special_notes": "KEYNOTE-394 confirmed benefit in Asian HCC population; max 35 cycles (2 years)"
    },

    # ===== LOCOREGIONAL / INTERVENTIONAL =====
    "doxorubicin_tace": {
        "names": ["doxorubicin", "adriamycin", "lipiodol", "dc-bead", "drug eluting bead"],
        "brand": "Adriamycin (doxorubicin), DC Bead (DEB-TACE)",
        "generic": "Doxorubicin / Drug-Eluting Beads",
        "category": "Locoregional Therapy - TACE (Transarterial Chemoembolization)",
        "cancer_types": ["Intermediate stage HCC (BCLC B)", "Unresectable HCC"],
        "intro": "TACE delivers chemotherapy (doxorubicin) directly into hepatic artery feeding the tumor, combined with embolization, for intermediate HCC.",
        "mechanism": "Intra-arterial drug delivery + ischemia from embolization; high local drug concentration with hepatic arterial embolization",
        "dose": "Doxorubicin 50-75 mg in lipiodol (conventional TACE) or loaded into DC beads (DEB-TACE)",
        "frequency": "On demand (every 6-8 weeks until no viable tumor)",
        "duration": "Repeated until complete response or progression",
        "route": "Intra-arterial (hepatic artery catheterization)",
        "usage_contexts": [
            "Standard of care for BCLC B (intermediate stage) HCC",
            "Unresectable multinodular HCC with preserved liver function",
            "Bridge therapy to liver transplant",
            "Downstaging to curative intent",
            "Can be combined with systemic therapy (atezolizumab + bevacizumab after TACE)",
            "Child-Pugh A or B7/8"
        ],
        "indications": ["Intermediate HCC (BCLC B)", "Unresectable multinodular HCC", "Bridge to transplant"],
        "common_side_effects": ["Post-embolization syndrome (fever, pain, nausea)", "Elevated liver enzymes", "Fatigue"],
        "serious_side_effects": ["Hepatic failure", "Biloma", "Hepatic abscess", "Non-target embolization"],
        "monitoring": ["Liver function before each session", "AFP response", "MRI/CT at 4-6 weeks post-TACE"],
        "special_notes": "Standard care for BCLC B; DEB-TACE reduces systemic doxorubicin exposure; response assessed by mRECIST"
    },

    "cisplatin_tace": {
        "names": ["cisplatin", "platinol"],
        "brand": "Platinol",
        "generic": "Cisplatin",
        "category": "Locoregional Therapy - TACE agent (Alternative)",
        "cancer_types": ["HCC (via TACE)"],
        "intro": "Cisplatin is used as an alternative to doxorubicin in TACE for HCC in some centers.",
        "mechanism": "Platinum-based alkylating agent; DNA cross-linking; combined with hepatic artery embolization",
        "dose": "Varies: 10-100 mg in lipiodol emulsion",
        "frequency": "Per TACE session",
        "duration": "Repeated sessions",
        "route": "Intra-arterial",
        "usage_contexts": [
            "TACE alternative to doxorubicin in some Asian centers",
            "Used in combination TACE protocols",
            "Comparable efficacy to doxorubicin in some trials"
        ],
        "indications": ["HCC TACE alternative agent"]
    },

    # ===== CHOLANGIOCARCINOMA THERAPY =====
    "gemcitabine_cisplatin": {
        "names": ["gemcitabine", "gemzar", "cisplatin"],
        "brand": "Gemzar (gemcitabine) + Platinol (cisplatin)",
        "generic": "Gemcitabine + Cisplatin",
        "category": "Chemotherapy - Standard First-Line (Biliary Tract / Cholangiocarcinoma)",
        "cancer_types": ["Cholangiocarcinoma (CCA)", "Biliary tract cancer (BTC)", "Gallbladder cancer"],
        "intro": "Gemcitabine + cisplatin is the standard first-line chemotherapy for advanced cholangiocarcinoma and biliary tract cancers.",
        "mechanism": "Gemcitabine: nucleoside analogue (DNA synthesis inhibitor); Cisplatin: platinum-based alkylating agent",
        "dose": "Gemcitabine 1000 mg/m² + Cisplatin 25 mg/m² IV",
        "frequency": "Days 1 and 8 of each 21-day cycle",
        "duration": "8 cycles (6 months) or until progression",
        "route": "IV infusion",
        "usage_contexts": [
            "First-line advanced/metastatic cholangiocarcinoma (ABC-02 trial)",
            "Improves OS vs gemcitabine alone",
            "Standard backbone for biliary tract cancers",
            "Combined with durvalumab in TOPAZ-1 trial (improved OS)",
            "Intrahepatic, extrahepatic, and perihilar CCA",
            "Gallbladder cancer first-line"
        ],
        "indications": ["Advanced cholangiocarcinoma", "Biliary tract cancer", "Gallbladder cancer"],
        "common_side_effects": ["Nausea/vomiting", "Myelosuppression", "Fatigue", "Peripheral neuropathy (cisplatin)", "Alopecia"],
        "serious_side_effects": ["Severe neutropenia", "Cisplatin nephrotoxicity", "Ototoxicity", "Peripheral neuropathy"],
        "monitoring": ["CBC before each cycle", "Renal function (cisplatin)", "Audiometry (cisplatin)", "Liver function"],
        "special_notes": "ABC-02 trial established standard; TOPAZ-1 added durvalumab improving OS; hydration mandatory with cisplatin"
    },

    "durvalumab_gem_cis": {
        "names": ["durvalumab", "imfinzi"],
        "brand": "Imfinzi",
        "generic": "Durvalumab",
        "category": "Immunotherapy + Chemotherapy (Cholangiocarcinoma First-Line)",
        "cancer_types": ["Biliary tract cancer", "Cholangiocarcinoma"],
        "intro": "Durvalumab added to gemcitabine + cisplatin is now a standard first-line option for advanced biliary tract cancers (TOPAZ-1 trial).",
        "mechanism": "PD-L1 inhibitor; enhances immune response against tumor when combined with chemotherapy",
        "dose": "Durvalumab 1500 mg IV Day 1 + Gemcitabine 1000 mg/m² + Cisplatin 25 mg/m² IV Days 1 and 8",
        "frequency": "Every 3 weeks (chemotherapy 8 cycles; durvalumab continues)",
        "duration": "Chemotherapy 8 cycles; durvalumab until progression",
        "route": "IV infusion",
        "usage_contexts": [
            "First-line advanced biliary tract cancer (TOPAZ-1 trial)",
            "Improved OS vs gem-cis alone",
            "Intrahepatic, extrahepatic CCA and gallbladder cancer",
            "PD-L1 expression not required for benefit",
            "Now a preferred first-line regimen"
        ],
        "indications": ["Advanced biliary tract cancer first-line"],
        "special_notes": "TOPAZ-1: OS benefit; durvalumab continues beyond chemotherapy completion"
    },

    "pemigatinib": {
        "names": ["pemigatinib", "pemazyre"],
        "brand": "Pemazyre",
        "generic": "Pemigatinib",
        "category": "Targeted Therapy - FGFR Inhibitor (Cholangiocarcinoma)",
        "cancer_types": ["Cholangiocarcinoma with FGFR2 fusion/rearrangement"],
        "intro": "Pemigatinib is an FGFR inhibitor specifically approved for cholangiocarcinoma with FGFR2 fusions or rearrangements.",
        "mechanism": "FGFR1/2/3 inhibitor; targets FGFR2 fusions common in intrahepatic CCA",
        "dose": "13.5 mg oral daily",
        "frequency": "Once daily for 14 days on / 7 days off (21-day cycle)",
        "duration": "Until progression",
        "route": "Oral tablet",
        "usage_contexts": [
            "Previously treated CCA with FGFR2 fusion/rearrangement (FIGHT-202 trial)",
            "FGFR2 testing mandatory before use",
            "ORR ~36% in FGFR2 fusion positive patients",
            "Durable responses seen",
            "Intrahepatic CCA most common FGFR2 fusion site"
        ],
        "indications": ["CCA with FGFR2 fusion/rearrangement (post-chemotherapy)"],
        "common_side_effects": ["Hyperphosphatemia", "Alopecia", "Diarrhea", "Dry eye/mouth", "Fatigue"],
        "serious_side_effects": ["Hyperphosphatemia (requires dietary restriction)", "Retinal detachment", "Severe skin toxicity"],
        "monitoring": ["Phosphate levels (weekly × 8 weeks)", "Ophthalmologic exams", "Liver function"],
        "special_notes": "FIGHT-202 trial; FGFR2 testing via NGS required; hyperphosphatemia managed with dietary restriction and phosphate binders; 14/7 schedule"
    },

    "ivosidenib": {
        "names": ["ivosidenib", "tibsovo"],
        "brand": "Tibsovo",
        "generic": "Ivosidenib",
        "category": "Targeted Therapy - IDH1 Inhibitor (Cholangiocarcinoma)",
        "cancer_types": ["Cholangiocarcinoma with IDH1 mutation"],
        "intro": "Ivosidenib is an IDH1 inhibitor approved for previously treated CCA with IDH1 mutations.",
        "mechanism": "Mutant IDH1 inhibitor; reduces 2-hydroxyglutarate oncometabolite; restores normal differentiation",
        "dose": "500 mg oral daily",
        "frequency": "Once daily (with or without food)",
        "duration": "Until progression",
        "route": "Oral tablet",
        "usage_contexts": [
            "Previously treated CCA with IDH1 mutation (ClarIDHy trial)",
            "IDH1 mutation testing mandatory (NGS or PCR)",
            "PFS benefit vs placebo",
            "Intrahepatic CCA most common IDH1 mutation site",
            "Second or third-line setting"
        ],
        "indications": ["Previously treated CCA with IDH1 mutation"],
        "common_side_effects": ["Fatigue", "Nausea", "Diarrhea", "Abdominal pain", "Ascites"],
        "serious_side_effects": ["QT prolongation", "Differentiation syndrome (rare)", "Guillain-Barré (rare)"],
        "monitoring": ["ECG (QT interval)", "IDH1 mutation status", "Liver function"],
        "special_notes": "ClarIDHy trial; IDH1 mutation NGS testing required; PFS benefit vs placebo; well tolerated"
    },

    # ===== HEPATOBLASTOMA CHEMOTHERAPY =====
    "cisplatin_hepatoblastoma": {
        "names": ["cisplatin", "carboplatin", "platin"],
        "brand": "Platinol, Paraplatin",
        "generic": "Cisplatin / Carboplatin",
        "category": "Chemotherapy - Hepatoblastoma (Pediatric)",
        "cancer_types": ["Hepatoblastoma (pediatric liver cancer)"],
        "intro": "Cisplatin-based chemotherapy is the backbone of hepatoblastoma treatment, used neoadjuvantly and adjuvantly.",
        "mechanism": "Platinum-based alkylating agent; DNA cross-linking; highly effective in hepatoblastoma",
        "dose": "Cisplatin 80-100 mg/m² IV; Carboplatin AUC-based dosing",
        "frequency": "Every 14-21 days per protocol",
        "duration": "Per SIOPEL/COG protocol (4-6 cycles neoadjuvant)",
        "route": "IV infusion",
        "usage_contexts": [
            "Neoadjuvant therapy before hepatoblastoma resection (SIOPEL protocol)",
            "Adjuvant therapy post-resection",
            "PRETEXT staging guides chemotherapy approach",
            "Carboplatin used when cisplatin ototoxicity is concern",
            "Combined with doxorubicin (PLADO regimen)"
        ],
        "indications": ["Hepatoblastoma (all stages)"],
        "serious_side_effects": ["Ototoxicity (cisplatin)", "Nephrotoxicity", "Myelosuppression", "Peripheral neuropathy"],
        "monitoring": ["Audiogram before each cycle", "Renal function", "CBC"],
        "special_notes": "SIOPEL protocol standard; ototoxicity requires audiogram monitoring; carboplatin as alternative; PRETEXT staging critical"
    },

    # ===== SUPPORTIVE / ADJUVANT =====
    "sorafenib_adjuvant": {
        "names": ["sorafenib adjuvant", "nexavar adjuvant"],
        "brand": "Nexavar (adjuvant setting)",
        "generic": "Sorafenib (adjuvant)",
        "category": "Targeted Therapy - Adjuvant (Post-Resection/Ablation HCC)",
        "cancer_types": ["Resected/ablated HCC"],
        "intro": "Sorafenib studied in adjuvant setting after HCC resection/ablation — STORM trial showed no benefit; not recommended.",
        "mechanism": "Same as sorafenib (VEGFR/RAF inhibitor)",
        "dose": "400 mg oral twice daily",
        "frequency": "Twice daily",
        "duration": "Up to 4 years (STORM trial)",
        "route": "Oral",
        "usage_contexts": [
            "STORM trial: adjuvant sorafenib post-resection/ablation — NO OS benefit",
            "Not recommended as standard adjuvant therapy",
            "Historical context — negative landmark trial",
            "Ongoing trials with immunotherapy in adjuvant setting"
        ],
        "indications": ["NOT recommended — negative STORM trial"],
        "special_notes": "STORM trial NEGATIVE — no DFS/OS benefit; adjuvant immunotherapy trials ongoing (CheckMate 9DX, KEYNOTE-937)"
    },

    "bevacizumab": {
        "names": ["bevacizumab", "avastin"],
        "brand": "Avastin",
        "generic": "Bevacizumab",
        "category": "Targeted Therapy - Anti-VEGF (Used with Atezolizumab)",
        "cancer_types": ["Advanced HCC (with atezolizumab)"],
        "intro": "Bevacizumab combined with atezolizumab is the preferred first-line regimen for advanced HCC.",
        "mechanism": "VEGF-A monoclonal antibody; anti-angiogenic; reduces immunosuppressive tumor microenvironment",
        "dose": "15 mg/kg IV",
        "frequency": "Every 3 weeks (with atezolizumab)",
        "duration": "Until progression",
        "route": "IV infusion",
        "usage_contexts": [
            "Always used WITH atezolizumab (not alone) for HCC",
            "IMbrave150 trial — improved OS vs sorafenib",
            "Requires endoscopy screening for varices first",
            "Contraindicated with active bleeding or untreated varices"
        ],
        "indications": ["Advanced HCC (always with atezolizumab)"],
        "serious_side_effects": ["GI bleeding from varices", "Arterial thrombosis", "Hypertension", "Wound healing impairment"],
        "special_notes": "Never used alone for HCC; always with atezolizumab; mandatory pre-treatment endoscopy"
    },
}

# Flatten drug keywords
ALL_DRUG_KEYWORDS = []
for drug, details in LIVER_CANCER_DRUG_DATABASE.items():
    ALL_DRUG_KEYWORDS.extend(details.get("names", []))
    if details.get("brand"):
        for b in details["brand"].split(","):
            ALL_DRUG_KEYWORDS.extend(b.strip().lower().split())

ALL_DRUG_KEYWORDS = list(set(ALL_DRUG_KEYWORDS))

# ================= FUNCTIONS =================

def search_pmc():
    try:
        print("🔍 Connecting to PubMed Central...")
        handle = Entrez.esearch(db="pmc", term=SEARCH_QUERY,
                                retmax=MAX_ARTICLES, sort="relevance")
        record = Entrez.read(handle)
        return record["IdList"]
    except Exception as e:
        print(f"❌ Search error: {e}")
        return []

def fetch_xml(pmc_id):
    handle = Entrez.efetch(db="pmc", id=pmc_id,
                           rettype="full", retmode="xml")
    return handle.read()

def extract_drug_paragraphs(xml_data):
    root = ET.fromstring(xml_data)
    drug_paragraphs = []
    try:
        for elem in root.iter():
            tag = elem.tag.lower()
            if tag.endswith("p") and elem.text:
                text = elem.text.strip().lower()
                if any(kw in text for kw in ALL_DRUG_KEYWORDS):
                    clean = re.sub(r'\s+', ' ', elem.text.strip())
                    if len(clean) > 100:
                        drug_paragraphs.append(clean)
    except:
        pass
    return drug_paragraphs

def identify_drugs(text):
    text_lower = text.lower()
    found = []
    for drug, details in LIVER_CANCER_DRUG_DATABASE.items():
        for name in details.get("names", []):
            if name in text_lower:
                found.append(drug)
                break
    return list(set(found))

def get_drug_details(drug_name):
    return LIVER_CANCER_DRUG_DATABASE.get(drug_name)

# ================= MAIN =================

def main():
    print("="*120)
    print("LIVER CANCER MEDICINES EXTRACTOR - PUBMED CENTRAL")
    print("Processing up to 1000 articles — HCC, CCA, Hepatoblastoma")
    print("="*120 + "\n")

    pmc_ids = search_pmc()
    print(f"✅ Found {len(pmc_ids)} relevant articles\n")

    if not pmc_ids:
        print("❌ No articles found.")
        return

    all_contexts = []
    processed = errors = 0

    for i, pmc_id in enumerate(pmc_ids[:MAX_ARTICLES], 1):
        if i % 50 == 0:
            print(f"📊 Progress: {i}/{min(len(pmc_ids), MAX_ARTICLES)} articles...")
        try:
            xml        = fetch_xml(pmc_id)
            paragraphs = extract_drug_paragraphs(xml)
            for para in paragraphs:
                drugs = identify_drugs(para)
                all_contexts.append({
                    "pmc_id":      pmc_id,
                    "article_url": f"https://pubmed.ncbi.nlm.nih.gov/pmc/{pmc_id}",
                    "paragraph":   para,
                    "drugs_found": drugs,
                    "drug_count":  len(drugs),
                    "drug_details":[get_drug_details(d) for d in drugs if get_drug_details(d)]
                })
            processed += 1
            time.sleep(0.2)
        except:
            errors += 1
            continue

    print(f"\n✅ Processed: {processed} | ❌ Errors: {errors}")
    print(f"📝 Drug contexts extracted: {len(all_contexts)}\n")

    # Save TEXT
    print("💾 Saving TEXT...")
    with open(OUTPUT_TEXT, "w", encoding="utf-8") as f:
        f.write("="*120 + "\n")
        f.write("LIVER CANCER MEDICINES — COMPREHENSIVE EXTRACTION\n")
        f.write(f"Contexts: {len(all_contexts)} | Articles: {processed}\n")
        f.write("="*120 + "\n\n")
        for i, ctx in enumerate(all_contexts, 1):
            f.write(f"\n{'='*120}\nCONTEXT #{i}\n{'='*120}\n")
            f.write(f"PMC: {ctx['pmc_id']} | Drugs: {ctx['drug_count']}\n\n")
            for drug in ctx['drug_details']:
                if drug:
                    f.write(f"\n{'─'*80}\n💊 {drug['names'][0].upper()} ({drug['brand']})\n{'─'*80}\n")
                    f.write(f"Category  : {drug['category']}\n")
                    f.write(f"Generic   : {drug['generic']}\n")
                    f.write(f"Cancers   : {', '.join(drug.get('cancer_types', ['N/A']))}\n\n")
                    f.write(f"INTRO:\n{drug['intro']}\n\n")
                    f.write(f"MECHANISM:\n{drug['mechanism']}\n\n")
                    f.write(f"DOSING:\n  Dose      : {drug.get('dose','N/A')}\n")
                    f.write(f"  Frequency : {drug.get('frequency','N/A')}\n")
                    f.write(f"  Duration  : {drug.get('duration','N/A')}\n")
                    f.write(f"  Route     : {drug.get('route','N/A')}\n\n")
                    if drug.get('usage_contexts'):
                        f.write("USAGE:\n")
                        for u in drug['usage_contexts'][:4]:
                            f.write(f"  • {u}\n")
                        f.write("\n")
                    if drug.get('common_side_effects'):
                        f.write("SIDE EFFECTS:\n")
                        for se in drug['common_side_effects'][:5]:
                            f.write(f"  • {se}\n")
                        f.write("\n")
                    if drug.get('special_notes'):
                        f.write(f"NOTES:\n{drug['special_notes']}\n\n")
            f.write(f"\nPARAGRAPH:\n{ctx['paragraph']}\n")

    # Save JSON
    print("💾 Saving JSON...")
    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(all_contexts, f, indent=2, ensure_ascii=False)

    # Stats
    all_drugs  = [d for ctx in all_contexts for d in ctx['drugs_found']]
    drug_freq  = Counter(all_drugs)
    stats = {
        "date": str(time.ctime()),
        "articles_processed": processed,
        "contexts": len(all_contexts),
        "unique_drugs": len(drug_freq),
        "top_drugs": dict(drug_freq.most_common(20)),
        "database_size": len(LIVER_CANCER_DRUG_DATABASE)
    }
    with open(OUTPUT_STATS, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)

    # Summary
    print("💾 Saving SUMMARY...")
    with open(OUTPUT_SUMMARY, "w", encoding="utf-8") as f:
        f.write("="*120 + "\n")
        f.write("LIVER CANCER MEDICINES — SUMMARY REPORT\n")
        f.write("="*120 + "\n\n")
        f.write(f"Articles: {processed} | Contexts: {len(all_contexts)} | Unique Drugs: {len(drug_freq)}\n\n")
        f.write("TOP 20 DRUGS:\n")
        for rank, (drug, count) in enumerate(drug_freq.most_common(20), 1):
            d   = get_drug_details(drug)
            cat = d.get('category','Unknown') if d else 'Unknown'
            f.write(f"{rank:2d}. {drug:30s} {count:4d} mentions — {cat}\n")
        f.write("\nKEY TRIALS:\n")
        trials = [
            "SHARP:    Sorafenib vs placebo in advanced HCC — OS benefit",
            "REFLECT:  Lenvatinib non-inferior to sorafenib in HCC",
            "IMbrave150: Atezo+bev superior to sorafenib — now preferred 1st line",
            "HIMALAYA: Durvalumab+tremelimumab (STRIDE) vs sorafenib — non-inferior",
            "RESORCE:  Regorafenib post-sorafenib — OS benefit",
            "CELESTIAL: Cabozantinib post-sorafenib — OS benefit",
            "REACH-2:  Ramucirumab in AFP≥400 post-sorafenib — OS benefit",
            "ABC-02:   Gem+cis standard for biliary tract cancers",
            "TOPAZ-1:  Durvalumab+gem+cis superior to gem+cis in BTC",
            "FIGHT-202: Pemigatinib for FGFR2 fusion CCA",
            "ClarIDHy: Ivosidenib for IDH1-mutant CCA",
            "STORM:    Adjuvant sorafenib — NEGATIVE trial",
        ]
        for t in trials:
            f.write(f"  • {t}\n")

    print("\n" + "🎉 " + "="*118)
    print("EXTRACTION COMPLETE — LIVER CANCER MEDICINES")
    print("="*120)
    print(f"📁 Files: {OUTPUT_TEXT}, {OUTPUT_JSON}, {OUTPUT_STATS}, {OUTPUT_SUMMARY}")
    print(f"📊 Contexts: {len(all_contexts)} | Unique Drugs: {len(drug_freq)}")
    print(f"\n🏆 TOP 10 DRUGS:")
    for rank, (drug, count) in enumerate(drug_freq.most_common(10), 1):
        print(f"   {rank}. {drug}: {count}")
    print("="*120)

if __name__ == "__main__":
    main()