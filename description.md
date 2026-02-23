Here is a professional **README.md** for your project:

---

# 🧬 Prostate Cancer Prevention Extractor

### Enhanced PubMed Central (PMC) Data Extraction for LLM Training

## 📌 Overview

The **Prostate Cancer Prevention Extractor** is an advanced biomedical data collection pipeline designed to extract prevention-related information on:

* Prostate Adenocarcinoma
* Advanced Prostate Cancer
* Metastatic Prostate Cancer

The system automatically searches **PubMed Central (PMC)** using the NCBI Entrez API, extracts prevention-focused paragraphs, maps them to a structured knowledge base, and generates:

* 📄 Detailed Text Dataset
* 📦 Structured JSON Dataset
* 🧠 Enriched JSON for LLM Training

This tool is built for **Cancer LLM development**, research automation, and structured biomedical knowledge engineering.

---

## 🎯 Project Purpose

This system is designed to:

* Build a high-quality **Prostate Cancer Prevention dataset**
* Support **LLM fine-tuning and RAG systems**
* Extract clinically relevant screening and prevention data
* Structure prevention knowledge into categorized topics
* Generate machine-learning-ready training datasets

---

## 🔍 Prevention Categories Covered

The knowledge base integrates structured medical evidence across:

### 🔎 1. Screening & Early Detection

* PSA Screening
* Digital Rectal Examination (DRE)
* Shared Decision-Making
* Active Surveillance
* Risk Calculators
* Multiparametric MRI (PI-RADS)

### 💊 2. Chemoprevention

* 5-Alpha Reductase Inhibitors (Finasteride, Dutasteride)
* Vitamin E & Selenium (SELECT Trial findings)
* FDA recommendations

### 🥗 3. Diet & Nutrition

* Mediterranean Diet
* Plant-Based Diet
* Lycopene (Tomatoes)
* Cruciferous Vegetables
* Red & Processed Meat Risks
* Calcium Intake

### 🏃 4. Physical Activity & Obesity

* Exercise recommendations (150+ min/week)
* BMI & Waist circumference guidelines
* Advanced disease risk reduction
* Hormonal & metabolic mechanisms

### 🧬 5. Genetic & Familial Risk

* BRCA1 / BRCA2
* HOXB13
* Lynch Syndrome
* Germline testing
* PARP inhibitors implications
* Cascade testing

### ⚖️ 6. Racial & Ethnic Disparities

* African American risk statistics
* Screening age adjustments
* Social determinants of health
* Health equity interventions

### ❤️ 7. Lifestyle & Emerging Observations

* Ejaculation frequency epidemiological findings
* Observational limitations

---

## 🏗️ System Architecture

```
PubMed Central (PMC)
        ↓
Entrez Search (Multiple Queries)
        ↓
XML Retrieval (Full Articles)
        ↓
Keyword-Based Paragraph Extraction
        ↓
Topic Identification Engine
        ↓
Prevention Knowledge Mapping
        ↓
Structured Output Generation
```

---

## 🔬 Search Queries Used

The system searches PMC using prevention-focused queries such as:

* prostate cancer prevention
* PSA screening prostate cancer
* finasteride prostate cancer prevention
* diet prostate cancer prevention
* BRCA prostate cancer risk
* racial disparities prostate cancer
* multiparametric MRI prostate

Each query retrieves up to 100 articles sorted by relevance.

---

## 📂 Output Files Generated

| File                                       | Description                                       |
| ------------------------------------------ | ------------------------------------------------- |
| `prostate_cancer_prevention_detailed.txt`  | Human-readable enriched prevention dataset        |
| `prostate_cancer_prevention_detailed.json` | Structured JSON with topic mapping                |
| `prostate_cancer_prevention_DETAILED.json` | Fully enriched JSON including prevention metadata |

Each context contains:

```json
{
  "pmc_id": "1234567",
  "article_url": "PMC URL",
  "paragraph": "Extracted prevention paragraph...",
  "topics_found": ["psa_screening", "genetic_familial_prostate"],
  "topic_details": [...]
}
```

---

## ⚙️ Configuration

Edit these parameters inside the script:

```python
Entrez.email = "your_email@example.com"
MAX_ARTICLES = 100
MAX_RETRIES = 3
RETRY_DELAY = 2
```

⚠️ NCBI requires a valid email address.

---

## 🚀 How to Run

### 1️⃣ Install Dependencies

```bash
pip install biopython
```

### 2️⃣ Run Script

```bash
python prostate_prevention_extractor.py
```

---

## 📊 Key Evidence Embedded in Database

* PSA Screening → 20–30% mortality reduction (with overdiagnosis risk)
* African American men → 70% higher incidence, 2x mortality
* BRCA2 mutation → 5–8x increased risk
* Vitamin E → 17% increased prostate cancer risk (SELECT Trial)
* 5-ARIs → 23–25% risk reduction but not recommended
* Physical activity → 30–50% reduction in advanced disease risk
* Active Surveillance → 50% avoid treatment long-term

---

## 🧠 Designed For

* Cancer LLM development
* Biomedical NLP
* Clinical AI assistants
* RAG-based oncology systems
* Research automation
* Evidence synthesis pipelines

---

## 🛡️ Ethical & Scientific Considerations

* Observational findings are labeled appropriately.
* Randomized Controlled Trial (RCT) evidence prioritized.
* Supplements with proven harm (Vitamin E) flagged.
* Shared decision-making emphasized.
* Health equity explicitly integrated.

---

## 📌 Limitations

* Keyword-based extraction (not semantic embedding-based)
* Observational associations cannot prove causation
* Dependent on PMC article availability
* Requires internet access
* Not a substitute for medical advice

---

## 🔮 Future Improvements

* Add embedding-based semantic filtering
* Add citation extraction
* Integrate UMLS ontology mapping
* Add RAG chunking for LLM ingestion
* Deduplicate near-identical paragraphs
* Add PubMed (non-PMC) abstract integration
* Add vector database support (FAISS / Chroma)

---

## 👨‍💻 Author

**Nishan Kharel**
Machine Learning Engineer
Biomedical AI & Cancer LLM Research

GitHub: [https://github.com/nishankhareln](https://github.com/nishankhareln)

---

## 📜 License

For research and educational purposes.

---

If you want, I can also generate:

* 🔥 A more research-paper-style README
* 📘 A version optimized for GitHub visibility
* 🧠 A README tailored for your Cancer LLM main repository
* 📊 A system design diagram version
* 📄 A citation-ready academic description

Just tell me your target audience (research, open-source, recruiter, or publication).
