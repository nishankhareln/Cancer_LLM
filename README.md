# 🧬 Cancer_LLM

**Cancer_LLM** is a research-oriented project for building high-quality text datasets focused on cancer outcomes, prognosis, and clinical research.  
It provides scripts to scrape, clean, and organize oncology articles into structured text that can be used for **medical large language model (LLM) training** and **AI research**.

> *Cancer based LLM project.* :contentReference[oaicite:1]{index=1}

---

## 📌 Project Overview

Modern research in oncology and medical AI often requires large, domain-specific corpora. This repository includes Python scripts that collect and clean biomedical literature related to cancer prognosis and symptoms.

The core idea is to gather research articles on various cancer types from public sources (like PubMed Central), extract the critical sections of the text (title, abstract, body), remove noise, and save them in a unified format suitable for LLM training.

---

## 📂 Repository Contents

| File | Description |
|------|-------------|
| `readme for prognosis` | (Draft) prognosis dataset README |
| `prevention_cleaning.py` | Script to clean prevention-related text |
| `symptoms_cleaned.py` | Script to clean symptom-related text |
| `.gitignore` | Standard ignore file for Python projects |

> This project currently focuses on **data preparation and cleaning for cancer literature**.

---

## 🧠 Key Features

✔ Extracts full-text biomedical articles  
✔ Cleans and normalizes text (removes references, figures, URLs)  
✔ Supports structured clinical NLP dataset creation  
✔ Prepares data for downstream AI/ML workflows  
✔ Modular and extensible Python code

---

## 🛠 How It Works

1. **Search & Download**  
   Uses NCBI’s E-utilities (`Bio.Entrez`) for searching and fetching PMC articles.

2. **XML Parsing**  
   Extracts meaningful sections (title, abstract, body) from PMC XML.

3. **Text Cleaning**  
   Removes reference markers, figure/table mentions, URLs, and hides noise.

4. **Data Saving**  
   Saves cleaned text into organized `.txt` files.

---

## 🚀 Getting Started

### Prerequisites

Make sure you have Python installed (≥3.8) and install dependencies:

```bash
pip install biopython
