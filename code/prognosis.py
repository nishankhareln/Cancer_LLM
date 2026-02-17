"""
Cancer LLM – Blood Cancer PROGNOSIS Dataset Builder
Collects prognosis, survival & outcome articles for ALL blood cancer types
Saves everything into ONE unified folder
"""

from Bio import Entrez
import xml.etree.ElementTree as ET
import time
import os
import re

Entrez.email = "nkharel57@gmail.com"

ARTICLES_PER_TYPE = 500   # 12 x 500 = 6,000 papers (safe for PMC)

# ALL MAJOR BLOOD CANCERS
BLOOD_CANCERS = {
    "AML": "acute myeloid leukemia",
    "ALL": "acute lymphoblastic leukemia",
    "CML": "chronic myeloid leukemia",
    "CLL": "chronic lymphocytic leukemia",
    "Hodgkin": "hodgkin lymphoma",
    "NonHodgkin": "non hodgkin lymphoma",
    "DLBCL": "diffuse large B cell lymphoma",
    "Follicular": "follicular lymphoma",
    "Myeloma": "multiple myeloma",
    "Plasmacytoma": "plasmacytoma",
    "MDS": "myelodysplastic syndrome",
    "MPN": "myeloproliferative neoplasm"
}

OUTPUT_FOLDER = "Cancer_LLM_Data/Blood_Cancer/Prognosis"
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# ---------------- XML Processing ---------------- #

def clean_text(text):
    """Remove references, figures, tables, URLs and normalize whitespace"""
    # Remove reference markers
    text = re.sub(r'\[\d+.*?\]', '', text)
    text = re.sub(r'\(\d+.*?\)', '', text)
    
    # Remove figure/table mentions
    text = re.sub(r'Figure\s+\d+', '', text, flags=re.I)
    text = re.sub(r'Table\s+\d+', '', text, flags=re.I)
    
    # Remove URLs
    text = re.sub(r'http\S+', '', text)
    
    # Normalize whitespace within lines but preserve paragraph breaks
    lines = text.split('\n')
    cleaned_lines = []
    for line in lines:
        line = re.sub(r'[ \t]+', ' ', line)  # Only collapse spaces/tabs
        line = line.strip()
        if line:  # Keep non-empty lines
            cleaned_lines.append(line)
    
    return '\n\n'.join(cleaned_lines)

def extract_text(xml):
    """Extract title, abstract, and body from PMC XML"""
    try:
        root = ET.fromstring(xml)
        parts = []

        # Extract title
        title = root.find(".//article-title")
        if title is not None and title.text:
            parts.append("TITLE: " + ''.join(title.itertext()))

        # Extract abstract
        abstract = root.find(".//abstract")
        if abstract is not None:
            abstract_text = ' '.join(abstract.itertext()).strip()
            if abstract_text:
                parts.append("\nABSTRACT:\n" + abstract_text)

        # Extract body paragraphs
        body = root.find(".//body")
        if body is not None:
            for p in body.findall(".//p"):
                txt = ' '.join(p.itertext()).strip()
                if len(txt) > 50:  # Only include substantial paragraphs
                    parts.append(txt)

        full_text = "\n\n".join(parts)
        return clean_text(full_text) if full_text else None
        
    except ET.ParseError as e:
        print(f"      ⚠️  XML Parse Error: {e}")
        return None
    except Exception as e:
        print(f"      ⚠️  Extraction Error: {e}")
        return None

# ---------------- PMC Search ---------------- #

def search_articles(cancer):
    """Search PMC for prognosis-related articles"""
    queries = [
        f'"{cancer}" AND (prognosis OR prognostic)',
        f'"{cancer}" AND (survival OR mortality)',
        f'"{cancer}" AND (outcome OR outcomes)',
        f'"{cancer}" AND (recurrence OR relapse)',
        f'"{cancer}" AND (staging OR risk)'
    ]

    ids = set()
    for q in queries:
        try:
            handle = Entrez.esearch(
                db="pmc", 
                term=q, 
                retmax=500,
                usehistory="y"
            )
            record = Entrez.read(handle)
            ids.update(record["IdList"])
            handle.close()
            time.sleep(0.4)  # Respect NCBI rate limits (3 requests/sec)
        except Exception as e:
            print(f"      ⚠️  Search error for '{q}': {e}")
            time.sleep(2)  # Back off on error

    return list(ids)

# ---------------- Download ---------------- #

def download(pmc_id, cancer_tag):
    """Download and save a single PMC article"""
    try:
        handle = Entrez.efetch(
            db="pmc", 
            id=pmc_id, 
            rettype="full", 
            retmode="xml"
        )
        xml = handle.read()
        handle.close()
        
        # Decode with error handling
        if isinstance(xml, bytes):
            xml = xml.decode("utf-8", errors="ignore")
        
        text = extract_text(xml)

        if text and len(text) > 500:
            filename = f"{OUTPUT_FOLDER}/{cancer_tag}_PMC{pmc_id}.txt"
            with open(filename, "w", encoding="utf-8") as f:
                f.write(text)
            return True
        else:
            return False
            
    except Exception as e:
        print(f"      ⚠️  Download error for PMC{pmc_id}: {e}")
        return False

# ---------------- Main ---------------- #

def main():
    print("\n🩸 BLOOD CANCER PROGNOSIS DATASET BUILDER\n")
    print(f"📧 Email: {Entrez.email}")
    print(f"📁 Output: {OUTPUT_FOLDER}\n")

    total = 0
    failed = 0

    for tag, name in BLOOD_CANCERS.items():
        print(f"\n🔬 {name.upper()} ({tag})")
        
        # Search for articles
        ids = search_articles(name)
        print(f"   Found {len(ids)} papers")

        if not ids:
            print(f"   ⚠️  No articles found, skipping...")
            continue

        # Download articles
        count = 0
        for i, pmc in enumerate(ids):
            if download(pmc, tag):
                count += 1
                total += 1
                if count % 50 == 0:
                    print(f"   Progress: {count}/{min(len(ids), ARTICLES_PER_TYPE)}")
            else:
                failed += 1
                
            if count >= ARTICLES_PER_TYPE:
                break
                
            time.sleep(1.0)  # Be conservative with rate limiting

        print(f"   ✅ Successfully collected: {count}")

    print("\n" + "="*60)
    print(f"🎯 TOTAL BLOOD CANCER PROGNOSIS PAPERS: {total}")
    print(f"❌ Failed downloads: {failed}")
    print(f"📁 Folder: {OUTPUT_FOLDER}")
    print("\n✨ Dataset ready for Medical LLM Training!")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()