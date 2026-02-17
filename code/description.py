"""
Cancer LLM - DIAGNOSTICS Articles Collector
Extracts clean text about cancer diagnosis methods
"""

from Bio import Entrez
import xml.etree.ElementTree as ET
import time
import os
import re

# ========== CONFIGURATION ==========
Entrez.email = "nkharel57@gmail.com"

ARTICLES_PER_CANCER = 1000  # Target 300 diagnostics articles

# Choose: Just Lung Cancer OR All 10 Cancers
CANCER_TYPES = {
    "Prostate_Cancer": "prostate cancer"
}

# ========== FUNCTIONS ==========

def extract_text_from_xml(xml_content):
    """
    Extract clean text from PMC XML
    """
    try:
        root = ET.fromstring(xml_content)
        extracted_parts = []
        
        # 1. EXTRACT TITLE
        title = root.find(".//article-title")
        if title is not None:
            title_text = ''.join(title.itertext()).strip()
            title_text = re.sub(r'<[^>]+>', '', title_text)
            extracted_parts.append("=" * 80)
            extracted_parts.append("TITLE")
            extracted_parts.append("=" * 80)
            extracted_parts.append(title_text)
            extracted_parts.append("")
        
        # 2. EXTRACT ABSTRACT
        abstract = root.find(".//abstract")
        if abstract is not None:
            abstract_text = ' '.join(abstract.itertext()).strip()
            abstract_text = re.sub(r'\s+', ' ', abstract_text)
            extracted_parts.append("=" * 80)
            extracted_parts.append("ABSTRACT")
            extracted_parts.append("=" * 80)
            extracted_parts.append(abstract_text)
            extracted_parts.append("")
        
        # 3. EXTRACT BODY SECTIONS
        body = root.find(".//body")
        if body is not None:
            for section in body.findall(".//sec"):
                sec_title = section.find("title")
                if sec_title is not None:
                    sec_title_text = ''.join(sec_title.itertext()).strip()
                    
                    # Skip unwanted sections
                    skip_sections = [
                        'reference', 'acknowledgment', 'conflict', 
                        'author contribution', 'funding', 'supplementary',
                        'abbreviation'
                    ]
                    if any(skip in sec_title_text.lower() for skip in skip_sections):
                        continue
                    
                    extracted_parts.append("=" * 80)
                    extracted_parts.append(sec_title_text.upper())
                    extracted_parts.append("=" * 80)
                
                # Get paragraphs
                paragraphs = []
                for para in section.findall(".//p"):
                    para_text = ' '.join(para.itertext()).strip()
                    para_text = re.sub(r'\s+', ' ', para_text)
                    if len(para_text) > 30:
                        paragraphs.append(para_text)
                
                if paragraphs:
                    extracted_parts.append('\n\n'.join(paragraphs))
                    extracted_parts.append("")
        
        full_text = '\n'.join(extracted_parts)
        full_text = clean_text(full_text)
        
        return full_text
    
    except Exception as e:
        print(f"      ❌ XML parsing error: {e}")
        return None

def clean_text(text):
    """
    Clean extracted text
    """
    # Remove citations
    text = re.sub(r'\[\d+(?:[-–,]\s*\d+)*\]', '', text)
    text = re.sub(r'\(\d+(?:[-–,]\s*\d+)*\)', '', text)
    
    # Remove figure/table references
    text = re.sub(r'Figure\s+\d+[A-Z]?', '', text, flags=re.IGNORECASE)
    text = re.sub(r'Table\s+\d+', '', text, flags=re.IGNORECASE)
    text = re.sub(r'Fig\.\s*\d+', '', text, flags=re.IGNORECASE)
    text = re.sub(r'as\s+shown\s+in\s*\.', '', text, flags=re.IGNORECASE)
text = re.sub(r'as\s+shown\s+in\s+\w+\s*\.', 'as shown.', text, flags=re.IGNORECASE)

# Fix remaining author citations like "Sun et al., 2004"
text = re.sub(r'\b[A-Z][a-z]+\s+et\s+al\.,\s+\d{4}\b', '', text)

# Fix scientific notation spacing (critical!)
text = re.sub(r'(\d+)\s*×\s*10\s+(\d+)', r'\1 × 10^\2', text)
text = re.sub(r'(\d+)\s*x\s*10\s+(\d+)', r'\1 × 10^\2', text)
    
    # Remove URLs and emails
    text = re.sub(r'http[s]?://\S+', '', text)
    text = re.sub(r'\S+@\S+', '', text)
    
    # Clean whitespace
    text = re.sub(r'\n\s*\n\s*\n+', '\n\n', text)
    text = re.sub(r' +', ' ', text)
    
    return text.strip()

def search_diagnostics_articles(cancer_name, max_results=1000):
    """
    Search PMC for DIAGNOSTICS articles
    """
    # DIAGNOSTICS-SPECIFIC SEARCH QUERIES
    queries = [
        f'"{cancer_name}" AND (diagnosis OR diagnostic) AND "open access"[filter]',
        f'"{cancer_name}" AND (biomarkers OR screening) AND "open access"[filter]',
        f'"{cancer_name}" AND (imaging OR CT OR MRI OR PET) AND "open access"[filter]',
        f'"{cancer_name}" AND (detection OR early detection) AND "open access"[filter]',
        f'"{cancer_name}" AND (genomic testing OR molecular testing) AND "open access"[filter]'
    ]
    
    all_pmc_ids = set()
    
    for i, query in enumerate(queries):
        print(f"   🔍 Diagnostics search {i+1}/5...")
        
        try:
            handle = Entrez.esearch(
                db="pmc",
                term=query,
                retmax=max_results,
                sort="relevance"
            )
            record = Entrez.read(handle)
            pmc_ids = record["IdList"]
            all_pmc_ids.update(pmc_ids)
            
            print(f"      Found {len(pmc_ids)} articles")
            time.sleep(0.5)
            
        except Exception as e:
            print(f"      ⚠️ Search error: {e}")
            continue
    
    print(f"   ✅ Total unique articles: {len(all_pmc_ids)}")
    return list(all_pmc_ids)

def download_and_extract_article(pmc_id, output_folder):
    """
    Download and extract article
    """
    try:
        handle = Entrez.efetch(
            db="pmc",
            id=pmc_id,
            rettype="full",
            retmode="xml"
        )
        
        xml_content = handle.read()
        
        if isinstance(xml_content, bytes):
            xml_content = xml_content.decode('utf-8', errors='ignore')
        
        clean_text_content = extract_text_from_xml(xml_content)
        
        if not clean_text_content or len(clean_text_content) < 500:
            return False
        
        filename = f"{output_folder}/PMC{pmc_id}_diagnostics.txt"
        with open(filename, "w", encoding="utf-8") as f:
            f.write(clean_text_content)
        
        return True
    
    except Exception as e:
        return False

def collect_diagnostics_for_cancer(cancer_folder_name, cancer_search_term, target_articles):
    """
    Collect DIAGNOSTICS articles for one cancer type
    """
    print(f"\n{'='*70}")
    print(f"🔬 {cancer_folder_name.upper()} - DIAGNOSTICS")
    print(f"   Target: {target_articles} articles")
    print(f"{'='*70}")
    
    # Create folder - DIAGNOSTICS subfolder
    diagnostics_folder = f"Cancer_LLM_Data/{cancer_folder_name}/Diagnostics"
    os.makedirs(diagnostics_folder, exist_ok=True)
    
    # Search
    print(f"   Step 1: Searching PMC for diagnostics articles...")
    pmc_ids = search_diagnostics_articles(cancer_search_term, max_results=4000)
    
    if not pmc_ids:
        print(f"   ❌ No articles found!")
        return 0
    
    print(f"   ✅ Found {len(pmc_ids)} diagnostics articles")
    
    # Download
    print(f"   Step 2: Downloading & extracting...")
    print(f"   ⏱️ Estimated time: {target_articles * 1.5 / 60:.0f} minutes")
    
    successful = 0
    failed = 0
    
    for i, pmc_id in enumerate(pmc_ids):
        if (i + 1) % 20 == 0:
            print(f"      Progress: {i+1} | Downloaded: {successful} | Failed: {failed}")
        
        if download_and_extract_article(pmc_id, diagnostics_folder):
            successful += 1
        else:
            failed += 1
        
        time.sleep(1)
        
        if successful >= target_articles:
            print(f"      🎯 Target reached! Stopping at {successful}")
            break
    
    print(f"\n   ✅ Downloaded: {successful}")
    print(f"   ⚠️ Failed: {failed}")
    print(f"   📁 Saved to: {diagnostics_folder}/")
    
    return successful

# ========== MAIN ==========

def main():
    print("\n" + "="*70)
    print("🔬 CANCER LLM - DIAGNOSTICS ARTICLES COLLECTOR")
    print("="*70)
    print(f"Target: {ARTICLES_PER_CANCER} diagnostics articles per cancer")
    print(f"Focus: Diagnosis methods, biomarkers, imaging, screening")
    print("="*70)
    
    print("\n📋 What this collects:")
    print("   - Diagnostic imaging (CT, MRI, PET)")
    print("   - Biomarkers and tumor markers")
    print("   - Screening and early detection")
    print("   - Genomic/molecular testing")
    print("   - Biopsy and histopathology")
    
    input("\n⚠️ Press ENTER to start collection...")
    
    total_collected = 0
    results = {}
    start_time = time.time()
    
    for idx, (cancer_folder, cancer_search) in enumerate(CANCER_TYPES.items(), 1):
        print(f"\n🎯 CANCER {idx}/{len(CANCER_TYPES)}")
        
        count = collect_diagnostics_for_cancer(
            cancer_folder,
            cancer_search,
            ARTICLES_PER_CANCER
        )
        
        results[cancer_folder] = count
        total_collected += count
    
    # Summary
    elapsed_hours = (time.time() - start_time) / 3600
    
    print("\n" + "="*70)
    print("🎉 DIAGNOSTICS COLLECTION COMPLETE!")
    print("="*70)
    
    for cancer, count in results.items():
        status = "✅" if count >= 250 else "⚠️"
        print(f"{status} {cancer}: {count} diagnostics articles")
    
    print(f"\n📊 TOTAL: {total_collected} diagnostics articles")
    print(f"⏱️ Time: {elapsed_hours:.1f} hours")
    print(f"📁 Location: Cancer_LLM_Data/")
    
    print("\n📂 Your folder structure now:")
    print("Cancer_LLM_Data/")
    print("└── prostate_Cancer/")
    print("    ├── Biology/        (1000 articles) ✅")
    print("    └── Diagnostics/    (1000 articles) ✅")
    
    print("="*70)

if __name__ == "__main__":
    main()


