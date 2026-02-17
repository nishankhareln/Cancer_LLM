"""
Cancer LLM - Biology Articles Collector with CLEAN TEXT EXTRACTION
Downloads articles and extracts only important text (no XML tags)
"""

from Bio import Entrez
import xml.etree.ElementTree as ET
import time
import os
import re

# ========== CONFIGURATION ==========
Entrez.email = "nkharel57@gmail.com"  # ← CHANGE THIS!

ARTICLES_PER_CANCER = 1000  # Per cancer

CANCER_TYPES = {
    "Prostate_Cancer": "prostate cancer"

}

# ========== FUNCTIONS ==========

def extract_text_from_xml(xml_content):
    """
    Extract ONLY the important text from PMC XML
    Returns clean text without any XML tags or metadata
    """
    try:
        # Parse XML
        root = ET.fromstring(xml_content)
        
        extracted_parts = []
        
        # 1. EXTRACT TITLE
        title = root.find(".//article-title")
        if title is not None:
            title_text = ''.join(title.itertext()).strip()
            extracted_parts.append(f"TITLE:\n{title_text}\n")
        
        # 2. EXTRACT ABSTRACT
        abstract = root.find(".//abstract")
        if abstract is not None:
            abstract_text = ' '.join(abstract.itertext()).strip()
            # Clean up extra spaces
            abstract_text = re.sub(r'\s+', ' ', abstract_text)
            extracted_parts.append(f"\nABSTRACT:\n{abstract_text}\n")
        
        # 3. EXTRACT BODY (Main content)
        body = root.find(".//body")
        if body is not None:
            # Get all sections
            for section in body.findall(".//sec"):
                # Get section title
                sec_title = section.find("title")
                if sec_title is not None:
                    sec_title_text = ''.join(sec_title.itertext()).strip()
                    extracted_parts.append(f"\n{sec_title_text.upper()}:")
                
                # Get all paragraphs in this section
                for para in section.findall(".//p"):
                    para_text = ' '.join(para.itertext()).strip()
                    # Clean up whitespace
                    para_text = re.sub(r'\s+', ' ', para_text)
                    if len(para_text) > 20:  # Only include substantial paragraphs
                        extracted_parts.append(para_text)
        
        # 4. EXTRACT CONCLUSION (if separate)
        conclusion = root.find(".//sec[@sec-type='conclusions']")
        if conclusion is not None:
            conclusion_text = ' '.join(conclusion.itertext()).strip()
            extracted_parts.append(f"\nCONCLUSION:\n{conclusion_text}")
        
        # Combine all parts
        full_text = '\n\n'.join(extracted_parts)
        
        # Final cleanup
        full_text = clean_text(full_text)
        
        return full_text
    
    except Exception as e:
        print(f"      ❌ XML parsing error: {e}")
        return None

def clean_text(text):
    """
    Final text cleaning - remove unwanted elements
    """
    # Remove citation markers [1], [2,3], etc.
    text = re.sub(r'\[\d+(?:[-–,]\s*\d+)*\]', '', text)
    text = re.sub(r'\(\d+(?:[-–,]\s*\d+)*\)', '', text)
    
    # Remove figure/table references
    text = re.sub(r'Figure\s+\d+[A-Z]?', '', text, flags=re.IGNORECASE)
    text = re.sub(r'Table\s+\d+', '', text, flags=re.IGNORECASE)
    text = re.sub(r'Fig\.\s*\d+', '', text, flags=re.IGNORECASE)
    
    # Remove URLs
    text = re.sub(r'http[s]?://\S+', '', text)
    
    # Remove email addresses
    text = re.sub(r'\S+@\S+', '', text)
    
    # Remove excessive whitespace
    text = re.sub(r'\n\s*\n\s*\n+', '\n\n', text)
    text = re.sub(r' +', ' ', text)
    
    # Remove very short lines (likely artifacts)
    lines = text.split('\n')
    lines = [line.strip() for line in lines if len(line.strip()) > 15 or line.strip().endswith(':')]
    text = '\n'.join(lines)
    
    return text.strip()

def search_biology_articles(cancer_name, max_results=4000):
    """
    Search PMC for biology articles
    """
    queries = [
        f'"{cancer_name}" AND (biology OR pathogenesis OR molecular mechanisms) AND "open access"[filter]',
        f'"{cancer_name}" AND (genetics OR genomics OR mutations) AND "open access"[filter]',
        f'"{cancer_name}" AND (cell biology OR tumor biology OR oncogenes) AND "open access"[filter]',
        f'"{cancer_name}" AND (molecular biology OR signaling pathways) AND "open access"[filter]'
    ]
    
    all_pmc_ids = set()
    
    for i, query in enumerate(queries):
        print(f"   🔍 Search {i+1}/4...")
        
        try:
            handle = Entrez.esearch(
                db="pmc",
                term=query,
                retmax=max_results,
                sort="relevance"
            )
            record = Entrez.read(handle)
            all_pmc_ids.update(record["IdList"])
            
            time.sleep(0.5)
            
        except Exception as e:
            print(f"      ⚠️ Search error: {e}")
            continue
    
    return list(all_pmc_ids)

def download_and_extract_article(pmc_id, output_folder):
    """
    Download article XML and extract clean text
    """
    try:
        # Download XML
        handle = Entrez.efetch(
            db="pmc",
            id=pmc_id,
            rettype="full",
            retmode="xml"
        )
        
        xml_content = handle.read()
        
        # Convert bytes to string if needed
        if isinstance(xml_content, bytes):
            xml_content = xml_content.decode('utf-8', errors='ignore')
        
        # Extract clean text from XML
        clean_text_content = extract_text_from_xml(xml_content)
        
        if not clean_text_content or len(clean_text_content) < 500:
            return False
        
        # Save clean text
        filename = f"{output_folder}/PMC{pmc_id}_biology.txt"
        with open(filename, "w", encoding="utf-8") as f:
            f.write(clean_text_content)
        
        return True
    
    except Exception as e:
        return False

def collect_biology_for_cancer(cancer_folder_name, cancer_search_term, target_articles):
    """
    Collect biology articles for one cancer type
    """
    print(f"\n{'='*70}")
    print(f"📚 {cancer_folder_name.upper()}")
    print(f"   Target: {target_articles} articles")
    print(f"{'='*70}")
    
    # Create folder
    biology_folder = f"Cancer_LLM_Data/{cancer_folder_name}/Biology"
    os.makedirs(biology_folder, exist_ok=True)
    
    # Search
    print(f"   Step 1: Searching PMC...")
    pmc_ids = search_biology_articles(cancer_search_term, max_results=4000)
    
    if not pmc_ids:
        print(f"   ❌ No articles found!")
        return 0
    
    print(f"   ✅ Found {len(pmc_ids)} unique articles")
    
    # Download
    print(f"   Step 2: Downloading & extracting clean text...")
    
    successful = 0
    failed = 0
    
    for i, pmc_id in enumerate(pmc_ids):
        if (i + 1) % 20 == 0:
            print(f"      Progress: {i+1} | Success: {successful} | Failed: {failed}")
        
        if download_and_extract_article(pmc_id, biology_folder):
            successful += 1
        else:
            failed += 1
        
        time.sleep(1)
        
        # Stop when target reached
        if successful >= target_articles:
            break
    
    print(f"\n   ✅ Downloaded: {successful}")
    print(f"   ⚠️ Failed: {failed}")
    print(f"   📁 Location: {biology_folder}/")
    
    return successful

# ========== MAIN ==========

def main():
    print("\n" + "="*70)
    print("🚀 CANCER LLM - CLEAN TEXT BIOLOGY COLLECTOR")
    print("="*70)
    print(f"Target: {ARTICLES_PER_CANCER} articles per cancer")
    print(f"Output: PURE TEXT (no XML, no HTML tags)")
    print("="*70)
    
    input("\n⚠️ This will take several hours. Press ENTER to start...")
    
    total_collected = 0
    results = {}
    start_time = time.time()
    
    for idx, (cancer_folder, cancer_search) in enumerate(CANCER_TYPES.items(), 1):
        print(f"\n🎯 CANCER {idx}/10")
        
        count = collect_biology_for_cancer(
            cancer_folder,
            cancer_search,
            ARTICLES_PER_CANCER
        )
        
        results[cancer_folder] = count
        total_collected += count
    
    # Summary
    elapsed_hours = (time.time() - start_time) / 3600
    
    print("\n" + "="*70)
    print("🎉 COLLECTION COMPLETE!")
    print("="*70)
    
    for cancer, count in results.items():
        status = "✅" if count >= 250 else "⚠️"
        print(f"{status} {cancer}: {count} articles")
    
    print(f"\n📊 TOTAL: {total_collected} clean text articles")
    print(f"⏱️ Time: {elapsed_hours:.1f} hours")
    print(f"📁 Location: Cancer_LLM_Data/")
    print("="*70)

if __name__ == "__main__":
    main()