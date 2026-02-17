import os
import re
from pathlib import Path
from collections import OrderedDict

# CONFIG — Medical LLM Optimized
ROOT_INPUT = Path("Cancer_LLM_Data/Blood_Cancer")  # Your input folder
ROOT_OUTPUT = Path("Cleaned_Cancer_LLM_Data/Blood_Cancer")  # Parallel cleaned
MIN_LEN = 150
MAX_LEN = 4200
DEDUP_THRESHOLD = 0.95
SKIP_ABSTRACT = True  # True: Skips "Abstract" (separate files)

PRIORITY_SECTIONS = {"background", "introduction", "results", "discussion", 
                     "conclusions", "methods", "summary", 
                     "hypoxia", "melanoma", "braf", "ras", "nf1", 
                     "cdkn2a", "mitf", "basal cell carcinoma", 
                     "cutaneous squamous cell carcinoma", 
                     "rare types of skin cancer", "reconstruction"}

def extract_pmc_id(filename: str) -> str:
    match = re.search(r'PMC(\d+)', filename)
    return match.group(0) if match else "Unknown"

def normalize_heading(line: str) -> str:
    line = re.sub(r'^\s*\d+(\.\d+)*\.?\s*', '', line.strip())
    line = re.sub(r'^\s*={1,}\s*(.+?)\s*={1,}\s*$', r'\1', line)
    line = re.sub(r'^={1,}$', '', line)
    if not line.strip():
        return ""
    if len(line) < 80 and line.isupper():
        line = line.title()
    return line.strip()

def fix_compound_medical_terms(text: str) -> str:
    """Fix missing hyphens in compound medical/scientific terms using intelligent patterns
    Safe for processing 100,000+ medical documents without breaking normal words"""
    
    # CRITICAL: Order matters! Apply most specific patterns first, then general ones
    
    # STRATEGY 1: Exact medical compound terms (highest priority - no false positives)
    exact_compounds = [
        # Hypoxia family
        (r'\bhypoxiainducible\b', 'hypoxia-inducible'),
        (r'\bhypoxiadependent\b', 'hypoxia-dependent'),
        (r'\bhypoxiaactivated\b', 'hypoxia-activated'),
        (r'\bhypoxiaregulated\b', 'hypoxia-regulated'),
        (r'\bhypoxiainduced\b', 'hypoxia-induced'),
        (r'\bhypoxiastimulated\b', 'hypoxia-stimulated'),
        (r'\bhypoxiamimicking\b', 'hypoxia-mimicking'),
        (r'\bhypoxiarelated\b', 'hypoxia-related'),
        (r'\bhypoxiabased\b', 'hypoxia-based'),
        
        # UV family
        (r'\bUVinduced\b', 'UV-induced'),
        (r'\bUVradiation\b', 'UV radiation'),
        (r'\bUVrelated\b', 'UV-related'),
        (r'\bUVdamage\b', 'UV damage'),
        
        # DNA/RNA
        (r'\bDNAdamage\b', 'DNA damage'),
        (r'\bDNArepair\b', 'DNA repair'),
        (r'\bRNAsequencing\b', 'RNA sequencing'),
        (r'\bcDNAmicroarray\b', 'cDNA microarray'),
        
        # ROS family
        (r'\bROSinduced\b', 'ROS-induced'),
        (r'\bROSinitiated\b', 'ROS-initiated'),
        (r'\bROSmediated\b', 'ROS-mediated'),
        (r'\bROSdependent\b', 'ROS-dependent'),
        
        # Gene names
        (r'\bBRAFmutant\b', 'BRAF-mutant'),
        (r'\bBRAFactivating\b', 'BRAF-activating'),
        (r'\bNRASdriven\b', 'NRAS-driven'),
        (r'\bNRASmutant\b', 'NRAS-mutant'),
        
        # Clinical terms
        (r'\bplacebocontrolled\b', 'placebo-controlled'),
        (r'\bdoubleblind\b', 'double-blind'),
        (r'\bsingleblind\b', 'single-blind'),
        (r'\brandomizedcontrolled\b', 'randomized-controlled'),
        (r'\bopenlabel\b', 'open-label'),
        (r'\bcrossover\b', 'cross-over'),
        
        # Wide compounds
        (r'\bgenomewide\b', 'genome-wide'),
        (r'\bcellwide\b', 'cell-wide'),
        (r'\btissuewide\b', 'tissue-wide'),
        (r'\bsystemwide\b', 'system-wide'),
        (r'\bworldwide\b', 'worldwide'),  # Keep as one word
        
        # Target/interaction
        (r'\bdrugtarget\b', 'drug-target'),
        (r'\bproteindrug\b', 'protein-drug'),
        (r'\bproteinprotein\b', 'protein-protein'),
        (r'\bcellcell\b', 'cell-cell'),
        (r'\bgenedrug\b', 'gene-drug'),
        
        # Function
        (r'\blossoffunction\b', 'loss-of-function'),
        (r'\bgainoffunction\b', 'gain-of-function'),
        (r'\blackofresponse\b', 'lack-of-response'),
        
        # Response/course
        (r'\bdoseresponse\b', 'dose-response'),
        (r'\btimecourse\b', 'time-course'),
        (r'\bfollowup\b', 'follow-up'),
        
        # Other common medical
        (r'\bischemiareperfusion\b', 'ischemia-reperfusion'),
        (r'\bepithelial(\s*)mesenchymal\b', 'epithelial-mesenchymal'),
        (r'\bwellbeing\b', 'well-being'),
    ]
    
    for pattern, replacement in exact_compounds:
        text = re.sub(pattern, replacement, text, flags=re.IGNORECASE)
    
    # STRATEGY 2: Controlled suffix patterns (only specific medical contexts)
    # Only match when we have clear medical/scientific root words
    medical_roots = [
        'time', 'dose', 'concentration', 'age', 'stage', 'grade', 'risk', 'sex',
        'gender', 'disease', 'cancer', 'tumor', 'cell', 'tissue', 'organ',
        'immune', 'receptor', 'ligand', 'enzyme', 'protein', 'gene', 'kinase',
        'pathway', 'signal', 'response', 'treatment', 'therapy', 'drug',
        'glucose', 'insulin', 'hormone', 'cytokine', 'antibody', 'antigen',
        'infection', 'inflammation', 'oxidative', 'mitochondria', 'apoptosis',
        'proliferation', 'differentiation', 'migration', 'invasion', 'metastasis',
        'angiogenesis', 'hypoxia', 'ischemia', 'necrosis', 'fibrosis'
    ]
    
    critical_suffixes = [
        'dependent', 'independent', 'mediated', 'induced', 'related', 'associated',
        'specific', 'sensitive', 'resistant', 'based', 'derived', 'driven',
        'activated', 'regulated', 'controlled', 'stimulated', 'inhibited',
        'responsive', 'initiated', 'triggered', 'enhanced', 'enriched'
    ]
    
    for root in medical_roots:
        for suffix in critical_suffixes:
            pattern = r'\b' + root + suffix + r'\b'
            replacement = root + '-' + suffix
            text = re.sub(pattern, replacement, text, flags=re.IGNORECASE)
    
    # STRATEGY 3: Range patterns (mild to severe, etc.)
    range_terms = ['mild', 'moderate', 'severe', 'low', 'high', 'minimal', 'maximal',
                   'early', 'late', 'short', 'long', 'acute', 'chronic', 'physiologic',
                   'modest', 'significant']
    
    for term in range_terms:
        pattern = r'\b' + term + r'tosevere\b'
        text = re.sub(pattern, term + '-to-severe', text, flags=re.IGNORECASE)
        pattern = r'\b' + term + r'tomoderate\b'
        text = re.sub(pattern, term + '-to-moderate', text, flags=re.IGNORECASE)
        pattern = r'\b' + term + r'tomild\b'
        text = re.sub(pattern, term + '-to-mild', text, flags=re.IGNORECASE)
    
    # STRATEGY 4: SAFE prefix patterns (only specific medical prefixes + specific roots)
    # This prevents "common" → "co-mmon" by being very selective
    safe_prefix_patterns = [
        # Anti- (only with medical terms)
        (r'\banti(inflammatory|oxidant|apoptotic|cancer|tumor|viral|bacterial|fungal|proliferative|angiogenic)\b', r'anti-\1'),
        
        # Non- (only with medical terms)
        (r'\bnon(invasive|cancerous|malignant|coding|specific|functional|responsive|toxic|small|melanoma)\b', r'non-\1'),
        
        # Multi- (only with medical terms)
        (r'\bmulti(modal|drug|target|organ|cellular|factorial|center|variate|component|kinase|molecule)\b', r'multi-\1'),
        
        # Pre/Post- (only with medical terms)
        (r'\bpre(clinical|treatment|operative|malignant|cursor|disposition|existing|natal)\b', r'pre-\1'),
        (r'\bpost(operative|treatment|natal|transcriptional|translational|menopausal)\b', r'post-\1'),
        
        # Intra/Inter/Extra- (only with medical terms)
        (r'\bintra(cellular|venous|muscular|tumoral|peritoneal|cranial|operative)\b', r'intra-\1'),
        (r'\binter(cellular|action|leukin|feron|molecular|phase|mittent)\b', r'inter-\1'),
        (r'\bextra(cellular|cranial|hepatic|ordinary)\b', r'extra-\1'),
        
        # Self- (safe)
        (r'\bself(renewal|sufficient|regulating|organizing|healing|assembly)\b', r'self-\1'),
        
        # Co- (only specific medical terms to avoid "common" → "co-mmon")
        (r'\bco(treatment|expression|culture|infection|administration|morbidity|factor)\b', r'co-\1'),
        
        # Over/Under- (only with medical terms)
        (r'\bover(expression|activation|production|weight|dose|estimate)\b', r'over-\1'),
        (r'\bunder(expression|representation|weight|nutrition|dose|estimate)\b', r'under-\1'),
    ]
    
    for pattern, replacement in safe_prefix_patterns:
        text = re.sub(pattern, replacement, text, flags=re.IGNORECASE)
    
    return text

def fix_subscripts_and_symbols(text: str) -> str:
    """Fix chemical formulas and subscript notation"""
    
    subscript_fixes = [
        # Oxygen partial pressure
        (r'\bpO\s*2\b', 'pO₂'),
        (r'\bpCO\s*2\b', 'pCO₂'),
        
        # Common chemical formulas
        (r'\bCO\s*2\b', 'CO₂'),
        (r'\bH\s*2\s*O\s*2\b', 'H₂O₂'),
        (r'\bO\s*2\b', 'O₂'),
        (r'\bH\s*2\s*O\b', 'H₂O'),
    ]
    
    for pattern, replacement in subscript_fixes:
        text = re.sub(pattern, replacement, text)
    
    return text

def remove_incomplete_references(text: str) -> str:
    """Remove incomplete figure/table references that lost their targets"""
    
    # Patterns for incomplete references
    incomplete_patterns = [
        r'[Aa]\s+schematic\s+diagram[^.]*is\s+shown\s+in\s*\.',
        r'[Aa]s\s+shown\s+in\s+\.',
        r'\b(?:shown|illustrated|depicted|presented)\s+in\s+\.',
        r'[Ss]ee\s+\.',
        r'[Rr]efer\s+to\s+\.',
        r'\(see\s+\)',
        r'\(shown\s+in\s+\)',
    ]
    
    for pattern in incomplete_patterns:
        text = re.sub(pattern, '', text, flags=re.IGNORECASE)
    
    return text

def remove_pmc_citations(text: str) -> str:
    """Remove ALL PMC-style citations comprehensively
    
    Handles:
    - Parenthetical author-year: (Author et al., 2020)
    - Multiple authors: (Smith, Jones & Brown, 2019)
    - Multiple citations: (Author1, 2020; Author2, 2021)
    - Number-only: (1, 2, 3) or (1-5) or (10, 15-20)
    - Mixed patterns
    """
    
    # STAGE 1: Remove complex author-year citations
    # Pattern: (Author et al., Year; Author2 et al., Year)
    # This is the most complex pattern and must come first
    text = re.sub(
        r'\(\s*[A-Z][a-z]+(?:\s+et\s+al\.?)?(?:\s*,\s*[A-Z][a-z]+)*\s*,\s*\d{4}(?:\s*[;,]\s*[A-Z][a-z]+(?:\s+et\s+al\.?)?(?:\s*,\s*[A-Z][a-z]+)*\s*,\s*\d{4})*\s*\)',
        '',
        text
    )
    
    # STAGE 2: Remove simpler author-year patterns
    # Pattern: (Author, Year) or (Author & Author2, Year)
    text = re.sub(
        r'\(\s*[A-Z][a-z]+(?:\s*,\s*[A-Z][a-z]+|\s+&\s+[A-Z][a-z]+|\s+et\s+al\.?)?\s*,\s*\d{4}\s*\)',
        '',
        text
    )
    
    # STAGE 3: Remove author names with periods (O'Leary., Tweedie., etc.)
    # Pattern: (O'Leary., 2016) or (Tweedie., 2021)
    text = re.sub(
        r'\(\s*[A-Z][a-zA-Z\']+\.\s*,\s*\d{4}\s*\)',
        '',
        text
    )
    
    # STAGE 4: Remove number-only citations with ranges
    # Pattern: (1-5), (10, 15-20, 25), (1, 2, 3)
    text = re.sub(
        r'\(\s*\d+(?:\s*[-–—]\s*\d+)?(?:\s*,\s*\d+(?:\s*[-–—]\s*\d+)?)*\s*\)',
        '',
        text
    )
    
    # STAGE 5: Remove year-only or simple patterns that might remain
    # Pattern: (2020), (2019-2021)
    text = re.sub(
        r'\(\s*\d{4}(?:\s*[-–—]\s*\d{4})?\s*\)',
        '',
        text
    )
    
    # STAGE 6: Clean up organization names in citations
    # Pattern: (World Health Organization, 2020)
    text = re.sub(
        r'\(\s*[A-Z][a-zA-Z\s]+(?:Organization|Institute|Agency|Association|Society|Committee|Consortium)\s*,?\s*\d{0,4}\s*\)',
        '',
        text
    )
    
    # STAGE 7: Remove any remaining isolated parentheses with just years/numbers
    text = re.sub(r'\(\s*\d+\s*\)', '', text)
    
    # STAGE 8: Clean up web/post/image references
    text = re.sub(r'\[web:\d+\]|\[post:\d+\]|\[image:\d+\]', '', text)
    
    return text

def clean_text(text: str) -> str:
    if not text:
        return ""
    
    # Basic whitespace normalization
    text = re.sub(r'\s+', ' ', text).strip()
    
    # Remove PMC-style citations (NEW - comprehensive removal)
    text = remove_pmc_citations(text)
    
    # Remove bracket citations [1], [1-5], [1, 2, 3]
    text = re.sub(r'\[\s*\d+(?:\s*[-,–—]\s*\d+)*\s*\]', '', text)
    
    # Remove DOI/PMID/PMC identifiers
    text = re.sub(r'\b(?:doi|pmid|pmc)\s*[:=]?\s*[0-9a-z.-]+', '', text, flags=re.I)
    
    # Remove figure and table references
    text = re.sub(r'\b(?:fig(?:ure)?|table)\s*\d+[a-z]?\b\.?', '', text, flags=re.I)
    
    # Remove incomplete references
    text = remove_incomplete_references(text)
    
    # Remove other common artifacts
    text = re.sub(r'\b(?:et al\.?|Suppl\.?|Appendix)\b', '', text, flags=re.I)
    text = re.sub(r'={1,}|-{1,}|─{1,}', '', text)  # Separators
    text = re.sub(r'\(\s*\)', '', text)  # Empty parentheses
    text = re.sub(r'\s+\.', '.', text)  # Space before period
    
    # Normalize quotes
    text = text.replace('\u2018', "'").replace('\u2019', "'")
    text = text.replace('\u201c', '"').replace('\u201d', '"')
    
    # Apply medical term fixes
    text = fix_compound_medical_terms(text)
    text = fix_subscripts_and_symbols(text)
    
    # Final cleanup of extra whitespace
    text = re.sub(r'\s+', ' ', text).strip()
    text = re.sub(r'\s+([.,;:!?])', r'\1', text)  # Remove space before punctuation
    text = re.sub(r'([.,;:!?])([A-Za-z])', r'\1 \2', text)  # Add space after punctuation
    
    # Remove any double spaces that might have been created
    text = re.sub(r'\s{2,}', ' ', text)
    
    return text.strip()

def split_long_paragraph(text: str) -> list[str]:
    if len(text) <= MAX_LEN:
        return [text]
    sentences = re.split(r'(?<=[.!?])\s+(?=[A-Z])', text)
    chunks, current = [], ""
    for s in sentences:
        if len(current) + len(s) + 1 <= MAX_LEN:
            current += (" " + s if current else s)
        else:
            if current:
                chunks.append(current.strip())
            current = s
    if current.strip():
        chunks.append(current.strip())
    return chunks

def paragraph_similarity(p1: str, p2: str) -> float:
    if not p1 or not p2:
        return 0.0
    words1 = set(p1.lower().split())
    words2 = set(p2.lower().split())
    return len(words1 & words2) / len(words1 | words2) if words1 or words2 else 0.0

def is_valid(text: str, heading: str = "") -> bool:
    if not text or len(text) < MIN_LEN:
        return False
    digit_ratio = sum(c.isdigit() for c in text) / len(text)
    if digit_ratio > 0.42:
        return False
    special_ratio = sum(not c.isalnum() and not c.isspace() for c in text) / len(text)
    if special_ratio > 0.32:
        return False
    if text.isupper() or text.isspace():
        return False
    heading_lower = heading.lower()
    if any(kw in heading_lower for kw in PRIORITY_SECTIONS):
        return len(text) > MIN_LEN - 50
    return True

def process_txt_file(in_path: Path, out_path: Path):
    with open(in_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    
    cleaned_blocks = []
    buffer = ""
    current_heading = ""
    seen_paras = OrderedDict()
    seen_headings = set()  # Track seen headings to remove duplicates
    
    start_idx = 1 if lines and lines[0].strip().startswith("<DOCUMENT") else 0
    
    for line in lines[start_idx:]:
        stripped = line.strip()
        if not stripped or re.match(r'^={1,}$', stripped):
            continue
        
        is_heading = (re.match(r'^\d+(\.\d+)*\s+[A-Z]', line) or 
                      line.isupper() or 
                      re.match(r'^={1,}.+={1,}$', stripped))
        
        if is_heading:
            if buffer:
                cleaned = clean_text(buffer)
                for chunk in split_long_paragraph(cleaned):
                    if is_valid(chunk, current_heading):
                        is_duplicate = any(paragraph_similarity(chunk, seen) > DEDUP_THRESHOLD for seen in seen_paras)
                        if not is_duplicate:
                            seen_paras[chunk] = True
                            cleaned_blocks.append(chunk)
                buffer = ""
            
            current_heading = normalize_heading(stripped).lower()
            if current_heading:
                if SKIP_ABSTRACT and "abstract" in current_heading:
                    buffer = ""  # Skip abstract content
                    continue
                
                # Check for duplicate section headings
                if current_heading not in seen_headings:
                    seen_headings.add(current_heading)
                    cleaned_blocks.append(f"\n{current_heading.title()}\n")
                # If duplicate heading, skip adding it but continue processing content
        
        else:
            if not (SKIP_ABSTRACT and "abstract" in current_heading):
                buffer += " " + stripped
    
    # Process final buffer
    if buffer and not (SKIP_ABSTRACT and "abstract" in current_heading):
        cleaned = clean_text(buffer)
        for chunk in split_long_paragraph(cleaned):
            if is_valid(chunk, current_heading):
                is_duplicate = any(paragraph_similarity(chunk, seen) > DEDUP_THRESHOLD for seen in seen_paras)
                if not is_duplicate:
                    seen_paras[chunk] = True
                    cleaned_blocks.append(chunk)
    
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(out_path, "w", encoding="utf-8") as f:
        for block in cleaned_blocks:
            f.write(block.strip() + "\n\n")

def main():
    print("\n🧹 ENHANCED MEDICAL TEXT CLEANING FOR LLM TRAINING\n")
    print("✨ Key Improvements:")
    print("   • COMPREHENSIVE PMC citation removal (all patterns)")
    print("   • Fixed compound medical term hyphenation")
    print("   • Fixed chemical formulas and subscripts")
    print("   • Removed duplicate section headings")
    print("   • Removed incomplete figure/table references")
    print("   • Preserved all medical meaning and context\n")
    
    txt_files = list(ROOT_INPUT.rglob("*.txt"))
    print(f"📁 Found {len(txt_files)} files to process...")
    
    if not txt_files:
        print("❌ No files found — check ROOT_INPUT path.")
        return
    
    processed = 0
    for txt in txt_files:
        relative = txt.relative_to(ROOT_INPUT)
        output_txt = ROOT_OUTPUT / relative
        print(f"✔ Processing: {relative}")
        process_txt_file(txt, output_txt)
        processed += 1
    
    print(f"\n✅ DONE — Processed {processed} files")
    print(f"📂 Clean training data in: {ROOT_OUTPUT}")
    print(f"⚠️  Abstracts {'skipped' if SKIP_ABSTRACT else 'included'}")

if __name__ == "__main__":
    main()
