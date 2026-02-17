import os
import re
from pathlib import Path
from collections import OrderedDict
import unicodedata
from concurrent.futures import ProcessPoolExecutor
# ============================================================================
# UNIVERSAL MEDICAL TEXT CLEANER FOR LLM TRAINING
# Optimized for 100,000+ PMC documents across ALL diseases
# ============================================================================

# CONFIG
ROOT_INPUT = Path("Cancer_LLM_Data/Prostate_Cancer")  # Your input folder
ROOT_OUTPUT = Path("Cleaned_Cancer_LLM_Data/Protaste_Cancer")  # Parallel cleaned
MIN_LEN = 150
MAX_LEN = 4200
DEDUP_THRESHOLD = 0.95
SKIP_ABSTRACT = True  # True: Skips "Abstract" sections

PRIORITY_SECTIONS = {
    "background", "introduction", "results", "discussion", 
    "conclusions", "methods", "summary", "findings",
    # Add any disease-specific sections you want to prioritize
}

def extract_pmc_id(filename: str) -> str:
    """Extract PMC ID from filename"""
    match = re.search(r'PMC(\d+)', filename)
    return match.group(0) if match else "Unknown"

def normalize_heading(line: str) -> str:
    """Normalize section headings"""
    # REMOVE ALL NUMBERING PATTERNS: 1., 1.2., 1.2.2., etc.
    line = re.sub(r'^\s*\d+(\.\d+)*\.?\s+', '', line.strip())
    line = re.sub(r'^={1,}\s*(.+?)\s*={1,}\s*$', r'\1', line)
    line = re.sub(r'^={1,}$', '', line)
    if not line.strip():
        return ""
    if len(line) < 80 and line.isupper():
        line = line.title()
    return line.strip()

def remove_pmc_citations_universal(text: str) -> str:
    """UNIVERSAL citation remover for ALL PMC documents
    
    Handles ALL patterns found across diverse medical literature:
    - Author-year citations (all formats)
    - Number citations (all formats) 
    - Figure/Table/Supplementary references
    - RRIDs, catalog numbers, strain IDs
    - Sample sizes, statistical notations
    - And 30+ other PMC-specific patterns
    """
    
    # ==========================================================================
    # STAGE 1: PMC-SPECIFIC IDENTIFIERS (Highest Priority)
    # ==========================================================================
    
    # Remove RRID identifiers: (RRID: CVCL_VI66), RRID: AB_123456, etc.
    text = re.sub(r'\(?\s*RRID:\s*[A-Z_]+[A-Z0-9_:]+\s*[,)]?', '', text, flags=re.IGNORECASE)
    
    # Remove catalog/clone numbers: (catalog no. 16-1471-82), (clone #5B8)
    text = re.sub(r'\(\s*(?:catalog|cat\.?|clone)\s+(?:no\.?|#)?\s*[0-9A-Z-]+(?:\s*,\s*RRID[^)]+)?\s*\)', '', text, flags=re.IGNORECASE)
    text = re.sub(r'(?:catalog|cat\.?|clone)\s+(?:no\.?|#)?\s*[0-9A-Z-]+', '', text, flags=re.IGNORECASE)
    
    # Remove strain/stock identifiers: (strain # 029015), (stock #12345)
    text = re.sub(r'\(\s*(?:strain|stock)\s*#?\s*[0-9A-Z-]+\s*\)', '', text, flags=re.IGNORECASE)
    
    # Remove accession numbers: (accession no. GSE12345), (GenBank: AB123456)
    text = re.sub(r'\(\s*(?:accession|GenBank|NCBI|GEO)\s*(?:no\.?|ID)?:?\s*[A-Z0-9_]+\s*\)', '', text, flags=re.IGNORECASE)
    text = re.sub(r'(?:accession|GenBank|NCBI|GEO)\s*(?:no\.?|ID)?:?\s*[A-Z0-9_]+', '', text, flags=re.IGNORECASE)
    
    # ADDED: GEO/SRA dataset identifiers: (GSE80632), GSE90784, SRA200820, (: GSE80632)
    text = re.sub(r'\(\s*:?\s*[GS][SE][EO]\d{4,}\s*\)', '', text)  # (GSE12345), (: GSE12345)
    text = re.sub(r'\b[GS][SE][EO]\d{4,}\b', '', text)  # GSE12345, SRA200820
    
    # ADDED: Experiment IDs: (ReMap2022, Experiment ID: ENCSR226NRS)
    text = re.sub(r'\([^)]*Experiment\s+ID:\s*[A-Z0-9]+[^)]*\)', '', text, flags=re.IGNORECASE)
    text = re.sub(r'Experiment\s+ID:\s*[A-Z0-9]+', '', text, flags=re.IGNORECASE)
    
    # ADDED: Clinical Trial IDs: NCT03906331, (NCT01644890)
    text = re.sub(r'\(?\s*NCT\d{8,}\s*\)?', '', text, flags=re.IGNORECASE)
    text = re.sub(r'\bNCT\d{8,}\b', '', text, flags=re.IGNORECASE)
    
    # ==========================================================================
    # STAGE 2: REFERENCE CITATIONS (Multiple Formats)
    # ==========================================================================
    
    # Reference citations: (refs. 1-5), (ref. 10), (references 1, 2), (see ref 5)
    text = re.sub(r'\(\s*(?:see\s+)?refs?\.?\s*\d+(?:\s*[-–—,]\s*\d+)*\s*\)', '', text, flags=re.IGNORECASE)
    text = re.sub(r'\b(?:see\s+)?refs?\.?\s*\d+(?:\s*[-–—,]\s*\d+)*\b', '', text, flags=re.IGNORECASE)
    
    # ==========================================================================
    # STAGE 3: FIGURE, TABLE, AND SUPPLEMENTARY REFERENCES
    # ==========================================================================
    
    # Single letter panel references: ( A ), ( B ), ( C-E ), (A, B)
    text = re.sub(r'\(\s*[A-Z](?:\s*[-–—,]\s*[A-Z])*\s*\)', '', text)
    
    # Figure references with all variations:
    # (Fig. 1), (Figure 1A), (Figs 1-3), (figure 1a-c)
    text = re.sub(
        r'\(\s*Figs?\.?\s*[S]?[0-9]+[A-Za-z]?(?:\s*[-–—,]\s*[S]?[0-9]+[A-Za-z]?)*\s*\)',
        '',
        text,
        flags=re.IGNORECASE
    )
    
    # Table references: (Table 1), (Tables 1-3), (table S1)
    text = re.sub(
        r'\(\s*Tables?\.?\s*[S]?[0-9]+[A-Za-z]?(?:\s*[-–—,]\s*[S]?[0-9]+[A-Za-z]?)*\s*\)',
        '',
        text,
        flags=re.IGNORECASE
    )
    
    # Supplementary references: (Supplementary Fig. S1), (Suppl. Table 1)
    text = re.sub(
        r'\(\s*(?:Supplementary|Suppl?\.?)\s+(?:Fig(?:ure)?|Table|Data|Material|Movie|Video)\.?\s*[S]?[0-9]+[A-Za-z]?(?:\s*[-–—,]\s*[S]?[0-9]+[A-Za-z]?)?\s*\)',
        '',
        text,
        flags=re.IGNORECASE
    )
    
    # ADDED: Non-parenthetical supplemental/supplementary references
    # "supplemental Table S1", "supplementary Figure S2", etc.
    text = re.sub(
        r'\bsupplemental(?:ary)?\s+(?:Fig(?:ure)?|Table|Data|Material|Movie|Video)\.?\s*[S]?[0-9]+[A-Za-z]?',
        '',
        text,
        flags=re.IGNORECASE
    )
    
    # Appendix references: (Appendix 1), (Appendix A)
    text = re.sub(r'\(\s*Appendix\s+[0-9A-Z]+\s*\)', '', text, flags=re.IGNORECASE)
    
    # Box/Scheme references: (Box 1), (Scheme 2)
    text = re.sub(r'\(\s*(?:Box|Scheme)\s+[0-9]+\s*\)', '', text, flags=re.IGNORECASE)
    
    # ==========================================================================
    # STAGE 3.5: COMPLEX FIGURE PANEL REFERENCES (MISSED PATTERNS)
    # ==========================================================================
    
    # ADDED: Figure references with mixed formats: ( A, 2B), ( E and supplemental Figure S 2 )
    # Pattern: (letter, number+letter), (letter and supplemental...)
    text = re.sub(
        r'\(\s*[A-Z]\s*,\s*\d+[A-Z]\s*\)',
        '',
        text
    )
    
    # ADDED: Complex references mixing panels and supplemental: ( B, supplemental Table S 5 )
    text = re.sub(
        r'\(\s*[A-Z]\s*,\s*supplemental(?:ary)?\s+(?:Fig(?:ure)?|Table)[^)]*\)',
        '',
        text,
        flags=re.IGNORECASE
    )
    
    # ADDED: References with "and": ( E and supplemental Figure S 2 )
    text = re.sub(
        r'\(\s*[A-Z]\s+and\s+supplemental(?:ary)?[^)]*\)',
        '',
        text,
        flags=re.IGNORECASE
    )
    
    # ==========================================================================
    # STAGE 4: AUTHOR-YEAR CITATIONS (All Formats)
    # ==========================================================================
    
    # Complex multi-author citations with semicolons:
    # (Smith et al., 2020; Jones & Brown, 2021; White, 2019)
    text = re.sub(
        r'\(\s*[A-Z][a-zA-Z\'-]+(?:\s+(?:et\s+al\.?|and|&)\s+[A-Z][a-zA-Z\'-]+)*\s*,\s*\d{4}[a-z]?(?:\s*[;,]\s*[A-Z][a-zA-Z\'-]+(?:\s+(?:et\s+al\.?|and|&)\s+[A-Z][a-zA-Z\'-]+)*\s*,\s*\d{4}[a-z]?)*\s*\)',
        '',
        text
    )
    
    # Simple author-year: (Smith, 2020), (Smith & Jones, 2020), (Smith et al., 2020)
    text = re.sub(
        r'\(\s*[A-Z][a-zA-Z\'-]+(?:\s+(?:et\s+al\.?|and|&)\s+[A-Z][a-zA-Z\'-]+)*\s*,\s*\d{4}[a-z]?\s*\)',
        '',
        text
    )
    
    # ADDED: et al. citations without parentheses (common in text)
    # "Smith et al. 2020", "Jones et al., 2021", "Brown et al. (2019)"
    text = re.sub(
        r'\b[A-Z][a-zA-Z\'-]+\s+et\s+al\.?\s*,?\s*\(?\d{4}[a-z]?\)?',
        '',
        text
    )
    
    # ADDED: Standalone "et al." with comma or period
    text = re.sub(r'\bet\s+al\.?\s*[,;]', '', text)
    
    # ADDED: Author name + "and colleagues" pattern: "Xiong J and colleagues"
    text = re.sub(r'\b[A-Z][a-z]+\s+[A-Z]\s+and\s+colleagues\b', '', text)
    text = re.sub(r'\b[A-Z][a-z]+\s+and\s+colleagues\b', '', text)
    
    # ==========================================================================
    # STAGE 4.5: SMART AUTHOR CITATION REMOVAL (MEDICAL TEXT SPECIFIC)
    # ==========================================================================
    
    # CRITICAL: Remove author citations WITHOUT removing medical terms
    # Pattern: "Name et al." - Most common in medical literature
    text = re.sub(r'\b[A-Z][a-z]+(?:\s+[A-Z][a-z]+)?\s+et\s+al\.?', '', text)
    
    # Pattern: "LastName, FirstInitial et al." - e.g., "Smith, J et al."
    text = re.sub(r'\b[A-Z][a-z]+\s*,\s*[A-Z]\.?\s+et\s+al\.?', '', text)
    
    # Pattern: "Name and Name" citations - e.g., "Smith and Jones", "Wong and Chiu"
   # Only remove "Smith and Jones" when near a year or citation context
    text = re.sub(
    r'\b[A-Z][a-z]+\s+and\s+[A-Z][a-z]+(?=\s*,?\s*\d{4})',
    '',
    text
)

    
    # Pattern: "Name, Name and Name" - e.g., "Smith, Brown and Jones"
    text = re.sub(r'\b[A-Z][a-z]+\s*,\s*[A-Z][a-z]+\s+and\s+[A-Z][a-z]+\b', '', text)
    
    # Pattern: Single author citations at end of sentence - e.g., "as shown by Smith."
    text = re.sub(r'\bby\s+[A-Z][a-z]+(?:\s+[A-Z][a-z]+)?\s*\.', 'by researchers.', text)
    
    # Pattern: "According to Name" or "as Name showed"
    text = re.sub(r'\baccording to\s+[A-Z][a-z]+(?:\s+[A-Z][a-z]+)?\b', 'according to researchers', text, flags=re.IGNORECASE)
    
    # Pattern: "Name reported", "Name demonstrated", "Name showed"
    text = re.sub(r'\b[A-Z][a-z]+(?:\s+[A-Z][a-z]+)?\s+(?:reported|demonstrated|showed|found|observed|noted|described)\b', 'researchers reported', text)
    
    # Author names with dots: (O'Leary., 2016), (van der Berg., 2020)
    text = re.sub(
        r'\(\s*(?:van\s+der\s+|van\s+|de\s+)?[A-Z][a-zA-Z\']+\.+\s*,\s*\d{4}[a-z]?\s*\)',
        '',
        text,
        flags=re.IGNORECASE
    )
    
    # Organization as author: (World Health Organization, 2020), (CDC, 2021)
    text = re.sub(
        r'\(\s*(?:[A-Z][a-zA-Z]+\s+)*(?:Organization|Institute|Agency|Association|Society|Committee|Consortium|Foundation|Administration|Department|CDC|WHO|NIH|FDA)\s*,?\s*\d{0,4}[a-z]?\s*\)',
        '',
        text
    )
    
    # ==========================================================================
    # STAGE 5: NUMBER-ONLY CITATIONS (All Formats) - ENHANCED
    # ==========================================================================
    
    # CRITICAL: Square bracket citations - ALL VARIATIONS
    # [23], [ 23 ], [1-5], [1,2,3], [10, 15-20], [23,24,25]
    text = re.sub(
        r'\[\s*\d+(?:\s*[-–—]\s*\d+)?(?:\s*,\s*\d+(?:\s*[-–—]\s*\d+)?)*\s*\]',
        '',
        text
    )
    
    # ADDED: Multiple sequential square brackets: [1][2][3]
    text = re.sub(r'(?:\[\s*\d+\s*\])+', '', text)
    
    # Number citations with ranges and lists: (1-5), (10, 15-20, 25), (1,2,3)
    text = re.sub(
        r'\(\s*\d+(?:\s*[-–—]\s*\d+)?(?:\s*[,;]\s*\d+(?:\s*[-–—]\s*\d+)?)*\s*\)',
        '',
        text
    )
    
    # ADDED: Multiple sequential parenthetical numbers: (1)(2)(3)
    text = re.sub(r'(?:\(\s*\d+\s*\))+', '', text)
    
    # Superscript-style citations: ¹, ², ³⁻⁵
    text = re.sub(r'[¹²³⁴⁵⁶⁷⁸⁹⁰]+(?:[-–—][¹²³⁴⁵⁶⁷⁸⁹⁰]+)?', '', text)
    
    # Year-only citations: (2020), (2019-2021), (2020a,b)
    text = re.sub(r'\(\s*\d{4}[a-z]?(?:\s*[-–—,]\s*\d{4}[a-z]?)*\s*\)', '', text)
    
    # ==========================================================================
    # STAGE 6: STATISTICAL AND SAMPLE SIZE NOTATIONS
    # ==========================================================================
    
    # Sample sizes: (n = 32), (n=5), n = 100, (N=50)
    text = re.sub(r'\(\s*[Nn]\s*=\s*\d+(?:\s+per\s+group)?\s*\)', '', text)
    text = re.sub(r'\b[Nn]\s*=\s*\d+\b', '', text)
    
    # P-values in isolation: (p < 0.05), (P=0.001)
    text = re.sub(r'\(\s*[Pp]\s*[<>=≤≥]\s*0?\.\d+\s*\)', '', text)
    
    # Percentage citations: (50%), (10-20%)
    # Note: Keep these ONLY if they appear to be data, not citations
    # text = re.sub(r'\(\s*\d+(?:\.\d+)?%\s*\)', '', text)
    
    # Confidence intervals: (95% CI), (CI: 1.2-3.4)
    text = re.sub(r'\(\s*\d+%\s*(?:CI|confidence interval)[^)]*\)', '', text, flags=re.IGNORECASE)
    
    # ==========================================================================
    # STAGE 7: URL, DOI, PMID, AND ELECTRONIC REFERENCES
    # ==========================================================================
    
    # DOI identifiers: doi:10.1234/xyz, (DOI: 10.1234/xyz)
    text = re.sub(r'\(?\s*DOI:?\s*10\.\d{4,}/[^\s)]+\s*\)?', '', text, flags=re.IGNORECASE)
    
    # PMID/PMC identifiers: PMID: 12345678, PMC12345678
    text = re.sub(r'\b(?:PMID|PMC)\s*:?\s*\d{6,}\b', '', text, flags=re.IGNORECASE)
    
    # URLs (but preserve domain names that are important)
    text = re.sub(r'https?://[^\s)]+', '', text)
    
    # Email addresses in parentheses (author correspondence)
    text = re.sub(r'\([a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}\)', '', text)
    
    # ==========================================================================
    # STAGE 8: INCOMPLETE REFERENCES AND PHRASES
    # ==========================================================================
    
    # "as shown in", "see", "refer to" without target
    incomplete_patterns = [
        r'[Aa]s\s+(?:shown|illustrated|depicted|presented|described)\s+in\s*[,.]',
        r'[Ss]ee\s+(?:also\s+)?[,.]',
        r'[Rr]efer\s+to\s+[,.]',
        r'\(\s*see\s+\)',
        r'\(\s*shown\s+in\s+\)',
        r'[Aa]\s+schematic\s+(?:diagram|representation)[^.]*is\s+shown\s+in\s*[,.]',
        # ADDED: More incomplete phrases
        r'\bas\s+previously\s+described\b',
        r'\bas\s+described\s+(?:previously|above|earlier|in\s+Methods)\b',
        r'\baccording\s+to\s+(?:the\s+)?manufacturer\'?s?\s+(?:instructions|protocol)\b',
    ]
    
    for pattern in incomplete_patterns:
        text = re.sub(pattern, '', text, flags=re.IGNORECASE)
    
    # ==========================================================================
    # STAGE 9: COMMON ARTIFACTS
    # ==========================================================================
    
    # "et al.", "i.e.", "e.g." without context (cleanup)
    text = re.sub(r'\bet\s+al\.?\s*[,;]?\s*$', '', text)  # Only at end of sentences
    
    # "Suppl.", "Appendix" standalone
    text = re.sub(r'\b(?:Suppl?\.?|Appendix)\s*[,;]?\s*$', '', text, flags=re.IGNORECASE)
    
    # Empty parentheses: ( ), (  )
    text = re.sub(r'\(\s*\)', '', text)
    
    # ADDED: Empty square brackets: [ ], [  ]
    text = re.sub(r'\[\s*\]', '', text)
    
    # Multiple punctuation marks: .. , ;;
    text = re.sub(r'\.{2,}', '.', text)
    text = re.sub(r',{2,}', ',', text)
    text = re.sub(r';{2,}', ';', text)
    
    # Separators: ===, ---, ─────
    text = re.sub(r'[=\-─]{3,}', '', text)
    
    # Space before punctuation: "word .", "word ,"
    text = re.sub(r'\s+([.,;:!?])', r'\1', text)
    
    # ==========================================================================
    # STAGE 10: CHEMICAL FORMULAS AND SUBSCRIPTS
    # ==========================================================================
    
    # Fix common chemical formulas (keep subscripts)
    subscript_map = {
        r'\bpO\s*2\b': 'pO₂',
        r'\bpCO\s*2\b': 'pCO₂',
        r'\bCO\s*2\b': 'CO₂',
        r'\bH\s*2\s*O\s*2\b': 'H₂O₂',
        r'\bO\s*2\b': 'O₂',
        r'\bH\s*2\s*O\b': 'H₂O',
        r'\bCa\s*2\s*\+': 'Ca²⁺',
        r'\bZn\s*2\s*\+': 'Zn²⁺',
    }
    
    for pattern, replacement in subscript_map.items():
        text = re.sub(pattern, replacement, text)
    text = text.replace('®', '').replace('™', '')
    text = re.sub(r'\bresearchers\'?s?\s+(?:instructions|protocol|reported)\b', '', text, flags=re.IGNORECASE)
    text = re.sub(r'\baccording to researchers\'?s?\b', '', text, flags=re.IGNORECASE)

    # ==========================================================================
    # STAGE 1: Dataset/Plasmid/IRB/Clinical Trial IDs
    # ==========================================================================
    text = re.sub(r'\b(?:GSE|SRA|EGA|E-MTAB|MSV)\d{4,}\b', '', text)
    text = re.sub(r'\bAddgene\s+plasmid\s*#?\s*[\d,]+\b', '', text, flags=re.IGNORECASE)
    text = re.sub(r'\b(?:NCT|CIRB\s+Ref|reference\s+number)[\d/:]+\b', '', text, flags=re.IGNORECASE)

    # ==========================================================================
    # STAGE 2: Antibody/Catalog/Company/Location blocks (aggressive & safe)
    # ==========================================================================
    text = re.sub(
        r'\(\s*[A-Za-z0-9\s&,.;#\-AP™®]+(?:Biotechnology|Technology|Biosciences|Bioscience|Proteintech|Pharmingen|DAKO|Akoya)\s*,?\s*[A-Za-z0-9#\-AP™®]+\s*(?:,\s*1:\d+)?\s*(?:,\s*[A-Za-z\s]+,\s*[A-Z]{2}(?:,\s*USA)?)?\s*\)',
        '',
        text,
        flags=re.IGNORECASE
    )
    text = re.sub(r'\(\s*[A-Za-z0-9#\-AP™®]+\s*\)', '', text)  # Standalone catalog

    # ==========================================================================
    # STAGE 3: Figure/Table/Supplemental/Appendix/Box/Scheme
    # ==========================================================================
    text = re.sub(r'\(\s*Figs?\.?\s*[S]?\d+[A-Za-z]?(?:\s*[-–,]\s*\d+[A-Za-z]?)*\s*\)', '', text, flags=re.IGNORECASE)
    text = re.sub(r'\(\s*Tables?\.?\s*[S]?\d+[A-Za-z]?(?:\s*[-–,]\s*\d+[A-Za-z]?)*\s*\)', '', text, flags=re.IGNORECASE)
    text = re.sub(r'\(\s*(?:Supplementary|Supplemental|Suppl\.?)\s+(?:Fig\.?|Figure|Table|Data|Material|Movie|Video)[^)]*\)', '', text, flags=re.IGNORECASE)
    text = re.sub(r'\bsupplement(?:al|ary)?\s+(?:Fig\.?|Figure|Table|Data)\s*[S]?\d+[A-Za-z]?', 'supplementary data', text, flags=re.IGNORECASE)
    text = re.sub(r'\(\s*Appendix\s+[A-Z0-9]+\s*\)', '', text, flags=re.IGNORECASE)
    text = re.sub(r'\(\s*(?:Box|Scheme)\s+\d+\s*\)', '', text, flags=re.IGNORECASE)

    # ==========================================================================
    # STAGE 4: Author-Year Citations (all formats)
    # ==========================================================================
    text = re.sub(r'\(\s*[A-Z][a-zA-Z\'-]+(?:\s+(?:et\s+al\.?|&|and)\s*[A-Z][a-zA-Z\'-]+)*\s*,\s*\d{4}[a-z]?(?:\s*[;,][^)]*)*\s*\)', '', text)
    text = re.sub(r'\b[A-Z][a-zA-Z\'-]+\s+et\s+al\.?\s*,?\s*\d{4}[a-z]?\b', '', text)
    text = re.sub(r'\bet\s+al\.?\b', '', text)

    # ==========================================================================
    # STAGE 5: Numbered Citations (bracketed, parenthetical, superscripts)
    # ==========================================================================
    text = re.sub(r'\[\s*\d+(?:\s*[-–,]\s*\d+)?(?:\s*,\s*\d+(?:\s*[-–,]\s*\d+)?)*\s*\]', '', text)
    text = re.sub(r'\(\s*\d+(?:\s*[-–,]\s*\d+)?(?:\s*[,;]\s*\d+(?:\s*[-–,]\s*\d+)?)*\s*\)', '', text)
    text = re.sub(r'[¹²³⁴⁵⁶⁷⁸⁹⁰]+(?:[-–][¹²³⁴⁵⁶⁷⁸⁹⁰]+)?', '', text)

    # ==========================================================================
    # STAGE 6: SMART Inline Number Removal (preserves real stats)
    # ==========================================================================
    # Only removes bare citation-like numbers (e.g., "1, 2.", "E-09") 
    # NOT percentages, p-values, ages, grades, etc.
    for _ in range(6):
        text = re.sub(r'(?<![\d.%<=>])\s*\d+(?:\s*(?:,\s*|-|–|—|and)\s*\d+)*\s*[\.,;]?(?![\d%])', ' ', text)

    # ==========================================================================
    # STAGE 7: DOI, PMID, PMC, URLs
    # ==========================================================================
    text = re.sub(r'\(?\s*DOI:?\s*10\.\d{4,}/[^\s)]+\s*\)?', '', text, flags=re.IGNORECASE)
    text = re.sub(r'\b(?:PMID|PMC)\s*:\s*\d{6,}\b', '', text, flags=re.IGNORECASE)
    text = re.sub(r'https?://[^\s]+', '', text)

    # ==========================================================================
    # STAGE 8: Final cleanup
    # ==========================================================================
    text = re.sub(r'\s{2,}', ' ', text)
    text = re.sub(r'\s+([.,;:!?])', r'\1', text)
    text = re.sub(r'\(\s*\)', '', text)
    text = re.sub(r'\[\s*\]', '', text)
    text = text.strip()

    # ==========================================================================
    # STAGE 9: Chemical subscripts (preserve chemistry)
    # ==========================================================================
    subscript_map = {
        r'\bpO\s*2\b': 'pO₂', r'\bpCO\s*2\b': 'pCO₂', r'\bCO\s*2\b': 'CO₂',
        r'\bH\s*2\s*O\s*2\b': 'H₂O₂', r'\bO\s*2\b': 'O₂', r'\bH\s*2\s*O\b': 'H₂O',
        r'\bCa\s*2\s*\+': 'Ca²⁺', r'\bZn\s*2\s*\+': 'Zn²⁺',
    }
    for pattern, replacement in subscript_map.items():
        text = re.sub(pattern, replacement, text)
    
    return text

def fix_compound_medical_terms(text: str) -> str:
    """Fix missing hyphens in compound medical/scientific terms
    
    Safe for 100K+ documents across all diseases
    """
    
    # EXACT medical compound terms (zero false positives)
    exact_compounds = [
        # Common across all diseases
        (r'\bhypoxiainducible\b', 'hypoxia-inducible'),
        (r'\bhypoxiadependent\b', 'hypoxia-dependent'),
        (r'\bplacebocontrolled\b', 'placebo-controlled'),
        (r'\bdoubleblind\b', 'double-blind'),
        (r'\bsingleblind\b', 'single-blind'),
        (r'\bcrossover\b', 'cross-over'),
        (r'\bfollowup\b', 'follow-up'),
        (r'\bdoseresponse\b', 'dose-response'),
        (r'\btimecourse\b', 'time-course'),
        (r'\bgenomewide\b', 'genome-wide'),
        (r'\bwellbeing\b', 'well-being'),
        (r'\blossoffunction\b', 'loss-of-function'),
        (r'\bgainoffunction\b', 'gain-of-function'),
        (r'\bproteinprotein\b', 'protein-protein'),
        (r'\bcellcell\b', 'cell-cell'),
        
        # DNA/RNA/UV
        (r'\bDNAdamage\b', 'DNA damage'),
        (r'\bDNArepair\b', 'DNA repair'),
        (r'\bRNAsequencing\b', 'RNA sequencing'),
        (r'\bUVinduced\b', 'UV-induced'),
        
        # Clinical/pathological
        (r'\bepithelial\s*mesenchymal\b', 'epithelial-mesenchymal'),
        (r'\bischemiareperfusion\b', 'ischemia-reperfusion'),
    ]
    
    for pattern, replacement in exact_compounds:
        text = re.sub(pattern, replacement, text, flags=re.IGNORECASE)
    
    # Medical roots + critical suffixes pattern
    medical_roots = [
        'time', 'dose', 'age', 'stage', 'grade', 'risk',
        'disease', 'cancer', 'tumor', 'cell', 'tissue',
        'immune', 'receptor', 'protein', 'gene',
        'treatment', 'therapy', 'drug', 'infection'
    ]
    
    suffixes = [
        'dependent', 'independent', 'mediated', 'induced', 
        'related', 'associated', 'specific', 'based',
        'derived', 'driven', 'activated', 'regulated'
    ]
    
    for root in medical_roots:
        for suffix in suffixes:
            pattern = r'\b' + root + suffix + r'\b'
            replacement = root + '-' + suffix
            text = re.sub(pattern, replacement, text, flags=re.IGNORECASE)
    
    # Safe prefix patterns (highly selective)
    safe_prefixes = [
        (r'\banti(inflammatory|oxidant|cancer|tumor|viral|bacterial)\b', r'anti-\1'),
        (r'\bnon(invasive|cancerous|specific|functional)\b', r'non-\1'),
        (r'\bmulti(modal|drug|center)\b', r'multi-\1'),
        (r'\bpre(clinical|treatment|operative|natal)\b', r'pre-\1'),
        (r'\bpost(operative|treatment|natal)\b', r'post-\1'),
        (r'\bintra(cellular|venous|operative)\b', r'intra-\1'),
        (r'\binter(cellular|action)\b', r'inter-\1'),
        (r'\bextra(cellular)\b', r'extra-\1'),
        (r'\bco(treatment|expression|infection)\b', r'co-\1'),
    ]
    
    for pattern, replacement in safe_prefixes:
        text = re.sub(pattern, replacement, text, flags=re.IGNORECASE)
    
    return text

def clean_text(text: str) -> str:
    if not text:
        return ""
    
    # 1. Normalize Unicode (Fixes greek letters, math symbols, weird spaces)
    text = unicodedata.normalize('NFKC', text)
    
    # Basic whitespace normalization
    text = re.sub(r'\s+', ' ', text).strip()
    
    # CRITICAL: Remove ALL citations (universal)
    text = remove_pmc_citations_universal(text)
    
    # Normalize quotes
    text = text.replace('\u2018', "'").replace('\u2019', "'")
    text = text.replace('\u201c', '"').replace('\u201d', '"')
    text = text.replace('\u2013', '-').replace('\u2014', '-')
    
    # Fix medical compound terms
    text = fix_compound_medical_terms(text)
    
    # Final cleanup
    text = re.sub(r'\s+', ' ', text).strip()
    text = re.sub(r'\s+([.,;:!?])', r'\1', text)
    text = re.sub(r'([.,;:!?])([A-Za-z])', r'\1 \2', text)
    text = re.sub(r'\s{2,}', ' ', text)
    
    # Remove orphaned punctuation at start
    text = re.sub(r'^[,;:]+\s*', '', text)
    text = re.sub(r'\s+', ' ', text)
    
    # 2. Fix dangling punctuation like " ( , ) " or " [ , ] "
    text = re.sub(r'\(\s*[,;.]?\s*\)', '', text)
    text = re.sub(r'\[\s*[,;.]?\s*\]', '', text)
    
    # 3. Fix cases where numbers were removed but the range dash remains
    text = re.sub(r'\s-\s', '-', text)

    text = re.sub(r'\(\s*Figures?\s*\)', '', text, flags=re.IGNORECASE)

    # 2. Remove In-text Figure References
    # This catches (Fig 1), (Figure 2), (Figs. 1-3)
    text = re.sub(r'\(\s*Figs?\.?\s*\d+[^)]*\)', '', text, flags=re.IGNORECASE)

    # 3. Handle standalone "Figures" on their own line (often from page breaks)
    text = re.sub(r'^\s*Figures?\s*$', '', text, flags=re.MULTILINE | re.IGNORECASE)
    
    # 4. Final cleanup of double spaces
    text = re.sub(r'\s+', ' ', text)
    text = re.sub(r'\[\s*[\d\s,.-]+\s*\]', '', text)

    # 2. UNIVERSAL LOCATION REGEX
    # Matches: (naning,china), (Chicago, IL, USA), (nanjing, china)
    # flags=re.IGNORECASE makes it work for lowercase names
    text = re.sub(r'\(\s*[A-Za-z\s\.]+,\s*[A-Za-z\s\.]+(?:,\s*[A-Za-z\s\.]+)?\s*\)', '', text, flags=re.IGNORECASE)
    location_pattern = r'\(\s*[A-Za-z\s\.]+(?:,\s*[A-Za-z\s\.]+)*\s*\)'
    text = re.sub(location_pattern, '', text, flags=re.IGNORECASE)
    text = re.sub(r'\(\s*[A-Za-z\d\s\.]+(?:,\s*[A-Za-z\d\s\.-]+)*\s*\)', '', text, flags=re.IGNORECASE)
    # 3. SUPPLEMENTAL & TABLE REGEX
    # Catches "supplemental Table S", "Table 1", etc.
    table_figure_pattern = r'\b(supplemental|supplementary|Table|Figure|Fig|S?Table|S?Fig)\.?\s+[A-Z\d]+\b'
    text = re.sub(table_figure_pattern, '', text, flags=re.IGNORECASE)

    # 4. FINAL PUNCTUATION CLEANUP
    # Removes empty brackets left behind: () or []
    text = re.sub(r'\(\s*[,;.\s]*\)', '', text)
    text = re.sub(r'\[\s*[,;.\s]*\]', '', text)
    text = re.sub(r'\(\s*[A-Za-z\d\s\.]+(?:,\s*[A-Za-z\d\s\.-]+)*\s*\)', '', text, flags=re.IGNORECASE)

    # 2. SPACED CITATIONS (Common in Diagnostics)
    # Matches: [ 30 ], [ 12 , 15 ], [1-3]
    text = re.sub(r'\[\s*[\d\s,.-]+\s*\]', '', text)

    # 3. SUPPLEMENTAL & FIGURE/TABLE LABELS
    # Matches: ( Supplementary ), ( Supplementary Methods ), Figure C, Table S1
    text = re.sub(r'\(\s*supplemental|supplementary[^)]*\)', '', text, flags=re.IGNORECASE)
    text = re.sub(r'\b(Table|Figure|Fig|S?Table|S?Fig)\.?\s+[A-Z\d]+\b', '', text, flags=re.IGNORECASE)

    # 4. GARBAGE COLLECTION
    # Clean up empty brackets/parentheses left behind
    text = re.sub(r'\(\s*[,;.\s]*\)', '', text)
    text = re.sub(r'\[\s*[,;.\s]*\]', '', text)
    
    # 5. WHITESPACE NORMALIZATION
    text = re.sub(r'\s+', ' ', text)

    
    # Standardize spaces
    text = re.sub(r'\s+', ' ', text)

    return text.strip()

def remove_supplemental_references(text: str) -> str:
    # 1. Matches: "supplemental Table S1", "supplementary Figure S2", "supplemental Methods"
    # Includes variations in spelling and potential missing numbers
    supp_pattern = r'\b(supplemental|supplementary)\s+(Table|Figure|Methods|Data|File|Note|sTable|sFig)\s+[A-Z]?\d*[a-z]?\b'
    text = re.sub(supp_pattern, '', text, flags=re.IGNORECASE)

    # 2. Matches: "(see supplemental Table S1)" or "(supplementary Data)"
    text = re.sub(r'\(\s*see\s+supplemental[^)]*\)', '', text, flags=re.IGNORECASE)
    
    return text

def final_gold_standard_cleaning(text: str) -> str:
    # 1. Remove Institutional Locations/Addresses in parentheses
    # Catches: (Chicago, IL, USA), (Basel, Switzerland), (Tokyo, Japan)
    # Logic: Look for parentheses containing a city/state and a country name
    location_pattern = r'\(\s*[A-Za-z\s\.]+,\s*[A-Za-z\s\.]+(?:,\s*[A-Za-z\s\.]+)?\s*\)'
    text = re.sub(location_pattern, '', text)

    # 2. Remove Remaining "Table" and "Supplemental" mentions
    # Catches: "supplemental Table S1", "Table S1", "supplementary Methods"
    # Found in your NKTL file: "described in supplemental Table S"
    table_pattern = r'\b(supplemental|supplementary|additional file|see)\s+(Table|Figure|Methods|Data|S?Table|S?Fig)\s+[A-Z\d]*\b'
    text = re.sub(table_pattern, '', text, flags=re.IGNORECASE)

    # 3. Garbage Collection (The most important final step)
    # Remove any empty parentheses () or brackets [] left behind
    text = re.sub(r'\(\s*[,;.\s]*\)', '', text)
    text = re.sub(r'\[\s*[,;.\s]*\]', '', text)
    
    # Standardize spaces
    text = re.sub(r'\s+', ' ', text)
    return text.strip()

def split_long_paragraph(text: str) -> list[str]:
    """Split paragraphs that exceed MAX_LEN"""
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
    """Calculate similarity between two paragraphs"""
    if not p1 or not p2:
        return 0.0
    words1 = set(p1.lower().split())
    words2 = set(p2.lower().split())
    return len(words1 & words2) / len(words1 | words2) if words1 or words2 else 0.0

def is_valid(text: str, heading: str = "") -> bool:
    """Validate if text chunk is suitable for training"""
    if not text or len(text) < MIN_LEN:
        return False
    
    # Reject if too many digits (likely tables)
    digit_ratio = sum(c.isdigit() for c in text) / len(text)
    if digit_ratio > 0.42:
        return False
    
    # Reject if too many special characters
    special_ratio = sum(not c.isalnum() and not c.isspace() for c in text) / len(text)
    if special_ratio > 0.32:
        return False
    
    # Reject all-caps or all-whitespace
    if text.isupper() or text.isspace():
        return False
    
    # Lower threshold for priority sections
    heading_lower = heading.lower()
    if any(kw in heading_lower for kw in PRIORITY_SECTIONS):
        return len(text) > MIN_LEN - 50
    
    text_lower = text.lower()
    table_keywords = ["confidence interval", "hazard ratio", "p value", "p-value", "vs.", "n =", "study design"]
    keyword_count = sum(1 for kw in table_keywords if kw in text_lower)
    
    # If a paragraph has 3+ of these keywords, it's likely a table row, not a sentence
    if keyword_count >= 3:
        return False

    return True

import unicodedata

def get_fingerprint(text: str) -> str:
    """Creates a unique ID for a paragraph by removing all non-alphanumeric characters."""
    return re.sub(r'\W+', '', text.lower())

def process_txt_file(in_path: Path, out_path: Path):
    """Process a single text file with caption filtering and fast fingerprinting"""
    try:
        # Errors='ignore' is good for batch processing 100k files to prevent crashes
        with open(in_path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
    except Exception as e:
        print(f"❌ Error reading {in_path}: {e}")
        return
    
    cleaned_blocks = []
    buffer = ""
    current_heading = ""
    
    # NEW: Using a set for O(1) lookup speed - essential for 100k articles
    seen_paras = set()  
    seen_headings = set()
    
    # Skip document tag if present
    start_idx = 1 if lines and lines[0].strip().startswith("<DOCUMENT") else 0
    
    for line in lines[start_idx:]:
        # 1. Unicode Normalization (Fixes Greek letters, math symbols, and ligatures)
        line = unicodedata.normalize('NFKC', line)
        stripped = line.strip()
        
        # Skip empty lines and separator lines
        if not stripped or re.match(r'^[=\-─]{3,}$', stripped):
            continue

        # ====================================================================
        # STEP 3: CAPTION FILTERING
        # ====================================================================
        # Skip lines that are clearly Figure or Table captions
        if re.match(r'^(?:Figure|Fig\.?|Table|Scheme|Box|Plate|Video|S\d+)\s*\d+[:.]', stripped, re.IGNORECASE):
            continue
            
        # Detect headings
        is_heading = (
            re.match(r'^\d+(\.\d+)*\.?\s+[A-Z]', line) or 
            line.isupper() or 
            re.match(r'^={1,}.+={1,}$', stripped)
        )
        
        if is_heading:
            # Process accumulated buffer before starting new section
            if buffer:
                cleaned = clean_text(buffer)
                for chunk in split_long_paragraph(cleaned):
                    if is_valid(chunk, current_heading):
                        # ====================================================
                        # STEP 4: FINGERPRINT DEDUPLICATION
                        # ====================================================
                        chunk_fp = get_fingerprint(chunk)
                        if chunk_fp not in seen_paras:
                            seen_paras.add(chunk_fp)
                            cleaned_blocks.append(chunk)
                buffer = ""
            
            # Process heading
            current_heading = normalize_heading(stripped).lower()
            if current_heading:
                if SKIP_ABSTRACT and "abstract" in current_heading:
                    buffer = ""
                    continue
                
                if current_heading not in seen_headings:
                    seen_headings.add(current_heading)
                    cleaned_blocks.append(f"\n{current_heading.title()}\n")
        
        else:
            # Accumulate text (skip if currently in an abstract section)
            if not (SKIP_ABSTRACT and "abstract" in current_heading):
                buffer += " " + stripped
    
    # Process final buffer at end of file
    if buffer and not (SKIP_ABSTRACT and "abstract" in current_heading):
        cleaned = clean_text(buffer)
        for chunk in split_long_paragraph(cleaned):
            if is_valid(chunk, current_heading):
                chunk_fp = get_fingerprint(chunk)
                if chunk_fp not in seen_paras:
                    seen_paras.add(chunk_fp)
                    cleaned_blocks.append(chunk)
    
    # Write output
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(out_path, "w", encoding="utf-8") as f:
            for block in cleaned_blocks:
                f.write(block.strip() + "\n\n")
    except Exception as e:
        print(f"❌ Error writing {out_path}: {e}")

def process_wrapper(txt_path):
    """
    This must be at the top level of the script to be 'picklable' 
    on Windows for multiprocessing.
    """
    try:
        # Note: You can define ROOT_INPUT/OUTPUT as global variables 
        # or use them directly if they are defined at the top of your script.
        relative = txt_path.relative_to(ROOT_INPUT)
        output_txt = ROOT_OUTPUT / relative
        process_txt_file(txt_path, output_txt)
        return True
    except Exception as e:
        print(f"❌ Failed: {txt_path.name} - {e}")
        return False
        

def main():
    print("\n" + "="*80)
    print("🧹 UNIVERSAL MEDICAL TEXT CLEANER FOR LLM TRAINING")
    print("="*80)
    print("\n✨ Optimized for 100,000+ PMC documents across ALL diseases")
    print("\n📋 Features:")
    print("   • Removes 30+ types of citations and references")
    print("   • Handles ALL PMC-specific patterns universally")
    print("   • Preserves medical terminology and meaning")
    print("   • Fixed compound medical term hyphenation")
    print("   • Deduplication and quality filtering")
    print("   • Chemical formula preservation")
    print()
    
    # Find all text files
    txt_files = list(ROOT_INPUT.rglob("*.txt"))
    total_files = len(txt_files)
    
    print(f"📁 Found {total_files} files to process...")
    print(f"🚀 Starting parallel processing using all CPU cores...")

    # 2. CALL THE TOP-LEVEL FUNCTION
    with ProcessPoolExecutor() as executor:
        # map() now points to the top-level process_wrapper
        results = list(executor.map(process_wrapper, txt_files))

    processed = sum(1 for r in results if r)
    errors = total_files - processed

    print("\n" + "="*80)
    print(f"✅ COMPLETE — Processed {processed} files")
    if errors > 0:
        print(f"⚠️  Errors: {errors} files failed")
    print("="*80 + "\n")

if __name__ == "__main__":
    # 3. THIS BLOCK IS MANDATORY ON WINDOWS
    main()