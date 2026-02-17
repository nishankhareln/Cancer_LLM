"""
Enhanced PubMed Central Prevention Extractor - PROSTATE CANCER
Extracts prevention strategies, risk reduction, screening, lifestyle, and preventive measures
Generates TEXT, JSON, and DETAILED formats for LLM training
Covers: Prostate Adenocarcinoma, Advanced Prostate Cancer, Metastatic Prostate Cancer
"""

from Bio import Entrez
import xml.etree.ElementTree as ET
import os
import re
import time
import json

# ================= CONFIG =================
Entrez.email = "nkharel57@gmail.com"  # REQUIRED - Your email

# Prostate cancer specific search queries
SEARCH_QUERIES = [
    "prostate cancer prevention",
    "prostate cancer risk factors",
    "prostate cancer screening",
    "PSA screening prostate cancer",
    "digital rectal examination prostate",
    "prostate cancer chemoprevention",
    "finasteride prostate cancer prevention",
    "dutasteride prostate cancer prevention",
    "5-alpha reductase inhibitors prostate",
    "diet prostate cancer prevention",
    "lycopene prostate cancer",
    "vitamin E selenium prostate cancer",
    "physical activity prostate cancer",
    "obesity prostate cancer prevention",
    "familial prostate cancer",
    "BRCA prostate cancer risk",
    "racial disparities prostate cancer"
]

MAX_ARTICLES = 100  # Per query
MAX_RETRIES = 3
RETRY_DELAY = 2  # seconds
OUTPUT_TEXT = "prostate_cancer_prevention_detailed.txt"
OUTPUT_JSON = "prostate_cancer_prevention_detailed.json"
OUTPUT_DETAILED = "prostate_cancer_prevention_DETAILED.json"

# ================= COMPREHENSIVE PROSTATE CANCER PREVENTION KNOWLEDGE BASE =================

PREVENTION_DATABASE = {
    
    # ===== SCREENING & EARLY DETECTION =====
    "psa_screening": {
        "category": "Secondary Prevention",
        "type": "Early Detection",
        "name": "PSA Screening for Prostate Cancer",
        "description": "Prostate-specific antigen (PSA) screening can detect early prostate cancer but involves complex tradeoffs",
        "current_guidelines": {
            "USPSTF_2018": {
                "ages_55_69": "Shared decision-making; individualized approach based on values and preferences",
                "age_70_plus": "Do not screen (grade D recommendation)",
                "rationale": "Small potential benefit; significant harms from overdiagnosis and overtreatment"
            },
            "ACS_2023": {
                "age_50_average_risk": "Discuss screening; make informed decision",
                "age_45_higher_risk": "African American men, family history",
                "age_40_highest_risk": "Multiple family members with early prostate cancer, BRCA mutations",
                "life_expectancy": "At least 10 years to benefit from screening"
            },
            "AUA_2023": {
                "ages_55_69": "Shared decision-making recommended",
                "age_40_45": "Baseline PSA for high-risk men",
                "screening_interval": "2+ years for most men if continuing screening"
            }
        },
        "screening_protocol": {
            "PSA_test": {
                "normal_range": "<4.0 ng/mL traditionally; newer approaches more nuanced",
                "PSA_velocity": "Rate of change over time",
                "PSA_density": "PSA level relative to prostate size",
                "free_PSA": "Ratio of free to total PSA; lower ratio suggests cancer"
            },
            "digital_rectal_exam": {
                "role": "Adjunct to PSA; can detect PSA-negative cancers",
                "limitations": "Only detects posterior/lateral prostate; subjective"
            },
            "risk_calculators": {
                "PCPT_calculator": "Prostate Cancer Prevention Trial risk calculator",
                "ERSPC_calculator": "European risk calculator",
                "use": "Incorporate age, race, family history, PSA, DRE"
            }
        },
        "benefits_of_screening": {
            "mortality_reduction": "20-30% relative reduction in prostate cancer mortality (ERSPC trial)",
            "early_detection": "Detects cancers at earlier, more treatable stages",
            "interval": "Number needed to screen: ~1000 men to prevent 1 death over 10-13 years"
        },
        "harms_of_screening": {
            "overdiagnosis": "23-50% of screen-detected cancers would never cause symptoms",
            "overtreatment": "Many indolent cancers treated unnecessarily",
            "biopsy_complications": "Infection (1-3%), bleeding, pain, anxiety",
            "treatment_complications": "Erectile dysfunction (20-70%), urinary incontinence (5-20%)",
            "false_positives": "Elevated PSA often not cancer; causes anxiety, unnecessary biopsies",
            "psychological_distress": "Anxiety from diagnosis and treatment decisions"
        },
        "shared_decision_making": {
            "discuss": [
                "Benefits: mortality reduction, early detection",
                "Harms: overdiagnosis, overtreatment, biopsy/treatment complications",
                "Individual risk factors (age, race, family history)",
                "Personal values and preferences",
                "Active surveillance as management option for low-risk cancer"
            ],
            "decision_aids": "Use validated tools to facilitate discussion",
            "documentation": "Document discussion and decision in medical record"
        },
        "active_surveillance": {
            "definition": "Close monitoring of low-risk prostate cancer without immediate treatment",
            "eligibility": "Low-grade (Gleason ≤6), low PSA, limited tumor volume",
            "protocol": "Regular PSA, DRE, repeat biopsies, MRI",
            "rationale": "Avoid/delay treatment complications while monitoring for progression",
            "outcomes": "Excellent cancer-specific survival; 50% avoid treatment long-term"
        },
        "evidence": "Large RCTs (ERSPC, PLCO) show mortality benefit but also substantial harms",
        "recommendations": [
            "Engage in shared decision-making for men 55-69 years",
            "Do not screen men <55 years (average risk) or >70 years",
            "Earlier discussion (45-50) for higher-risk men",
            "Use risk calculators to personalize approach",
            "Offer active surveillance for low-risk cancers",
            "Screening interval 2+ years if PSA low and stable"
        ]
    },

    "multiparametric_mri": {
        "category": "Secondary Prevention",
        "type": "Advanced Detection",
        "name": "Multiparametric MRI in Prostate Cancer Detection",
        "description": "mpMRI can improve detection of clinically significant prostate cancer and reduce overdiagnosis",
        "role_in_detection": {
            "pre_biopsy_mri": "Identify suspicious lesions before biopsy; target biopsy to lesions",
            "MRI_targeted_biopsy": "Fusion or cognitive targeting of suspicious areas",
            "benefits": "Detects more high-grade cancers, fewer low-grade cancers vs systematic biopsy"
        },
        "PI_RADS_scoring": {
            "scale": "1-5 scoring system for likelihood of clinically significant cancer",
            "PI_RADS_1_2": "Low risk; may avoid biopsy",
            "PI_RADS_3": "Equivocal; consider biopsy based on other factors",
            "PI_RADS_4_5": "High risk; biopsy recommended"
        },
        "advantages": {
            "detect_significant_cancer": "Better detection of clinically significant (Gleason ≥7) cancers",
            "reduce_overdetection": "May avoid detection of insignificant cancers",
            "guide_biopsy": "Target suspicious areas; improve diagnostic yield",
            "monitoring": "Useful for active surveillance"
        },
        "limitations": {
            "cost": "Expensive compared to PSA/DRE",
            "availability": "Requires expertise; not universally available",
            "false_positives": "Can miss some cancers; false positives occur",
            "standardization": "Quality varies by center and radiologist"
        },
        "evidence": "Multiple trials show improved detection of significant cancers; reducing unnecessary biopsies",
        "recommendations": [
            "Consider pre-biopsy MRI in biopsy-naive men with elevated PSA",
            "Use MRI-targeted biopsy when lesions identified",
            "Helpful for active surveillance monitoring",
            "Ensure high-quality MRI and experienced interpretation",
            "Not yet standard of care everywhere; evolving role"
        ]
    },

    # ===== CHEMOPREVENTION =====
    "five_alpha_reductase_inhibitors": {
        "category": "Medical Prevention",
        "type": "Chemoprevention",
        "name": "5-Alpha Reductase Inhibitors (Finasteride, Dutasteride)",
        "description": "5-ARIs reduce prostate cancer risk but increase high-grade cancer detection; not recommended for prevention",
        "mechanism": {
            "action": "Block conversion of testosterone to dihydrotestosterone (DHT)",
            "prostate_effect": "Reduce prostate size, lower PSA by ~50%",
            "drugs": "Finasteride (Proscar) 5mg daily, Dutasteride (Avodart) 0.5mg daily"
        },
        "clinical_trials": {
            "PCPT_finasteride": {
                "design": "18,882 men randomized to finasteride vs placebo for 7 years",
                "results": "25% reduction in overall prostate cancer risk",
                "concern": "Increased high-grade cancers (Gleason 7-10) in finasteride group",
                "debate": "Unclear if true increase or detection bias"
            },
            "REDUCE_dutasteride": {
                "design": "8,231 men randomized to dutasteride vs placebo for 4 years",
                "results": "23% reduction in overall prostate cancer risk",
                "concern": "Small increase in high-grade cancers",
                "consistency": "Similar findings to PCPT"
            },
            "long_term_followup": {
                "PCPT_18_years": "No difference in survival; possible slight reduction in high-grade cancers",
                "interpretation": "Initial concerns about high-grade cancers likely detection artifact"
            }
        },
        "FDA_stance": {
            "status": "Not approved for prostate cancer prevention",
            "approved_use": "Benign prostatic hyperplasia (BPH), male pattern baldness",
            "reason": "Unclear benefit-risk profile for prevention"
        },
        "who_might_benefit": {
            "high_risk_men": "Strong family history, BRCA mutations, African American",
            "BPH_patients": "Already taking for urinary symptoms",
            "informed_decision": "After discussion of risks/benefits"
        },
        "side_effects": {
            "sexual": "Decreased libido (5-10%), erectile dysfunction (5-10%), ejaculation disorders",
            "other": "Gynecomastia (breast enlargement), mood changes",
            "PSA_effect": "Lowers PSA ~50%; adjust interpretation"
        },
        "current_recommendations": {
            "USPSTF": "Recommends against routine use for prevention (Grade D)",
            "AUA": "Not recommended for prevention in general population",
            "individual_cases": "May discuss with high-risk men on case-by-case basis"
        },
        "evidence": "Strong evidence of reduced cancer risk; concerns about high-grade cancers; not recommended for prevention",
        "recommendations": [
            "Do not use 5-ARIs solely for prostate cancer prevention",
            "May be reasonable for high-risk men already on medication for BPH",
            "Discuss risks/benefits if considering for high-risk individuals",
            "Adjust PSA interpretation (multiply by 2) if on 5-ARI",
            "Monitor for sexual side effects"
        ]
    },

    "vitamin_e_selenium": {
        "category": "Medical Prevention",
        "type": "Chemoprevention",
        "name": "Vitamin E and Selenium Supplementation",
        "description": "SELECT trial showed no benefit and possible harm from vitamin E and selenium for prostate cancer prevention",
        "SELECT_trial": {
            "design": "35,533 men randomized to selenium, vitamin E, both, or placebo",
            "duration": "Planned 7-12 years; stopped early due to lack of benefit",
            "groups": [
                "Selenium 200 mcg/day from L-selenomethionine",
                "Vitamin E 400 IU/day (all-rac-alpha-tocopheryl acetate)",
                "Both selenium + vitamin E",
                "Placebo"
            ]
        },
        "results": {
            "selenium": "No reduction in prostate cancer; slight increase in diabetes risk",
            "vitamin_E": "17% increased risk of prostate cancer (statistically significant)",
            "combination": "No benefit; possible increased cancer risk",
            "mortality": "No impact on overall mortality"
        },
        "previous_optimism": {
            "observational_studies": "Suggested benefit from selenium, vitamin E",
            "NPC_trial": "Selenium reduced prostate cancer in secondary analysis",
            "ATBC_trial": "Vitamin E reduced prostate cancer in smokers (secondary analysis)",
            "lesson": "Observational studies and secondary analyses can mislead"
        },
        "mechanisms_explored": {
            "selenium": "Antioxidant; may affect androgen metabolism",
            "vitamin_E": "Antioxidant; anti-inflammatory properties",
            "paradox": "Despite plausible mechanisms, no benefit in RCT"
        },
        "current_status": {
            "recommendation": "Do NOT take vitamin E or selenium for prostate cancer prevention",
            "harm": "Vitamin E may increase prostate cancer risk",
            "dietary_sources": "Adequate nutrition from food; supplements not beneficial"
        },
        "evidence": "Large definitive RCT (SELECT) showed no benefit and harm from vitamin E",
        "recommendations": [
            "Do NOT recommend vitamin E supplements for prostate cancer prevention",
            "Do NOT recommend selenium supplements for prevention",
            "May increase prostate cancer risk (vitamin E)",
            "Adequate nutrition from balanced diet sufficient",
            "Counsel patients against these supplements for prevention"
        ]
    },

    # ===== DIET & NUTRITION =====
    "diet_nutrition_prostate": {
        "category": "Lifestyle Prevention",
        "type": "Primary Prevention",
        "name": "Diet and Nutrition for Prostate Cancer Prevention",
        "description": "Dietary patterns and specific nutrients may influence prostate cancer risk",
        "dietary_patterns": {
            "plant_based_diet": {
                "evidence": "Associated with lower prostate cancer risk and progression",
                "components": "High fruits, vegetables, whole grains, legumes",
                "mechanisms": "Antioxidants, phytochemicals, fiber, lower IGF-1"
            },
            "mediterranean_diet": {
                "evidence": "Associated with lower risk of advanced/fatal prostate cancer",
                "components": "Olive oil, fish, vegetables, fruits, nuts, whole grains",
                "mechanisms": "Anti-inflammatory, antioxidant effects"
            },
            "western_diet": {
                "evidence": "Associated with higher risk, especially advanced disease",
                "components": "High red/processed meat, refined grains, sugar, saturated fat",
                "mechanisms": "Pro-inflammatory, obesity, IGF-1 elevation"
            }
        },
        "specific_nutrients": {
            "lycopene": {
                "sources": "Tomatoes (especially cooked), watermelon, pink grapefruit",
                "evidence": "Modest protective effect in some studies; inconsistent",
                "bioavailability": "Better absorbed from cooked tomatoes with fat",
                "supplements": "Not proven effective; food sources preferred"
            },
            "cruciferous_vegetables": {
                "sources": "Broccoli, cauliflower, Brussels sprouts, kale",
                "evidence": "Possible protective effect",
                "mechanisms": "Sulforaphane, indole-3-carbinol; anti-cancer properties",
                "recommendation": "Include regularly in diet"
            },
            "soy_foods": {
                "sources": "Tofu, tempeh, edamame, soy milk",
                "evidence": "Possible protective effect, especially in Asian populations",
                "mechanisms": "Phytoestrogens (isoflavones); anti-androgenic",
                "recommendation": "Whole soy foods may be beneficial"
            },
            "omega_3_fatty_acids": {
                "sources": "Fatty fish (salmon, sardines), flaxseed, walnuts",
                "evidence": "Mixed results; possible benefit for advanced disease prevention",
                "recommendation": "2+ servings fish per week"
            },
            "vitamin_D": {
                "sources": "Sunlight, fatty fish, fortified foods",
                "evidence": "Low levels associated with aggressive disease; supplementation unproven",
                "recommendation": "Adequate levels important; avoid deficiency"
            },
            "calcium": {
                "concern": "Very high intake (>2000mg/day) may increase risk",
                "recommendation": "Moderate intake; avoid excessive supplementation",
                "balance": "Adequate but not excessive"
            }
        },
        "foods_to_limit": {
            "red_meat": {
                "evidence": "High intake associated with advanced prostate cancer",
                "mechanisms": "Heme iron, heterocyclic amines from cooking, saturated fat",
                "recommendation": "Limit to <2 servings/week; choose lean cuts"
            },
            "processed_meat": {
                "evidence": "Associated with increased risk",
                "examples": "Bacon, sausage, hot dogs, deli meats",
                "recommendation": "Minimize consumption"
            },
            "high_fat_dairy": {
                "evidence": "Possible increased risk with very high intake",
                "recommendation": "Moderate intake; choose low-fat options",
                "nuance": "Not as clear as red meat association"
            },
            "charred_grilled_meat": {
                "concern": "Heterocyclic amines (HCAs), polycyclic aromatic hydrocarbons (PAHs)",
                "recommendation": "Avoid heavily charred meat; marinate before grilling"
            }
        },
        "supplements_not_recommended": {
            "vitamin_E": "May increase risk (SELECT trial)",
            "selenium": "No benefit; possible harm (SELECT trial)",
            "beta_carotene": "No benefit for prostate cancer",
            "multivitamins": "No proven benefit for prevention",
            "high_dose_calcium": "May increase risk at very high doses"
        },
        "evidence": "Consistent observational evidence for dietary patterns; mixed for specific nutrients; supplements not beneficial",
        "recommendations": [
            "Plant-based or Mediterranean dietary pattern",
            "5+ servings fruits and vegetables daily",
            "Include tomatoes, cruciferous vegetables",
            "Limit red and processed meat",
            "Choose whole grains over refined",
            "Fatty fish 2+ times per week",
            "Avoid excessive calcium supplementation",
            "Do NOT use vitamin E or selenium supplements"
        ]
    },

    # ===== PHYSICAL ACTIVITY & OBESITY =====
    "physical_activity_obesity": {
        "category": "Lifestyle Prevention",
        "type": "Primary Prevention",
        "name": "Physical Activity and Obesity Prevention",
        "description": "Physical activity and healthy weight may reduce risk of advanced prostate cancer and improve outcomes",
        "physical_activity": {
            "evidence": {
                "overall_risk": "Modest reduction in overall prostate cancer risk (10-30%)",
                "advanced_disease": "Stronger protection against advanced/fatal prostate cancer (30-50% reduction)",
                "post_diagnosis": "Improves survival after diagnosis; reduces recurrence",
                "dose_response": "More activity associated with greater benefit"
            },
            "mechanisms": {
                "hormonal": "Reduces testosterone, IGF-1; improves insulin sensitivity",
                "obesity_prevention": "Helps maintain healthy weight",
                "inflammation": "Anti-inflammatory effects",
                "immune_function": "Enhanced immune surveillance"
            },
            "recommendations": {
                "moderate_intensity": "150+ minutes per week (brisk walking, cycling, swimming)",
                "vigorous_intensity": "75+ minutes per week (running, sports)",
                "strength_training": "2+ days per week",
                "reduce_sedentary": "Limit sitting time; break up prolonged sitting"
            }
        },
        "obesity": {
            "risk_relationship": {
                "overall_risk": "Inverse association with overall prostate cancer (obesity protective)",
                "advanced_disease": "Positive association with advanced/aggressive disease",
                "mortality": "Obesity increases prostate cancer mortality",
                "paradox": "Lower PSA in obese men may delay detection"
            },
            "mechanisms": {
                "hormonal_changes": "Lower testosterone, altered estrogen, higher insulin/IGF-1",
                "inflammation": "Chronic low-grade inflammation",
                "detection_bias": "Larger prostate, lower PSA may delay diagnosis",
                "treatment_effects": "Worse outcomes after treatment in obese men"
            },
            "body_composition": {
                "visceral_fat": "Abdominal obesity particularly harmful",
                "muscle_mass": "Sarcopenia associated with worse outcomes",
                "waist_circumference": ">40 inches (men) associated with higher risk"
            },
            "weight_management": {
                "healthy_BMI": "18.5-24.9 kg/m²",
                "waist_circumference": "<40 inches for men",
                "weight_loss": "5-10% loss improves metabolic health",
                "maintain": "Prevent weight gain with age"
            }
        },
        "evidence": "Strong evidence for physical activity benefits; complex obesity relationship",
        "recommendations": [
            "Regular physical activity: 150+ min/week moderate or 75+ min/week vigorous",
            "Strength training 2+ days per week",
            "Maintain healthy weight (BMI 18.5-24.9)",
            "Avoid weight gain with age",
            "Reduce sedentary time",
            "Even modest activity beneficial; some activity better than none",
            "Post-diagnosis: Continue/start exercise to improve outcomes"
        ]
    },

    # ===== GENETIC & FAMILIAL RISK =====
    "genetic_familial_prostate": {
        "category": "Genetic Risk Assessment",
        "type": "Risk Stratification",
        "name": "Genetic and Familial Risk of Prostate Cancer",
        "description": "Family history and genetic mutations significantly increase prostate cancer risk",
        "family_history": {
            "one_first_degree": "2-3 times increased risk with one affected first-degree relative",
            "two_or_more": "5-11 times increased risk with 2+ affected relatives",
            "early_onset": "Higher risk if relative diagnosed <60 years",
            "hereditary_prostate": "5-10% of prostate cancers are hereditary",
            "implications": "Earlier screening, genetic testing consideration"
        },
        "high_risk_genes": {
            "BRCA2": {
                "risk": "5-8 times increased risk; younger age at diagnosis",
                "aggressiveness": "More aggressive disease; worse prognosis",
                "prevalence": "~2% of prostate cancers",
                "treatment": "May respond to PARP inhibitors",
                "screening": "Earlier screening (age 40-45); more intensive"
            },
            "BRCA1": {
                "risk": "2-3 times increased risk",
                "aggressiveness": "Possible more aggressive disease",
                "implications": "Earlier screening recommended"
            },
            "HOXB13_G84E": {
                "risk": "3-5 times increased risk",
                "prevalence": "~1-2% in men of European ancestry",
                "early_onset": "Associated with early-onset disease",
                "screening": "Earlier screening if mutation present"
            },
            "Lynch_syndrome": {
                "genes": "MLH1, MSH2, MSH6, PMS2",
                "risk": "2-5 times increased risk",
                "other_cancers": "Colorectal, endometrial, other cancers",
                "immunotherapy": "May respond to checkpoint inhibitors"
            },
            "ATM": {
                "risk": "Moderate increased risk (2-3x)",
                "DNA_repair": "DNA damage response gene"
            }
        },
        "genetic_testing": {
            "indications": [
                "Metastatic prostate cancer (all patients)",
                "High-risk localized cancer (consider testing)",
                "Strong family history (3+ relatives or early onset)",
                "Ashkenazi Jewish ancestry",
                "Personal/family history of BRCA-associated cancers"
            ],
            "panel_testing": "Multi-gene panels test multiple genes simultaneously",
            "germline_vs_somatic": "Germline testing (inherited); somatic testing (tumor-only)",
            "cascade_testing": "Test family members if pathogenic mutation found"
        },
        "clinical_implications": {
            "screening": [
                "Earlier screening (age 40-45 for BRCA2)",
                "More frequent screening",
                "Consider MRI in addition to PSA",
                "Lower threshold for biopsy"
            ],
            "treatment": [
                "PARP inhibitors for BRCA1/2 metastatic cancer",
                "Immunotherapy for Lynch syndrome",
                "May influence treatment decisions",
                "More aggressive treatment consideration"
            ],
            "family_counseling": [
                "Genetic counseling for patient and family",
                "Inform relatives about inherited risk",
                "Cascade testing for family members",
                "Female relatives: breast/ovarian cancer risk"
            ]
        },
        "polygenic_risk_scores": {
            "concept": "Combined effect of 100+ common genetic variants",
            "potential": "Identify men at higher/lower risk",
            "status": "Research stage; not yet clinical standard",
            "future": "May personalize screening approach"
        },
        "evidence": "Strong evidence for BRCA2, BRCA1, HOXB13; moderate for others; genetic testing increasingly important",
        "recommendations": [
            "Assess family history at age 40",
            "Genetic testing for metastatic prostate cancer",
            "Consider testing for high-risk localized disease or strong family history",
            "Earlier screening (40-45) for BRCA2 carriers",
            "Genetic counseling when testing performed",
            "Inform family members if mutation found",
            "Use genetic information to guide treatment (PARP inhibitors, immunotherapy)"
        ]
    },

    # ===== RACIAL DISPARITIES =====
    "racial_ethnic_disparities": {
        "category": "Health Equity",
        "type": "Risk Stratification",
        "name": "Racial and Ethnic Disparities in Prostate Cancer",
        "description": "African American and Caribbean men have highest prostate cancer incidence and mortality worldwide",
        "african_american_men": {
            "incidence": "70% higher incidence than white men",
            "mortality": "2-2.5 times higher mortality than white men",
            "age_at_diagnosis": "Younger age at diagnosis (median 2-3 years earlier)",
            "aggressiveness": "Higher rates of advanced disease at diagnosis",
            "biology": "May have more aggressive tumor biology",
            "outcomes": "Worse outcomes even after controlling for access and treatment"
        },
        "contributing_factors": {
            "genetic": {
                "ancestry": "African ancestry associated with higher risk",
                "genetic_variants": "Different distribution of risk alleles",
                "biology": "Possible biological differences in tumors"
            },
            "social_determinants": {
                "access_to_care": "Barriers to screening, diagnosis, treatment",
                "socioeconomic": "Poverty, lack of insurance, healthcare deserts",
                "structural_racism": "Systemic inequities in healthcare",
                "trust": "Historical mistrust of medical system"
            },
            "lifestyle_factors": {
                "diet": "Higher red meat, lower vegetable intake in some populations",
                "obesity": "Higher obesity rates",
                "vitamin_D": "Lower vitamin D levels (darker skin)",
                "comorbidities": "Higher rates of diabetes, hypertension"
            }
        },
        "other_populations": {
            "hispanic_latino": "Lower incidence than white men; improving but still disparities in access",
            "asian_american": "Lower incidence overall; varies by country of origin",
            "native_american": "Limited data; possible underdiagnosis",
            "caribbean": "Among highest incidence rates globally"
        },
        "screening_recommendations": {
            "earlier_discussion": "Begin screening discussion at age 40-45 for African American men",
            "family_history": "Age 40 if family history in addition to African American race",
            "informed_decision": "Discuss higher risk and potential benefits/harms",
            "access": "Ensure access to screening and follow-up"
        },
        "interventions_to_reduce_disparities": {
            "policy_level": [
                "Expand healthcare coverage",
                "Address social determinants of health",
                "Increase diversity in clinical trials",
                "Community-based interventions",
                "Culturally tailored education"
            ],
            "healthcare_system": [
                "Patient navigation programs",
                "Reduce structural barriers to care",
                "Implicit bias training for providers",
                "Culturally competent care",
                "Community health workers"
            ],
            "research": [
                "Include diverse populations in research",
                "Study biological differences",
                "Understand social determinants",
                "Develop tailored interventions"
            ],
            "community_engagement": [
                "Partner with community organizations",
                "Barbershop health programs",
                "Faith-based interventions",
                "Trusted messengers",
                "Address mistrust through transparency"
            ]
        },
        "evidence": "Overwhelming evidence of disparities; multiple contributing factors; interventions can reduce gaps",
        "recommendations": [
            "Earlier screening discussion (age 40-45) for African American men",
            "Aggressive outreach and navigation for underserved populations",
            "Address social determinants and structural barriers",
            "Culturally tailored education and interventions",
            "Increase diversity in clinical trials and research",
            "Community-based participatory approaches",
            "Address implicit bias in healthcare",
            "Ensure equal access to high-quality treatment"
        ]
    },

    # ===== SEXUAL ACTIVITY & EJACULATION =====
    "ejaculation_frequency": {
        "category": "Lifestyle Prevention",
        "type": "Primary Prevention",
        "name": "Ejaculation Frequency and Prostate Cancer Risk",
        "description": "Higher ejaculation frequency may be associated with lower prostate cancer risk",
        "epidemiological_evidence": {
            "health_professionals_study": "21+ ejaculations/month associated with 20-30% lower risk vs ≤7/month",
            "age_relationship": "Association strongest for ejaculations in 20s and 40s-50s",
            "dose_response": "Higher frequency associated with lower risk",
            "consistency": "Multiple studies show similar associations"
        },
        "proposed_mechanisms": {
            "clearance_hypothesis": "Regular ejaculation clears potentially harmful substances from prostate",
            "cell_turnover": "Increases cell turnover, reduces time for mutations",
            "hormonal": "May affect hormonal environment",
            "psychological": "Stress reduction, immune effects",
            "speculative": "Mechanisms remain theoretical"
        },
        "limitations_of_evidence": {
            "observational": "All evidence observational; can't prove causation",
            "recall_bias": "Self-reported ejaculation frequency",
            "confounding": "Healthier men may have higher ejaculation frequency",
            "reverse_causation": "Possible that early disease reduces ejaculation",
            "no_rct": "Randomized trial not feasible"
        },
        "clinical_implications": {
            "not_recommendation": "Insufficient evidence to make clinical recommendation",
            "patient_inquiry": "Can discuss if patient asks; explain observational nature",
            "no_harm": "Regular sexual activity not harmful; may have other health benefits",
            "context": "Part of overall healthy lifestyle"
        },
        "evidence": "Consistent observational associations; mechanisms plausible but unproven; cannot prove causation",
        "recommendations": [
            "Not sufficient evidence for clinical recommendation",
            "Can discuss observational findings if patient inquires",
            "Emphasize healthy lifestyle overall",
            "Regular sexual activity part of overall health and well-being"
        ]
    }
}

# Flatten keywords for searching
ALL_PREVENTION_KEYWORDS = [
    "prevention", "risk reduction", "screening", "surveillance",
    "PSA", "prostate-specific antigen", "digital rectal exam", "DRE",
    "biopsy", "active surveillance", "watchful waiting",
    "finasteride", "dutasteride", "5-alpha reductase",
    "chemoprevention", "vitamin E", "selenium", "SELECT",
    "diet", "nutrition", "lycopene", "tomato",
    "physical activity", "exercise", "obesity", "weight",
    "genetic", "familial", "BRCA", "BRCA1", "BRCA2", "family history",
    "African American", "racial disparities", "health equity",
    "MRI", "multiparametric", "PI-RADS",
    "shared decision making", "overdiagnosis", "overtreatment"
]

SECTION_ALLOW = [
    "prevention", "screening", "surveillance", "risk", "risk factors",
    "PSA", "prostate-specific antigen", "digital rectal",
    "active surveillance", "watchful waiting",
    "chemoprevention", "finasteride", "dutasteride",
    "vitamin", "selenium", "supplement",
    "diet", "nutrition", "lifestyle",
    "physical activity", "exercise", "obesity",
    "genetic", "familial", "hereditary", "family history",
    "BRCA", "mutation",
    "disparity", "African American", "race", "ethnicity",
    "MRI", "imaging", "detection",
    "early detection", "shared decision"
]

# ================= FUNCTIONS =================

def search_pmc():
    """Search PubMed Central with multiple queries"""
    all_ids = set()  # Use set to avoid duplicates
    
    for query_num, query in enumerate(SEARCH_QUERIES, 1):
        print(f"\n🔍 Query {query_num}/{len(SEARCH_QUERIES)}: '{query}'")
        
        for attempt in range(MAX_RETRIES):
            try:
                handle = Entrez.esearch(
                    db="pmc", 
                    term=query, 
                    retmax=MAX_ARTICLES, 
                    sort="relevance",
                    usehistory="y"
                )
                record = Entrez.read(handle)
                handle.close()
                
                new_ids = record.get("IdList", [])
                all_ids.update(new_ids)
                print(f"  ✅ Found {len(new_ids)} articles ({len(all_ids)} total unique)")
                time.sleep(0.5)  # Be nice to NCBI servers
                break
                
            except Exception as e:
                print(f"  ⚠️ Attempt {attempt + 1}/{MAX_RETRIES} failed: {str(e)[:60]}")
                if attempt < MAX_RETRIES - 1:
                    time.sleep(RETRY_DELAY * (attempt + 1))
                else:
                    print(f"  ❌ Skipping query after {MAX_RETRIES} attempts")
    
    return list(all_ids)

def fetch_xml(pmc_id):
    """Fetch article XML with retry logic"""
    for attempt in range(MAX_RETRIES):
        try:
            handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="full", retmode="xml")
            xml_data = handle.read()
            handle.close()
            return xml_data
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
            else:
                raise e

def is_prevention_section(title):
    """Check if section is prevention-related"""
    if not title:
        return False
    title = title.lower()
    return any(key in title for key in SECTION_ALLOW)

def get_prevention_details(topic_name):
    """Retrieve prevention details from database"""
    topic_name_lower = topic_name.lower()
    for topic, details in PREVENTION_DATABASE.items():
        if topic_name_lower == topic or topic_name_lower in topic.lower():
            return details
    return None

def extract_prevention_paragraphs(xml_data):
    """Extract prevention-related paragraphs from XML"""
    root = ET.fromstring(xml_data)
    prevention_paragraphs = []
    
    try:
        for elem in root.iter():
            tag = elem.tag.lower()
            if tag.endswith("p") and elem.text:
                text = elem.text.strip().lower()
                # Check if paragraph contains prevention keywords
                if any(keyword in text for keyword in ALL_PREVENTION_KEYWORDS):
                    clean = re.sub(r'\s+', ' ', elem.text.strip())
                    if len(clean) > 100:
                        prevention_paragraphs.append(clean)
    except:
        pass
    
    return prevention_paragraphs

def identify_topics(text):
    """Identify prevention topics in text"""
    text_lower = text.lower()
    found_topics = []
    
    for topic in PREVENTION_DATABASE.keys():
        # Check for topic keywords
        topic_lower = topic.replace("_", " ")
        
        if topic_lower in text_lower:
            found_topics.append(topic)
        # Check specific keywords for each topic
        elif topic == "psa_screening" and any(keyword in text_lower for keyword in ["psa", "prostate-specific antigen", "screening", "digital rectal"]):
            found_topics.append(topic)
        elif topic == "multiparametric_mri" and any(keyword in text_lower for keyword in ["mri", "multiparametric", "pi-rads", "magnetic resonance"]):
            found_topics.append(topic)
        elif topic == "five_alpha_reductase_inhibitors" and any(keyword in text_lower for keyword in ["finasteride", "dutasteride", "5-alpha reductase", "5-ari"]):
            found_topics.append(topic)
        elif topic == "vitamin_e_selenium" and any(keyword in text_lower for keyword in ["vitamin e", "selenium", "select trial", "tocopherol"]):
            found_topics.append(topic)
        elif topic == "diet_nutrition_prostate" and any(keyword in text_lower for keyword in ["diet", "nutrition", "lycopene", "tomato", "dietary"]):
            found_topics.append(topic)
        elif topic == "physical_activity_obesity" and any(keyword in text_lower for keyword in ["physical activity", "exercise", "obesity", "body mass index", "bmi"]):
            found_topics.append(topic)
        elif topic == "genetic_familial_prostate" and any(keyword in text_lower for keyword in ["genetic", "familial", "family history", "brca", "hereditary"]):
            found_topics.append(topic)
        elif topic == "racial_ethnic_disparities" and any(keyword in text_lower for keyword in ["african american", "racial", "ethnic", "disparity", "disparities"]):
            found_topics.append(topic)
        elif topic == "ejaculation_frequency" and any(keyword in text_lower for keyword in ["ejaculation", "sexual activity", "ejaculatory"]):
            found_topics.append(topic)
    
    return list(set(found_topics))

# ================= MAIN =================

def main():
    print("="*130)
    print("🔬 PROSTATE CANCER PREVENTION EXTRACTOR - PubMed Central")
    print("="*130)
    print("📚 Searching for: Prostate Adenocarcinoma, Advanced & Metastatic Prostate Cancer Prevention")
    print(f"📊 Target: ~{MAX_ARTICLES * len(SEARCH_QUERIES)} articles across {len(SEARCH_QUERIES)} queries\n")
    
    pmc_ids = search_pmc()
    print(f"\n✅ Found {len(pmc_ids)} unique articles total")
    
    if not pmc_ids:
        print("❌ No articles found. Check your internet connection and try again.")
        return
    
    all_contexts = []
    processed = 0
    errors = 0
    
    for i, pmc_id in enumerate(pmc_ids, 1):
        try:
            print(f"\r📄 Processing {i}/{len(pmc_ids)}: PMC{pmc_id}... ", end='', flush=True)
            
            xml = fetch_xml(pmc_id)
            paragraphs = extract_prevention_paragraphs(xml)
            
            if paragraphs:
                for para in paragraphs:
                    topics = identify_topics(para)
                    if topics:  # Only save if topics found
                        all_contexts.append({
                            "pmc_id": pmc_id,
                            "article_url": f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}",
                            "paragraph": para,
                            "topics_found": topics,
                            "topic_count": len(topics),
                            "topic_details": [get_prevention_details(t) for t in topics if get_prevention_details(t)]
                        })
                processed += 1
            
            # Rate limiting
            time.sleep(0.34)  # ~3 requests per second (NCBI limit)
            
        except Exception as e:
            errors += 1
            if errors <= 5:  # Only show first 5 errors
                print(f"\n⚠️ Error PMC{pmc_id}: {str(e)[:50]}")
    
    print(f"\n\n{'='*130}")
    print(f"✅ Successfully processed: {processed}/{len(pmc_ids)} articles")
    print(f"⚠️ Errors encountered: {errors}")
    print(f"📊 Total prevention contexts extracted: {len(all_contexts)}")
    print(f"{'='*130}\n")
    
    # ================= SAVE DETAILED TEXT FORMAT =================
    print(f"\n💾 Saving DETAILED TEXT format...")
    with open(OUTPUT_TEXT, "w", encoding="utf-8") as f:
        f.write("=" * 130 + "\n")
        f.write("PROSTATE CANCER PREVENTION - DETAILED EXTRACTION FROM PUBMED\n")
        f.write("Includes: Prostate Adenocarcinoma, Advanced & Metastatic Prostate Cancer Prevention\n")
        f.write("=" * 130 + "\n\n")
        
        for i, ctx in enumerate(all_contexts, 1):
            f.write(f"\n{'='*130}\n")
            f.write(f"CONTEXT #{i}\n")
            f.write(f"{'='*130}\n")
            f.write(f"PMC ID: {ctx['pmc_id']}\n")
            f.write(f"Topics Found: {len(ctx['topics_found'])}\n")
            f.write(f"URL: {ctx['article_url']}\n\n")
            
            # Write prevention topic details
            for topic_detail in ctx['topic_details']:
                if topic_detail:
                    f.write(f"\n{'─'*130}\n")
                    f.write(f"📋 TOPIC: {topic_detail.get('name', 'N/A').upper()}\n")
                    f.write(f"{'─'*130}\n")
                    f.write(f"Category: {topic_detail.get('category', 'N/A')}\n")
                    f.write(f"Type: {topic_detail.get('type', 'N/A')}\n\n")
                    
                    if topic_detail.get('description'):
                        f.write(f"📝 DESCRIPTION:\n{topic_detail['description']}\n\n")
                    
                    # Handle nested dictionaries for various fields
                    for key in ['current_guidelines', 'screening_protocol', 'benefits_of_screening', 
                                'harms_of_screening', 'mechanism', 'clinical_trials', 'dietary_patterns',
                                'specific_nutrients', 'foods_to_limit', 'physical_activity', 'obesity',
                                'family_history', 'high_risk_genes', 'african_american_men']:
                        if topic_detail.get(key):
                            f.write(f"📊 {key.replace('_', ' ').upper()}:\n")
                            data = topic_detail[key]
                            if isinstance(data, dict):
                                for subkey, subvalue in data.items():
                                    f.write(f"\n  {subkey.replace('_', ' ').title()}:\n")
                                    if isinstance(subvalue, dict):
                                        for k, v in subvalue.items():
                                            f.write(f"    • {k}: {v}\n")
                                    elif isinstance(subvalue, list):
                                        for item in subvalue:
                                            f.write(f"    • {item}\n")
                                    else:
                                        f.write(f"    {subvalue}\n")
                            else:
                                f.write(f"  {data}\n")
                            f.write("\n")
                    
                    if topic_detail.get('evidence'):
                        f.write(f"📊 EVIDENCE:\n{topic_detail['evidence']}\n\n")
                    
                    if topic_detail.get('recommendations'):
                        f.write(f"✅ RECOMMENDATIONS:\n")
                        if isinstance(topic_detail['recommendations'], list):
                            for rec in topic_detail['recommendations']:
                                f.write(f"  • {rec}\n")
                        else:
                            f.write(f"  {topic_detail['recommendations']}\n")
                        f.write("\n")
            
            f.write(f"\n📌 CONTEXT PARAGRAPH FROM ARTICLE:\n{ctx['paragraph']}\n")
    
    # ================= SAVE DETAILED JSON =================
    print(f"💾 Saving DETAILED JSON format...")
    with open(OUTPUT_DETAILED, "w", encoding="utf-8") as f:
        json.dump(all_contexts, f, indent=2, ensure_ascii=False)
    
    # ================= SAVE SIMPLE JSON =================
    print(f"💾 Saving simple JSON format...")
    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(all_contexts, f, indent=2, ensure_ascii=False)
    
    # ================= STATISTICS =================
    print("\n" + "🎉 " + "="*128)
    print("EXTRACTION COMPLETE")
    print("=" * 130)
    print(f"📝 Text file: {OUTPUT_TEXT}")
    print(f"💾 JSON file: {OUTPUT_JSON}")
    print(f"💾 Detailed JSON: {OUTPUT_DETAILED}")
    print(f"📊 Total contexts: {len(all_contexts)}")
    
    all_topics = []
    for ctx in all_contexts:
        all_topics.extend(ctx['topics_found'])
    
    from collections import Counter
    topic_freq = Counter(all_topics)
    
    print(f"\n🏆 Top 10 Most Mentioned Prevention Topics:")
    for topic, count in topic_freq.most_common(10):
        print(f"  {topic}: {count} times")
    
    print("\n💡 Key Prostate Cancer Prevention Categories:")
    print("  🔍 Screening: PSA, DRE, MRI, shared decision-making")
    print("  💊 Chemoprevention: 5-ARIs (not recommended), Vitamin E/Selenium (harmful)")
    print("  🥗 Diet: Plant-based, Mediterranean diet, limit red meat")
    print("  🏃 Activity: 150+ min/week, maintain healthy weight")
    print("  🧬 Genetic: BRCA1/2, HOXB13, family history assessment")
    print("  ⚖️ Equity: African American men at higher risk")
    print("  👁️ Active Surveillance: Monitor low-risk cancers without treatment")
    
    print("\n⚠️ KEY INSIGHTS:")
    print("  • PSA screening: 20-30% mortality reduction but high overdiagnosis")
    print("  • Shared decision-making essential for screening (age 55-69)")
    print("  • African American men: 70% higher incidence, start screening 40-45")
    print("  • Vitamin E supplements INCREASE prostate cancer risk (17%)")
    print("  • 5-ARIs reduce risk 23-25% but not recommended for prevention")
    print("  • BRCA2 carriers: 5-8x increased risk, screen at age 40-45")
    print("  • Active surveillance appropriate for 50% of low-risk cancers")
    print("  • Physical activity reduces advanced disease risk by 30-50%")
    
    print("\n" + "=" * 130)

if __name__ == "__main__":
    main()