"""
OncoDx — Rule-Based Clinical Chatbot
Asks fixed clinical questions → BioBERT classifies → Medical summary generated
Zero API. Zero internet. Fully offline.

Run: streamlit run app.py
"""

import streamlit as st
import torch
import torch.nn.functional as F
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from pathlib import Path
from datetime import datetime

# ─────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────

MODEL_PATH = Path("D:/bilogy/biobert_output/best_model")
DEVICE     = torch.device("cuda" if torch.cuda.is_available() else "cpu")

CANCER_INFO = {
    "Blood_Cancer":    {"icon": "🩸", "color": "#ef4444", "light": "#fef2f2", "desc": "Leukemia · Lymphoma · Myeloma"},
    "Brain_Cancer":    {"icon": "🧠", "color": "#7c3aed", "light": "#f5f3ff", "desc": "Glioma · Meningioma · Astrocytoma"},
    "Liver_Cancer":    {"icon": "🫀", "color": "#d97706", "light": "#fffbeb", "desc": "HCC · Cholangiocarcinoma"},
    "Lung_Cancer":     {"icon": "🫁", "color": "#2563eb", "light": "#eff6ff", "desc": "NSCLC · SCLC · Adenocarcinoma"},
    "Mouth_Cancer":    {"icon": "🦷", "color": "#db2777", "light": "#fdf2f8", "desc": "OSCC · Oropharyngeal Cancer"},
    "Ovarian_Cancer":  {"icon": "🔬", "color": "#0891b2", "light": "#ecfeff", "desc": "Epithelial · Germ Cell · Stromal"},
    "Protaste_Cancer": {"icon": "⚕️",  "color": "#059669", "light": "#ecfdf5", "desc": "Adenocarcinoma · CRPC · mHSPC"},
    "Skin_Cancer":     {"icon": "☀️",  "color": "#b45309", "light": "#fef9c3", "desc": "Melanoma · BCC · SCC"},
}

# ─────────────────────────────────────────────────────────────
# CLINICAL QUESTIONS — rule-based flow
# Each question has: key, question text, hint, type
# ─────────────────────────────────────────────────────────────

QUESTIONS = [
    {
        "key":   "greeting",
        "bot":   "👋 Hello! I am OncoDx Assistant.\n\nI will ask you a few clinical questions to help identify possible cancer type. Your answers will be analysed by our BioBERT model.\n\nLet's begin. What is the **patient's age and gender**?",
        "hint":  "e.g. 58 year old male",
        "field": "Age / Gender",
    },
    {
        "key":   "chief_complaint",
        "bot":   "Thank you. What are the **main symptoms** the patient is experiencing?",
        "hint":  "e.g. persistent cough, weight loss, fatigue",
        "field": "Chief Complaint / Main Symptoms",
    },
    {
        "key":   "duration",
        "bot":   "How **long** have these symptoms been present?",
        "hint":  "e.g. 3 weeks, 2 months, 6 months",
        "field": "Duration",
    },
    {
        "key":   "severity",
        "bot":   "On a scale of **1 to 10**, how severe are the symptoms? (1 = mild, 10 = severe)\n\nAlso describe the character — is the pain sharp, dull, burning, or pressure-like?",
        "hint":  "e.g. 7/10, sharp stabbing pain in chest",
        "field": "Severity",
    },
    {
        "key":   "associated",
        "bot":   "Are there any **associated symptoms**?\n\nFor example: fever, night sweats, nausea, vomiting, bleeding, swelling, skin changes, urinary changes, bowel changes?",
        "hint":  "e.g. fever, night sweats, swollen lymph nodes",
        "field": "Associated Symptoms",
    },
    {
        "key":   "constitutional",
        "bot":   "Has the patient experienced any **constitutional symptoms**?\n\n• Unintentional weight loss\n• Persistent fatigue or weakness\n• Loss of appetite\n• Recurrent fever",
        "hint":  "e.g. 8kg weight loss over 2 months, severe fatigue",
        "field": "Constitutional Symptoms",
    },
    {
        "key":   "medical_history",
        "bot":   "Does the patient have any **past medical history**?\n\nPrevious cancers, chronic illnesses, surgeries, or hospitalizations?",
        "hint":  "e.g. Type 2 diabetes, hypertension, no previous cancer",
        "field": "Medical History",
    },
    {
        "key":   "family_history",
        "bot":   "Is there any **family history of cancer**?\n\nWhich relative and what type of cancer?",
        "hint":  "e.g. Father had prostate cancer, mother had breast cancer",
        "field": "Family History",
    },
    {
        "key":   "medications",
        "bot":   "What **medications** is the patient currently taking?\n\nInclude any supplements or over-the-counter drugs.",
        "hint":  "e.g. Metformin 500mg, Aspirin 75mg, or none",
        "field": "Current Medications",
    },
    {
        "key":   "lifestyle",
        "bot":   "What is the patient's **lifestyle history**?\n\n• Smoking (how many years, pack/day)\n• Alcohol consumption\n• Occupation\n• Diet",
        "hint":  "e.g. Smoker 20 pack-years, occasional alcohol, office worker",
        "field": "Lifestyle",
    },
    {
        "key":   "labs",
        "bot":   "Are there any **lab results or investigations** available?\n\nBlood tests, imaging, biopsy, PSA, tumour markers, CBC, etc.?",
        "hint":  "e.g. PSA 22 ng/mL elevated, CT shows 4cm mass in right lung",
        "field": "Lab Results / Investigations",
    },
    {
        "key":   "additional",
        "bot":   "Is there **anything else** clinically relevant you would like to add?\n\nAny other symptoms, recent travel, allergies, or concerns?",
        "hint":  "e.g. Recent weight loss, skin rash, no known allergies",
        "field": "Additional Information",
    },
]

DONE_MSG = "✅ Thank you — I have all the information needed.\n\nClick **🔬 Analyse & Generate Report** to run BioBERT classification and generate the medical summary."

# ─────────────────────────────────────────────────────────────
# RED FLAG DETECTOR — pure rule-based, no API
# ─────────────────────────────────────────────────────────────

RED_FLAG_KEYWORDS = [
    "weight loss","hemoptysis","blood in","bleeding","haematuria","hematuria",
    "dysphagia","difficulty swallowing","seizure","vision loss","paralysis",
    "bone pain","night sweats","persistent fever","lymph node","swollen gland",
    "mass","lump","lesion","ulcer","non-healing","jaundice","ascites",
    "psa elevated","psa high","elevated psa","ca-125","afp elevated",
    "biopsy","malignant","carcinoma","metastasis","metastatic",
]

URGENCY_KEYWORDS = {
    "emergency": ["severe bleeding","cord compression","seizure","paralysis",
                  "airway obstruction","acute","emergency"],
    "urgent":    ["hemoptysis","bone pain","weight loss","elevated psa","night sweats",
                  "persistent","swollen lymph","jaundice","mass","lump"],
}

def detect_red_flags(text):
    tl = text.lower()
    return [kw for kw in RED_FLAG_KEYWORDS if kw in tl]

def detect_urgency(text):
    tl = text.lower()
    for kw in URGENCY_KEYWORDS["emergency"]:
        if kw in tl: return "EMERGENCY"
    for kw in URGENCY_KEYWORDS["urgent"]:
        if kw in tl: return "URGENT"
    return "ROUTINE"

def detect_investigations(answers):
    """Rule-based investigation suggestions based on symptoms."""
    combined = " ".join(answers.values()).lower()
    inv = []
    if any(w in combined for w in ["psa","prostate","urinary","urine"]): inv.append("PSA blood test")
    if any(w in combined for w in ["cough","lung","chest","hemoptysis","breathing"]): inv.append("CT Chest")
    if any(w in combined for w in ["liver","jaundice","abdominal","hepatic"]): inv.append("LFTs + AFP + Liver USS")
    if any(w in combined for w in ["blood","fatigue","pallor","bruising","lymph"]): inv.append("Full Blood Count + Blood film")
    if any(w in combined for w in ["brain","headache","seizure","neurological","vision"]): inv.append("MRI Brain with contrast")
    if any(w in combined for w in ["skin","mole","lesion","melanoma","pigmented"]): inv.append("Dermatology referral + Dermoscopy")
    if any(w in combined for w in ["mouth","oral","tongue","swallowing","throat"]): inv.append("OPG + Panendoscopy")
    if any(w in combined for w in ["ovarian","pelvic","bloating","ca-125"]): inv.append("CA-125 + Transvaginal USS")
    if any(w in combined for w in ["weight loss","fatigue","appetite"]): inv.append("Full metabolic panel")
    if any(w in combined for w in ["biopsy","mass","lump","lesion"]): inv.append("Tissue biopsy + Histopathology")
    if not inv: inv.append("Full Blood Count", )
    return inv if inv else ["Clinical review + Full Blood Count"]

# ─────────────────────────────────────────────────────────────
# GENERATE MEDICAL SUMMARY — pure Python, no API
# ─────────────────────────────────────────────────────────────

def generate_summary(answers, results):
    top      = results[0]
    pct      = top["prob"] * 100
    info     = CANCER_INFO.get(top["cancer"], {"icon":"🔬","desc":""})
    combined = " ".join(str(v) for v in answers.values())
    flags    = detect_red_flags(combined)
    urgency  = detect_urgency(combined)
    inv      = detect_investigations(answers)

    # Build overview sentence
    age_gender = answers.get("greeting", "Patient")
    complaint  = answers.get("chief_complaint", "multiple symptoms")
    duration   = answers.get("duration", "unspecified duration")
    overview   = f"{age_gender} presenting with {complaint} for {duration}."

    return {
        "overview":              overview,
        "chief_complaint":       answers.get("chief_complaint", "—"),
        "duration":              answers.get("duration", "—"),
        "severity":              answers.get("severity", "—"),
        "associated_symptoms":   answers.get("associated", "—"),
        "constitutional":        answers.get("constitutional", "—"),
        "medical_history":       answers.get("medical_history", "—"),
        "family_history":        answers.get("family_history", "—"),
        "medications":           answers.get("medications", "—"),
        "lifestyle":             answers.get("lifestyle", "—"),
        "labs":                  answers.get("labs", "—"),
        "additional":            answers.get("additional", "—"),
        "red_flags":             flags,
        "investigations":        inv,
        "urgency":               urgency,
        "classification":        top["cancer"].replace("_", " "),
        "confidence":            f"{pct:.1f}%",
        "cancer_desc":           info["desc"],
        "generated_at":          datetime.now().strftime("%Y-%m-%d %H:%M"),
    }

# ─────────────────────────────────────────────────────────────
# MODEL LOADER
# ─────────────────────────────────────────────────────────────

@st.cache_resource(show_spinner=False)
def load_model():
    if not MODEL_PATH.exists():
        return None, None
    tok = AutoTokenizer.from_pretrained(str(MODEL_PATH))
    mdl = AutoModelForSequenceClassification.from_pretrained(str(MODEL_PATH))
    mdl.to(DEVICE).eval()
    return tok, mdl

def classify(text, tok, mdl):
    enc = tok(text, max_length=512, padding="max_length",
               truncation=True, return_tensors="pt")
    enc = {k: v.to(DEVICE) for k, v in enc.items()}
    with torch.no_grad():
        probs = F.softmax(mdl(**enc).logits, dim=-1)[0]
    id2l = mdl.config.id2label
    return sorted(
        [{"cancer": id2l[i], "prob": float(probs[i])} for i in range(len(probs))],
        key=lambda x: x["prob"], reverse=True
    )

# ─────────────────────────────────────────────────────────────
# PAGE CONFIG
# ─────────────────────────────────────────────────────────────

st.set_page_config(
    page_title = "OncoDx Clinical Chat",
    page_icon  = "🔬",
    layout     = "wide",
)

# ─────────────────────────────────────────────────────────────
# CSS
# ─────────────────────────────────────────────────────────────

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Plus+Jakarta+Sans:wght@300;400;500;600;700&family=Lora:ital,wght@0,400;0,600&display=swap');

html,body,[class*="css"]{ font-family:'Plus Jakarta Sans',sans-serif !important; }
.stApp{ background:#f0f4f8 !important; }
#MainMenu,footer,header{ visibility:hidden; }
.block-container{ padding:1.5rem 2rem 2rem !important; max-width:1300px !important; }

/* Navbar */
.navbar{
    display:flex; align-items:center; justify-content:space-between;
    background:white; border-radius:18px; padding:1rem 1.75rem;
    margin-bottom:1.5rem;
    box-shadow:0 2px 8px rgba(0,0,0,0.07),0 0 0 1px rgba(0,0,0,0.04);
}
.nav-logo{ font-family:'Lora',serif; font-size:1.7rem; font-weight:600; color:#1e293b; }
.nav-logo span{ color:#2563eb; }
.pills{ display:flex; gap:0.5rem; flex-wrap:wrap; }
.pill{ padding:0.28rem 0.85rem; border-radius:100px; font-size:0.7rem; font-weight:700; }
.pb{ background:#eff6ff; color:#1d4ed8; }
.pg{ background:#f0fdf4; color:#15803d; }
.pp{ background:#faf5ff; color:#6d28d9; }
.po{ background:#fff7ed; color:#c2410c; }

/* Progress steps */
.steps{
    display:flex; gap:0; margin-bottom:1.25rem;
    background:white; border-radius:14px; overflow:hidden;
    box-shadow:0 2px 8px rgba(0,0,0,0.06),0 0 0 1px rgba(0,0,0,0.04);
}
.step{
    flex:1; padding:0.5rem 0.25rem; text-align:center;
    font-size:0.65rem; font-weight:600; color:#94a3b8;
    border-right:1px solid #f1f5f9; letter-spacing:0.03em;
}
.step:last-child{ border-right:none; }
.step.done{ background:#f0fdf4; color:#15803d; }
.step.active{ background:#eff6ff; color:#1d4ed8; font-weight:700; }

/* Chat */
.chat-box{
    background:white; border-radius:18px; padding:1.25rem 1.1rem;
    min-height:400px; max-height:480px; overflow-y:auto;
    display:flex; flex-direction:column; gap:0.85rem;
    box-shadow:0 2px 8px rgba(0,0,0,0.06),0 0 0 1px rgba(0,0,0,0.04);
    margin-bottom:0.9rem;
}
.bwrap{ display:flex; align-items:flex-end; gap:0.55rem; }
.bwrap.u{ flex-direction:row-reverse; }
.av{
    width:34px; height:34px; border-radius:50%; flex-shrink:0;
    display:flex; align-items:center; justify-content:center; font-size:0.95rem;
}
.av-ai  { background:linear-gradient(135deg,#eff6ff,#dbeafe); }
.av-usr { background:linear-gradient(135deg,#f0fdf4,#bbf7d0); }
.bub{
    max-width:82%; padding:0.75rem 1rem; border-radius:16px;
    font-size:0.86rem; line-height:1.7; color:#1e293b;
}
.bub-ai { background:#f8fafc; border:1px solid #e2e8f0; border-bottom-left-radius:4px; }
.bub-usr{ background:linear-gradient(135deg,#2563eb,#3b82f6); color:white; border-bottom-right-radius:4px; }
.bt{ font-size:0.62rem; color:#cbd5e1; margin-top:3px; padding:0 3px; }
.bwrap.u .bt{ text-align:right; }

/* Input */
.stTextInput input{
    background:white !important; border:2px solid #e2e8f0 !important;
    border-radius:12px !important; color:#1e293b !important;
    font-size:0.9rem !important; padding:0.65rem 1rem !important;
    font-family:'Plus Jakarta Sans',sans-serif !important;
}
.stTextInput input:focus{
    border-color:#2563eb !important;
    box-shadow:0 0 0 3px rgba(37,99,235,0.1) !important;
}

/* Buttons */
.stButton>button{
    border-radius:12px !important; border:none !important;
    font-family:'Plus Jakarta Sans',sans-serif !important;
    font-weight:700 !important; font-size:0.88rem !important;
    padding:0.62rem 1.4rem !important; transition:all 0.2s !important;
    color:white !important;
}
.stButton>button:hover{
    transform:translateY(-1px) !important;
    box-shadow:0 6px 20px rgba(37,99,235,0.25) !important;
}

/* Result card */
.rcard{
    border-radius:18px; padding:1.4rem; margin-bottom:1rem;
    box-shadow:0 2px 8px rgba(0,0,0,0.07),0 0 0 1px rgba(0,0,0,0.04);
}
.rname{ font-family:'Lora',serif; font-size:1.65rem; font-weight:600; color:#1e293b; margin:0.2rem 0 0.1rem; }
.rdesc{ font-size:0.78rem; color:#64748b; margin-bottom:0.75rem; }
.rconf{ font-size:2.8rem; font-weight:700; line-height:1; }
.rcl  { font-size:0.67rem; color:#94a3b8; text-transform:uppercase; letter-spacing:0.07em; }
.slbl { font-size:0.67rem; font-weight:700; letter-spacing:0.1em; text-transform:uppercase; color:#94a3b8; margin-bottom:0.55rem; }

/* Prob bars */
.prow{ display:flex; align-items:center; gap:0.65rem; margin-bottom:0.45rem; }
.plbl{ font-size:0.77rem; color:#475569; min-width:155px; white-space:nowrap; font-weight:500; }
.pbg { flex:1; height:7px; background:#f1f5f9; border-radius:100px; overflow:hidden; }
.pf  { height:100%; border-radius:100px; }
.ppct{ font-size:0.73rem; font-weight:600; min-width:40px; text-align:right; }

/* Summary */
.scard{
    background:white; border-radius:18px; padding:1.4rem;
    box-shadow:0 2px 8px rgba(0,0,0,0.06),0 0 0 1px rgba(0,0,0,0.04);
    margin-bottom:1rem;
}
.stitle{
    font-family:'Lora',serif; font-size:1.05rem; font-weight:600; color:#1e293b;
    margin-bottom:0.9rem; padding-bottom:0.7rem; border-bottom:1px solid #f1f5f9;
    display:flex; align-items:center; gap:0.5rem;
}
.srow{ display:flex; gap:0.75rem; margin-bottom:0.6rem; font-size:0.8rem; line-height:1.55; }
.skey{ font-weight:700; color:#94a3b8; min-width:130px; flex-shrink:0; font-size:0.68rem; text-transform:uppercase; letter-spacing:0.05em; padding-top:2px; }
.sval{ color:#1e293b; }
.tag{ display:inline-block; padding:0.12rem 0.5rem; border-radius:6px; font-size:0.7rem; font-weight:600; margin:2px; }
.t-blue { background:#eff6ff; color:#1d4ed8; }
.t-red  { background:#fef2f2; color:#dc2626; }
.t-green{ background:#f0fdf4; color:#15803d; }
.ubadge{ padding:0.25rem 0.75rem; border-radius:100px; font-size:0.7rem; font-weight:700; letter-spacing:0.06em; text-transform:uppercase; }
.u-ROUTINE  { background:#f0fdf4; color:#15803d; }
.u-URGENT   { background:#fff7ed; color:#c2410c; }
.u-EMERGENCY{ background:#fef2f2; color:#dc2626; }

.disclaimer{
    background:#fff7ed; border:1px solid #fed7aa; border-radius:12px;
    padding:0.7rem 1rem; font-size:0.74rem; color:#92400e;
    line-height:1.5; margin-bottom:1rem;
}
</style>
""", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────
# SESSION STATE
# ─────────────────────────────────────────────────────────────

defaults = {
    "q_index":   0,          # current question index
    "answers":   {},         # collected answers
    "chat":      [],         # chat display history
    "started":   False,
    "finished":  False,
    "results":   None,
    "summary":   None,
}
for k, v in defaults.items():
    if k not in st.session_state:
        st.session_state[k] = v

# ─────────────────────────────────────────────────────────────
# LOAD MODEL
# ─────────────────────────────────────────────────────────────

with st.spinner("Loading BioBERT..."):
    tok, mdl = load_model()

if mdl is None:
    st.error(f"❌ Model not found at `{MODEL_PATH}` — run `biobert_finetune.py` first.")
    st.stop()

# ─────────────────────────────────────────────────────────────
# NAVBAR
# ─────────────────────────────────────────────────────────────

gpu = torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU"

st.markdown(f"""
<div class="navbar">
  <div class="nav-logo">Onco<span>Dx</span></div>
  <div class="pills">
    <span class="pill pb">🔬 BioBERT · 91% Accuracy</span>
    <span class="pill pg">● {gpu}</span>
    <span class="pill pp">8 Cancer Types</span>
    <span class="pill po">⚡ Fully Offline · No API</span>
  </div>
</div>
""", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────
# PROGRESS STEPS
# ─────────────────────────────────────────────────────────────

total_q   = len(QUESTIONS)
cur_q     = st.session_state.q_index
step_html = '<div class="steps">'
labels    = ["Age","Symptoms","Duration","Severity","Associated",
             "Constitutional","Med Hx","Family Hx","Meds","Lifestyle","Labs","Other"]
for i, lbl in enumerate(labels):
    if i < cur_q:
        cls = "done"
    elif i == cur_q:
        cls = "active"
    else:
        cls = ""
    icon = "✓ " if i < cur_q else ""
    step_html += f'<div class="step {cls}">{icon}{lbl}</div>'
step_html += "</div>"
st.markdown(step_html, unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────
# LAYOUT
# ─────────────────────────────────────────────────────────────

left, right = st.columns([1.05, 1], gap="large")

# ══════════════════════════════════
# LEFT — CHAT
# ══════════════════════════════════
with left:

    st.markdown('<div class="slbl">💬 Clinical Interview</div>', unsafe_allow_html=True)
    st.markdown('<div class="disclaimer">⚠️ <strong>Screening only.</strong> All findings must be verified by a qualified oncologist. This tool does not replace professional medical diagnosis.</div>', unsafe_allow_html=True)

    # ── Render chat bubbles ──
    if st.session_state.chat:
        html = '<div class="chat-box">'
        for m in st.session_state.chat:
            txt = m["text"]
            t   = m.get("time", "")
            if m["role"] == "bot":
                html += f'<div class="bwrap"><div class="av av-ai">🔬</div><div><div class="bub bub-ai">{txt}</div><div class="bt">{t}</div></div></div>'
            else:
                html += f'<div class="bwrap u"><div class="av av-usr">👤</div><div><div class="bub bub-usr">{txt}</div><div class="bt">{t}</div></div></div>'
        html += "</div>"
        st.markdown(html, unsafe_allow_html=True)
    else:
        st.markdown("""
        <div class="chat-box" style="align-items:center;justify-content:center;">
          <div style="text-align:center;color:#94a3b8;">
            <div style="font-size:2.5rem;margin-bottom:0.5rem;">🩺</div>
            <div style="font-family:'Lora',serif;font-size:1rem;color:#475569;">Clinical interview ready</div>
            <div style="font-size:0.78rem;margin-top:0.3rem;">Click Start to begin</div>
          </div>
        </div>""", unsafe_allow_html=True)

    # ── Controls ──
    now = datetime.now().strftime("%H:%M")

    if not st.session_state.started:
        if st.button("▶  Start Clinical Interview", use_container_width=True):
            st.session_state.started = True
            first_q = QUESTIONS[0]
            st.session_state.chat.append({"role":"bot","text":first_q["bot"],"time":now})
            st.rerun()

    elif not st.session_state.finished:
        # Show hint for current question
        if cur_q < total_q:
            hint = QUESTIONS[cur_q]["hint"]
            st.markdown(f'<div style="font-size:0.72rem;color:#94a3b8;margin-bottom:0.4rem;">💡 Hint: {hint}</div>', unsafe_allow_html=True)

        with st.form("reply_form", clear_on_submit=True):
            c1, c2 = st.columns([5, 1])
            with c1:
                user_input = st.text_input("", placeholder="Type your answer here...", label_visibility="collapsed")
            with c2:
                sent = st.form_submit_button("Send")

        if sent and user_input.strip():
            # Save answer
            q_key = QUESTIONS[cur_q]["key"]
            st.session_state.answers[q_key] = user_input
            st.session_state.chat.append({"role":"user","text":user_input,"time":now})

            next_q = cur_q + 1
            st.session_state.q_index = next_q

            if next_q < total_q:
                # Ask next question
                next_question = QUESTIONS[next_q]
                st.session_state.chat.append({"role":"bot","text":next_question["bot"],"time":now})
            else:
                # All questions done
                st.session_state.finished = True
                st.session_state.chat.append({"role":"bot","text":DONE_MSG,"time":now})

            st.rerun()

    else:
        # Finished — show action buttons
        st.success(f"✅ All {total_q} questions answered!")
        c1, c2 = st.columns(2)
        with c1:
            if st.button("🔬 Analyse & Generate Report", use_container_width=True):
                # Build full text from all answers
                full_text = " ".join([
                    f"{QUESTIONS[i]['field']}: {st.session_state.answers.get(QUESTIONS[i]['key'], '')}"
                    for i in range(total_q)
                    if st.session_state.answers.get(QUESTIONS[i]["key"])
                ])
                with st.spinner("Running BioBERT classification..."):
                    st.session_state.results = classify(full_text, tok, mdl)
                with st.spinner("Generating medical summary..."):
                    st.session_state.summary = generate_summary(
                        st.session_state.answers,
                        st.session_state.results
                    )
                st.rerun()

        with c2:
            if st.button("🔄 New Interview", use_container_width=True):
                for k in defaults:
                    st.session_state[k] = defaults[k]
                st.rerun()


# ══════════════════════════════════
# RIGHT — RESULTS + SUMMARY
# ══════════════════════════════════
with right:

    if st.session_state.results and st.session_state.summary:
        results = st.session_state.results
        s       = st.session_state.summary
        top     = results[0]
        info    = CANCER_INFO.get(top["cancer"], {"icon":"🔬","color":"#2563eb","light":"#eff6ff","desc":""})
        pct     = top["prob"] * 100

        # ── Classification ──
        st.markdown('<div class="slbl">🎯 BioBERT Classification</div>', unsafe_allow_html=True)
        st.markdown(f"""
        <div class="rcard" style="background:linear-gradient(135deg,white 55%,{info['light']});
             border:1px solid {info['color']}20; border-left:5px solid {info['color']};">
          <div style="display:flex;justify-content:space-between;align-items:flex-start;">
            <div>
              <div class="slbl">Primary Classification</div>
              <div class="rname">{info['icon']} {top['cancer'].replace('_',' ')}</div>
              <div class="rdesc">{info['desc']}</div>
            </div>
            <div style="text-align:right">
              <div class="rconf" style="color:{info['color']}">{pct:.1f}%</div>
              <div class="rcl">Confidence</div>
            </div>
          </div>
          {"<div style='margin-top:0.75rem;background:#fff7ed;border:1px solid #fed7aa;border-radius:8px;padding:0.45rem 0.75rem;font-size:0.74rem;color:#92400e;'>⚠️ Low confidence — review with clinician</div>" if pct < 60 else ""}
        </div>
        """, unsafe_allow_html=True)

        # ── Probability bars ──
        st.markdown('<div class="slbl">All Probabilities</div>', unsafe_allow_html=True)
        bars = ""
        for i, r in enumerate(results):
            ci   = CANCER_INFO.get(r["cancer"], {"icon":"🔬","color":"#cbd5e1"})
            p    = r["prob"] * 100
            col  = ci["color"] if i == 0 else "#cbd5e1"
            bold = "font-weight:700;" if i == 0 else ""
            bars += f"""<div class="prow">
              <div class="plbl" style="{bold}">{ci['icon']} {r['cancer'].replace('_',' ')}</div>
              <div class="pbg"><div class="pf" style="width:{p}%;background:{col};"></div></div>
              <div class="ppct" style="color:{col if i==0 else '#94a3b8'};{bold}">{p:.1f}%</div>
            </div>"""
        st.markdown(f'<div style="background:white;border-radius:16px;padding:1.1rem 1.25rem;box-shadow:0 2px 8px rgba(0,0,0,0.06),0 0 0 1px rgba(0,0,0,0.04);margin-bottom:1rem;">{bars}</div>', unsafe_allow_html=True)

        # ── Medical Summary ──
        urg = s.get("urgency", "ROUTINE")
        def tags(items, cls="t-blue"):
            if not items: return "<span style='color:#94a3b8;font-size:0.76rem;'>None reported</span>"
            if isinstance(items, str): return f'<span class="tag {cls}">{items}</span>'
            return " ".join(f'<span class="tag {cls}">{i}</span>' for i in items)

        st.markdown('<div class="slbl">📋 Medical Intake Summary</div>', unsafe_allow_html=True)
        st.markdown(f"""
        <div class="scard">
          <div class="stitle">
            📄 Clinical Report — {s['generated_at']}
            <span class="ubadge u-{urg}" style="margin-left:auto;">{urg}</span>
          </div>
          <div class="srow"><div class="skey">Overview</div><div class="sval">{s['overview']}</div></div>
          <div class="srow"><div class="skey">Chief Complaint</div><div class="sval">{s['chief_complaint']}</div></div>
          <div class="srow"><div class="skey">Duration</div><div class="sval">{s['duration']}</div></div>
          <div class="srow"><div class="skey">Severity</div><div class="sval">{s['severity']}</div></div>
          <div class="srow"><div class="skey">Associated</div><div class="sval">{s['associated_symptoms']}</div></div>
          <div class="srow"><div class="skey">Constitutional</div><div class="sval">{s['constitutional']}</div></div>
          <div class="srow"><div class="skey">Medical History</div><div class="sval">{s['medical_history']}</div></div>
          <div class="srow"><div class="skey">Family History</div><div class="sval">{s['family_history']}</div></div>
          <div class="srow"><div class="skey">Medications</div><div class="sval">{s['medications']}</div></div>
          <div class="srow"><div class="skey">Lifestyle</div><div class="sval">{s['lifestyle']}</div></div>
          <div class="srow"><div class="skey">Labs / Imaging</div><div class="sval">{s['labs']}</div></div>
          <div class="srow"><div class="skey">Red Flags</div><div class="sval">{tags(s['red_flags'], 't-red')}</div></div>
          <div class="srow"><div class="skey">Investigations</div><div class="sval">{tags(s['investigations'], 't-green')}</div></div>
        </div>
        """, unsafe_allow_html=True)

        # ── Download ──
        report_txt = f"""ONCODX CLINICAL INTAKE REPORT
Generated  : {s['generated_at']}
Model      : BioBERT Fine-tuned | Accuracy: 91%
{'='*55}
CLASSIFICATION : {s['classification']} ({s['confidence']} confidence)
URGENCY        : {urg}
{'='*55}
PATIENT OVERVIEW:
{s['overview']}

CHIEF COMPLAINT  : {s['chief_complaint']}
DURATION         : {s['duration']}
SEVERITY         : {s['severity']}
ASSOCIATED SX    : {s['associated_symptoms']}
CONSTITUTIONAL   : {s['constitutional']}
MEDICAL HISTORY  : {s['medical_history']}
FAMILY HISTORY   : {s['family_history']}
MEDICATIONS      : {s['medications']}
LIFESTYLE        : {s['lifestyle']}
LABS / IMAGING   : {s['labs']}
ADDITIONAL       : {s['additional']}

RED FLAGS        : {', '.join(s['red_flags']) if s['red_flags'] else 'None identified'}
INVESTIGATIONS   : {', '.join(s['investigations'])}
{'='*55}
ALL CANCER PROBABILITIES:
""" + "\n".join(
    f"  {r['cancer'].replace('_',' '):<25} {r['prob']*100:.1f}%"
    for r in results
) + f"""
{'='*55}
DISCLAIMER: For clinical screening only.
Not a substitute for professional medical diagnosis."""

        st.download_button(
            "⬇️  Download Full Report",
            data      = report_txt,
            file_name = f"oncodx_report_{datetime.now().strftime('%Y%m%d_%H%M')}.txt",
            mime      = "text/plain",
            use_container_width = True,
        )

    else:
        st.markdown('<div class="slbl">Result</div>', unsafe_allow_html=True)
        st.markdown("""
        <div style="background:white;border-radius:18px;padding:3.5rem 2rem;text-align:center;
             box-shadow:0 2px 8px rgba(0,0,0,0.06),0 0 0 1px rgba(0,0,0,0.04);">
          <div style="font-size:3rem;margin-bottom:0.75rem;">🩺</div>
          <div style="font-family:'Lora',serif;font-size:1.15rem;color:#475569;margin-bottom:0.5rem;">
            Awaiting consultation
          </div>
          <div style="font-size:0.82rem;color:#94a3b8;line-height:1.7;max-width:230px;margin:0 auto;">
            Complete the clinical interview on the left, then click Analyse
          </div>
          <div style="margin-top:1.5rem;display:flex;justify-content:center;gap:0.5rem;flex-wrap:wrap;">
            <span style="background:#eff6ff;color:#1d4ed8;padding:0.22rem 0.7rem;border-radius:100px;font-size:0.68rem;font-weight:700;">Zero API</span>
            <span style="background:#f0fdf4;color:#15803d;padding:0.22rem 0.7rem;border-radius:100px;font-size:0.68rem;font-weight:700;">Fully Offline</span>
            <span style="background:#faf5ff;color:#6d28d9;padding:0.22rem 0.7rem;border-radius:100px;font-size:0.68rem;font-weight:700;">12 Questions</span>
          </div>
        </div>
        """, unsafe_allow_html=True)