"""
═══════════════════════════════════════════════════════════════
  Cancer Classification API — FastAPI Backend
  
  Serves BioBERT predictions for the web UI
  
  Run: python api.py
  Docs: http://localhost:8000/docs
═══════════════════════════════════════════════════════════════
"""

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel
import torch
import torch.nn.functional as F
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from pathlib import Path
import time
import json

# ─────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────

MODEL_PATH = Path("D:/bilogy/biobert_output1/best_model")
DEVICE     = torch.device("cuda" if torch.cuda.is_available() else "cpu")

CANCER_INFO = {
    "Blood_Cancer":   {"icon": "🩸", "desc": "Leukemia, Lymphoma, Myeloma"},
    "Brain_Cancer":   {"icon": "🧠", "desc": "Glioma, Meningioma, Astrocytoma"},
    "Liver_Cancer":   {"icon": "🫀", "desc": "HCC, Cholangiocarcinoma"},
    "Lung_Cancer":    {"icon": "🫁", "desc": "NSCLC, SCLC, Adenocarcinoma"},
    "Mouth_Cancer":   {"icon": "🦷", "desc": "OSCC, Oropharyngeal Cancer"},
    "Ovarian_Cancer": {"icon": "🔬", "desc": "Epithelial, Germ Cell, Stromal"},
    "Protaste_Cancer":{"icon": "⚕️",  "desc": "Adenocarcinoma, CRPC, mHSPC"},
    "Skin_Cancer":    {"icon": "🌡️", "desc": "Melanoma, BCC, SCC"},
}

# ─────────────────────────────────────────────────────────────
# APP
# ─────────────────────────────────────────────────────────────

app = FastAPI(
    title       = "Cancer Classification API",
    description = "BioBERT-powered cancer type classification",
    version     = "1.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins  = ["*"],
    allow_methods  = ["*"],
    allow_headers  = ["*"],
)

# ─────────────────────────────────────────────────────────────
# MODEL LOADING (once at startup)
# ─────────────────────────────────────────────────────────────

tokenizer = None
model     = None

@app.on_event("startup")
async def load_model():
    global tokenizer, model

    if not MODEL_PATH.exists():
        print(f"❌ Model not found at {MODEL_PATH}")
        print("   Run biobert_finetune.py first!")
        return

    print(f"Loading BioBERT from {MODEL_PATH}...")
    tokenizer = AutoTokenizer.from_pretrained(str(MODEL_PATH))
    model     = AutoModelForSequenceClassification.from_pretrained(str(MODEL_PATH))
    model     = model.to(DEVICE)
    model.eval()

    gpu = torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU"
    print(f"✅ Model ready on {gpu}")


# ─────────────────────────────────────────────────────────────
# SCHEMAS
# ─────────────────────────────────────────────────────────────

class PredictRequest(BaseModel):
    text: str
    top_k: int = 3

class CancerResult(BaseModel):
    cancer:      str
    confidence:  float
    percentage:  float
    icon:        str
    description: str

class PredictResponse(BaseModel):
    prediction:   str
    confidence:   float
    percentage:   float
    all_results:  list
    inference_ms: float
    text_length:  int
    warning:      str = ""


# ─────────────────────────────────────────────────────────────
# ENDPOINTS
# ─────────────────────────────────────────────────────────────

@app.get("/")
async def root():
    return FileResponse("ui.html")


@app.get("/health")
async def health():
    return {
        "status":    "ready" if model else "model_not_loaded",
        "device":    str(DEVICE),
        "gpu":       torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU",
        "model":     str(MODEL_PATH),
    }


@app.post("/predict", response_model=PredictResponse)
async def predict(req: PredictRequest):
    if model is None or tokenizer is None:
        raise HTTPException(503, "Model not loaded — run biobert_finetune.py first")

    text = req.text.strip()
    if len(text.split()) < 5:
        raise HTTPException(400, "Text too short — provide at least 5 words")

    t0 = time.time()

    # Tokenize
    inputs = tokenizer(
        text,
        max_length     = 512,
        padding        = "max_length",
        truncation     = True,
        return_tensors = "pt",
    )
    inputs = {k: v.to(DEVICE) for k, v in inputs.items()}

    # Inference
    with torch.no_grad():
        outputs = model(**inputs)
        probs   = F.softmax(outputs.logits, dim=-1)[0]

    inference_ms = (time.time() - t0) * 1000
    id2label     = model.config.id2label

    # Build results sorted by confidence
    all_results = sorted([
        {
            "cancer":      id2label[i],
            "confidence":  round(probs[i].item(), 4),
            "percentage":  round(probs[i].item() * 100, 2),
            "icon":        CANCER_INFO.get(id2label[i], {}).get("icon", "🔬"),
            "description": CANCER_INFO.get(id2label[i], {}).get("desc", ""),
        }
        for i in range(len(probs))
    ], key=lambda x: x["confidence"], reverse=True)

    top = all_results[0]

    warning = ""
    if top["confidence"] < 0.60:
        warning = "Low confidence — text may not contain sufficient clinical information."

    return PredictResponse(
        prediction   = top["cancer"],
        confidence   = top["confidence"],
        percentage   = top["percentage"],
        all_results  = all_results,
        inference_ms = round(inference_ms, 1),
        text_length  = len(text.split()),
        warning      = warning,
    )


@app.get("/model-info")
async def model_info():
    if model is None:
        raise HTTPException(503, "Model not loaded")
    return {
        "model_name":   "dmis-lab/biobert-base-cased-v1.2 (fine-tuned)",
        "cancer_types": list(CANCER_INFO.keys()),
        "accuracy":     "89.77%",
        "parameters":   "108M",
        "device":       str(DEVICE),
        "cancer_info":  CANCER_INFO,
    }


# ─────────────────────────────────────────────────────────────
# RUN
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import uvicorn
    print("═" * 55)
    print("  Cancer Classification API")
    print("  http://localhost:8000       ← Web UI")
    print("  http://localhost:8000/docs  ← API Docs")
    print("═" * 55)
    uvicorn.run("api:app", host="0.0.0.0", port=8000, reload=False)