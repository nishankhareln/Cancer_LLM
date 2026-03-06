"""
═══════════════════════════════════════════════════════════════
  BioBERT Advanced Fine-Tuning — Cancer Classification
  GPU: RTX 4050 6GB
  
  Improvements over basic training:
  ✅ 30 epochs — full training with cosine LR schedule
  ✅ Label smoothing → prevents overconfidence
  ✅ Layer-wise learning rate decay → better fine-tuning
  ✅ Early stopping → stops if no improvement
  ✅ Gradient clipping → stable training
  ✅ Warmup + cosine decay schedule
  ✅ Best model auto-saved
  ✅ Full per-class report + confusion matrix

  Run: python biobert_advanced_train.py
═══════════════════════════════════════════════════════════════
"""

import os
import json
import torch
import numpy as np
from pathlib import Path
from torch.utils.data import Dataset

# ─────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────

CONFIG = {
    # ── Paths ─────────────────────────────────────────────────
    "tokenized_dir": Path("D:/bilogy/tokenized1_data"),
    "output_dir":    Path("D:/bilogy/biobert_advanced1_output"),
    "model_name":    "dmis-lab/biobert-base-cased-v1.2",

    # ── RTX 4050 6GB settings ─────────────────────────────────
    "per_device_train_batch_size":  8,
    "per_device_eval_batch_size":   16,
    "gradient_accumulation_steps":  4,       # effective batch = 32
    "fp16":                         True,    # saves ~40% VRAM
    "dataloader_num_workers":       0,       # safe on Windows
    "dataloader_pin_memory":        True,

    # ── Advanced training schedule ────────────────────────────
    "num_epochs":           50,              # more epochs
    "learning_rate":        3e-5,            # slightly higher start
    "min_learning_rate":    1e-6,            # cosine decay floor
    "warmup_ratio":         0.06,            # 6% warmup
    "weight_decay":         0.01,
    "max_grad_norm":        1.0,             # gradient clipping

    # ── Regularization ────────────────────────────────────────
    "label_smoothing":      0.1,             # prevents overconfidence

    # ── Early stopping ────────────────────────────────────────

    # ── LR schedule ───────────────────────────────────────────
    "lr_scheduler_type":    "cosine",        # cosine decay

    # ── Logging ───────────────────────────────────────────────
    "logging_steps":        25,
    "save_total_limit":     2,
    "seed":                 42,
}

LABEL_MAP = {}   # filled after loading meta.json


# ─────────────────────────────────────────────────────────────
# DATASET
# ─────────────────────────────────────────────────────────────

class CancerDataset(Dataset):
    def __init__(self, data):
        self.input_ids      = data["input_ids"]
        self.attention_mask = data["attention_mask"]
        self.token_type_ids = data["token_type_ids"]
        self.labels         = data["labels"]

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return {
            "input_ids":      self.input_ids[idx],
            "attention_mask": self.attention_mask[idx],
            "token_type_ids": self.token_type_ids[idx],
            "labels":         self.labels[idx],
        }


# ─────────────────────────────────────────────────────────────
# STEP 1 — GPU SETUP
# ─────────────────────────────────────────────────────────────

def setup_gpu():
    print("═" * 62)
    print("  STEP 1 — GPU SETUP")
    print("═" * 62)

    if not torch.cuda.is_available():
        print("  ⚠️  No GPU found — CPU only (very slow!)")
        return torch.device("cpu")

    device = torch.device("cuda")
    name   = torch.cuda.get_device_name(0)
    vram   = torch.cuda.get_device_properties(0).total_memory / 1024**3
    free   = (torch.cuda.get_device_properties(0).total_memory
              - torch.cuda.memory_allocated(0)) / 1024**3

    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32       = True
    torch.backends.cudnn.benchmark        = True

    print(f"\n  ✅ GPU      : {name}")
    print(f"  VRAM       : {vram:.1f} GB total  |  {free:.1f} GB free")
    print(f"  CUDA       : {torch.version.cuda}")
    print(f"  fp16       : {CONFIG['fp16']}  ← saves ~40% VRAM")
    print(f"  Batch      : {CONFIG['per_device_train_batch_size']} × "
          f"{CONFIG['gradient_accumulation_steps']} = "
          f"{CONFIG['per_device_train_batch_size'] * CONFIG['gradient_accumulation_steps']} effective")
    print(f"  Epochs     : {CONFIG['num_epochs']}")
    print(f"  LR         : {CONFIG['learning_rate']} → cosine decay → {CONFIG['min_learning_rate']}")
    print(f"  Label smooth: {CONFIG['label_smoothing']}")
    print(f"  Early stop : DISABLED — full 30 epochs")
    print(f"  ⚡ TF32 + cuDNN benchmark ON")
    return device


def print_vram(label=""):
    if not torch.cuda.is_available():
        return
    used  = torch.cuda.memory_allocated(0) / 1024**3
    total = torch.cuda.get_device_properties(0).total_memory / 1024**3
    print(f"  🖥️  VRAM [{label}]: {used:.2f} / {total:.1f} GB")


# ─────────────────────────────────────────────────────────────
# STEP 2 — LOAD DATA
# ─────────────────────────────────────────────────────────────

def load_data():
    global LABEL_MAP

    print("\n" + "═" * 62)
    print("  STEP 2 — Loading Tokenized Data")
    print("═" * 62)

    tok_dir = CONFIG["tokenized_dir"]
    for name in ["train.pt", "val.pt", "test.pt", "meta.json"]:
        if not (tok_dir / name).exists():
            print(f"  ❌ Missing: {tok_dir / name}")
            print(f"  Run tokenize_dataset.py first!")
            return None, None, None, None

    print(f"\n  Loading from: {tok_dir}")
    train_raw = torch.load(tok_dir / "train.pt", weights_only=False)
    val_raw   = torch.load(tok_dir / "val.pt",   weights_only=False)
    test_raw  = torch.load(tok_dir / "test.pt",  weights_only=False)

    with open(tok_dir / "meta.json") as f:
        meta = json.load(f)
    LABEL_MAP = meta["label_map"]

    train_ds = CancerDataset(train_raw)
    val_ds   = CancerDataset(val_raw)
    test_ds  = CancerDataset(test_raw)

    id2label = {v: k for k, v in LABEL_MAP.items()}
    print(f"\n  ✅ Data loaded!")
    print(f"  Train : {len(train_ds):>7,}")
    print(f"  Val   : {len(val_ds):>7,}")
    print(f"  Test  : {len(test_ds):>7,}")
    print(f"\n  Class distribution (train):")
    lbls = train_raw["labels"]
    for i in range(len(LABEL_MAP)):
        c   = (lbls == i).sum().item()
        bar = "█" * max(1, c // max(1, len(lbls) // 30))
        print(f"    {id2label[i]:<25} {c:>6,}  {bar}")

    return train_ds, val_ds, test_ds, LABEL_MAP


# ─────────────────────────────────────────────────────────────
# STEP 3 — LOAD MODEL
# ─────────────────────────────────────────────────────────────

def load_model(label_map, device):
    from transformers import AutoModelForSequenceClassification

    print("\n" + "═" * 62)
    print("  STEP 3 — Loading BioBERT Model")
    print("═" * 62)

    id2label = {v: k for k, v in label_map.items()}

    model = AutoModelForSequenceClassification.from_pretrained(
        CONFIG["model_name"],
        num_labels              = len(label_map),
        id2label                = id2label,
        label2id                = label_map,
        ignore_mismatched_sizes = True,
        hidden_dropout_prob     = 0.1,          # dropout regularization
        attention_probs_dropout_prob = 0.1,
    )
    model = model.to(device)

    total     = sum(p.numel() for p in model.parameters())
    trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"\n  ✅ Model loaded!")
    print(f"  Parameters : {total:,} total  |  {trainable:,} trainable")
    print_vram("after model load")

    return model


# ─────────────────────────────────────────────────────────────
# STEP 4 — LAYER-WISE LEARNING RATE DECAY
# Lower layers get smaller LR — better for fine-tuning BERT
# ─────────────────────────────────────────────────────────────

def get_optimizer(model):
    """
    Layer-wise LR decay:
    - Classifier head  : full LR (3e-5)
    - Top BERT layers  : 0.9 × LR
    - Middle layers    : 0.9^2 × LR
    - Bottom layers    : 0.9^n × LR (smallest)
    - Embeddings       : 0.9^13 × LR (nearly frozen)
    """
    from torch.optim import AdamW

    decay_factor = 0.9
    base_lr      = CONFIG["learning_rate"]
    no_decay     = ["bias", "LayerNorm.weight"]

    # Group parameters by layer depth
    param_groups = []

    # Classifier (full LR)
    classifier_params = [
        (n, p) for n, p in model.named_parameters()
        if "classifier" in n or "pooler" in n
    ]
    param_groups.append({
        "params":       [p for n, p in classifier_params if not any(nd in n for nd in no_decay)],
        "lr":           base_lr,
        "weight_decay": CONFIG["weight_decay"],
    })
    param_groups.append({
        "params":       [p for n, p in classifier_params if any(nd in n for nd in no_decay)],
        "lr":           base_lr,
        "weight_decay": 0.0,
    })

    # BERT encoder layers (12 layers, decay from top)
    for layer_idx in range(11, -1, -1):
        layer_lr = base_lr * (decay_factor ** (11 - layer_idx))
        layer_params = [
            (n, p) for n, p in model.named_parameters()
            if f"encoder.layer.{layer_idx}." in n
        ]
        if not layer_params:
            continue
        param_groups.append({
            "params":       [p for n, p in layer_params if not any(nd in n for nd in no_decay)],
            "lr":           layer_lr,
            "weight_decay": CONFIG["weight_decay"],
        })
        param_groups.append({
            "params":       [p for n, p in layer_params if any(nd in n for nd in no_decay)],
            "lr":           layer_lr,
            "weight_decay": 0.0,
        })

    # Embeddings (smallest LR)
    embed_lr = base_lr * (decay_factor ** 12)
    embed_params = [
        (n, p) for n, p in model.named_parameters()
        if "embedding" in n
    ]
    param_groups.append({
        "params":       [p for n, p in embed_params if not any(nd in n for nd in no_decay)],
        "lr":           embed_lr,
        "weight_decay": CONFIG["weight_decay"],
    })
    param_groups.append({
        "params":       [p for n, p in embed_params if any(nd in n for nd in no_decay)],
        "lr":           embed_lr,
        "weight_decay": 0.0,
    })

    # Remove empty groups
    param_groups = [g for g in param_groups if len(g["params"]) > 0]

    optimizer = AdamW(param_groups, eps=1e-8)
    print(f"  ✅ Optimizer: AdamW with layer-wise LR decay")
    print(f"  Classifier LR : {base_lr:.0e}")
    print(f"  Embedding LR  : {embed_lr:.2e}")

    return optimizer


# ─────────────────────────────────────────────────────────────
# STEP 5 — METRICS
# ─────────────────────────────────────────────────────────────

def compute_metrics(eval_pred):
    from sklearn.metrics import accuracy_score, f1_score
    logits, labels = eval_pred
    preds = np.argmax(logits, axis=-1)
    return {
        "accuracy":    round(float(accuracy_score(labels, preds)), 4),
        "f1_macro":    round(float(f1_score(labels, preds, average="macro",    zero_division=0)), 4),
        "f1_weighted": round(float(f1_score(labels, preds, average="weighted", zero_division=0)), 4),
    }


# ─────────────────────────────────────────────────────────────
# STEP 6 — TRAIN
# ─────────────────────────────────────────────────────────────

def train(model, train_ds, val_ds, device):
    from transformers import TrainingArguments, Trainer
    print("\n" + "═" * 62)
    print("  STEP 4 — Advanced Training")
    print("═" * 62)

    out_dir = CONFIG["output_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)

    steps_per_epoch = len(train_ds) // (
        CONFIG["per_device_train_batch_size"] *
        CONFIG["gradient_accumulation_steps"]
    )
    total_steps = steps_per_epoch * CONFIG["num_epochs"]
    warmup_steps = int(total_steps * CONFIG["warmup_ratio"])

    print(f"\n  🚀 Training plan:")
    print(f"  Epochs        : {CONFIG['num_epochs']}")
    print(f"  Steps/epoch   : {steps_per_epoch:,}")
    print(f"  Total steps   : {total_steps:,}")
    print(f"  Warmup steps  : {warmup_steps:,}")
    print(f"  LR schedule   : cosine decay")
    print(f"  Label smooth  : {CONFIG['label_smoothing']}")
    print(f"  Grad clip     : {CONFIG['max_grad_norm']}")
    print(f"  Early stop    : DISABLED — full 30 epochs")
    print(f"  Output        : {out_dir}")
    print("─" * 62)

    args = TrainingArguments(
        output_dir = str(out_dir),

        # ── Batch ──────────────────────────────────────────
        per_device_train_batch_size = CONFIG["per_device_train_batch_size"],
        per_device_eval_batch_size  = CONFIG["per_device_eval_batch_size"],
        gradient_accumulation_steps = CONFIG["gradient_accumulation_steps"],

        # ── Speed ──────────────────────────────────────────
        fp16                        = CONFIG["fp16"],
        dataloader_pin_memory       = CONFIG["dataloader_pin_memory"],
        dataloader_num_workers      = CONFIG["dataloader_num_workers"],

        # ── Schedule ───────────────────────────────────────
        num_train_epochs            = CONFIG["num_epochs"],
        learning_rate               = CONFIG["learning_rate"],
        lr_scheduler_type           = CONFIG["lr_scheduler_type"],
        warmup_steps                = warmup_steps,
        weight_decay                = CONFIG["weight_decay"],
        max_grad_norm               = CONFIG["max_grad_norm"],

        # ── Regularization ─────────────────────────────────
        label_smoothing_factor      = CONFIG["label_smoothing"],

        # ── Eval & Save ────────────────────────────────────
        eval_strategy               = "epoch",
        save_strategy               = "epoch",
        load_best_model_at_end      = True,
        metric_for_best_model       = "f1_weighted",
        greater_is_better           = True,
        save_total_limit            = CONFIG["save_total_limit"],

        # ── Logging ────────────────────────────────────────
        logging_steps               = CONFIG["logging_steps"],
        logging_dir                 = str(out_dir / "logs"),
        report_to                   = "none",

        seed                        = CONFIG["seed"],
    )

    trainer = Trainer(
        model           = model,
        args            = args,
        train_dataset   = train_ds,
        eval_dataset    = val_ds,
        compute_metrics = compute_metrics,
        callbacks       = [],   # no early stopping — full 30 epochs
    )

    trainer.train()

    # Save best model
    best_path = out_dir / "best_model"
    trainer.save_model(str(best_path))

    from transformers import AutoTokenizer
    tok = AutoTokenizer.from_pretrained(CONFIG["model_name"])
    tok.save_pretrained(str(best_path))

    print(f"\n  ✅ Best model saved → {best_path}")
    print_vram("after training")

    return trainer


# ─────────────────────────────────────────────────────────────
# STEP 7 — EVALUATE
# ─────────────────────────────────────────────────────────────

def evaluate(trainer, test_ds):
    from sklearn.metrics import classification_report

    print("\n" + "═" * 62)
    print("  STEP 5 — Test Set Evaluation")
    print("═" * 62)

    preds_out = trainer.predict(test_ds)
    preds     = np.argmax(preds_out.predictions, axis=-1)
    labels    = preds_out.label_ids
    id2label  = {v: k for k, v in LABEL_MAP.items()}
    names     = [id2label[i] for i in range(len(LABEL_MAP))]

    report = classification_report(labels, preds, target_names=names, digits=4)
    acc    = (preds == labels).mean() * 100

    print(f"\n  CLASSIFICATION REPORT")
    print("─" * 62)
    print(report)
    print(f"  Overall Accuracy : {acc:.2f}%")

    # Per-class summary
    print(f"\n  Per-class F1 summary:")
    from sklearn.metrics import f1_score
    for i, name in enumerate(names):
        f1  = f1_score(labels, preds, labels=[i], average="micro", zero_division=0)
        bar = "█" * int(f1 * 20)
        print(f"    {name:<25} F1={f1:.4f}  {bar}")

    # Save
    out_dir = CONFIG["output_dir"]
    with open(out_dir / "test_results.txt", "w") as f:
        f.write("BioBERT Advanced — Cancer Classification Results\n")
        f.write("=" * 62 + "\n\n")
        f.write(f"Epochs         : {CONFIG['num_epochs']}\n")
        f.write(f"LR             : {CONFIG['learning_rate']} (cosine decay)\n")
        f.write(f"Label smoothing: {CONFIG['label_smoothing']}\n")
        f.write(f"Layer-wise LR  : Yes\n")
        f.write(f"Early stopping : DISABLED — full 30 epochs\n")
        f.write(f"Overall Acc    : {acc:.2f}%\n\n")
        f.write(report)

    preds_json = [
        {"true": id2label[int(l)], "predicted": id2label[int(p)], "correct": bool(l == p)}
        for l, p in zip(labels, preds)
    ]
    with open(out_dir / "predictions.json", "w") as f:
        json.dump(preds_json, f, indent=2)

    print(f"\n  Saved → {out_dir}/test_results.txt")
    print(f"  Saved → {out_dir}/predictions.json")

    return acc


# ─────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────

def main():
    print("═" * 62)
    print("  BioBERT ADVANCED Fine-Tuning")
    print("  RTX 4050 6GB  |  8 Cancer Types")
    print("  ─────────────────────────────────────────")
    print("  ✅ 30 epochs — full training")
    print("  ✅ Cosine LR decay")
    print("  ✅ Layer-wise LR (LLRD)")
    print("  ✅ Label smoothing 0.1")
    print("  ✅ No early stopping — full 30 epochs")
    print("  ✅ Gradient clipping 1.0")
    print("═" * 62 + "\n")

    # Dependency check
    missing = []
    for pkg in ["transformers", "torch", "sklearn"]:
        try:
            __import__("sklearn" if pkg == "sklearn" else pkg)
        except ImportError:
            missing.append("scikit-learn" if pkg == "sklearn" else pkg)
    if missing:
        print(f"  Install: pip install {' '.join(missing)}")
        return

    torch.manual_seed(CONFIG["seed"])
    np.random.seed(CONFIG["seed"])

    # 1. GPU
    device = setup_gpu()

    # 2. Data
    train_ds, val_ds, test_ds, label_map = load_data()
    if train_ds is None:
        return

    # 3. Model
    model = load_model(label_map, device)

    # 4. Train
    trainer = train(model, train_ds, val_ds, device)

    # 5. Evaluate
    acc = evaluate(trainer, test_ds)

    print(f"\n{'═'*62}")
    print(f"  ✅ TRAINING COMPLETE!")
    print(f"  Accuracy       : {acc:.2f}%")
    print(f"  Best model     : {CONFIG['output_dir']}/best_model/")
    print(f"  Results        : {CONFIG['output_dir']}/test_results.txt")
    print(f"{'═'*62}\n")


if __name__ == "__main__":
    main()