import pickle, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")                        
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import (RepeatedStratifiedKFold, StratifiedKFold,
                                     cross_validate, cross_val_predict)
from sklearn.metrics import (make_scorer, matthews_corrcoef, roc_curve,
                             precision_recall_curve, roc_auc_score,
                             average_precision_score, brier_score_loss,
                             confusion_matrix)
from sklearn.calibration import calibration_curve

warnings.filterwarnings("ignore")
out = Path("results/ML_comparison"); out.mkdir(parents=True, exist_ok=True)

# Load data
with open("data/dataset.pkl", "rb") as f:
    data = pickle.load(f)

by_seq = {}                                   # sequence -> set of labels seen
for s, l in zip(data["sequences"], data["labels"]):
    by_seq.setdefault(s.upper().strip(), set()).add(l)

# keep a sequence only if it appears with ONE label 
pairs = [(s, next(iter(ls))) for s, ls in by_seq.items()
         if len(ls) == 1 and 100 <= len(s) <= 1000]
print(f"Raw: {len(data['sequences'])}  ->  after dedup/length filter: {len(pairs)}")

# Features
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
FEATURE_NAMES = list(AMINO_ACIDS) + ["Length", "Aromatic", "Hydrophobic", "Charged"]

def extract_features(seq):
    n = len(seq)
    f = [seq.count(aa) / n for aa in AMINO_ACIDS]
    f.append(n)
    f.append(sum(seq.count(a) for a in "FYW") / n)
    f.append(sum(seq.count(a) for a in "AILMFVPWG") / n)
    f.append(sum(seq.count(a) for a in "DEKR") / n)
    return f

X = np.array([extract_features(s) for s, _ in pairs])
y = np.array([l for _, l in pairs])
print(f"X shape {X.shape} | PETase={y.sum()} Non-PETase={(y==0).sum()}")

# Models
try:
    from xgboost import XGBClassifier
    boost_name = "XGBoost"
    boost = XGBClassifier(n_estimators=200, max_depth=3, learning_rate=0.05,
                          subsample=0.8, eval_metric="logloss",
                          random_state=42, verbosity=0)
except ImportError:                            # conda install -c conda-forge xgboost
    boost_name = "GradBoost"
    boost = GradientBoostingClassifier(n_estimators=200, max_depth=3,
                                       learning_rate=0.05, random_state=42)

MODELS = {
    "RandomForest": Pipeline([("clf", RandomForestClassifier(
        n_estimators=100, class_weight="balanced", random_state=42))]),
    "SVM_RBF": Pipeline([("sc", StandardScaler()),
        ("clf", SVC(kernel="rbf", C=10, gamma="scale", class_weight="balanced",
                    probability=True, random_state=42))]),
    boost_name: Pipeline([("clf", boost)]),
    "LogisticReg": Pipeline([("sc", StandardScaler()),
        ("clf", LogisticRegression(max_iter=2000, class_weight="balanced",
                                   random_state=42))]),
    "MLP": Pipeline([("sc", StandardScaler()),
        ("clf", MLPClassifier(hidden_layer_sizes=(64, 32), alpha=1e-3,
                              max_iter=1000, early_stopping=True,
                              random_state=42))]),
}
COLORS = dict(zip(MODELS, ["#2196F3", "#FF5722", "#4CAF50", "#9C27B0", "#FF9800"]))

# Repeated Stratified CV
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10,
                             random_state=42)
scoring = {"accuracy": "accuracy", "balanced_acc": "balanced_accuracy",
           "precision": "precision", "recall": "recall", "f1": "f1",
           "mcc": make_scorer(matthews_corrcoef), "roc_auc": "roc_auc",
           "pr_auc": "average_precision", "neg_brier": "neg_brier_score"}

raw, rows = {}, []
for name, pipe in MODELS.items():
    print(f"CV: {name}")
    res = cross_validate(pipe, X, y, cv=cv, scoring=scoring, n_jobs=-1)
    raw[name] = res
    row = {"model": name}
    for m in scoring:
        v = res[f"test_{m}"]
        if m == "neg_brier": v = -v           # flip so lower = better
        row[m] = f"{v.mean():.3f} ± {v.std():.3f}"
    rows.append(row)

table = pd.DataFrame(rows).set_index("model")
table.rename(columns={"neg_brier": "brier(lower=better)"}, inplace=True)
table.to_csv(out / "cv_metrics.csv")
print("\n", table.to_string())

# -Out of Fold Probabilities
oof = {}
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
for name, pipe in MODELS.items():
    oof[name] = cross_val_predict(pipe, X, y, cv=skf, method="predict_proba")[:, 1]

# Plots
#  ROC
plt.figure(figsize=(6.5, 6))
for n, p in oof.items():
    fpr, tpr, _ = roc_curve(y, p)
    plt.plot(fpr, tpr, lw=2, color=COLORS[n], label=f"{n} (AUC={roc_auc_score(y,p):.3f})")
plt.plot([0, 1], [0, 1], "k--", alpha=.5)
plt.xlabel("False positive rate"); plt.ylabel("True positive rate")
plt.title("ROC (out-of-fold predictions)"); plt.legend(loc="lower right"); plt.grid(alpha=.3)
plt.tight_layout(); plt.savefig(out / "roc.png", dpi=200); plt.close()

#  Precision-Recall
plt.figure(figsize=(6.5, 6))
for n, p in oof.items():
    pr, rc, _ = precision_recall_curve(y, p)
    plt.plot(rc, pr, lw=2, color=COLORS[n], label=f"{n} (AP={average_precision_score(y,p):.3f})")
plt.axhline(y.mean(), color="gray", ls="--", label=f"baseline={y.mean():.2f}")
plt.xlabel("Recall"); plt.ylabel("Precision")
plt.title("Precision-Recall (out-of-fold)"); plt.legend(loc="lower left"); plt.grid(alpha=.3)
plt.tight_layout(); plt.savefig(out / "pr.png", dpi=200); plt.close()

#  Metric bars with SD error bars
show = ["accuracy", "balanced_acc", "f1", "mcc", "roc_auc"]
w = 0.16; xs = np.arange(len(show))
plt.figure(figsize=(11, 5))
for i, n in enumerate(MODELS):
    means = [raw[n][f"test_{m}"].mean() for m in show]
    sds   = [raw[n][f"test_{m}"].std() for m in show]
    plt.bar(xs + i * w, means, w, yerr=sds, capsize=2, label=n, color=COLORS[n])
plt.xticks(xs + 2 * w, show); plt.ylim(0, 1.05); plt.ylabel("Score (mean ± SD)")
plt.title(f"{5}-fold x {10} repeated CV"); plt.legend(ncol=5, fontsize=8)
plt.grid(axis="y", alpha=.3); plt.tight_layout()
plt.savefig(out / "metric_bars.png", dpi=200); plt.close()

#  MCC distribution across all folds (shows whether differences are real)
plt.figure(figsize=(7, 5))
plt.boxplot([raw[n]["test_mcc"] for n in MODELS], showmeans=True)
plt.ylabel("MCC per fold"); plt.title("MCC spread across CV folds")
plt.xticks(range(1, len(MODELS) + 1), list(MODELS), rotation=20)
plt.savefig(out / "mcc_boxplot.png", dpi=200); plt.close()

#  Confusion matrices (out-of-fold, threshold 0.5)
fig, axes = plt.subplots(1, len(MODELS), figsize=(4 * len(MODELS), 3.8))
for ax, (n, p) in zip(axes, oof.items()):
    cm = confusion_matrix(y, (p >= 0.5).astype(int))
    ax.imshow(cm, cmap="Blues")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, cm[i, j], ha="center", va="center", fontsize=13)
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
    ax.set_xticklabels(["Non", "PET"]); ax.set_yticklabels(["Non", "PET"])
    ax.set_xlabel("Predicted"); ax.set_ylabel("True"); ax.set_title(n)
plt.tight_layout(); plt.savefig(out / "confusion.png", dpi=200); plt.close()

#  Calibration (are probabilities trustworthy?)
plt.figure(figsize=(6, 6))
for n, p in oof.items():
    fr, mp = calibration_curve(y, p, n_bins=8, strategy="quantile")
    plt.plot(mp, fr, "o-", color=COLORS[n], label=f"{n} (Brier={brier_score_loss(y,p):.3f})")
plt.plot([0, 1], [0, 1], "k--", alpha=.5)
plt.xlabel("Mean predicted probability"); plt.ylabel("Observed PETase fraction")
plt.title("Calibration"); plt.legend(fontsize=8); plt.grid(alpha=.3)
plt.tight_layout(); plt.savefig(out / "calibration.png", dpi=200); plt.close()

#  RF feature importance (fit on all data, for interpretation only)
rf = MODELS["RandomForest"].fit(X, y).named_steps["clf"]
idx = np.argsort(rf.feature_importances_)[::-1][:15][::-1]
plt.figure(figsize=(6, 6))
plt.barh([FEATURE_NAMES[i] for i in idx], rf.feature_importances_[idx], color=COLORS["RandomForest"])
plt.xlabel("Importance (MDI)"); plt.title("RF top 15 features")
plt.tight_layout(); plt.savefig(out / "rf_importance.png", dpi=200); plt.close()

print(f"\nDone. Everything is in {out}/")