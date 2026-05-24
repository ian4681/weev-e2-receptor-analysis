"""
分析 AlphaFold3 预测的 WEEV E2 - 受体复合物结果（两轮共6个job）
输入：data/af3_results/<job_name>/ 目录下的 JSON 文件
"""
import json, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

DATA_DIR = "F:/CC/projects/cancer_perturb/文书/套磁信/章新政/02_定制层/data/af3_results"
FIG_DIR  = "F:/CC/projects/cancer_perturb/文书/套磁信/章新政/02_定制层/figures"
os.makedirs(FIG_DIR, exist_ok=True)

JOBS = [
    "WEEV_WT_VLDLR",
    "WEEV_V265F_VLDLR",
    "WEEV_WT_PCDH10",
    "WEEV_V265F_PCDH10",
    "WEEV_WT_VLDLR_LA3",
    "WEEV_V265F_VLDLR_LA3",
]

WEEV_LEN = 416
V265_IDX = 264  # 0-indexed

# ── 1. 读取所有job，选最优model ──────────────────────────────────────────────
def load_best_model(job_dir):
    files = os.listdir(job_dir)
    summary_files = [f for f in files if "summary_confidences" in f and f.endswith(".json")]
    best_summary, best_full, best_idx, best_score = None, None, 0, -1
    for sf in summary_files:
        try:
            txt = open(f"{job_dir}/{sf}", encoding="utf-8").read().strip()
            if not txt:
                continue
            s = json.loads(txt)
        except Exception:
            continue
        score = s.get("ranking_score") or 0
        if score > best_score:
            best_score = score
            best_summary = s
            m = re.search(r"_(\d+)\.json$", sf)
            best_idx = int(m.group(1)) if m else 0
            full_name = sf.replace("summary_confidences", "full_data")
            try:
                txt2 = open(f"{job_dir}/{full_name}", encoding="utf-8").read().strip()
                best_full = json.loads(txt2) if txt2 else None
            except Exception:
                best_full = None
    return best_summary, best_full, best_idx

print("=" * 72)
print(f"{'Job':<28} {'mdl':>3} {'iptm':>6} {'ptm':>6} {'pae_AB':>8} {'iptm_AB':>8}")
print("=" * 72)

results = {}
for job in JOBS:
    job_dir = f"{DATA_DIR}/{job}"
    if not os.path.isdir(job_dir):
        print(f"  {job:<28} -- 目录不存在")
        continue
    summary, full, idx = load_best_model(job_dir)
    if summary is None:
        print(f"  {job:<28} -- 无有效summary")
        continue

    iptm    = summary.get("iptm")
    ptm     = summary.get("ptm")
    pae_min = summary.get("chain_pair_pae_min")
    pae_AB  = pae_min[0][1] if pae_min else None
    cp_iptm = summary.get("chain_pair_iptm")
    iptm_AB = cp_iptm[0][1] if cp_iptm else None

    results[job] = dict(iptm=iptm, ptm=ptm, pae_AB=pae_AB, iptm_AB=iptm_AB,
                        full=full, model_idx=idx)

    def f(v): return f"{v:.3f}" if v is not None else " N/A"
    print(f"  {job:<28} {idx:>3} {f(iptm):>6} {f(ptm):>6} {f(pae_AB):>8} {f(iptm_AB):>8}")

print("=" * 72)

# ── 2. PAE heatmap（6张）────────────────────────────────────────────────────
for job in JOBS:
    if job not in results or results[job]["full"] is None:
        print(f"  PAE skip: {job} (no full data)")
        continue
    pae = np.array(results[job]["full"].get("pae", []))
    if pae.size == 0:
        print(f"  PAE skip: {job} (empty full_pae)")
        continue

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(pae, cmap="bwr_r", vmin=0, vmax=30, aspect="auto")
    plt.colorbar(im, ax=ax, label="PAE (Angstrom)")

    ax.axhline(WEEV_LEN - 0.5, color="lime", lw=1.5, linestyle="--")
    ax.axvline(WEEV_LEN - 0.5, color="lime", lw=1.5, linestyle="--")
    ax.axhline(V265_IDX - 0.5, color="yellow", lw=1, linestyle=":")
    ax.axvline(V265_IDX - 0.5, color="yellow", lw=1, linestyle=":")

    iptm_v = results[job]["iptm"] or 0
    pae_v  = results[job]["pae_AB"] or 0
    ax.set_title(f"{job}\nipTM={iptm_v:.3f}  pae_AB={pae_v:.2f} Ang",
                 fontsize=9, fontweight="bold")
    ax.set_xlabel("Token index")
    ax.set_ylabel("Token index")

    n = pae.shape[1]
    ax.text(WEEV_LEN / 2, -12, "WEEV E2",
            ha="center", fontsize=8, color="lime", clip_on=False)
    ax.text(WEEV_LEN + (n - WEEV_LEN) / 2, -12, "Receptor",
            ha="center", fontsize=8, color="lime", clip_on=False)

    fig.tight_layout()
    out = f"{FIG_DIR}/pae_{job}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: pae_{job}.png")

# ── 3. V265 接触概率（3对比较）──────────────────────────────────────────────
for receptor, wt_job, mut_job in [
    ("VLDLR",     "WEEV_WT_VLDLR",     "WEEV_V265F_VLDLR"),
    ("PCDH10",    "WEEV_WT_PCDH10",    "WEEV_V265F_PCDH10"),
    ("VLDLR_LA3", "WEEV_WT_VLDLR_LA3", "WEEV_V265F_VLDLR_LA3"),
]:
    if wt_job not in results or mut_job not in results:
        continue
    if results[wt_job]["full"] is None or results[mut_job]["full"] is None:
        print(f"  Contact skip: {receptor} (no full data)")
        continue

    cp_wt  = np.array(results[wt_job]["full"].get("contact_probs", []))
    cp_mut = np.array(results[mut_job]["full"].get("contact_probs", []))
    if cp_wt.size == 0 or cp_mut.size == 0:
        print(f"  Contact skip: {receptor} (empty contact_probs)")
        continue

    v_wt  = cp_wt[V265_IDX, WEEV_LEN:]
    v_mut = cp_mut[V265_IDX, WEEV_LEN:]

    fig, axes = plt.subplots(2, 1, figsize=(10, 5), sharex=True)
    x = np.arange(len(v_wt))
    axes[0].bar(x, v_wt,  color="#2166ac", alpha=0.8, width=1)
    axes[0].set_title(f"V265 (WT) contact prob with {receptor}", fontsize=9)
    axes[0].set_ylabel("Contact prob"); axes[0].set_ylim(0, 1)

    axes[1].bar(x, v_mut, color="#d73027", alpha=0.8, width=1)
    axes[1].set_title(f"V265F (mut) contact prob with {receptor}", fontsize=9)
    axes[1].set_ylabel("Contact prob"); axes[1].set_ylim(0, 1)
    axes[1].set_xlabel(f"{receptor} residue index (0-based from chain start)")

    fig.suptitle(f"Position 265 contact with {receptor}  [blue=WT, red=V265F]",
                 fontweight="bold")
    fig.tight_layout()
    out = f"{FIG_DIR}/contact_V265_{receptor}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: contact_V265_{receptor}.png")

# ── 4. 综合对比：iptm + pae_AB（6 jobs，分两轮着色）────────────────────────
labels, iptm_vals, pae_vals, bar_colors = [], [], [], []
round_colors = {
    "WEEV_WT_VLDLR":       ("#1a6faf", "R1"),
    "WEEV_V265F_VLDLR":    ("#e05c2c", "R1"),
    "WEEV_WT_PCDH10":      ("#1a6faf", "R1"),
    "WEEV_V265F_PCDH10":   ("#e05c2c", "R1"),
    "WEEV_WT_VLDLR_LA3":   ("#0a3d6b", "R2"),
    "WEEV_V265F_VLDLR_LA3":("#8b1a00", "R2"),
}
for job in JOBS:
    if job not in results:
        continue
    short = job.replace("WEEV_", "").replace("_VLDLR_LA3", "\nLA3").replace("_VLDLR", "\nVLDLR").replace("_PCDH10", "\nPCDH10")
    labels.append(short)
    iptm_vals.append(results[job]["iptm"] or 0)
    pae_vals.append(results[job]["pae_AB"] or 0)
    bar_colors.append(round_colors.get(job, ("#888", ""))[0])

x = np.arange(len(labels))
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

axes[0].bar(x, iptm_vals, color=bar_colors, alpha=0.88, width=0.6, edgecolor="white")
axes[0].axhline(0.7, color="gray", lw=1, linestyle="--", alpha=0.6)
axes[0].text(len(labels) - 0.5, 0.725, "0.7 confidence threshold",
             fontsize=7, color="gray", ha="right")
axes[0].set_xticks(x); axes[0].set_xticklabels(labels, fontsize=7.5)
axes[0].set_ylabel("ipTM"); axes[0].set_ylim(0, 0.8)
axes[0].set_title("Interface ipTM (higher = better)", fontweight="bold")
for i, v in enumerate(iptm_vals):
    axes[0].text(i, v + 0.008, f"{v:.3f}", ha="center", fontsize=7.5)

axes[1].bar(x, pae_vals, color=bar_colors, alpha=0.88, width=0.6, edgecolor="white")
axes[1].set_xticks(x); axes[1].set_xticklabels(labels, fontsize=7.5)
axes[1].set_ylabel("pae_AB (Angstrom)"); axes[1].set_ylim(0, 35)
axes[1].set_title("Cross-chain PAE min (lower = better)", fontweight="bold")
for i, v in enumerate(pae_vals):
    axes[1].text(i, v + 0.4, f"{v:.1f}", ha="center", fontsize=7.5)

legend_elems = [
    Patch(facecolor="#1a6faf", label="WT — Round 1 (full-length receptor)"),
    Patch(facecolor="#e05c2c", label="V265F — Round 1"),
    Patch(facecolor="#0a3d6b", label="WT — Round 2 (LA3 truncation)"),
    Patch(facecolor="#8b1a00", label="V265F — Round 2"),
]
axes[0].legend(handles=legend_elems, fontsize=7, loc="upper right")

fig.suptitle("WEEV E2 Receptor Binding — AF3 Confidence Summary (2 rounds)",
             fontweight="bold", fontsize=11)
fig.tight_layout()
out = f"{FIG_DIR}/af3_summary_all.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: af3_summary_all.png")

print(f"\nDone. Figures in: {FIG_DIR}")
