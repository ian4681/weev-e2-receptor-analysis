"""生成最终双面板热图，用于套磁信附件"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm

DATA_DIR = "F:/CC/projects/cancer_perturb/文书/套磁信/章新政/02_定制层/data"
FIG_DIR  = "F:/CC/projects/cancer_perturb/文书/套磁信/章新政/02_定制层/figures"
os.makedirs(FIG_DIR, exist_ok=True)

AAS = list("ACDEFGHIKLMNPQRSTVWY")

with open(f"{DATA_DIR}/all_scores.json") as f:
    scores = json.load(f)

def make_matrix(label, positions):
    """返回 (n_pos x 20) score matrix 和 WT 列表"""
    ps = scores[label]
    mat = np.zeros((len(positions), len(AAS)))
    wts = []
    for i, p in enumerate(positions):
        entry = ps.get(str(p), ps.get(p, {}))
        wt = entry.get("WT", "?")
        wts.append(wt)
        for j, aa in enumerate(AAS):
            mat[i, j] = entry.get(aa, 0.0)
        if wt in AAS:
            mat[i, AAS.index(wt)] = 0.0  # WT 定义为 0
    return mat, wts

fig, axes = plt.subplots(1, 2, figsize=(15, 4.5),
                         gridspec_kw={"width_ratios": [1, 2.5], "wspace": 0.4})
fig.patch.set_facecolor("white")

# ─── 左图：WEEV E2 V265 ────────────────────────────────────────────────
mat_w, wts_w = make_matrix("WEEV_E2_71V", [265])
vmax_w = max(abs(mat_w).max(), 0.5)
norm_w = TwoSlopeNorm(vmin=-vmax_w, vcenter=0, vmax=vmax_w)

ax = axes[0]
im = ax.imshow(mat_w, aspect="auto", cmap="RdBu_r", norm=norm_w)
ax.set_xticks(range(len(AAS))); ax.set_xticklabels(AAS, fontsize=8.5)
ax.set_yticks([0]); ax.set_yticklabels([f"V265\n(WT=V)"], fontsize=9, fontweight="bold")
ax.set_title("WEEV E2 (71V strain)\nV265 mutational landscape", fontsize=10, fontweight="bold", pad=8)
ax.set_xlabel("Substitution amino acid", fontsize=9)

# 标注 V265F（关键受体偏好切换）
j_F = AAS.index("F")
score_F = mat_w[0, j_F]
ax.add_patch(mpatches.FancyBboxPatch((j_F-0.48, -0.48), 0.96, 0.96,
    boxstyle="round,pad=0.05", linewidth=2.2, edgecolor="black", facecolor="none"))
ax.text(j_F, 0, f"{score_F:.2f}", ha="center", va="center",
        fontsize=7.5, fontweight="bold", color="black")

cbar = plt.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
cbar.set_label("log p(mut|context) − log p(WT|context)", fontsize=7.5)

ax.text(0.5, -0.18,
        "■ V265F: score=−0.77 (switches from PCDH10→VLDLR binding, Cell Reports 2025)",
        transform=ax.transAxes, ha="center", fontsize=7.5, color="#333333")

# ─── 右图：VLDLR LA3 结合界面 5 位点 ─────────────────────────────────
positions_v = [129, 132, 135, 137, 139]
mat_v, wts_v = make_matrix("VLDLR_human", positions_v)
ylabels = [f"{wt}{p}" for wt, p in zip(wts_v, positions_v)]

vmax_v = max(abs(mat_v).max(), 1)
norm_v = TwoSlopeNorm(vmin=-vmax_v, vcenter=0, vmax=vmax_v)

ax2 = axes[1]
im2 = ax2.imshow(mat_v, aspect="auto", cmap="RdBu_r", norm=norm_v)
ax2.set_xticks(range(len(AAS))); ax2.set_xticklabels(AAS, fontsize=8.5)
ax2.set_yticks(range(len(positions_v))); ax2.set_yticklabels(ylabels, fontsize=9, fontweight="bold")
ax2.set_title("Human VLDLR\nSFV-binding interface residues (LA domain)", fontsize=10, fontweight="bold", pad=8)
ax2.set_xlabel("Substitution amino acid", fontsize=9)

# 标注 P129H
i_P129 = positions_v.index(129)
j_H    = AAS.index("H")
score_P129H = mat_v[i_P129, j_H]
ax2.add_patch(mpatches.FancyBboxPatch((j_H-0.48, i_P129-0.48), 0.96, 0.96,
    boxstyle="round,pad=0.05", linewidth=2.2, edgecolor="black", facecolor="none"))
ax2.text(j_H, i_P129, f"{score_P129H:.2f}", ha="center", va="center",
         fontsize=7, fontweight="bold", color="black")

# 标注 W132A
i_W132 = positions_v.index(132)
j_A    = AAS.index("A")
score_W132A = mat_v[i_W132, j_A]
ax2.add_patch(mpatches.FancyBboxPatch((j_A-0.48, i_W132-0.48), 0.96, 0.96,
    boxstyle="round,pad=0.05", linewidth=2.0, edgecolor="dimgray", facecolor="none",
    linestyle="--"))
ax2.text(j_A, i_W132, f"{score_W132A:.1f}", ha="center", va="center",
         fontsize=6.5, fontweight="bold", color="dimgray")

cbar2 = plt.colorbar(im2, ax=ax2, shrink=0.75, pad=0.02)
cbar2.set_label("log p(mut|context) − log p(WT|context)", fontsize=7.5)

ax2.text(0.5, -0.18,
         "■ P129H: score=−2.83 (disrupts VLDLR-SFV binding, Cell 2023)   "
         "□ W132A: score=−7.4 (essential contact residue)",
         transform=ax2.transAxes, ha="center", fontsize=7.5, color="#333333")

fig.suptitle(
    "ESM-2 650M Zero-Shot Mutation Effect Scores — Alphavirus Receptor Interface",
    fontsize=11.5, fontweight="bold", y=1.01
)

for fmt in ["png", "pdf"]:
    out = f"{FIG_DIR}/esm2_mutation_heatmap.{fmt}"
    fig.savefig(out, bbox_inches="tight", dpi=200)
    print(f"Saved: {out}")

plt.close()
print("Done.")
