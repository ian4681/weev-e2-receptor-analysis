# WEEV E2 受体界面计算探索

**西方马脑炎病毒 E2 V265F 受体向性切换的结构基础 — AlphaFold3 + ESM-2 分析**

> 独立计算探索项目 · 2026年5月  
> 参考文献：Byrd et al., *Cell Reports* 2025 | Chen et al., *Cell* 2023

---

## 项目背景

*Cell Reports* 2025 报道，WEEV E2 蛋白第265位的单点突变（V265F）足以将病毒受体偏好从 PCDH10 切换至 VLDLR。本项目尝试从纯计算角度回答：**这一界面差异能否被 AlphaFold3 结构预测捕捉到？ESM-2 对该位点的零样本打分与实验结果是否一致？**

---

## 📁 项目结构

```
.
├── index.html          ← 分析报告（即 GitHub Pages 展示页）
├── scripts/
├── data/
└── figures/
```

---

## 📊 主要结果

### ESM-2 零样本突变打分

| 位点 | WT | 关键突变 | 打分 | 解读 |
|---|---|---|---|---|
| WEEV E2 pos265 | V | **V265F** | **−0.77** | 功能改变而非丧失，与受体切换实验一致 |
| VLDLR LA3 P129 | P | P129H | −2.83 | 高度功能性保守，替换破坏结合 |
| VLDLR LA3 W132 | W | W132A | −7.37 | 必需接触残基，极度不耐突变 |

### AlphaFold3 结构预测（两轮迭代）

| Job | 轮次 | ipTM | pae_AB (Å) | 备注 |
|---|---|---|---|---|
| WEEV_WT_VLDLR | R1 全长 | 0.210 | 16.35 | 信号被8个LA重复稀释 |
| WEEV_V265F_VLDLR | R1 全长 | 0.160 | 20.18 | |
| WEEV_WT_PCDH10 | R1 全长 | 0.150 | 27.04 | 无模板，信号极低 |
| WEEV_V265F_PCDH10 | R1 全长 | 0.140 | 27.07 | |
| **WEEV_WT_VLDLR_LA3** | **R2 截短** | **0.290** | **13.33** | **最优** |
| WEEV_V265F_VLDLR_LA3 | R2 截短 | 0.180 | 17.29 | |

**R2 设计动机**：全长 VLDLR 含8个 LA 重复，AF3 无法定位真实结合区（LA3）。截短至 LA2-LA4（125aa）后界面置信度显著提升。

### 关键发现

1. **V265 不在直接接触界面**：265位残基的接触概率在所有条件下均接近零，提示 V265F 通过影响 E2 整体构象（变构机制）而非直接接触受体来切换受体偏好
2. **AF3 存在模板偏差**：PDB 9L1N 是 WEEV E2 WT 晶体结构，AF3 对 WT 预测有完整模板支撑，导致 WT 系统性高于 V265F，不能直接解读为结合亲和力差异
3. **整体置信度偏低**：所有 ipTM < 0.5，属于低置信预测；PDB 中无 E2-VLDLR 复合物实验结构是根本限制

---

## 🔁 迭代过程亮点

本项目经历了几个有代表性的设计迭代，记录在 `scripts/README.md` 中：

- **AF3 JSON 上传失败** → 改用手动网页 UI（了解工具边界）
- **遗漏受体侧打分** → `02b_score_vldlr.py` 补丁脚本（发现问题即补救）
- **全长 VLDLR 信号稀释** → R2 截短设计（分析 R1 结果→诊断问题→重新实验）
- **`full_pae` 键名错误** → inspect 实际 JSON 键名修正（调试流程）

---

## ⚙️ 运行环境

```
Python 3.x
PyTorch 2.6.0 + CUDA 12.4
fair-esm 2.0.0（ESM-2 650M，本地权重约 1.3 GB，不含于此仓库）
matplotlib / numpy / Pillow / biopython / requests
```

ESM-2 模型权重需自行下载或通过 `torch.hub` 获取：
```python
import torch, esm
torch.hub.set_dir('/path/to/cache')
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
```

AF3 预测通过 [alphafoldserver.com](https://alphafoldserver.com) 网页界面手动提交，不需要本地 GPU。

---

## 完整报告

完整分析报告（含图表与代码解读）：[在线阅读](https://ian4681.github.io/weev-e2-receptor-analysis/) | 或本地打开 [`index.html`](./index.html)

---

## 📚 参考文献

1. Byrd et al. *Cell Reports* 2025 — WEEV dual-receptor recognition (V265F)
2. Chen et al. *Cell* 2023 — SFV-VLDLR complex structure
3. Lin et al. *Science* 2023 — ESM-2 protein language model
4. Meier et al. *NeurIPS* 2021 — Masked marginal scoring for mutation effect prediction
