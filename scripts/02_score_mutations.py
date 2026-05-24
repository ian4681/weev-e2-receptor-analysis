"""
02_score_mutations.py
用 ESM-2 masked marginal log-ratio 对关键位点的全部 20 种氨基酸替换打分。

打分公式（Meier et al. 2021 / Notin et al. 2022）：
  score(WT→MUT, pos) = log p(MUT | x_masked_at_pos) - log p(WT | x_masked_at_pos)
  负分 = 相对 WT 的序列 plausibility 下降（倾向于功能破坏）
  正分 = 序列 plausibility 上升

目标：
  - WEEV E2（WEEV_71V_E2）V265 位点：期待 V265F 显示大负分（原本不能结合 VLDLR）
    但在 WEEV_McM_E2 中 V265 本来就是 F，所以分析 71V 株
  - SFV E1 W132/D135/E137/D139 位点：期待这些保守残基替换显示负分

输出：
  data/scores_WEEV71V_E2_V265.json
  data/scores_SFV_E1_interface.json
  data/all_scores.json

运行环境：F:/CC/conda_envs/esm2
"""

import os, json, torch, esm
import numpy as np

SCRIPT_DIR = os.path.dirname(__file__)
DATA_DIR   = os.path.join(SCRIPT_DIR, "..", "data")
FASTA_PATH = os.path.join(DATA_DIR, "alphavirus_sequences.fasta")
MODEL_CACHE = "F:/CC/.cache/torch/hub"

AAS = list("ACDEFGHIKLMNPQRSTVWY")

# ── 关键位点定义（UniProt 全长编号，1-indexed）──────────────────────────
# 需在拿到序列后手动确认编号；这里先用文献中报告的成熟蛋白编号
# 如果 UniProt 序列包含信号肽，实际位置 = 成熟蛋白位置 + signal_peptide_length
SITES = {
    "WEEV_71V_E2": {
        "positions_1idx": [265],           # V265，切换 VLDLR 结合能力的关键位点
        "wt_context": "V265 (WT=V，McM株此位为F)",
    },
    "SFV_E1": {
        "positions_1idx": [132, 135, 137, 139],  # W132,D135,E137,D139
        "wt_context": "E1-DIII VLDLR 结合界面",
    },
}

def parse_fasta(path):
    seqs = {}
    cur, seq_lines = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur: seqs[cur] = "".join(seq_lines)
                cur = line[1:].split("|")[0]
                seq_lines = []
            else:
                seq_lines.append(line)
    if cur: seqs[cur] = "".join(seq_lines)
    return seqs

def load_model():
    torch.hub.set_dir(MODEL_CACHE)
    print("Loading ESM-2 650M from local cache...", flush=True)
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device).eval()
    print(f"  Model on: {device}")
    return model, alphabet, device

def masked_marginal_scores(model, alphabet, device, sequence, positions_1idx):
    """
    对 sequence 中的每个目标位点（1-indexed），
    返回 dict: {pos: {aa: score}} 其中 score = log p(aa|masked) - log p(WT|masked)
    """
    batch_converter = alphabet.get_batch_converter()
    all_scores = {}

    for pos1 in positions_1idx:
        pos0 = pos1 - 1                     # 转为 0-indexed
        if pos0 >= len(sequence):
            print(f"  WARNING: position {pos1} out of range (seq len={len(sequence)})")
            continue

        wt_aa = sequence[pos0]

        # 先 tokenize 原始序列，再在 tensor 上直接 mask 该位点
        _, _, tokens = batch_converter([("seq", sequence)])
        tokens = tokens.clone().to(device)
        tokens[0, pos0 + 1] = alphabet.mask_idx   # +1 因为 BOS 占位 0

        with torch.no_grad():
            results = model(tokens, repr_layers=[], return_contacts=False)

        # tokens 中 +1 是因为第 0 位是 <cls>
        logits = results["logits"][0, pos0 + 1]
        log_probs = torch.log_softmax(logits, dim=-1).cpu()

        wt_idx  = alphabet.get_idx(wt_aa)
        wt_lp   = log_probs[wt_idx].item()

        pos_scores = {"WT": wt_aa}
        for aa in AAS:
            idx = alphabet.get_idx(aa)
            pos_scores[aa] = round(log_probs[idx].item() - wt_lp, 4)

        all_scores[pos1] = pos_scores
        print(f"  pos {pos1} (WT={wt_aa}): scored {len(AAS)} mutations")

    return all_scores

def main():
    seqs = parse_fasta(FASTA_PATH)
    print(f"Loaded {len(seqs)} sequences: {list(seqs.keys())}")

    # 快速检查关键序列存在
    for label in SITES:
        if label not in seqs:
            print(f"ERROR: '{label}' not found in FASTA. Available: {list(seqs.keys())}")
            return

    model, alphabet, device = load_model()

    all_results = {}
    for label, cfg in SITES.items():
        seq = seqs[label]
        print(f"\n── {label} (len={len(seq)}) ──")
        print(f"   Context: {cfg['wt_context']}")

        scores = masked_marginal_scores(model, alphabet, device, seq, cfg["positions_1idx"])
        all_results[label] = scores

        # 打印 V265F 这个关键突变的分数
        for pos1, pos_scores in scores.items():
            wt = pos_scores["WT"]
            print(f"   pos {pos1} WT={wt}")
            # 打印最极端的 top5 negative 和 positive
            mut_scores = {aa: pos_scores[aa] for aa in AAS}
            sorted_muts = sorted(mut_scores.items(), key=lambda x: x[1])
            print(f"   最负（功能破坏风险）: {sorted_muts[:3]}")
            print(f"   最正（较为保守）:     {sorted_muts[-3:]}")
            if "F" in mut_scores and label == "WEEV_71V_E2":
                print(f"   *** V265F score = {mut_scores['F']:.4f}  (论文中 V265F 使 71V 株获得 VLDLR 结合能力)")

    # 保存
    out_path = os.path.join(DATA_DIR, "all_scores.json")
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n打分结果 → {out_path}")

if __name__ == "__main__":
    main()
