"""对 VLDLR LA3 关键结合位点打 ESM-2 masked marginal scores"""
import torch, esm, json, os

torch.hub.set_dir("F:/CC/.cache/torch/hub")
DATA_DIR = "F:/CC/projects/cancer_perturb/文书/套磁信/章新政/02_定制层/data"
AAS = list("ACDEFGHIKLMNPQRSTVWY")

# 读取 VLDLR 序列
with open(f"{DATA_DIR}/vldlr_human.fasta") as f:
    lines = f.read().strip().split("\n")
seq = "".join(lines[1:])
print(f"VLDLR len={len(seq)}, checking positions:")
for p in [129, 132, 135, 137, 139]:
    print(f"  pos {p}: {seq[p-1]}")

# 加载模型
print("\nLoading ESM-2 650M...", flush=True)
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = model.to(device).eval()
print(f"  Device: {device}")

batch_converter = alphabet.get_batch_converter()
_, _, tokens_orig = batch_converter([("seq", seq)])

positions = [129, 132, 135, 137, 139]
scores = {}
for pos1 in positions:
    pos0 = pos1 - 1
    wt_aa = seq[pos0]
    tokens = tokens_orig.clone().to(device)
    tokens[0, pos0 + 1] = alphabet.mask_idx
    with torch.no_grad():
        res = model(tokens, repr_layers=[], return_contacts=False)
    logits = res["logits"][0, pos0 + 1]
    lp = torch.log_softmax(logits, dim=-1).cpu()
    wt_lp = lp[alphabet.get_idx(wt_aa)].item()
    pos_scores = {"WT": wt_aa}
    for aa in AAS:
        pos_scores[aa] = round(lp[alphabet.get_idx(aa)].item() - wt_lp, 4)
    scores[pos1] = pos_scores
    muts = [(aa, pos_scores[aa]) for aa in AAS]
    muts.sort(key=lambda x: x[1])
    print(f"  pos {pos1} WT={wt_aa}: most negative={muts[:2]}, V265-equiv check:")
    if pos1 == 129:
        print(f"    P129H score = {pos_scores.get('H', 'N/A'):.4f}  (should be large negative)")
    if pos1 == 132:
        print(f"    W132A score = {pos_scores.get('A', 'N/A'):.4f}  (W is essential, A should be very negative)")

# 加载已有 WEEV 结果并合并
existing = json.load(open(f"{DATA_DIR}/all_scores.json"))
existing["VLDLR_human"] = scores

with open(f"{DATA_DIR}/all_scores.json", "w") as f:
    json.dump(existing, f, indent=2)
print(f"\nUpdated all_scores.json with VLDLR data")
