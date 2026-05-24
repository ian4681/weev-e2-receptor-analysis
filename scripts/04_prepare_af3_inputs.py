"""准备 AlphaFold3 复合物预测的输入序列"""
import requests, os, textwrap

DATA_DIR = "F:/CC/projects/cancer_perturb/文书/套磁信/章新政/02_定制层/data"
OUT_DIR  = f"{DATA_DIR}/af3_inputs"
os.makedirs(OUT_DIR, exist_ok=True)

# ── 读取已有序列 ────────────────────────────────────────────────────────────
def read_fasta(path):
    seqs = {}
    with open(path) as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(seq)
                header, seq = line[1:].split("|")[0], []
            else:
                seq.append(line)
        if header:
            seqs[header] = "".join(seq)
    return seqs

alphav = read_fasta(f"{DATA_DIR}/alphavirus_sequences.fasta")
vldlr  = read_fasta(f"{DATA_DIR}/vldlr_human.fasta")

weev_wt = alphav["WEEV_E2_71V"]
vldlr_seq = list(vldlr.values())[0]

# 验证 V265
assert weev_wt[264] == "V", f"pos265={weev_wt[264]}"
print(f"WEEV E2 WT  len={len(weev_wt)},  pos265={weev_wt[264]} OK")
print(f"VLDLR       len={len(vldlr_seq)}")

# 生成 V265F 突变体
weev_f = weev_wt[:264] + "F" + weev_wt[265:]
assert weev_f[264] == "F"
print(f"WEEV E2 V265F: pos265={weev_f[264]} OK")

# ── 获取 PCDH10 (Q9Y5H8) ────────────────────────────────────────────────────
print("\n获取 PCDH10 (Q9Y5H8)...")
r = requests.get("https://rest.uniprot.org/uniprotkb/Q9Y5H8.fasta", timeout=30)
lines = r.text.strip().split("\n")
pcdh10_full = "".join(lines[1:])
print(f"PCDH10 full length: {len(pcdh10_full)}")

# 提取胞外域：去掉信号肽（~20aa）和跨膜域（~720aa之后）
# UniProt Q9Y5H8 注释：Signal 1-25, EC repeats 26-～706, TM ～707-727
# 保守起见取 25-710
pcdh10_ec = pcdh10_full[24:710]
print(f"PCDH10 EC domain (25-710): len={len(pcdh10_ec)}")
print(f"  起始: {pcdh10_ec[:10]}...")
print(f"  结尾: ...{pcdh10_ec[-10:]}")

# 保存 PCDH10 完整序列
with open(f"{DATA_DIR}/pcdh10_human.fasta", "w") as f:
    f.write(f">PCDH10_human|Q9Y5H8|PCDH10|Homo sapiens\n")
    f.write("\n".join(textwrap.wrap(pcdh10_full, 60)) + "\n")
print(f"已保存 pcdh10_human.fasta")

# ── 写 AF3 输入 FASTA（每个job一个文件）────────────────────────────────────
def wrap60(s):
    return "\n".join(textwrap.wrap(s, 60))

jobs = {
    "job1_WEEV_WT_VLDLR": [
        (">WEEV_E2_WT",  weev_wt),
        (">VLDLR_human", vldlr_seq),
    ],
    "job2_WEEV_V265F_VLDLR": [
        (">WEEV_E2_V265F", weev_f),
        (">VLDLR_human",   vldlr_seq),
    ],
    "job3_WEEV_WT_PCDH10": [
        (">WEEV_E2_WT",    weev_wt),
        (">PCDH10_EC",     pcdh10_ec),
    ],
    "job4_WEEV_V265F_PCDH10": [
        (">WEEV_E2_V265F", weev_f),
        (">PCDH10_EC",     pcdh10_ec),
    ],
}

print("\n── AF3 输入文件 ──")
for name, chains in jobs.items():
    path = f"{OUT_DIR}/{name}.fasta"
    total = sum(len(s) for _, s in chains)
    with open(path, "w") as f:
        for header, seq in chains:
            f.write(f"{header}\n{wrap60(seq)}\n")
    chains_summary = " + ".join(f"{h[1:]}({len(s)}aa)" for h, s in chains)
    print(f"  {name}.fasta  [{chains_summary}  total={total}aa]")

print("\n全部准备完成。")
print(f"输入文件在: {OUT_DIR}/")
print("\n下一步: 登录 alphafoldserver.com，依次提交4个job（每个上传对应fasta）")
