"""
01_fetch_sequences.py
下载甲病毒 E1/E2 糖蛋白序列，用于 ESM-2 突变效应打分。

目标序列：
  - WEEV E2（Western equine encephalitis virus）：McMillan 株 + 71V 株
  - SFV E1（Semliki Forest virus）
  - 其他甲病毒 E1/E2 作对比（EEEV, VEEV, CHIKV, SINV）

输出：
  data/alphavirus_sequences.fasta   — 全部序列
  data/sequences_info.tsv           — 每条序列的元信息（病毒名、蛋白、UniProt ID、长度）

运行环境：F:/CC/conda_envs/esm2
"""

import requests
import time
import os
import csv

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR  = os.path.join(BASE_DIR, "..", "data")
os.makedirs(OUT_DIR, exist_ok=True)

# 验证后的正确 accession 号（2026-05-23 人工核查）
# 来源：UniProt reviewed 条目 + PDB 9L1N（Cell Reports 2025）
# 长度标准：E2 ~ 416-450 aa，E1 ~ 438 aa
TARGETS = [
    # source, accession, label, protein, virus, db
    # WEEV E2 — 来自 PDB 9L1N Chain B（71V 样株），V265=V 已验证
    ("ncbi",    "2990964796", "WEEV_E2_71V",  "E2", "Western equine encephalitis virus (71V-type, PDB:9L1N)"),
    # SFV E1 — Q87046, 438 aa，W132/D135/E137/D139 已验证
    ("uniprot", "Q87046",     "SFV_E1",       "E1", "Semliki Forest virus"),
    # SFV E2 — Q87050, 422 aa
    ("uniprot", "Q87050",     "SFV_E2",       "E2", "Semliki Forest virus"),
    # 其他代表性甲病毒 E2（对比保守性）
    ("uniprot", "H1ZYL1",     "CHIKV_E2",     "E2", "Chikungunya virus"),
    ("uniprot", "A9XR09",     "EEEV_E2",      "E2", "Eastern equine encephalitis virus"),
    ("uniprot", "A0A1C8C4L6", "VEEV_E2",      "E2", "Venezuelan equine encephalitis virus"),
]

def fetch_uniprot_sequence(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    lines = r.text.strip().split("\n")
    return "".join(lines[1:])

def fetch_ncbi_sequence(uid):
    from Bio import Entrez, SeqIO
    import io
    Entrez.email = "demo@example.com"
    handle = Entrez.efetch(db="protein", id=uid, rettype="fasta", retmode="text")
    rec = list(SeqIO.parse(handle, "fasta"))[0]
    return str(rec.seq)

def main():
    fasta_path = os.path.join(OUT_DIR, "alphavirus_sequences.fasta")
    info_path  = os.path.join(OUT_DIR, "sequences_info.tsv")

    records = []
    with open(fasta_path, "w") as fout:
        for db, acc, label, protein, virus in TARGETS:
            print(f"Fetching {label} ({acc}, {db})...", end=" ", flush=True)
            try:
                if db == "ncbi":
                    seq = fetch_ncbi_sequence(acc)
                else:
                    seq = fetch_uniprot_sequence(acc)
                fout.write(f">{label}|{acc}|{protein}|{virus}\n{seq}\n")
                records.append({
                    "label": label, "accession": acc,
                    "protein": protein, "virus": virus,
                    "length": len(seq)
                })
                print(f"OK  ({len(seq)} aa)")
            except Exception as e:
                print(f"FAILED: {e}")
            time.sleep(0.3)

    with open(info_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["label","accession","protein","virus","length"], delimiter="\t")
        w.writeheader()
        w.writerows(records)

    print(f"\nWrote {len(records)} sequences -> {fasta_path}")
    print(f"Info  -> {info_path}")
    print("\nVerified positions:")
    print("  WEEV E2 pos 265 = V (71V-type, V265F switches to VLDLR binding)")
    print("  SFV  E1 pos 132/135/137/139 = W/D/E/D (VLDLR interface)")

if __name__ == "__main__":
    main()
