"""从 pdbfixer 输出的 PDB 文件提取结构摘要（不依赖生成进程完成）。

用法:
    python pdb_summary.py <pdb_dir> <out_json>
"""
import json
import os
import sys
from collections import Counter


def parse_pdb(path):
    atoms = []
    residues = Counter()
    chains = []
    for line in open(path):
        rec = line[:6]
        if rec == "ATOM  " or rec == "HETATM":
            name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21]
            resseq = line[22:26].strip()
            element = line[76:78].strip()
            if not element:
                # 从原子名推断
                element = name[0] if name[0].isalpha() else name[1]
            atoms.append(element)
            residues[(chain, resseq, resname)] += 1
            if chain not in chains:
                chains.append(chain)
    elems = Counter(atoms)
    return {
        "n_atoms": len(atoms),
        "n_residues": len(residues),
        "n_chains": len(chains),
        "chain_ids": chains,
        "elements": dict(sorted(elems.items())),
        "residue_name_counts": dict(sorted(Counter(k[2] for k in residues).items())),
    }


def main():
    pdb_dir = sys.argv[1]
    out_json = sys.argv[2]
    summary = {}
    for fn in sorted(os.listdir(pdb_dir)):
        if not fn.endswith(".pdb"):
            continue
        pdbid = fn[:-4]
        path = os.path.join(pdb_dir, fn)
        try:
            s = parse_pdb(path)
            s["ok"] = True
            s["out_pdb"] = path
        except Exception as e:
            s = {"ok": False, "error": str(e)}
        summary[pdbid] = s
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
    print("wrote %s (%d cases)" % (out_json, len(summary)))


if __name__ == "__main__":
    main()
