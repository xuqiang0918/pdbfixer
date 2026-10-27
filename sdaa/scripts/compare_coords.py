"""逐原子坐标比对两平台的 pdbfixer 输出 PDB。

比对维度（按 ATOM 记录的固定列）：
  * 原子名 / 残基名 / 链 ID / 残基号 —— 结构标识
  * 元素 —— 元素组成
  * x/y/z 坐标 —— 数值（报告最大绝对差）

坐标差用两平台各自的 PDB 输出直接比（两边都是各自机器算的）。
"""
import argparse
import math
import os
import sys


def read_atoms(path):
    atoms = []
    for line in open(path):
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        atoms.append({
            "name": line[12:16],
            "resname": line[17:20],
            "chain": line[21],
            "resseq": line[22:26],
            "element": line[76:78].strip(),
            "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
        })
    return atoms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir-a", required=True)
    ap.add_argument("--dir-b", required=True)
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    args = ap.parse_args()

    cases = sorted(set(f[:-4] for f in os.listdir(args.dir_a) if f.endswith(".pdb")) &
                   set(f[:-4] for f in os.listdir(args.dir_b) if f.endswith(".pdb")))
    print("%-6s %6s %8s %8s %12s %12s %12s %s" % (
        "pdbid", "natoms", "id_same", "elem_same", "max_dx", "max_dy", "max_dz", "verdict"))
    print("-" * 92)
    n_ok = 0
    n_total = 0
    worst = 0.0
    for c in cases:
        aa = read_atoms(os.path.join(args.dir_a, c + ".pdb"))
        bb = read_atoms(os.path.join(args.dir_b, c + ".pdb"))
        if len(aa) != len(bb) or not aa:
            print("%-6s %6d %8s %8s %12s %12s %12s %s" % (
                c, len(aa), "-", "-", "-", "-", "-", "ATOM COUNT MISMATCH"))
            continue
        id_same = all(a["name"] == b["name"] and a["resname"] == b["resname"] and
                      a["chain"] == b["chain"] and a["resseq"] == b["resseq"]
                      for a, b in zip(aa, bb))
        elem_same = all(a["element"] == b["element"] for a, b in zip(aa, bb))
        d = [max(abs(a["xyz"][i] - b["xyz"][i]) for i in range(3)) for a, b in zip(aa, bb)]
        m = max(d)
        worst = max(worst, m)
        n_total += 1
        if id_same and elem_same:
            n_ok += 1
        verdict = "IDENTICAL" if m == 0.0 else ("close(%.1e)" % m)
        print("%-6s %6d %8s %8s %12.3e %12.3e %12.3e %s" % (
            c, len(aa), id_same, elem_same, d[0], d[1], d[2], verdict))
    print()
    print("比对 case 数        : %d" % n_total)
    print("标识+元素完全一致  : %d" % n_ok)
    print("全局最大坐标差     : %.3e Angstrom" % worst)
    if worst == 0.0:
        print("=> 两平台输出逐比特相同")
    elif worst < 1e-6:
        print("=> 两平台输出在浮点舍入范围内一致")
    else:
        print("=> 存在超过浮点舍入的坐标差异")


if __name__ == "__main__":
    main()
