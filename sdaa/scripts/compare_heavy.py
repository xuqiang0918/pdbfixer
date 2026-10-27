"""细分坐标差异：重原子 vs 氢原子。

pdbfixer 不移动已有重原子，只新增缺失原子和加氢。如果差异全部落在
「新增原子 + 氢的位置」上，说明是加氢优化在两个平台的浮点结果不同
导致的分叉，而不是结构重建出错。
"""
import argparse
import os


def read_atoms(path):
    atoms = []
    for line in open(path):
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        atoms.append({
            "name": line[12:16].strip(),
            "element": (line[76:78].strip() or line[12:16].strip()[0]).upper(),
            "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
        })
    return atoms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir-a", required=True)
    ap.add_argument("--dir-b", required=True)
    args = ap.parse_args()
    cases = sorted(set(f[:-4] for f in os.listdir(args.dir_a) if f.endswith(".pdb")) &
                   set(f[:-4] for f in os.listdir(args.dir_b) if f.endswith(".pdb")))
    print("%-6s %8s %14s %14s %14s %14s" % (
        "pdbid", "n_heavy", "heavy_max_d", "H_max_d", "nH", "heavy_same_atoms"))
    print("-" * 78)
    for c in cases:
        aa = read_atoms(os.path.join(args.dir_a, c + ".pdb"))
        bb = read_atoms(os.path.join(args.dir_b, c + ".pdb"))
        if len(aa) != len(bb):
            print("%-6s %8s  (atom count %d vs %d)" % (c, "-", len(aa), len(bb)))
            continue
        hd, hh = [], []
        nH = 0
        same_heavy = 0
        n_heavy = 0
        for a, b in zip(aa, bb):
            d = max(abs(a["xyz"][i] - b["xyz"][i]) for i in range(3))
            if a["element"] == "H":
                hh.append(d)
                nH += 1
            else:
                hd.append(d)
                n_heavy += 1
                if a["name"] == b["name"] and d == 0.0:
                    same_heavy += 1
        print("%-6s %8d %14.3e %14.3e %14d %14d" % (
            c, n_heavy,
            max(hd) if hd else 0.0,
            max(hh) if hh else 0.0,
            nH, same_heavy))


if __name__ == "__main__":
    main()
