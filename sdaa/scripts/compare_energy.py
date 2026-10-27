"""最小化后能量对照 —— SDAA vs x86，以及各平台自身重复的噪声地板。

判据：平台间能量差应当与「同平台两次运行的差」同量级。
"""
import json
import sys

a = json.load(open(sys.argv[1]))  # SDAA
b = json.load(open(sys.argv[2]))  # x86
la = sys.argv[3] if len(sys.argv) > 3 else "SDAA"
lb = sys.argv[4] if len(sys.argv) > 4 else "x86"

keys = sorted(set(a) & set(b))
print("%-6s %14s %14s %12s %8s %10s" % (
    "pdbid", la + "_min", lb + "_min", "dE", "dE/natom", "dE_pct"))
print("-" * 76)
diffs = []
for k in keys:
    if not (a[k].get("ok") and b[k].get("ok")):
        continue
    ea = a[k]["e_min"]
    eb = b[k]["e_min"]
    d = ea - eb
    n = a[k].get("n_atoms") or 1
    diffs.append(d)
    print("%-6s %14.2f %14.2f %12.2f %8.4f %9.2f%%" % (
        k, ea, eb, d, d / n, 100.0 * d / max(abs(eb), 1.0)))

if diffs:
    print()
    print("能量差 |dE| 均值      : %.2f kcal/mol" % (sum(abs(d) for d in diffs) / len(diffs)))
    print("能量差 |dE| 最大      : %.2f kcal/mol" % max(abs(d) for d in diffs))
    print("相对差均值            : %.3f%%" %
          (sum(100.0 * abs(d) / max(abs(b[k]["e_min"]), 1.0)
               for k, d in zip([k for k in keys if a[k].get("ok") and b[k].get("ok")], diffs))
           / len(diffs)))
