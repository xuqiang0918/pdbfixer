"""比对两平台的 pdbfixer 输出摘要（summary.json）。

用法:
    python compare_xcheck.py a/summary.json b/summary.json
"""
import json
import sys


def main():
    a = json.load(open(sys.argv[1]))
    b = json.load(open(sys.argv[2]))
    la = sys.argv[3] if len(sys.argv) > 3 else "A"
    lb = sys.argv[4] if len(sys.argv) > 4 else "B"

    keys = sorted(set(a) | set(b))
    print("%-6s %-5s %6s %6s | %7s %7s %7s | %7s %7s %7s | %6s %6s" % (
        "pdbid", "ok", "atomA", "atomB", "nresA", "nresB", "nH_A", "nH_B",
        "nH_pct_A", "nH_pct_B", "wallA", "wallB"))
    print("-" * 96)
    n_same_atoms = 0
    n_same_h = 0
    n_total = 0
    n_ok = 0
    for k in keys:
        ra = a.get(k, {})
        rb = b.get(k, {})
        oka = ra.get("ok", False)
        okb = rb.get("ok", False)
        if oka and okb:
            n_ok += 1
        ea = ra.get("elements", {})
        eb = rb.get("elements", {})
        ha = ea.get("H", 0)
        hb = eb.get("H", 0)
        ta = ra.get("n_atoms", 0)
        tb = rb.get("n_atoms", 0)
        if oka and okb:
            n_total += 1
            if ta == tb:
                n_same_atoms += 1
            if ha == hb:
                n_same_h += 1
        print("%-6s %-5s %6d %6d | %7d %7d %7d %7d | %7s %7s | %6s %6s" % (
            k, "%s/%s" % ("T" if oka else "F", "T" if okb else "F"),
            ta, tb,
            ra.get("n_residues", 0), rb.get("n_residues", 0), ha, hb,
            "%.2f%%" % (100.0 * ha / ta) if ta else "-",
            "%.2f%%" % (100.0 * hb / tb) if tb else "-",
            ra.get("wall_s", "-"), rb.get("wall_s", "-")))

    print()
    print("两边都成功        : %d / %d" % (n_ok, len(keys)))
    if n_total:
        print("原子数完全一致    : %d / %d" % (n_same_atoms, n_total))
        print("氢原子数完全一致  : %d / %d" % (n_same_h, n_total))
    # 逐 case 详细差异
    diffs = []
    for k in keys:
        ra, rb = a.get(k, {}), b.get(k, {})
        if not (ra.get("ok") and rb.get("ok")):
            continue
        d = {}
        for f in ("n_atoms", "n_residues", "n_chains", "n_bonds"):
            if ra.get(f) != rb.get(f):
                d[f] = (ra.get(f), rb.get(f))
        ea, eb = ra.get("elements", {}), rb.get("elements", {})
        for sym in set(ea) | set(eb):
            if ea.get(sym, 0) != eb.get(sym, 0):
                d["elem_" + sym] = (ea.get(sym, 0), eb.get(sym, 0))
        if ra.get("residue_name_counts") != rb.get("residue_name_counts"):
            d["residue_counts"] = "differ"
        if d:
            diffs.append((k, d))
    print("\n有差异的 case: %d" % len(diffs))
    for k, d in diffs:
        print("  %-6s %s" % (k, d))
    if not diffs:
        print("  （全部 case 的原子/残基/链/键/元素计数完全一致）")


if __name__ == "__main__":
    main()
