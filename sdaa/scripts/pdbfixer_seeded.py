"""固定 Python random 种子后验证 pdbfixer 加氢的可重复性。

归因假设：OpenMM Modeller.addHydrogens 把氢先放在随机位置再做能量最小化，
而 Python random 未设种子，导致同一输入两次运行的氢位置略有不同。

做法：设 random.seed(fixed) 后跑两遍，比对输出。
"""
import argparse
import json
import os
import random
import time
from collections import Counter

from pdbfixer import PDBFixer
from openmm import app


def summarize(fixer):
    elems = Counter()
    for atom in fixer.topology.atoms():
        elems[atom.element.symbol if atom.element else "?"] += 1
    residues = [r.name for r in fixer.topology.residues()]
    return {
        "n_atoms": fixer.topology.getNumAtoms(),
        "n_residues": fixer.topology.getNumResidues(),
        "n_chains": fixer.topology.getNumChains(),
        "n_bonds": fixer.topology.getNumBonds(),
        "elements": dict(sorted(elems.items())),
        "residue_name_counts": dict(sorted(Counter(residues).items())),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--pdbids", required=True)
    ap.add_argument("--ph", type=float, default=7.0)
    ap.add_argument("--seed", type=int, default=-1, help="-1 = 不设种子")
    ap.add_argument("--platform", default="")
    args = ap.parse_args()

    platform = None
    if args.platform:
        import openmm
        platform = openmm.Platform.getPlatformByName(args.platform)

    os.makedirs(args.outdir, exist_ok=True)
    summary = {}
    for pdbid in [p.strip() for p in args.pdbids.split(",") if p.strip()]:
        t0 = time.time()
        rec = {"pdbid": pdbid}
        try:
            if args.seed >= 0:
                random.seed(args.seed)
            fixer = PDBFixer(pdbid=pdbid, platform=platform)
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.removeHeterogens(False)
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(args.ph)
            rec.update(summarize(fixer))
            outpdb = os.path.join(args.outdir, pdbid + ".pdb")
            with open(outpdb, "w") as f:
                app.PDBFile.writeFile(fixer.topology, fixer.positions, f)
            rec["out_pdb"] = outpdb
            rec["ok"] = True
        except Exception as e:
            rec["ok"] = False
            rec["error"] = "%s: %s" % (type(e).__name__, e)
        rec["wall_s"] = round(time.time() - t0, 2)
        summary[pdbid] = rec
        print("%-6s ok=%s  %s" % (pdbid, rec["ok"], rec.get("error", "")), flush=True)

    with open(os.path.join(args.outdir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
    print("summary -> %s" % os.path.join(args.outdir, "summary.json"))


if __name__ == "__main__":
    main()
