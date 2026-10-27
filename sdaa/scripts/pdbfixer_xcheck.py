"""pdbfixer 跨平台输出对照脚本（在两边各跑一次，比对结果文件）。

对每个 PDB ID 走一遍 pdbfixer 的标准流程：
    findMissingResidues -> findNonstandardResidues -> replaceNonstandardResidues
    -> removeHeterogens -> findMissingAtoms -> addMissingAtoms -> addMissingHydrogens
    -> save

输出一个 PDB 文件 + 一份 JSON 摘要（原子数/残基数/氢数/元素计数），
两个平台比对这两个产物即可。

用法:
    python pdbfixer_xcheck.py --outdir /tmp/pdbx_sdaa --pdbids 1LE1,1BNA
"""
import argparse
import json
import os
import time
from collections import Counter

from pdbfixer import PDBFixer
from openmm import app


def summarize(fixer):
    """从 topology 提取与平台无关的结构摘要。"""
    elems = Counter()
    for atom in fixer.topology.atoms():
        elems[atom.element.symbol if atom.element else "?"] += 1
    residues = [r.name for r in fixer.topology.residues()]
    chains = [c.id for c in fixer.topology.chains()]
    return {
        "n_atoms": fixer.topology.getNumAtoms(),
        "n_residues": fixer.topology.getNumResidues(),
        "n_chains": fixer.topology.getNumChains(),
        "n_bonds": fixer.topology.getNumBonds(),
        "elements": dict(sorted(elems.items())),
        "residue_names": sorted(set(residues)),
        "residue_name_counts": dict(sorted(Counter(residues).items())),
        "chain_ids": chains,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--pdbids", required=True)
    ap.add_argument("--ph", type=float, default=7.0)
    ap.add_argument("--platform", default="",
                    help="OpenMM platform name, e.g. Reference / CPU. "
                         "空 = OpenMM 自动选择")
    args = ap.parse_args()

    platform = None
    if args.platform:
        import openmm
        platform = openmm.Platform.getPlatformByName(args.platform)
        print("platform = %s" % platform.getName(), flush=True)

    os.makedirs(args.outdir, exist_ok=True)
    summary = {}
    for pdbid in [p.strip() for p in args.pdbids.split(",") if p.strip()]:
        t0 = time.time()
        rec = {"pdbid": pdbid}
        try:
            fixer = PDBFixer(pdbid=pdbid, platform=platform)
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            rec["nonstandard_found"] = len(fixer.nonstandardResidues)
            fixer.replaceNonstandardResidues()
            fixer.removeHeterogens(False)
            fixer.findMissingAtoms()
            rec["missing_atoms"] = len(fixer.missingAtoms)
            rec["missing_terminals"] = len(fixer.missingTerminals)
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
