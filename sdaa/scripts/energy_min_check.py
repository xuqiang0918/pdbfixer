"""能量最小化后的能量对照 —— 物理等价性的硬判据。

原始产物的能量对氢位置敏感（坐标差 0.1 Å 就能让能量差几百 kcal/mol），
所以直接比初始能量没有意义。做法：
  1. 对两平台的产物各做一次能量最小化（相同力场、相同收敛判据）
  2. 比最小化后的能量

最小化后仍存在的能量差，反映的是「两个局部极小点的能量差」，
若在同量级（相对于体系大小），就说明产物是等价的构象。
"""
import argparse
import json
import os

from openmm import app, unit
import openmm as mm


def minimized_energy(pdb_path, maxiter=2000):
    pdb = app.PDBFile(pdb_path)
    ff = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds)
    integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)
    ctx = mm.Context(system, integrator)
    ctx.setPositions(pdb.positions)
    e0 = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilocalories_per_mole)
    mm.LocalEnergyMinimizer.minimize(ctx, 1.0 * unit.kilocalories_per_mole /
                                     unit.angstroms, maxiter)
    e1 = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilocalories_per_mole)
    del ctx, integrator, system
    return e0, e1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--maxiter", type=int, default=2000)
    args = ap.parse_args()
    out = {}
    for fn in sorted(os.listdir(args.dir)):
        if not fn.endswith(".pdb"):
            continue
        p = os.path.join(args.dir, fn)
        try:
            e0, e1 = minimized_energy(p, args.maxiter)
            out[fn[:-4]] = {"e_init": e0, "e_min": e1, "ok": True}
            print("%-6s init %14.2f -> min %14.2f kcal/mol" % (fn[:-4], e0, e1),
                  flush=True)
        except Exception as e:
            out[fn[:-4]] = {"ok": False, "error": str(e)}
            print("%-6s FAILED %s" % (fn[:-4], e), flush=True)
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
