"""数 1R8J addMissingAtoms/addMissingHydrogens 里 L-BFGS 的迭代次数。

Monkey-patch LocalEnergyMinimizer.minimize，挂 MinimizationReporter 计数每次调用的迭代数，
然后走 pdbfixer_seeded 同样的流程（seed=12345，只跑 1R8J）。
"""
import json
import random
import sys
import time

import openmm
from openmm.openmm import LocalEnergyMinimizer, MinimizationReporter

CALLS = []


class IterCounter(MinimizationReporter):
    def __init__(self):
        super().__init__()
        self.n = 0

    def report(self, *args):
        self.n += 1
        return False


_orig = LocalEnergyMinimizer.minimize


def _patched(context, tolerance=10, maxIterations=0, reporter=None):
    c = IterCounter()
    t0 = time.time()
    r = _orig(context, tolerance, maxIterations, c)
    CALLS.append({"iters": c.n, "wall_s": round(time.time() - t0, 2)})
    return r


LocalEnergyMinimizer.minimize = staticmethod(_patched)

from pdbfixer import PDBFixer  # noqa: E402
from openmm import app  # noqa: E402

random.seed(12345)
t0 = time.time()
fixer = PDBFixer(pdbid="1R8J")
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(False)
t_add_atoms = time.time()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
t_atoms_done = time.time()
CALLS_atoms = list(CALLS)
CALLS.clear()
fixer.addMissingHydrogens(7.0)
t_h_done = time.time()
CALLS_h = list(CALLS)

out = {
    "platforms": [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())],
    "minimize_calls_during_addMissingAtoms": CALLS_atoms,
    "minimize_calls_during_addMissingHydrogens": CALLS_h,
    "total_iters_atoms": sum(c["iters"] for c in CALLS_atoms),
    "total_iters_h": sum(c["iters"] for c in CALLS_h),
    "wall_addMissingAtoms_s": round(t_atoms_done - t_add_atoms, 2),
    "wall_addMissingHydrogens_s": round(t_h_done - t_atoms_done, 2),
    "wall_total_s": round(time.time() - t0, 2),
}
print(json.dumps(out, indent=2))
with open(sys.argv[1], "w") as f:
    json.dump(out, f, indent=2)
