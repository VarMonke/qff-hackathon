"""Local benchmark: python bench.py [noise] [reps]"""

import sys
import time

import solver
from silentgroup import MockGroup, Oracle, load_public
from silentgroup.check import score

noise = float(sys.argv[1]) if len(sys.argv) > 1 else 0.002
reps = int(sys.argv[2]) if len(sys.argv) > 2 else 1

cases = [
    ("public_1", lambda s: load_public("public_1")),
    ("public_2", lambda s: load_public("public_2")),
    ("cyclic60", lambda s: MockGroup([(1 << 60) - 93], seed=s)),
    ("smooth", lambda s: MockGroup([2**20 * 3**10 * 5**6, 2**10 * 3**4], seed=s)),
    ("forgedC", lambda s: MockGroup([2**20 * 3**5, 2**8 * 3**3], forged="C", seed=s)),
    ("forgedAB", lambda s: MockGroup([1000003 * 999983, 1000003 * 7], forged="AB", seed=s)),
    ("tiny", lambda s: MockGroup([6, 4], seed=s)),
]

only = sys.argv[3].split(",") if len(sys.argv) > 3 else None
cases = [c for c in cases if not only or c[0] in only]
total = 0
for name, make in cases:
    case_pts, case_q = 0, 0
    for s in range(reps):
        inst = make(s)
        o = Oracle(inst, budget=80_000, noise=noise)
        t0 = time.monotonic()
        try:
            res = solver.solve(o)
        except Exception as e:
            res = repr(e)
        pts = score(res, inst.answer, 80_000, o.queries_used)
        total += pts
        case_pts += pts
        case_q += o.queries_used
        if pts < 1:
            print(f"BAD {name:9s} pts={pts:.3f} q={o.queries_used:6d} t={time.monotonic()-t0:5.1f}s"
                  f"  got={res} want={inst.answer} log={getattr(solver.solve, 'log', None)}")
    print(f"{name:9s} avg pts {case_pts / reps:.3f}  avg queries {case_q / reps:8.0f}", flush=True)
print(f"avg pts {total / (len(cases) * reps):.3f}")
