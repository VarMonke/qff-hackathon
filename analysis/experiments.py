"""Experiments behind the cost / noise analysis. Writes analysis/results.json.

    uv run python analysis/experiments.py [A|B|C ...]     (default: all)

A  shot level:   one 132-bit order shot, classified clean / rescued by list decoding / bad,
                 versus noise; compared with model.order_shot.
B  solver level: success rate, queries and score of solve() versus noise.
C  size:         queries of solve() versus log2(n) at the competition noise, prime orders;
                 compared with model.expected_queries and model.floor_queries.
"""

import json
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import solver  # noqa: E402
from silentgroup import MockGroup, Oracle, load_public  # noqa: E402
from silentgroup.check import score  # noqa: E402
from silentgroup.program import compile_circuit, validate  # noqa: E402
from silentgroup.simulator import simulate  # noqa: E402

OUT = Path(__file__).with_name("results.json")
NOISES = [0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.008]


def random_prime(bits, rng):
    while True:
        p = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if solver._is_prime(p):
            return p


def exp_a(shots=300, seed=1):
    rng = random.Random(seed)
    n = random_prime(59, rng)
    inst = MockGroup([n], seed=seed)  # cyclic of prime order: chi(g) = j/n, j uniform
    m = solver.ORDER_BITS if hasattr(solver, "ORDER_BITS") else 2 * solver.NMAX_BITS + 12
    qc, _ = solver.ipe_circuit([("g", m)])
    nq, nc, ops = validate(compile_circuit(qc), ["A", "B", "C", "g"], 16)
    rows = []
    for sigma in NOISES:
        clean = rescued = bad = 0
        for _ in range(shots):
            bits = next(iter(simulate(inst, nq, nc, ops, 1, rng, sigma)))[::-1]
            est = sum(int(bits[k]) << (m - 1 - k) for k in range(m))
            c, dec = solver.decode(est, m, solver.NMAX_BITS)
            if c == n:
                clean += 1
            elif n in dec[:6]:  # the solver confirms up to 6 decodings
                rescued += 1
            else:
                bad += 1
        rows.append({"sigma": sigma, "shots": shots, "clean": clean, "rescued": rescued, "bad": bad})
        print("A", rows[-1], flush=True)
    return {"n": n, "m": m, "rows": rows}


def run_solver(inst, sigma, seed):
    o = Oracle(inst, budget=80_000, seed=seed, noise=sigma)
    t = time.monotonic()
    try:
        res = solver.solve(o)
    except Exception as exc:  # solve() already guards; this is belt and braces
        res = repr(exc)
    return score(res, inst.answer, 80_000, o.queries_used), o.queries_used, time.monotonic() - t


def exp_b(reps=15, seed=2):
    rng = random.Random(seed)
    cases = {"public_1": load_public("public_1"), "public_2": load_public("public_2"),
             "prime60": MockGroup([random_prime(60, rng)], seed=seed)}
    rows = []
    for name, inst in cases.items():
        for sigma in NOISES:
            runs = [run_solver(inst, sigma, rng.randrange(1 << 30)) for _ in range(reps)]
            rows.append({"case": name, "sigma": sigma, "reps": reps,
                         "success": sum(p >= 1 for p, _, _ in runs) / reps,
                         "order_only": sum(0 < p < 1 for p, _, _ in runs) / reps,
                         "queries": sum(q for _, q, _ in runs) / reps,
                         "score": sum(p for p, _, _ in runs) / reps,
                         "seconds": max(t for _, _, t in runs)})
            print("B", rows[-1], flush=True)
    return {"rows": rows}


def exp_c(reps=12, seed=3, sigma=0.002):
    rng = random.Random(seed)
    rows = []
    for L in range(12, 61, 4):
        qs = []
        for _ in range(reps):
            inst = MockGroup([random_prime(L, rng)], seed=rng.randrange(1 << 30))
            pts, q, _ = run_solver(inst, sigma, rng.randrange(1 << 30))
            assert pts >= 1, (L, inst.answer)
            qs.append(q)
        qs.sort()
        rows.append({"L": L, "reps": reps, "mean": sum(qs) / reps, "median": qs[reps // 2],
                     "min": qs[0], "max": qs[-1]})
        print("C", rows[-1], flush=True)
    return {"sigma": sigma, "rows": rows}


if __name__ == "__main__":
    which = sys.argv[1:] or ["A", "B", "C"]
    results = json.loads(OUT.read_text()) if OUT.exists() else {}
    for w in which:
        results[w] = {"A": exp_a, "B": exp_b, "C": exp_c}[w]()
        OUT.write_text(json.dumps(results, indent=1))
    print("wrote", OUT)
