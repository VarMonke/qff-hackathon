"""Compare solver parameter settings on identical instances.

    python sweep.py RUNS SEED NAME=VALUE[,NAME=VALUE...] [more settings ...]

Each setting is a comma-separated list of solver module overrides (e.g. PROBE_BITS=0).
Runs the random stress mix plus the fixed bench cases, in parallel, and prints per-setting
average points, bad runs, and average points on the fixed cases.
"""

import multiprocessing as mp
import random
import sys

import solver
from silentgroup import MockGroup, Oracle, load_public
from silentgroup.check import score
from stress import random_instance

FIXED = {
    "public_1": lambda s: load_public("public_1"),
    "public_2": lambda s: load_public("public_2"),
    "cyclic60": lambda s: MockGroup([(1 << 60) - 93], seed=s),
    "smooth": lambda s: MockGroup([2**20 * 3**10 * 5**6, 2**10 * 3**4], seed=s),
    "forgedC": lambda s: MockGroup([2**20 * 3**5, 2**8 * 3**3], forged="C", seed=s),
    "forgedAB": lambda s: MockGroup([1000003 * 999983, 1000003 * 7], forged="AB", seed=s),
}


def parse(setting):
    out = {}
    for kv in filter(None, setting.split(",")):
        k, v = kv.split("=")
        out[k] = int(v)
    return out


def job(args):
    setting, kind, idx, seed = args
    for k, v in parse(setting).items():
        setattr(solver, k, v)
    if kind == "stress":
        inst = random_instance(random.Random(seed * 100003 + idx))[0]
    else:
        inst = FIXED[kind](idx)
    o = Oracle(inst, budget=80_000, seed=seed * 7919 + idx * 31 + hash(kind) % 1000)
    try:
        res = solver.solve(o)
    except Exception as e:
        res = repr(e)
    pts = score(res, inst.answer, 80_000, o.queries_used)
    info = None
    if pts < 1:
        info = f"{kind}#{idx} got={res} want={inst.answer} log={solver.solve.log}"
    return setting, kind, pts, o.queries_used, info


def main():
    runs, seed = int(sys.argv[1]), int(sys.argv[2])
    settings = sys.argv[3:] or [""]
    fixed_reps = max(1, runs // 20)
    jobs = [(s, "stress", i, seed) for s in settings for i in range(runs)]
    jobs += [(s, k, i, seed) for s in settings for k in FIXED for i in range(fixed_reps)]
    stats = {s: {"stress": [], "fixed": {}, "bad": []} for s in settings}
    with mp.Pool() as pool:
        for s, kind, pts, q, info in pool.imap_unordered(job, jobs, chunksize=4):
            st = stats[s]
            (st["stress"] if kind == "stress" else st["fixed"].setdefault(kind, [])).append(pts)
            if info:
                st["bad"].append(info)
    for s in settings:
        st = stats[s]
        fx = {k: sum(v) / len(v) for k, v in st["fixed"].items()}
        print(f"[{s or 'baseline'}] stress avg {sum(st['stress']) / len(st['stress']):.4f}  "
              f"fixed avg {sum(fx.values()) / len(fx):.4f}  bad {len(st['bad'])}")
        print("    " + "  ".join(f"{k}={v:.3f}" for k, v in fx.items()))
        for b in st["bad"]:
            print("    BAD", b)


if __name__ == "__main__":
    main()
