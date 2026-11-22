"""Randomized stress test: python stress.py [runs] [seed]

Random groups Z_d1 x ... x Z_dk (k = 1..3, |G| <= 2^60), random g, random forgeries,
scored exactly like the judge. Prints every failure with the solver log.
"""

import math
import random
import sys
import time

import solver
from silentgroup import MockGroup, Oracle
from silentgroup.check import score


def random_modulus(rng, max_bits):
    max_bits = max(max_bits, 2)
    kind = rng.random()
    if kind < 0.3:  # prime-ish
        return rng.randrange(2, 1 << max_bits) | 1
    if kind < 0.6:  # smooth
        d = 1
        while True:
            p = rng.choice([2, 3, 5, 7, 11, 13])
            if (d * p).bit_length() > max_bits:
                return max(d, 2)
            d *= p
    if kind < 0.8:  # near a power of two
        b = rng.randrange(2, max_bits + 1)
        return max(2, (1 << b) - rng.randrange(0, 200))
    return rng.randrange(2, 1 << rng.randrange(2, max_bits + 1))


def random_instance(rng):
    k = rng.choice([1, 1, 2, 2, 3])
    budget_bits = rng.randrange(8, 61)
    moduli, left = [], budget_bits
    for i in range(k):
        bits = left if i == k - 1 else rng.randrange(1, max(2, left - (k - 1 - i)))
        bits = max(1, min(bits, left))
        moduli.append(random_modulus(rng, bits))
        left -= moduli[-1].bit_length()
        if left < 2:
            break
    while math.prod(moduli).bit_length() > 60:
        moduli[-1] //= 2
    moduli = [max(2, d) for d in moduli]
    forged = ""
    if len(moduli) > 1 and rng.random() < 0.4:
        forged = rng.choice(["A", "B", "C", "AB", "ABC"])
    try:
        return MockGroup(moduli, forged=forged, seed=rng.randrange(1 << 30)), moduli, forged
    except ValueError:  # could not place a forgery (cyclic quotient); fall back to honest
        return MockGroup(moduli, seed=rng.randrange(1 << 30)), moduli, ""


def main():
    runs = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    total = bad = 0
    qs = []
    t0 = time.monotonic()
    for i in range(runs):
        inst, moduli, forged = random_instance(rng)
        o = Oracle(inst, budget=80_000)
        t = time.monotonic()
        try:
            res = solver.solve(o)
        except Exception as e:
            res = repr(e)
        dt = time.monotonic() - t
        pts = score(res, inst.answer, 80_000, o.queries_used)
        total += pts
        qs.append(o.queries_used)
        if pts < 1 or dt > 60 or o.queries_used > 50_000:
            bad += pts < 1
            print(f"{'BAD' if pts < 1 else 'SLOW'} #{i} moduli={moduli} forged={forged!r} pts={pts:.2f} q={o.queries_used} t={dt:.1f}s"
                  f"\n    got={res} want={inst.answer}\n    log={getattr(solver.solve, 'log', None)}", flush=True)
    qs.sort()
    print(f"{runs} runs, {bad} bad, avg pts {total / runs:.4f}, median q {qs[len(qs) // 2]}, "
          f"p90 q {qs[int(len(qs) * 0.9)]}, max q {qs[-1]}, {time.monotonic() - t0:.0f}s")


if __name__ == "__main__":
    main()
