"""Open track: the full structure of H = <g, A, B, C> as Z_n1 x Z_n2 x ...  (n1 | n2 | ...)

    uv run python structure.py            # demo + ground-truth check on mock and public groups

Why it works
------------
Every shot samples a uniformly random character chi of G, and restriction G^ -> H^ is a
surjective homomorphism, so chi|H is a uniformly random character of H. One shot of joint
phase estimation on the generators h = (g, A, B, C) reads

    chi(h) = r / E,   r in Z_E^k,   E = exponent of H (lcm of the generators' orders).

The vectors r that occur form a subgroup R of Z_E^k isomorphic to H^ ~ H. Lift it to the
lattice R~ = R + E Z^k. Its dual (scaled by E) is the relation lattice

    Lambda = { x in Z^k : h_1^x_1 ... h_k^x_k = e } = E (R~)^-1,

and H ~ Z^k / Lambda, whose invariant factors are the Smith normal form of Lambda.

Noise: a shot whose estimate is corrupted but still looks clean adds a vector outside R,
making R too big (relations go missing). A shot set that is too small (missing characters)
makes Lambda too big (claims false relations). Both are caught on the device:
  * every basis relation of Lambda is checked with the chi-independent identity test, so a
    false relation is rejected (R too small -> keep sampling);
  * leave-one-out: if dropping shot i shrinks R and all the relations that reappear verify,
    shot i contradicted true relations and is discarded.

Elements whose discrete log the solver already verified add nothing to H, so only g and the
unresolved elements are sampled. An honest meeting costs no extra queries: H = Z_n.
"""

import math
from fractions import Fraction
from itertools import combinations

import solver

NAMES = ["g", "A", "B", "C"]
STABLE = 3  # consecutive clean shots that leave R unchanged before we try to certify it


# ---------------------------------------------------------------- integer lattices

def _egcd(a, b):
    if b == 0:
        return (abs(a), 1 if a >= 0 else -1, 0)
    g, s, t = _egcd(b, a % b)
    return g, t, s - (a // b) * t


def hnf(rows, E, k):
    """Triangular basis (k rows) of span(rows) + E Z^k. Columns are eliminated from the last
    to the first, so row i has zeros after column i and pivot rows[i][i] > 0. For a relation
    lattice that makes row i read  h_i^(pivot) * prod_{j<i} h_j^(x_j) = e."""
    rows = [[x % E for x in r] for r in rows]
    basis = [None] * k
    for c in range(k - 1, -1, -1):
        piv, rest = None, []
        # E e_c is in the lattice; adding it here (not upfront) survives the mod-E reduction.
        for r in rows + [[E * (j == c) for j in range(k)]]:
            if r[c] == 0:
                rest.append(r)
            elif piv is None:
                piv = r
            else:
                g, s, t = _egcd(piv[c], r[c])
                a, b = piv[c] // g, r[c] // g
                piv, r2 = [s * x + t * y for x, y in zip(piv, r)], [b * x - a * y for x, y in zip(piv, r)]
                rest.append(r2)
        if piv[c] < 0:
            piv = [-x for x in piv]
        basis[c] = [x % E if j < c else x for j, x in enumerate(piv)]
        rows = [[x % E for x in r] for r in rest if any(r)]
    # Reduce entries left of each pivot by the rows below it (canonical form).
    for i in range(k):
        for j in range(i - 1, -1, -1):
            q = basis[i][j] // basis[j][j]
            if q:
                basis[i] = [x - q * y for x, y in zip(basis[i], basis[j])]
    return tuple(tuple(r) for r in basis)


def _det(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    return sum((-1) ** j * M[0][j] * _det([r[:j] + r[j + 1:] for r in M[1:]]) for j in range(n) if M[0][j])


def invariant_factors(M):
    """Invariant factors > 1 of Z^cols / rowspan(M) (M must have full column rank), via
    determinantal divisors d_i = gcd of the i x i minors: s_i = d_i / d_(i-1)."""
    rows, cols = len(M), len(M[0])
    d = [1]
    for i in range(1, cols + 1):
        g = 0
        for R in combinations(range(rows), i):
            for C in combinations(range(cols), i):
                g = math.gcd(g, _det([[M[r][c] for c in C] for r in R]))
                if g == d[-1]:
                    break
            if g == d[-1]:
                break
        d.append(g)
    return [s for s in (d[i] // d[i - 1] for i in range(1, cols + 1)) if s != 1]


def relations(Rb, E, k):
    """Basis of Lambda = E * R~^-1 (columns), from a triangular basis of R~, in HNF."""
    inv = [[Fraction(0)] * k for _ in range(k)]  # inverse of the lower-triangular Rb
    for i in range(k):
        inv[i][i] = Fraction(1, Rb[i][i])
        for j in range(i - 1, -1, -1):
            inv[i][j] = -sum(Rb[i][t] * inv[t][j] for t in range(j, i)) / Rb[i][i]
    cols = [[int(E * inv[i][j]) for i in range(k)] for j in range(k)]  # column j of E Rb^-1
    return hnf(cols, E, k)


# ---------------------------------------------------------------- device side

def element_order(oracle, X, n, log):
    """Exact order of X: cheap if X^n = e (strip primes of n), else full order finding."""
    if solver.is_identity(oracle, [(X, n)]):
        N = n
    else:
        N = solver.find_order(oracle, log, elem=X)
    return solver.strip_primes(oracle, N, solver.prime_factors(N), X, strict=True)


def sample(oracle, gens, E, m):
    """One shot: r with chi(gens) = r / E, or None if any estimate looks corrupted."""
    (est,) = solver.run_ipe(oracle, [(X, m) for X in gens])
    out = []
    for e in est:
        r, ok = solver.round_to_n(e, m, E)
        if not ok:
            return None
        out.append(r)
    return tuple(out)


class Verifier:
    def __init__(self, oracle, gens, orders):
        self.oracle, self.gens, self.orders, self.cache = oracle, gens, orders, {}

    def __call__(self, rel):
        terms = [(X, solver.short(x, self.orders[X])) for X, x in zip(self.gens, rel)]
        terms = [(X, x) for X, x in terms if x]
        key = tuple(terms)
        if key not in self.cache:
            self.cache[key] = solver.is_identity(self.oracle, terms)
        return self.cache[key]


def certify(samples, E, k, verify, log):
    """Drop shots that contradict verified relations; return (R basis, Lambda basis) once every
    basis relation verifies, or None if R is still missing characters."""
    while True:
        R = hnf(samples, E, k)
        lam = relations(R, E, k)
        if not all(verify(rel) for rel in lam):
            return None
        for i in range(len(samples)):
            Ri = hnf(samples[:i] + samples[i + 1:], E, k)
            if Ri != R and all(verify(rel) for rel in relations(Ri, E, k)):
                log.append(f"structure: dropped corrupted shot {samples[i]}")
                samples.pop(i)
                break
        else:
            return R, lam


def structure(oracle, log=None, max_shots=60):
    log = [] if log is None else log
    n, logs = solver.order_and_logs(oracle, log)  # n is exact: find_logs strips any excess
    known = {X: x for X, x in logs.items() if x is not None}
    gens = ["g"] + [X for X in "ABC" if X not in known]
    orders = {"g": n}
    for X in gens[1:]:
        orders[X] = element_order(oracle, X, n, log)
    E = math.lcm(*orders.values())
    k = len(gens)
    m = E.bit_length() + solver.LOG_MARGIN
    verify = Verifier(oracle, gens, orders)
    samples, stable, shots, cert = [], 0, 0, None
    if k == 1:
        cert = ((E,),), ((n,),)
    while cert is None and shots < max_shots:
        shots += 1
        s = sample(oracle, gens, E, m)
        if s is None:
            continue
        before = hnf(samples, E, k) if samples else None
        samples.append(s)
        stable = stable + 1 if hnf(samples, E, k) == before else 0
        if stable >= STABLE:
            cert = certify(samples, E, k, verify, log)
            stable = 0
    if cert is None:
        raise RuntimeError("structure not certified")
    _, lam = cert
    # Full relation lattice on (g, A, B, C): add the known logs X = g^x.
    rows = []
    for rel in lam:
        row = [0] * 4
        for X, x in zip(gens, rel):
            row[NAMES.index(X)] = x
        rows.append(row)
    for X, x in known.items():
        row = [0] * 4
        row[0], row[NAMES.index(X)] = -x, 1
        rows.append(row)
    lam4 = hnf(rows, E, 4)
    factors = invariant_factors(lam4)
    quotient = invariant_factors([r[1:] for r in lam4])  # H / <g> = Z^3 / proj(Lambda)
    honest = {}
    for i, X in enumerate(NAMES[1:], 1):
        r = lam4[i]
        if r[i] == 1 and not any(r[1:i]):
            honest[X] = (-r[0]) % n  # row reads X * g^r0 = e
    secret = "none" if len(honest) < 3 else honest["A"] * honest["B"] * honest["C"] % n
    return {
        "order": n,
        "secret": secret,
        "invariant_factors": factors,  # H = Z_f1 x Z_f2 x ...  with f1 | f2 | ...
        "size": math.prod(factors),
        "quotient_by_g": quotient,  # H / <g>: trivial iff no forgery
        "orders": {X: orders[X] if X in orders else n // math.gcd(n, known[X]) for X in NAMES},
        "relations": lam4,  # lower-triangular basis of {x : g^x0 A^x1 B^x2 C^x3 = e}
        "shots": shots,
    }


# ---------------------------------------------------------------- ground truth

def true_structure(inst):
    """Invariant factors of <g, A, B, C> and of H/<g>, computed from the group itself."""
    G, el = inst.group, inst.elements
    hs = [el[X] for X in NAMES]
    E = math.lcm(*(G.element_order(h) for h in hs))
    # Characters chi_j(h) = h_j / d_j generate G^; their restrictions generate H^.
    R = [[E * h[j] // d % E for h in hs] for j, d in enumerate(G.moduli)]
    lam = relations(hnf(R, E, 4), E, 4)
    return invariant_factors(lam), invariant_factors([r[1:] for r in lam])


def brute_structure(inst):
    """|H| by enumeration (small groups only), as an independent check of true_structure."""
    G, el = inst.group, inst.elements
    seen = {G.reduce((0,) * len(G.moduli))}
    frontier = list(seen)
    while frontier:
        nxt = []
        for x in frontier:
            for X in NAMES:
                y = G.add(x, el[X])
                if y not in seen:
                    seen.add(y)
                    nxt.append(y)
        frontier = nxt
    return len(seen)


def _demo():
    import random
    import sys
    import time

    from silentgroup import MockGroup, Oracle, load_public

    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    # Independent check of the lattice algebra against enumeration.
    for _ in range(200):
        mods = [rng.randrange(2, 40) for _ in range(rng.randrange(1, 4))]
        forged = rng.choice(["", "A", "BC", "ABC"]) if len(mods) > 1 else ""
        try:
            inst = MockGroup(mods, forged=forged, seed=rng.randrange(1 << 30))
        except ValueError:
            inst = MockGroup(mods, seed=rng.randrange(1 << 30))
        f, _ = true_structure(inst)
        assert math.prod(f) == brute_structure(inst), (mods, f)
    print("lattice algebra matches enumeration on 200 small groups")

    cases = [
        ("public_1", load_public("public_1")),
        ("public_2", load_public("public_2")),
        ("Z_2^20*3^5 x Z_2^8*3^3, C forged", MockGroup([2**20 * 3**5, 2**8 * 3**3], forged="C", seed=1)),
        ("Z_p*q x Z_p*7, AB forged", MockGroup([1000003 * 999983, 1000003 * 7], forged="AB", seed=1)),
        ("Z_N x Z_5, ABC forged (smooth N)", MockGroup([2171410315074000, 5], forged="ABC", seed=3)),
        ("Z_60 x Z_36 x Z_90, ABC forged", MockGroup([60, 36, 90], forged="ABC", seed=2)),
        ("Z_2^30 x Z_2^20 x Z_2^5, BC forged", MockGroup([2**30, 2**20, 2**5], forged="BC", seed=4)),
    ]
    for name, inst in cases:
        o = Oracle(inst)
        t = time.monotonic()
        log = []
        res = structure(o, log)
        want, want_q = true_structure(inst)
        ok = res["invariant_factors"] == want and res["quotient_by_g"] == want_q
        ok &= res["order"] == inst.answer["order"] and res["secret"] == inst.answer["secret"]
        print(f"{'OK ' if ok else 'BAD'} {name}: H = {' x '.join(f'Z_{f}' for f in res['invariant_factors']) or '1'}"
              f"  |H/<g>| = {math.prod(res['quotient_by_g'])} {res['quotient_by_g']}"
              f"  queries={o.queries_used} shots={res['shots']} {time.monotonic() - t:.1f}s")
        if not ok:
            print("    want", want, want_q, inst.answer, "got", res, log)


if __name__ == "__main__":
    _demo()
