import math
from fractions import Fraction
from itertools import combinations

from qiskit import QuantumCircuit

from silentgroup import oracle_gate

NMAX_BITS = 60
PROBE_BITS = 36
COMBO_BITS = 52
EVIDENCE_NEEDED = 0
LOG_MARGIN = 8
CAP = 12  # semiclassical qft


def ipe_circuit(specs, mult=1):
    """
    specs: list of (element, m bits); phases of element^mult

    Returns circuit, clbit offsets
    """
    total = sum(m for _, m in specs)
    qc = QuantumCircuit(1, total)
    offsets, c0 = [], 0
    for elem, m in specs:
        offsets.append(c0)
        for k in range(m - 1, -1, -1):
            qc.reset(0)
            qc.h(0)
            qc.append(oracle_gate(elem, mult << k), [0])
            for l in range(1, min(CAP, m - 1 - k) + 1):
                with qc.if_test((qc.clbits[c0 + k + l], 1)):
                    qc.p(-2 * math.pi / (1 << (l + 1)), 0)
            qc.h(0)
            qc.measure(0, c0 + k)
        c0 += m
    return qc, offsets


def run_ipe(oracle, specs, shots=1, mult=1):
    """
    Per shot, one integer estimate M per spec, with chi(elem^mult) ~ M / 2^m
    """
    qc, offsets = ipe_circuit(specs, mult)
    counts = oracle.run(qc, shots=shots)
    out = []
    for key, cnt in counts.items():
        bits = key[::-1]  # bits[c] = clbit c
        est = []
        for (elem, m), c0 in zip(specs, offsets):
            # clbit c0+k holds binary digit k+1 of chi, which has weight 2^(m-1-k) in M
            est.append(sum(int(bits[c0 + k]) << (m - 1 - k) for k in range(m)))
        out.extend([est] * cnt)
    return out


def is_identity(oracle, terms, max_shots=24, reject_at=3, batch=4, fast_accept=16, strict=False):
    """
    chi-independent test that prod h^k == e. terms: list of (elem, power).

    Non-identity gives outcome 1 with probability exactly 1/2 per shot (chi(t) is uniform
    over multiples of 1/ord(t)); the identity only through noise (~0.2%). Sequential:
    accept after 16 straight zeros (false accept 2^-16), otherwise after 24 shots with
    <= 2 ones (301/2^24); reject at 3 ones (false reject ~2e-5). strict (for candidates that
    are likely wrong, where false accepts add up): 18 zeros or <= 2 ones in 28, ~5e-6.
    """
    if strict:
        max_shots, fast_accept = 28, 18
    terms = [(e, k) for e, k in terms if k != 0]
    if not terms:
        return True
    qc = QuantumCircuit(1, 1)
    qc.h(0)
    for e, k in terms:
        qc.append(oracle_gate(e, k), [0])
    qc.h(0)
    qc.measure(0, 0)
    ones = used = 0
    while used < max_shots:
        counts = oracle.run(qc, shots=batch)
        ones += counts.get("1", 0)
        used += batch
        if ones >= reject_at:
            return False
        if ones == 0 and used >= fast_accept:
            return True
    return True


def crt(eqs):
    """
    Solve x = r_i (mod m_i).
    Returns (x0, M) or None if inconsistent
    """
    x0, M = 0, 1
    for r, m in eqs:
        g = math.gcd(M, m)
        if (r - x0) % g:
            return None
        l = M // g * m
        t = ((r - x0) // g) * pow(M // g, -1, m // g) % (m // g) if m // g > 1 else 0
        x0, M = (x0 + M * t) % l, l
    return x0, M


def solve_lin(j, r, n):
    """
    x * j = r (mod n) -> (x0, m) with x = x0 (mod m), or None.
    """
    g = math.gcd(j, n)
    if r % g:
        return None
    m = n // g
    return ((r // g) * pow(j // g, -1, m) % m if m > 1 else 0), m


def _is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _rho(n):
    if n % 2 == 0:
        return 2
    for c in range(1, 100):
        x = y = 2
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = math.gcd(abs(x - y), n)
        if d != n:
            return d
    raise ValueError("rho failed")


def prime_factors(n):
    """
    Distinct prime factors of n.
    """
    out, stack = set(), [n]
    while stack:
        x = stack.pop()
        if x == 1:
            continue
        if _is_prime(x):
            out.add(x)
            continue
        d = _rho(x)
        stack += [d, x // d]
    return sorted(out)


def short(x, n):
    """
    Representative of x mod n with the smallest absolute value (cheaper oracle power).
    """
    x %= n
    return x - n if x > n // 2 else x


def simplest_between(lo, hi):
    """
    The fraction with the smallest denominator in [lo, hi] (Stern-Brocot descent)
    """
    a = math.ceil(lo)
    if a <= hi:
        return Fraction(a)
    f = math.floor(lo)
    return f + 1 / simplest_between(1 / (hi - f), 1 / (lo - f))


def _fit(est, m, bound):
    """
    Smallest d <= 2^bound with est / 2^m within 8 LSB of some j/d, else None.
    """
    theta, tol = Fraction(est, 1 << m), Fraction(8, 1 << m)
    d = simplest_between(theta - tol, theta + tol).denominator
    return d if d.bit_length() <= bound else None


def estimate_denominators(oracle, bound, mult=1, elem="g"):
    """
    One shot of phase estimation on elem^mult. Returns (clean, decoded): the denominator d
    of chi(elem^mult) = j/d if the raw estimate fits one, and alternatives from list decoding.

    With m = 2 bound + 12 bits, fractions with denominator <= 2^bound cover only ~0.1% of
    [0, 1) at 8-LSB resolution, so a correct estimate fits and a corrupted one almost never
    does. Noise mostly flips a single semiclassical step, which shifts the estimate by
    exactly +-2^p LSB, so we also try undoing each such shift.
    """
    m = 2 * bound + 12
    (est,) = run_ipe(oracle, [(elem, m)], mult=mult)[0]
    estimate_denominators.raw = (est, m)
    return decode(est, m, bound)


def decode(est, m, bound):
    """
    (clean fit, sorted single-flip decodings smaller than it) for an m-bit estimate
    """
    clean = _fit(est, m, bound)
    decoded = set()
    for p in range(m):
        for sgn in (1, -1):
            d = _fit((est + sgn * (1 << p)) % (1 << m), m, bound)
            if d is not None and (clean is None or d < clean):
                decoded.add(d)
    return clean, sorted(decoded)


def _lll(B, delta=Fraction(3, 4)):
    """
    LLL-reduce a small list of integer row vectors (exact arithmetic)
    """
    B = [list(b) for b in B]
    n = len(B)

    def gso():
        Bs, mu = [], [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = [Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = sum(Fraction(a) * b for a, b in zip(B[i], Bs[j])) / sum(b * b for b in Bs[j])
                v = [a - mu[i][j] * b for a, b in zip(v, Bs[j])]
            Bs.append(v)
        return Bs, mu

    Bs, mu = gso()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [a - q * b for a, b in zip(B[k], B[j])]
                for t in range(j):
                    mu[k][t] -= q * mu[j][t]
                mu[k][j] -= q
        nk, nk1 = sum(x * x for x in Bs[k]), sum(x * x for x in Bs[k - 1])
        if nk >= (delta - mu[k][k - 1] ** 2) * nk1:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            Bs, mu = gso()
            k = max(k - 1, 1)
    return B


def combine(e1, m1, e2, m2, bound=NMAX_BITS):
    """
    Smallest q <= 2^bound with q * e_i / 2^m_i within 8q LSB of an integer for both shots
    (simultaneous Diophantine approximation, so needs m1 + m2 >~ 3 bound + 16), or None.
    Short vector of the lattice spanned by (8, e1, e2), (0, 2^m1, 0), (0, 0, 2^m2)
    """
    B = _lll([[8, e1, e2], [0, 1 << m1, 0], [0, 0, 1 << m2]])
    best = None
    for a in range(-2, 3):
        for b in range(-2, 3):
            for c in range(-1, 2):
                v = [a * x + b * y + c * z for x, y, z in zip(*B)]
                if v[0] < 0:
                    v = [-x for x in v]
                q = v[0] // 8
                if 0 < q and q.bit_length() <= bound and abs(v[1]) <= 8 * q and abs(v[2]) <= 8 * q:
                    best = q if best is None else min(best, q)
    return best


def combine_decoded(e1, m1, e2, m2):
    """
    combine() on the raw pair, else on single +-2^p flips of either estimate, most likely
    first (a flip at weight 2^p comes from the step with power 2^(m-1-p), whose flip
    probability grows with that step's cost), stopping at the first hit
    """
    q = combine(e1, m1, e2, m2)
    if q:
        return [q]
    tries = [(m - 1 - p, which, p) for which, m in ((0, m1), (1, m2)) for p in range(m)]
    for _, which, p in sorted(tries, reverse=True):
        for sgn in (1, -1):
            if which == 0:
                q = combine((e1 + sgn * (1 << p)) % (1 << m1), m1, e2, m2)
            else:
                q = combine(e1, m1, (e2 + sgn * (1 << p)) % (1 << m2), m2)
            if q:
                return [q]
    return []


SMALL_L = math.lcm(*range(1, 25))  # 5354228880, 33 bits


def strip_primes(oracle, N, primes, elem="g", strict=False):
    """
    Given elem^N = e, divide out every p^t (p in primes) with elem^(N / p^t) still = e
    """
    for p in primes:
        e, x = 0, N
        while x % p == 0:
            x, e = x // p, e + 1
        lo, hi = 0, e  # largest t with g^(N / p^t) = e; monotone in t, so binary search
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if is_identity(oracle, [(elem, N // p**mid)], strict=strict):
                lo = mid
            else:
                hi = mid - 1
        N //= p**lo
    return N


def find_order(oracle, log, elem="g"):
    """
    Order of elem (possibly a multiple of it) via phase estimation + continued fractions
    + verification.
    """
    # Tiny orders (all prime powers <= 24) cost a few hundred queries this way,
    # and they are where list decoding produces spurious multiples of n
    if is_identity(oracle, [(elem, SMALL_L)]):
        n = strip_primes(oracle, SMALL_L, prime_factors(SMALL_L), elem)
        log.append(f"order: {n} (divides lcm(1..24))")
        return n
    rejections, plausible = {}, []  # plausible: values believed to divide n

    def confirm(c):
        # A noise-induced false rejection must not ban the true order forever: allow a
        # retest when an independent shot proposes the same value again.
        if c.bit_length() > NMAX_BITS or rejections.get(c, 0) >= 2:
            return False
        if is_identity(oracle, [(elem, c)]):
            return True
        rejections[c] = rejections.get(c, 0) + 1
        return False

    def refine(c):
        """If c | n then g^c has small order; a cheap estimate on g^c recovers it.
        Returns the order (or None), and whether c looked like a divisor of n."""
        # Each pass sees a random character of <g^c>, so it may recover only part of the
        # cofactor, or nothing (chi(g^c) = 0 has probability 1/s): keep sampling.
        for _ in range(4):
            s, _ = estimate_denominators(oracle, 8, mult=c, elem=elem)
            if s is None:
                return None, False
            if s == 1:
                continue
            c *= s
            if confirm(c):
                return c, True
        return None, True

    bounds = [PROBE_BITS, COMBO_BITS] if PROBE_BITS else []
    bounds = [b for b in bounds if b]
    raw = []
    for shot in range(7):
        bound = bounds[shot] if shot < len(bounds) else NMAX_BITS
        clean, decoded = estimate_denominators(oracle, bound, elem=elem)
        raw.append(estimate_denominators.raw)
        trusted = [clean] if clean else decoded if len(decoded) == 1 else []
        if shot == 1 and len(bounds) == 2 and not clean:
            # n in (2^52, 2^60]: neither shot fits alone, but together they pin n down.
            trusted = list(dict.fromkeys(combine_decoded(*raw[0], *raw[1])[:3] + trusted))
        found = None
        # Decodings smaller than the clean fit go first: for small n a bit flip can shift
        # j/n onto a fraction with denominator n * 2^t that also fits cleanly.
        for d in [d for d in decoded if d not in trusted][:6]:
            if confirm(d):  # may be a multiple of n; find_logs strips any excess
                found = d
                break
        for d in trusted if not found else []:
            for c in dict.fromkeys([math.lcm(d, p) for p in plausible] + [d]):
                if c.bit_length() > NMAX_BITS:
                    continue  # not an order: |G| <= 2^60
                if confirm(c):
                    found = c
                    break
                found, ok = refine(c)
                if found:
                    break
                if ok and c not in plausible:
                    plausible.append(c)
            if found:
                break
        if found:
            log.append(f"order: {found} (shot {shot}, clean={clean is not None}, {len(decoded)} decodings)")
            return found
        log.append(f"order: clean={clean} decoded={decoded[:6]} not confirmed; plausible={plausible}")
    raise RuntimeError("order not found")


def round_to_n(e, m, n):
    """
    Nearest r with chi ~ r/n, and whether the estimate e / 2^m actually sits near r/n.

    m = log n + 8 bits puts a good estimate within 1/32 of r/n. A bit flip shifts it by
    n 2^(p-m), which is usually far from an integer so the 1/16 threshold flags it, but not
    when n is close to a power of two; so "clean" is a hint, not a proof. Every element of
    <g> (or of G[n]) has chi in (1/n)Z, so a forgery inside G[n] still looks clean here
    """
    x = Fraction(e * n, 1 << m)
    r = round(x)
    return r % n, abs(x - r) <= Fraction(1, 16)


def strip_order(oracle, N, j, checked, log):
    """
    N is a confirmed multiple of ord(g) and chi(g) = j/N exactly. If N = k * n then
    k | gcd(j, N), so only primes of gcd(j, N) can be excess; test those (usually none)
    """
    primes = [p for p in prime_factors(math.gcd(j, N)) if p not in checked]
    checked.update(primes)
    # Strict: a false strip leaves n a divisor of ord(g), and logs can still verify mod it.
    # Only accepting tests (real excess, rare) pay the extra shots.
    N2 = strip_primes(oracle, N, primes, strict=True)
    if N2 != N:
        log.append(f"order: stripped {N // N2} -> {N2}")
    return N2


class OrderWrong(RuntimeError):
    """
    g^n != e after all: an identity test upstream accepted falsely
    """


def p_clean_wrong(n, m, sigma=0.002):
    """
    Model: P(an m-bit estimate of chi in (1/n)Z passes round_to_n's check but is wrong).
    A flip at the step with power 2^k shifts chi * n by n / 2^(k+1); it goes unnoticed when
    that is within 1/16 of a nonzero integer (n near a power of two is the bad case)
    """
    tot = 0.0
    for k in range(m):
        s = n / 2 ** (k + 1)
        if s < 0.5:
            break
        if min(s % 1, 1 - s % 1) <= 1 / 16:
            tot += (1 - math.exp(-2 * math.pi**2 * sigma**2 * (k + 1))) / 2
    return tot


def evidence_needed_for(n, m):
    """
    Inconsistent clean rounds needed to call X forged: an honest X needs that many bad
    rounds (each needing a clean-but-wrong estimate of g or X) before its first good one.
    1.5x margin over the model, which matches simulation within around 20%.
    """
    p = 1.5 * 2 * p_clean_wrong(n, m)
    e = 3  # never 2: one unexplained honest double-failure seen at small n (~1 in 3000 runs)
    while e < 6 and p**e > 1e-5:
        e += 1
    return e


def find_logs(oracle, n, log, max_rounds=8, evidence_needed=None, checked=None):
    """
    Discrete logs of A, B, C base g mod n. Value None means forged.
    Returns (n, logs): n may shrink if the order stage handed us a multiple of ord(g)
    """
    checked = set() if checked is None else checked
    find_logs.current_n = n
    names = ["A", "B", "C"]
    m = n.bit_length() + LOG_MARGIN
    evidence_needed = evidence_needed or EVIDENCE_NEEDED or evidence_needed_for(n, m)
    eqs = {X: [] for X in names}
    result = {}
    tested = {X: set() for X in names}
    evidence = dict.fromkeys(names, 0)
    checked_pow = set()
    corrupted = 0
    for rnd in range(max_rounds):
        pending = [X for X in names if X not in result]
        if not pending:
            break
        # One forgery decides the answer: once something looks forged, chase only it.
        suspect = max(pending, key=lambda X: evidence[X])
        if evidence[suspect]:
            pending = [suspect]
        (est,) = run_ipe(oracle, [("g", m)] + [(X, m) for X in pending])
        j, ok = round_to_n(est[0], m, n)
        if not ok:
            log.append(f"logs: round {rnd} g estimate corrupted")
            corrupted += 1
            # With the right n this is around 1% likely; with a wrong n (a false accept upstream)
            # most rounds look corrupted.
            if corrupted == 2 and not is_identity(oracle, [("g", n)], strict=True):
                raise OrderWrong(n)
            continue
        corrupted = 0
        n2 = strip_order(oracle, n, j, checked, log)
        if n2 != n:
            return find_logs(oracle, n2, log, max_rounds, evidence_needed, checked)
        fresh = {}
        for X, e in zip(pending, est[1:]):
            r, ok = round_to_n(e, m, n)
            if not ok:
                continue
            sol = solve_lin(j, r, n)
            fresh[X] = sol
            if sol is None:
                evidence[X] += 1  # clean data, but r is not in j * Z_n: X is not in <g>
            else:
                eqs[X].append(sol)
        # (A joint test of prod X * g^-(sum x) would be cheaper, but two forgeries can cancel
        # in G/<g>, and ambiguous candidates off by multiples of n/gcd cancel too: unsound.)
        for X in pending:
            # Every wrong candidate is another chance of a false accept (and a forgery
            # only has wrong ones), so test only this round's class plus a few of the
            # best-supported combinations.
            own = []
            if fresh.get(X) and n // fresh[X][1] <= 8:
                x0, M = fresh[X]
                own = [(x0 + t * M) % n for t in range(n // M)]
            queue = list(dict.fromkeys(own + candidates(eqs[X], n)[:4]))
            unique = fresh[X][0] if fresh.get(X) and fresh[X][1] == n else None
            for cand in [c for c in queue if c not in tested[X]]:
                tested[X].add(cand)
                if is_identity(oracle, [(X, 1), ("g", -short(cand, n))], strict=cand != unique):
                    result[X] = cand
                    log.append(f"log {X}: {cand} (round {rnd})")
                    break
            if X not in result and fresh.get(X):
                # Clean round whose whole candidate class x0 + t * M failed: evidence
                # Large classes are tested at once: X in x0 + M Z  =>  (X g^-x0)^(n/M) = e
                x0, M = fresh[X]
                d = n // M
                if d <= 8:
                    bad = all((x0 + t * M) % n in tested[X] for t in range(d))
                else:
                    bad = not is_identity(oracle, [(X, d), ("g", -short(x0 * d, n))])
                if bad:
                    evidence[X] += 1
                    eqs[X].remove(fresh[X])
        # A forgery anywhere already decides the secret. X^n != e proves one outright
        for X in pending:
            if X in result:
                continue
            if X not in checked_pow:
                checked_pow.add(X)
                if not is_identity(oracle, [(X, n)]):
                    log.append(f"{X}^n != e: forged")
                    return n, {**result, X: None}
            if evidence[X] >= evidence_needed:
                log.append(f"{X}: {evidence[X]} clean rounds inconsistent with <g>: forged")
                return n, {**result, X: None}
    for X in names:
        if X not in result:
            log.append(f"log {X}: never verified, declaring forged")
            result[X] = None
    return n, result


def candidates(eqs, n, max_expand=8):
    """Candidate logs mod n from per-shot congruences, most supported first."""
    votes = {}
    for size in range(1, min(3, len(eqs)) + 1):
        for sub in combinations(eqs, size):
            sol = crt(sub)
            if sol is None:
                continue
            x0, M = sol
            if n % M or n // M > max_expand:
                continue
            for t in range(n // M):
                c = (x0 + t * M) % n
                votes[c] = votes.get(c, 0) + 1
    return sorted(votes, key=lambda c: -votes[c])


def order_and_logs(oracle, log, attempts=3):
    """find_order + find_logs, starting over if the logs stage proves the order wrong."""
    for _ in range(attempts):
        n = find_logs.current_n = find_order(oracle, log)
        try:
            return find_logs(oracle, n, log)
        except OrderWrong:
            log.append(f"order {n} disproved (g^n != e); starting over")
    raise RuntimeError("order keeps failing")


def solve(oracle) -> dict:
    log = solve.log = []
    find_logs.current_n = 1
    try:
        n, logs = order_and_logs(oracle, log)
    except Exception as exc:
        # Out of budget or shots. The order alone is worth 0.3, and "none" is right whenever
        # a forgery is why we could not finish.
        log.append(f"failed: {exc!r}")
        return {"order": find_logs.current_n, "secret": "none"}
    if any(v is None for v in logs.values()):
        return {"order": n, "secret": "none"}
    return {"order": n, "secret": logs["A"] * logs["B"] * logs["C"] % n}
