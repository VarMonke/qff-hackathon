"""Analytic cost / noise model of the solver.

Derivation (all from the device rules):

* One semiclassical phase-estimation step at power 2^k costs bit_length(2^k) = k + 1
  queries, so an m-bit estimate costs  C(m) = sum_{k<m} (k+1) = m (m+1) / 2.
  Resolving a phase to 2^-m needs some power >= 2^(m-1), so each bit at scale 2^-k costs
  >= k queries: m^2/2 is also (up to constants) the floor for a 1-qubit measurement.

* The step at power 2^k kicks back phase + eps, eps ~ N(0, s^2), s = sigma sqrt(k+1) turns.
  The measured bit flips with probability E[sin^2(pi eps)] = (1 - exp(-2 pi^2 s^2)) / 2:
      q_k = (1 - exp(-2 pi^2 sigma^2 (k+1))) / 2   ~=  pi^2 sigma^2 (k+1).
  So errors are concentrated on the expensive, high-power steps.

* A flip at step k shifts the estimate by +-2^-(k+1) (it re-converges afterwards). It only
  matters if that shift exceeds the tolerance of the decoder:
    - order shot (m = 132, fit tolerance 8 LSB = 2^-(m-3)): steps k <= m - 5 matter;
    - log shot   (round to the nearest multiple of 1/n): shift n 2^-(k+1) >= 1/2, i.e.
      steps k <= log2 n matter.
  Flips are independent, so the number of harmful flips is ~ Poisson(lambda), lambda = sum q_k.

* List decoding (undo one +-2^p shift) rescues shots with exactly one harmful flip:
      P(order shot usable) = e^-lambda (1 + lambda).
"""

import math

NMAX_BITS = 60
ORDER_BITS = 2 * NMAX_BITS + 12  # 132
TEST_SHOTS = 16  # an accepted identity test (16 straight zeros)
REJECT_SHOTS = 6  # typical shots until a non-identity test sees 3 ones


def ipe_cost(m: int) -> int:
    return m * (m + 1) // 2


def flip_prob(sigma: float, k: int) -> float:
    return (1 - math.exp(-2 * math.pi**2 * sigma**2 * (k + 1))) / 2


def lam(sigma: float, k_max: int) -> float:
    """Expected number of harmful flips over steps 0..k_max."""
    return sum(flip_prob(sigma, k) for k in range(k_max + 1))


def order_shot(sigma: float, m: int = ORDER_BITS) -> dict:
    l = lam(sigma, m - 5)
    return {"lambda": l, "p_clean": math.exp(-l), "p_usable": math.exp(-l) * (1 + l)}


PROBE_BITS, COMBO_BITS = 36, 52


def order_cost(L: int, sigma: float) -> float:
    """Order stage of the current solver: an 84-bit probe (n <= 2^36), then a 116-bit shot that
    alone covers n <= 2^52 and together with the probe (lattice) n <= 2^60, then 132-bit shots."""
    m1, m2 = 2 * PROBE_BITS + 12, 2 * COMBO_BITS + 12
    p1, p2 = order_shot(sigma, m1)["p_usable"], order_shot(sigma, m2)["p_usable"]
    p3 = order_shot(sigma)["p_usable"]
    tail = ipe_cost(ORDER_BITS) / p3
    if L <= PROBE_BITS:
        return ipe_cost(m1) + (1 - p1) * (ipe_cost(m2) + (1 - p2) * tail)
    # L > 52 needs both shots usable (about p1 * p2 with single-flip decoding of the pair).
    ok = p2 if L <= COMBO_BITS else p1 * p2
    return ipe_cost(m1) + ipe_cost(m2) + (1 - ok) * tail


def log_estimate_ok(sigma: float, L: int) -> float:
    return math.exp(-lam(sigma, L))


def expected_queries(L: int, sigma: float) -> dict:
    """Expected solver queries for an honest instance with prime order n ~ 2^L.

    Prime n avoids the gcd(j, n) > 1 refinements, so every term is modelled exactly up to
    the retry approximation. Logs: after a round, X is resolved iff both its own and g's
    estimate survived (prob p^2); unresolved X are redone together with g.
    """
    small_probe = REJECT_SHOTS * 33
    order = order_cost(L, sigma) + TEST_SHOTS * L
    m = L + 8
    C = ipe_cost(m)
    p = log_estimate_ok(sigma, L)
    p2 = p * p
    # Round 1 measures 4 elements; each failed X costs geometric extra rounds of (g + X),
    # and failed rounds share their g (counted once per extra round as an approximation).
    extra_rounds = (1 - p2) / p2
    logs = 4 * C + extra_rounds * 2 * C
    tests = 3 * TEST_SHOTS * (L + 1) + 3 * (1 - p2) / p2 * REJECT_SHOTS * (L + 1)
    total = small_probe + order + logs + tests
    return {"order": order, "logs": logs, "tests": tests + small_probe, "total": total}


def floor_queries(L: int) -> int:
    """No noise, no margins, no verification: C(2L) for the order plus 4 C(L) for logs."""
    return ipe_cost(2 * L) + 4 * ipe_cost(L)


def score(queries: float, budget: int = 80_000) -> float:
    return 1 + 0.5 * (budget - queries) / budget


if __name__ == "__main__":
    s = 0.002
    print("order shot at sigma=0.002:", order_shot(s))
    for L in (20, 35, 56, 60):
        e = expected_queries(L, s)
        print(f"L={L}: model {e['total']:.0f} (score {score(e['total']):.3f}), floor {floor_queries(L)}"
              f" (score {score(floor_queries(L)):.3f})")
