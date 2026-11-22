# Why the solver works, and what it costs

Reproduce: `uv run python analysis/experiments.py && uv run python analysis/plots.py`
(model: `analysis/model.py`, raw data: `analysis/results.json`).

## 1. Cost model

| Quantity | Formula | Why |
|---|---|---|
| m-bit phase estimate | C(m) = m(m+1)/2 queries | step at power 2^k costs bit_length(2^k) = k+1 |
| Floor for one qubit | ~m^2/2 | resolving scale 2^-k needs a power >= 2^(k-1), i.e. >= k queries, 1 bit per measurement |
| Bit flip at step k | q_k = (1 - e^{-2 pi^2 sigma^2 (k+1)})/2 ~ pi^2 sigma^2 (k+1) | kicked phase has N(0, sigma^2 (k+1)) error; P(flip) = E[sin^2(pi eps)] |
| Harmful flips | Poisson(lambda), lambda = sum of q_k over steps whose shift 2^-(k+1) exceeds the decoder tolerance | flips are independent and each shifts the estimate by exactly +-2^-(k+1) |
| Order shot usable | e^-lambda (1 + lambda) | list decoding undoes one flip |

At the competition noise (sigma = 0.002) the 132-bit order shot has lambda = 0.32:
**72% clean, 96% usable after list decoding** (model).

## 2. Model vs simulation: order shot (`shot_vs_noise.png`)

| sigma | clean: sim / model | usable: sim / model |
|---|---|---|
| 0.001 | 0.90 / 0.92 | 0.98 / 1.00 |
| 0.002 | 0.72 / 0.72 | 0.91 / 0.96 |
| 0.003 | 0.47 / 0.48 | 0.77 / 0.83 |
| 0.004 | 0.24 / 0.28 | 0.58 / 0.63 |
| 0.006 | 0.04 / 0.06 | 0.17 / 0.22 |

- The clean rate matches the model across the whole range.
- "Usable" is 5% optimistic. A flip does not always re-converge to an exact +-2^p shift, and two flips close together look like one.
- About 2% of noiseless shots are not clean. That is the inherent phase-estimation tail, P(error > c LSB) ~ 1/(pi^2 c), and list decoding rescues those too.

## 3. Cost vs group size (`queries_vs_size.png`, sigma = 0.002, prime n)

Current solver (probe + lattice order stage, section 4b), regenerated Nov 4. The v1 column is the
original single 132-bit order shot.

| log2 n | measured median | model | v1 measured | floor (no noise, no tests) |
|---|---|---|---|---|
| 12 | 5.3k | 5.5k | 10.5k | 0.6k |
| 36 | 10.2k | 10.3k | 15.3k | 5.3k |
| 40 | 18.2k | 18.3k | — | 6.5k |
| 56 | 24.0k | 23.5k | 21.0k | 12.7k |
| 60 | 25.5k | 25.1k | 22.2k | 14.6k |

- The model (`model.order_cost` + logs + tests) tracks the measurements within ~3% at every size.
- v1 paid a constant ~8.8k for the 132-bit order shot, because it had to assume n could be 2^60. At L = 12 that was 93% of the cost.
- **The trade-off, honestly:** the 84-bit probe halves the cost when n <= 2^36. It costs ~3k more when n > 2^36. Even with the lattice reuse, two shots that together have >= 196 bits cost more than one 132-bit shot, because cost is quadratic in bits.
  - The jump at L = 37 is exactly that.
  - On the two public meetings the net is positive: public_2 +0.05, public_1 -0.01. Random stress improves too.
- The rest of the gap to the floor is verification (~16 shots × L per check) and retries (the long upper error bars).

## 4. Robustness to noise (`solver_vs_noise.png`)

| sigma | public_1 correct | prime60 correct | public_2 correct |
|---|---|---|---|
| 0 to 0.003 | 100% | 100% | 100% |
| 0.004 (2x) | 93% | 93% | 100% |
| 0.005 | 73% | 60% | 100% |
| 0.006 | 27% | 27% | 93% |
| 0.008 | 0% | 0% | 33% |

(Regenerated Nov 4 with the current solver, 15 runs per point. public_2 now holds up better: its order is found by the short probe shot, which carries less noise. The 2^56–2^60 cases drop off earlier at 0.006, because their order needs two shots to both be usable.)

- **The solver is flat up to 1.5x the competition noise.** It breaks down where the model says the order shot stops being usable (lambda ~ 1, around sigma ~ 0.004).
- When it degrades, it usually still gets the order right: wrong answers at 0.005–0.006 are mostly the 0.3-point "order only" case. The fallback that returns the order when the logs fail is what does this.
- public_2 (a forgery, n ~ 2^35) degrades latest because its shots are shortest: errors grow with the power used, ~ sigma^2 log n.
- public_2 costs ~30k even at zero noise. Its forged B lies inside G[n], so B^n = e and only repeated failed rounds prove the forgery (4 rounds of evidence). That cost comes from deciding the forgery, not from noise.

## 4b. Solver changes after this analysis (Nov 4)

| Change | Why | Effect (1000-run sweeps, identical instances) |
|---|---|---|
| Probe order shot at n <= 2^36 first (C(84) = 3.6k), full 2^60 shot only if it fails | the 132-bit shot is 93% of cost for small n | stress 1.393 -> 1.410; public_2 +0.035 |
| Forgery evidence threshold derived from n | `p_clean_wrong(n)`: a flip at power 2^k shifts chi·n by n/2^(k+1), invisible when that is near an integer. Model vs simulation: random n 1.7% / 1.3%, n near 2^k 7.3% / 7.3%, smooth 3.9% / 3.5%. Need (1.5p)^e <= 1e-5 | e = 3 for most n, 4–5 near powers of two |
| Whole-class test (X g^-x0)^(n/M) = e when gcd(j, n) > 8 | smooth n used to give no evidence, ever: Z_N x Z_5 with smooth N and ABC forged cost 60–80k every run | that case 75k -> ~57k |
| Once something looks forged, chase only it | one forgery decides the answer | fewer elements per round on multi-forgery cases |
| Probe miss is reused: 2nd shot at 116 bits + probe's 84 bits -> LLL (3D simultaneous Diophantine approximation, m1 + m2 >= 3·60 + 16), with single-flip list decoding | the probe used to be pure waste for n > 2^36 | order stage for 2^52 < n <= 2^60: 12.3k -> 10.4k (public_1 24.5k -> 22.5k typical) |
| Strict identity test (18 zeros) for enumerated class members and for order strips | a false accept there gives a wrong log / an order that is a proper divisor (0 pts). Seen once per ~3000 runs each | removed both failure modes in later sweeps |
| `OrderWrong` recovery: two corrupted g rounds in a row -> strict g^n = e; if it fails, redo the order | backstop for any upstream false accept | fault injection: a forced false strip is recovered |
| Skip lcm(d, plausible) > 2^60 in the order stage | they were refined anyway, burning the budget (seen once) | — |
| Forgery threshold floor e >= 3 | measured honest bad-round rate at small n matches the model (16 vs 13 predicted over 1500 runs); floor is insurance | rounds at small n are cheap |

Known remaining cost sinks, inherent to per-round phase estimation:
- Smooth n (many small primes) makes gcd(j, n) large, so one round rarely pins a log. An honest cyclic group with n = 2^3·3^2·5·7^4·11^2·13^5 takes a median of 33k queries.
- Forgeries living in a tiny torsion factor (e.g. Z_2 x Z_N) are only visible when gcd(j, n) is coprime to that factor. That's ~half the rounds, so median ~47k.

## 5. Qubits

Every circuit uses **1 qubit**, reset and reused between steps. Mid-circuit measurement plus `if_test` phase corrections turn the inverse QFT into a classically controlled sequence (Griffiths–Niu semiclassical QFT). That is the minimum possible.

## 6. Open track: full structure of H = <g, A, B, C>  (`structure.py`)

**Idea.** Restricting characters G^ -> H^ is a surjective homomorphism, so each shot's chi is a
uniformly random character of H. One joint phase-estimation shot on the generators reads
chi(h) = r/E with r in Z_E^k (E = exponent of H). The r's form a subgroup R of Z_E^k with
R ~ H^ ~ H. Then:

- Relation lattice: Lambda = {x : g^x0 A^x1 B^x2 C^x3 = e} = E · (R + E Z^k)^-1 (the scaled dual lattice).
- H ~ Z^4 / Lambda. Its invariant factors come from the determinantal divisors (gcds of i x i minors).
- H/<g> ~ Z^3 / proj(Lambda) is trivial iff nothing is forged, and its size says how far the forgeries sit outside <g>.

**Certification, no trust in single shots.** Every basis relation of Lambda is checked on the
device with the chi-independent identity test, so missing characters (R too small) can't
pass. A corrupted-but-clean shot makes R too big. Leave-one-out catches it: if dropping
shot i shrinks R and every relation that comes back verifies, shot i contradicted true
relations and is discarded. Elements whose log the solver already verified add only the
relation X = g^x, so an honest meeting costs nothing extra (H = Z_n).

**Results** (competition noise, `uv run python structure.py`, ground truth from the group itself,
cross-checked against brute-force enumeration on 200 small groups):

| Group | Recovered H | H/<g> |
|---|---|---|
| public_1 | Z_66165815680692533 | trivial (honest) |
| public_2 | Z_585585 x Z_37199287125 | Z_585585 |
| Z_60 x Z_36 x Z_90, ABC forged | Z_6 x Z_180 x Z_180 | Z_6 x Z_180 |
| Z_2^30 x Z_2^20 x Z_2^5, BC forged | Z_16 x Z_2^16 x Z_2^30 | Z_16 x Z_2^16 |

Random stress groups (rank 1–3, |G| <= 2^60, random forgeries), all with noise: 1 wrong in the
first 900 runs. The cause was a false accept while stripping an element's order; strict tests
fixed it, and the next **800 runs had 0 wrong**. Median 11k queries, including the solve.
Limit: two corrupted shots that push R out in the same direction would survive leave-one-out.
That's ~p^2 per run and was never seen.
