"""Plots from analysis/results.json: uv run python analysis/plots.py"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from cycler import cycler  # noqa: E402

# Monochrome: C0..C3 map to black and greys; markers/line styles in each plot tell series apart.
plt.rcParams["axes.prop_cycle"] = cycler(color=["#000000", "#555555", "#888888", "#bbbbbb"])

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model  # noqa: E402

R = json.loads((HERE / "results.json").read_text())
COMP = 0.002  # competition noise


def mark_competition(ax):
    ax.axvline(COMP, color="grey", ls=":", lw=1)
    ax.text(COMP, ax.get_ylim()[1], " competition", va="top", fontsize=8, color="grey")


def plot_a():
    rows = R["A"]["rows"]
    s = [r["sigma"] for r in rows]
    fine = [x / 10000 for x in range(0, 81)]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(fine, [model.order_shot(x)["p_clean"] for x in fine], "C0-", label="model: clean  $e^{-\\lambda}$")
    ax.plot(fine, [model.order_shot(x)["p_usable"] for x in fine], "C1-",
            label="model: + list decoding  $e^{-\\lambda}(1+\\lambda)$")
    ax.plot(s, [r["clean"] / r["shots"] for r in rows], "C0o", label="simulated: clean")
    ax.plot(s, [(r["clean"] + r["rescued"]) / r["shots"] for r in rows], "C1s", label="simulated: usable")
    ax.set_xlabel("noise $\\sigma$ (turns per $\\sqrt{\\mathrm{query}}$)")
    ax.set_ylabel("P(order shot gives n)")
    ax.set_title(f"132-bit order shot, n prime ~ 2^{R['A']['n'].bit_length()}")
    ax.legend(fontsize=8)
    mark_competition(ax)
    fig.tight_layout()
    fig.savefig(HERE / "shot_vs_noise.png", dpi=150)


def plot_b():
    rows = R["B"]["rows"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
    for i, case in enumerate(dict.fromkeys(r["case"] for r in rows)):
        rr = [r for r in rows if r["case"] == case]
        s = [r["sigma"] for r in rr]
        style = f"C{i}{'os^'[i % 3]}{['-', '--', ':'][i % 3]}"
        a1.plot(s, [r["success"] for r in rr], style, label=case)
        a2.plot(s, [r["queries"] for r in rr], style, label=case)
    a1.set_ylabel("fraction fully correct")
    a2.set_ylabel("mean queries")
    a2.axhline(80_000, color="black", ls="--", lw=1, label="budget")
    for a in (a1, a2):
        a.set_xlabel("noise $\\sigma$")
        a.legend(fontsize=8)
        mark_competition(a)
    a1.set_title("solve(): correctness vs noise")
    a2.set_title("solve(): cost vs noise")
    fig.tight_layout()
    fig.savefig(HERE / "solver_vs_noise.png", dpi=150)


def plot_c():
    rows = R["C"]["rows"]
    L = [r["L"] for r in rows]
    fine = list(range(L[0], L[-1] + 1))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.errorbar(L, [r["mean"] for r in rows],
                yerr=[[r["mean"] - r["min"] for r in rows], [r["max"] - r["mean"] for r in rows]],
                fmt="ko", capsize=3, label="measured (mean, min-max)")
    ax.plot(fine, [model.expected_queries(x, R["C"]["sigma"])["total"] for x in fine], "C1-", label="model")
    ax.plot(fine, [model.floor_queries(x) for x in fine], "C2--", label="floor: $C(2L) + 4C(L)$, no noise/tests")
    ax.set_xlabel("$L = \\log_2 n$ (prime order)")
    ax.set_ylabel("queries")
    ax.set_title(f"cost vs group size at $\\sigma$ = {R['C']['sigma']}")
    sec = ax.secondary_yaxis("right", functions=(model.score, lambda s: 80_000 * (3 - 2 * s)))
    sec.set_ylabel("score per meeting")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "queries_vs_size.png", dpi=150)


if __name__ == "__main__":
    for name, f in (("A", plot_a), ("B", plot_b), ("C", plot_c)):
        if name in R:
            f()
            print("plotted", name)
