"""
Scratch-start comparison: build a network from zero stations.

Sequential methods (one station at a time, same as compare_methods.py):
  LEA    - Largest Empty Area. Needs a starting rule because Voronoi vertices
           do not exist yet:
             0 stations -> pole of inaccessibility (point farthest from the boundary)
             1 station  -> candidate centers are the boundary vertices; circles
                           clipped and ranked by clipped area, as in the method
             2+         -> the published method, unchanged
  BUF-r  - Station buffer, unchanged (first pick = centroid of the whole region).

All-at-once reference:
  CVT    - Centroidal Voronoi layout of exactly N stations (k-means on the
           region's grid). Near-hexagonal in the interior and adapted to the
           boundary: the layout you would choose if all N stations were
           installed at once and build order did not matter.

Each sequential network is scored after N stations, for several N, against a
CVT layout with the same N.
"""
import csv, json, os, sys, time
import numpy as np
from shapely.geometry import Point
from shapely.ops import polylabel
from scipy.cluster.vq import kmeans2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compare_methods import (Domain, load_kansas, random_polygon, pick_lea,
                             pick_buffer_centroid, metrics)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results_scratch")
os.makedirs(OUT, exist_ok=True)
BUFFER_RADII = (25, 50)


# ----------------------------------------------------------------------------
# Scratch-start pickers
# ----------------------------------------------------------------------------
def pick_lea_scratch(dom, st, r=None):
    if len(st) == 0:
        p = polylabel(dom.poly, tolerance=dom.res)
        return np.array([p.x, p.y])
    if len(st) == 1:
        C = np.array(dom.poly.exterior.coords)[:-1]
        R = np.hypot(*(C - st[0]).T)
        circ = [dom.poly.intersection(Point(c).buffer(rr, 32)) for c, rr in zip(C, R)]
        g = circ[int(np.argmax([c.area for c in circ]))].centroid
        return np.array([g.x, g.y])
    return pick_lea(dom, st)


def pick_buffer_scratch(dom, st, r):
    if len(st) == 0:
        cx, cy = dom.pts.mean(0)          # the whole region is the gap
        if dom.poly.contains(Point(cx, cy)):
            return np.array([cx, cy])
        return dom.pts[np.argmin(np.hypot(dom.pts[:, 0] - cx, dom.pts[:, 1] - cy))]
    return pick_buffer_centroid(dom, st, r)


def build(dom, picker, r, n_max):
    st = np.zeros((0, 2))
    for _ in range(n_max):
        p = picker(dom, st, r)
        if p is None:
            break
        st = np.vstack([st, p])
    return st


def cvt_layout(dom, n, seed=0, restarts=4):
    """Best of a few k-means runs on the grid points (minimises mean squared distance)."""
    best, best_cost = None, np.inf
    for k in range(restarts):
        c, lab = kmeans2(dom.pts, n, minit="++", iter=50, seed=seed + k)
        cost = ((dom.pts - c[lab]) ** 2).sum()
        if cost < best_cost:
            best, best_cost = c, cost
    return best


# ----------------------------------------------------------------------------
# 1) Kansas from scratch
# ----------------------------------------------------------------------------
def kansas_scratch(checkpoints=(5, 10, 20, 40)):
    poly, *_ = load_kansas()
    dom = Domain(poly, res_km=2.0)
    n_max = max(checkpoints)
    seqs = {"Largest Empty Area": build(dom, pick_lea_scratch, None, n_max)}
    for r in BUFFER_RADII:
        seqs[f"Buffer r={r} km"] = build(dom, pick_buffer_scratch, r, n_max)
    for k, v in seqs.items():
        print(f"  Kansas {k:22s} placed {len(v)}")

    curves = {k: {} for k in list(seqs) + ["All-at-once (CVT)"]}
    cvts = {}
    for n in range(1, n_max + 1):
        for k, st in seqs.items():
            if len(st) >= n:
                curves[k][n] = metrics(dom, st[:n], radii=(50,), Ls=(25, 50, 100))
        if n in checkpoints or n % 5 == 0:
            cvts[n] = cvt_layout(dom, n)
            curves["All-at-once (CVT)"][n] = metrics(dom, cvts[n], radii=(50,), Ls=(25, 50, 100))

    # maps at N = 10 and N = 40
    shows = [10, 40]
    names = list(seqs) + ["All-at-once (CVT)"]
    fig, axs = plt.subplots(len(names), len(shows), figsize=(11, 3.1 * len(names)), constrained_layout=True)
    bx, by = dom.poly.exterior.xy
    for i, name in enumerate(names):
        for j, n in enumerate(shows):
            ax = axs[i, j]
            st = cvts[n] if name.startswith("All") else seqs[name][:n]
            ax.plot(bx, by, "k-", lw=1)
            ax.plot(st[:, 0], st[:, 1], "o", ms=9 if n <= 10 else 6, mfc="#1d4ed8", mec="w")
            if not name.startswith("All"):
                for q, p in enumerate(st):
                    ax.text(p[0], p[1], str(q + 1), color="w", ha="center", va="center",
                            fontsize=6 if n <= 10 else 4.5, weight="bold")
            m = curves[name].get(n) if len(st) == n else None
            sub = (f"largest gap {m['max_d']:.0f} km, mean dist {m['mean_d']:.0f} km" if m
                   else f"stopped at {len(st)} stations (buffers cover everything)")
            ax.set_title(f"{name}, N = {n}\n{sub}", fontsize=8.5)
            ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Kansas from scratch (numbers = build order)", fontsize=10)
    fig.savefig(os.path.join(OUT, "kansas_scratch_maps.png"), dpi=130)
    plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(14, 3.8), constrained_layout=True)
    for name in names:
        ns = sorted(curves[name]); c = curves[name]
        style = dict(ls="none", marker="D", ms=5, color="k") if name.startswith("All") else dict(lw=1.6)
        for ax, key in zip(axs, ["max_d", "mean_d", "ksd50"]):
            ax.plot(ns, [c[n][key] for n in ns], label=name, **style)
    for ax, t in zip(axs, ["Largest gap (km)", "Mean distance to nearest station (km)", "Mean kriging SD, L = 50 km"]):
        ax.set_title(t, fontsize=9); ax.set_xlabel("number of stations"); ax.grid(alpha=.3)
    axs[0].legend(fontsize=8)
    fig.savefig(os.path.join(OUT, "kansas_scratch_curves.png"), dpi=130)
    plt.close(fig)

    return {name: {n: {k: round(float(v), 3) for k, v in curves[name][n].items()}
                   for n in checkpoints if n in curves[name]} for name in names} | \
           {"stopped_at": {k: len(v) for k, v in seqs.items()}}


# ----------------------------------------------------------------------------
# 2) Monte Carlo from scratch
# ----------------------------------------------------------------------------
def monte_carlo_scratch(n_trials=150, checkpoints=(5, 10, 20, 30), seed=2):
    rng = np.random.default_rng(seed)
    regions = [("Kansas-like", 38.5, -98.3), ("Sahel", 14, 2), ("Scandinavia", 64, 16),
               ("Patagonia", -45, -69), ("SE Asia", 15, 103), ("Tropics", -3, 25)]
    rows, t0 = [], time.time()
    n_max = max(checkpoints)
    for t in range(n_trials):
        reg, la, lo = regions[t % len(regions)]
        poly = random_polygon(rng, la + rng.uniform(-3, 3), lo + rng.uniform(-3, 3), rng.uniform(150, 400))
        if poly.geom_type != "Polygon":
            poly = max(poly.geoms, key=lambda g: g.area)
        dom = Domain(poly, res_km=4.0)
        nets = {"LEA": build(dom, pick_lea_scratch, None, n_max)}
        for r in BUFFER_RADII:
            nets[f"BUF{r}"] = build(dom, pick_buffer_scratch, r, n_max)
        rec = {"trial": t, "region": reg, "area": dom.area}
        for n in checkpoints:
            cvt = cvt_layout(dom, n, seed=t)
            for k, v in metrics(dom, cvt, radii=(50,), Ls=(25, 100)).items():
                rec[f"CVT_{n}_{k}"] = v
            for name, st in nets.items():
                ok = len(st) >= n
                rec[f"{name}_{n}_ok"] = int(ok)
                if ok:
                    for k, v in metrics(dom, st[:n], radii=(50,), Ls=(25, 100)).items():
                        rec[f"{name}_{n}_{k}"] = v
        rows.append(rec)
        if t % 15 == 0:
            print(f"  trial {t}/{n_trials}  {reg:12s} area {dom.area:8.0f}  ({time.time()-t0:.0f}s)")
    return rows


def summarize(rows, checkpoints=(5, 10, 20, 30)):
    """% above the all-at-once CVT layout (lower is better), and LEA vs buffer wins."""
    keys = ["max_d", "mean_d", "ksd25", "ksd100"]
    out = {"n_trials": len(rows)}
    for n in checkpoints:
        res = {}
        for name in ["LEA"] + [f"BUF{r}" for r in BUFFER_RADII]:
            sub = [r for r in rows if r[f"{name}_{n}_ok"]]
            d = {"no_answer_%": round(100 * (1 - len(sub) / len(rows)), 1)}
            for k in keys:
                x = np.array([r[f"{name}_{n}_{k}"] / r[f"CVT_{n}_{k}"] - 1 for r in sub]) * 100
                d[f"{k}_median_%_above_CVT"] = round(float(np.median(x)), 1) if len(x) else None
            if name != "LEA":
                both = [r for r in sub if r[f"LEA_{n}_ok"]]
                for k in ["max_d", "mean_d"]:
                    g = np.array([r[f"{name}_{n}_{k}"] / r[f"LEA_{n}_{k}"] - 1 for r in both]) * 100
                    d[f"{k}_LEA_better_%"] = round(float((g > 0.5).mean() * 100), 1)
                    d[f"{k}_BUF_better_%"] = round(float((g < -0.5).mean() * 100), 1)
            res[name] = d
        out[f"N={n}"] = res
    return out


if __name__ == "__main__":
    n_trials = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    ks = kansas_scratch()
    rows = monte_carlo_scratch(n_trials)
    summ = summarize(rows)
    json.dump({"kansas": ks, "monte_carlo": summ}, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k.split("_")[0] != "trial", k))
    with open(os.path.join(OUT, "monte_carlo_trials.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(json.dumps(summ, indent=1))
