"""
Compare two station-siting methods for mesoscale monitoring networks.

Methods (both greedy: pick one station, add it, repeat):
  LEA    - Largest Empty Area (Patrignani et al., 2020). Candidates = Voronoi
           vertices inside the boundary + Voronoi-edge/boundary intersections;
           empty circle radius = distance to nearest station; circles clipped by
           the boundary and ranked by clipped area; next station at the centroid
           of the largest clipped area.
  BUF-r  - Station buffer (what index.html does today). Buffer each station by
           r km, subtract the buffers from the boundary, take the largest
           remaining gap polygon, place the station at its centroid (or at a
           point inside it if the centroid falls outside).

Evaluation (on a grid inside the boundary, in an equal-area projection):
  max distance to nearest station (largest gap), mean and 95th-percentile
  distance, % area within 50 km of a station, and mean kriging (GP) standard
  deviation for exponential covariance with several correlation lengths. The
  kriging metric is the expected interpolation error of the network and does
  not favour either method's own objective.
"""
import csv, json, os, sys, time
import numpy as np
from pyproj import Transformer
from shapely.geometry import shape, Polygon, MultiPoint, Point, box
from shapely import contains_xy, voronoi_polygons
from shapely.ops import unary_union
from scipy.spatial import cKDTree
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)


# ----------------------------------------------------------------------------
# Domain: boundary + grid in a local Lambert azimuthal equal-area projection
# ----------------------------------------------------------------------------
class Domain:
    def __init__(self, poly_lonlat, res_km):
        c = poly_lonlat.centroid
        self.lon0, self.lat0 = c.x, c.y
        crs = f"+proj=laea +lat_0={c.y} +lon_0={c.x} +units=km +ellps=WGS84"
        self.fwd = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
        self.inv = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
        # densify so projected edges follow the geodesic shape
        dens = poly_lonlat.segmentize(0.05)
        xs, ys = self.fwd.transform(*np.array(dens.exterior.coords).T)
        self.poly = Polygon(np.c_[xs, ys])
        self.res = res_km
        minx, miny, maxx, maxy = self.poly.bounds
        self.gx = np.arange(minx + res_km / 2, maxx, res_km)
        self.gy = np.arange(miny + res_km / 2, maxy, res_km)
        X, Y = np.meshgrid(self.gx, self.gy)
        self.X, self.Y = X, Y
        self.mask = contains_xy(self.poly, X, Y)
        self.pts = np.c_[X[self.mask], Y[self.mask]]
        self.area = self.poly.area

    def to_xy(self, lon, lat):
        return np.c_[self.fwd.transform(np.asarray(lon), np.asarray(lat))]

    def to_lonlat(self, xy):
        xy = np.atleast_2d(xy)
        return np.c_[self.inv.transform(xy[:, 0], xy[:, 1])]

    def dist_grid(self, st):
        """Distance (km) from every grid cell to nearest station; NaN outside."""
        d = np.full(self.X.shape, np.nan)
        d[self.mask] = cKDTree(st).query(self.pts)[0]
        return d


# ----------------------------------------------------------------------------
# Siting methods: each returns the xy of the next station (or None)
# ----------------------------------------------------------------------------
def pick_buffer_centroid(dom, st, r):
    d = dom.dist_grid(st)
    gap = dom.mask & (d > r)
    if not gap.any():
        return None  # buffers already cover everything -> method has no answer
    lab, n = ndimage.label(gap, structure=np.ones((3, 3)))
    sizes = ndimage.sum(gap, lab, index=np.arange(1, n + 1))
    k = int(np.argmax(sizes)) + 1
    comp = lab == k
    cx, cy = dom.X[comp].mean(), dom.Y[comp].mean()
    # like turf.centerOfMass -> fallback turf.pointOnFeature if outside gap
    ci = np.argmin(np.abs(dom.gx - cx)); cj = np.argmin(np.abs(dom.gy - cy))
    if comp[cj, ci]:
        return np.array([cx, cy])
    cand = np.c_[dom.X[comp], dom.Y[comp]]
    return cand[np.argmin(np.hypot(cand[:, 0] - cx, cand[:, 1] - cy))]


def lea_candidates(dom, st):
    """Patrignani et al. (2020), steps i-vi.
    (i) Voronoi vertices inside the boundary, (ii) intersections of Voronoi
    edges with the boundary, (iii) empty circle radius = distance to nearest
    station, (iv) clip circle with the boundary, (v) clipped area,
    (vi) rank by clipped area. Returns centers, radii, clipped circles, areas."""
    pad = 10 * max(dom.poly.bounds[2] - dom.poly.bounds[0], dom.poly.bounds[3] - dom.poly.bounds[1])
    env = box(*dom.poly.buffer(pad).bounds)
    cells = voronoi_polygons(MultiPoint(st), extend_to=env)
    edges = unary_union([c.exterior for c in cells.geoms])
    verts = np.unique(np.round(np.vstack([np.array(c.exterior.coords) for c in cells.geoms]), 6), axis=0)
    verts = verts[contains_xy(dom.poly, verts[:, 0], verts[:, 1])]
    hits = edges.intersection(dom.poly.boundary)
    bpts = np.array([[p.x, p.y] for p in getattr(hits, "geoms", [hits]) if p.geom_type == "Point"]).reshape(-1, 2)
    C = np.vstack([verts, bpts])
    R = cKDTree(st).query(C)[0]
    circ = [dom.poly.intersection(Point(c).buffer(r, 32)) for c, r in zip(C, R)]
    A = np.array([g.area for g in circ])
    return C, R, circ, A


def pick_lea(dom, st, r=None):
    """Largest Empty Area: next station at the centroid of the largest clipped empty circle."""
    C, R, circ, A = lea_candidates(dom, st)
    g = circ[int(np.argmax(A))].centroid
    return np.array([g.x, g.y])


def greedy(dom, st0, picker, r, n_add):
    st = np.array(st0, float)
    added = []
    for _ in range(n_add):
        p = picker(dom, st, r)
        if p is None:
            break
        added.append(p)
        st = np.vstack([st, p])
    return st, np.array(added).reshape(-1, 2)


# ----------------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------------
def kriging_sd(dom, st, L, nugget=0.05, sub=2):
    """Mean simple-kriging std (sill=1, exp. covariance, corr. length L km)."""
    P = dom.pts[::sub]
    K = np.exp(-np.hypot(*(st[:, None, :] - st[None, :, :]).transpose(2, 0, 1)) / L)
    K += nugget * np.eye(len(st))
    k = np.exp(-np.hypot(P[:, None, 0] - st[None, :, 0], P[:, None, 1] - st[None, :, 1]) / L)
    sol = np.linalg.solve(K, k.T)
    var = 1.0 - np.einsum("ij,ji->i", k, sol)
    return float(np.sqrt(np.clip(var, 0, None)).mean())


def metrics(dom, st, radii=(25, 50, 75, 100), Ls=(25, 50, 100)):
    d = cKDTree(st).query(dom.pts)[0]
    m = {"max_d": d.max(), "mean_d": d.mean(), "p95_d": np.percentile(d, 95)}
    for r in radii:
        m[f"cov{r}"] = 100 * (d <= r).mean()
    for L in Ls:
        m[f"ksd{L}"] = kriging_sd(dom, st, L)
    return m


# ----------------------------------------------------------------------------
# 1) Kansas case study
# ----------------------------------------------------------------------------
def load_kansas():
    gj = json.load(open(os.path.join(ROOT, "kansas_boundary.geojson")))
    poly = shape(gj["features"][0]["geometry"])
    rows = list(csv.DictReader(open(os.path.join(ROOT, "kansas_stations.csv"))))
    lon = [float(r["lon"]) for r in rows]; lat = [float(r["lat"]) for r in rows]
    return poly, lon, lat, [r["station"] for r in rows]


BUFFER_RADII = (25, 50)


def kansas_study(N=10):
    poly, lon, lat, names = load_kansas()
    dom = Domain(poly, res_km=2.0)
    st0 = dom.to_xy(lon, lat)
    print(f"Kansas: area {dom.area:,.0f} km2, {len(st0)} stations, grid {dom.mask.sum()} cells")

    methods = [("Largest Empty Area", pick_lea, None)]
    methods += [(f"Buffer r={r} km", pick_buffer_centroid, r) for r in BUFFER_RADII]
    seqs, curves = {}, {}
    for name, f, r in methods:
        t = time.time()
        st, added = greedy(dom, st0, f, r, N)
        seqs[name] = added
        curves[name] = [metrics(dom, st0)] + [metrics(dom, np.vstack([st0, added[:k]])) for k in range(1, len(added) + 1)]
        print(f"  {name:24s} placed {len(added):2d}  ({time.time()-t:.1f}s)")

    # first-station picks + sequence table
    table = {}
    for name, added in seqs.items():
        ll = dom.to_lonlat(added) if len(added) else np.zeros((0, 2))
        table[name] = [(round(a[1], 3), round(a[0], 3)) for a in ll]

    # sensitivity of the *first* pick to r for the buffer-centroid method
    d0 = dom.dist_grid(st0)
    R = float(np.nanmax(d0))
    rs = np.arange(10, R + 5, 5)
    lea = pick_lea(dom, st0)
    sens = []
    for r in rs:
        pc = pick_buffer_centroid(dom, st0, r)
        sens.append((r,
                     np.nan if pc is None else np.hypot(*(pc - lea)),
                     np.nan if pc is None else cKDTree(st0).query(pc)[0]))
    sens = np.array(sens)

    # ---- figure: maps
    show = list(seqs)
    fig, axs = plt.subplots(len(show), 1, figsize=(7.5, 3.6 * len(show)), constrained_layout=True)
    bx, by = dom.poly.exterior.xy
    for ax, name in zip(axs.flat, show):
        ax.imshow(np.ma.masked_invalid(d0), origin="lower", cmap="YlOrRd",
                  extent=[dom.gx[0], dom.gx[-1], dom.gy[0], dom.gy[-1]], vmin=0, vmax=R)
        ax.plot(bx, by, "k-", lw=1)
        ax.plot(st0[:, 0], st0[:, 1], "^", ms=7, color="k")
        a = seqs[name]
        ax.plot(a[:, 0], a[:, 1], "o", ms=13, mfc="#1d4ed8", mec="w")
        for i, p in enumerate(a):
            ax.text(p[0], p[1], str(i + 1), color="w", ha="center", va="center", fontsize=7, weight="bold")
        m = curves[name][-1]
        ax.set_title(f"{name}\nafter {len(a)} added: max gap {m['max_d']:.0f} km, mean dist {m['mean_d']:.0f} km, "
                     f"kriging SD(L=50) {m['ksd50']:.3f}", fontsize=9)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Kansas (8 stations in the example file): 10 greedy additions. Background = distance to nearest existing station (km)", fontsize=10)
    fig.savefig(os.path.join(OUT, "kansas_sequences.png"), dpi=130)
    plt.close(fig)

    # ---- figure: convergence curves
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.6), constrained_layout=True)
    for name in curves:
        c = curves[name]; k = np.arange(len(c))
        ls = "-" if name.startswith("Largest") else "--"
        for ax, key in zip(axs, ["max_d", "mean_d", "cov50", "ksd50"]):
            ax.plot(k, [x[key] for x in c], ls, marker="o", ms=3, label=name)
    for ax, t in zip(axs, ["Largest gap: max distance to station (km)", "Mean distance to station (km)",
                           "% area within 50 km", "Mean kriging SD, L = 50 km"]):
        ax.set_title(t, fontsize=9); ax.set_xlabel("stations added"); ax.grid(alpha=.3)
    axs[0].legend(fontsize=7)
    fig.savefig(os.path.join(OUT, "kansas_convergence.png"), dpi=130)
    plt.close(fig)

    # ---- figure: sensitivity to r
    fig, ax = plt.subplots(figsize=(6.5, 3.6), constrained_layout=True)
    ax.plot(sens[:, 0], sens[:, 1], "o-", ms=3, label="distance between buffer pick and LEA pick")
    ax.plot(sens[:, 0], sens[:, 2], "s-", ms=3, label="distance from buffer pick to nearest station")
    ax.axvline(R, color="k", ls=":", lw=1); ax.text(R, ax.get_ylim()[1] * .9, f" covering\n radius {R:.0f} km", fontsize=8)
    ax.set_xlabel("buffer radius r (km)"); ax.set_ylabel("km"); ax.grid(alpha=.3); ax.legend(fontsize=8)
    ax.set_title("Where does the buffer method put station #1, as a function of r?", fontsize=9)
    fig.savefig(os.path.join(OUT, "kansas_sensitivity_r.png"), dpi=130)
    plt.close(fig)

    final = {n: {k: round(float(v), 3) for k, v in curves[n][-1].items()} | {"n_added": len(seqs[n])} for n in curves}
    start = {k: round(float(v), 3) for k, v in curves["Largest Empty Area"][0].items()}
    return {"start": start, "final": final, "sequences_latlon": table,
            "covering_radius_km": R, "lea_first_latlon": dom.to_lonlat(lea)[0][::-1].round(3).tolist(),
            "sensitivity": sens.round(1).tolist(),
            "curves": {n: [{k: round(float(v), 3) for k, v in m.items()} for m in c] for n, c in curves.items()}}


# ----------------------------------------------------------------------------
# 2) Monte-Carlo: random boundaries worldwide x random (clustered) networks
# ----------------------------------------------------------------------------
def random_polygon(rng, lat0, lon0, size_km):
    """Irregular, often concave polygon (radial noise) in lon/lat."""
    n = 24
    th = np.sort(rng.uniform(0, 2 * np.pi, n))
    base = rng.uniform(0.35, 1.0, n)
    rad = ndimage.gaussian_filter1d(base, 1.2, mode="wrap") * size_km
    rad *= np.where(rng.random() < 0.5, 1 + 0.5 * np.cos(th - rng.uniform(0, 6)), 1)  # elongation
    x, y = rad * np.cos(th), rad * np.sin(th)
    lat = lat0 + y / 111.0
    lon = lon0 + x / (111.0 * np.cos(np.radians(lat)))
    return Polygon(np.c_[lon, lat]).buffer(0)


def random_network(rng, dom, n):
    """Mix of uniform stations and stations clustered around 'towns'."""
    P = dom.pts
    n_cl = rng.integers(0, n // 2 + 1)
    idx = list(rng.choice(len(P), n - n_cl, replace=False))
    towns = P[rng.choice(len(P), max(1, n_cl // 3 + 1), replace=False)]
    tree = cKDTree(P)
    for _ in range(n_cl):
        t = towns[rng.integers(len(towns))]
        q = t + rng.normal(0, 25, 2)
        idx.append(tree.query(q)[1])
    return P[np.unique(idx)]


METHODS = {"LEA": (pick_lea, None), **{f"BUF{r}": (pick_buffer_centroid, float(r)) for r in BUFFER_RADII}}


def monte_carlo(n_trials=200, n_add=10, seed=1):
    rng = np.random.default_rng(seed)
    regions = [("Kansas-like", 38.5, -98.3), ("Sahel", 14, 2), ("Scandinavia", 64, 16),
               ("Patagonia", -45, -69), ("SE Asia", 15, 103), ("Tropics", -3, 25)]
    rows = []
    t0 = time.time()
    for t in range(n_trials):
        reg, la, lo = regions[t % len(regions)]
        size = rng.uniform(150, 400)
        poly = random_polygon(rng, la + rng.uniform(-3, 3), lo + rng.uniform(-3, 3), size)
        if poly.geom_type != "Polygon":
            poly = max(poly.geoms, key=lambda g: g.area)
        dom = Domain(poly, res_km=4.0)
        n0 = int(rng.integers(4, 40))
        st0 = random_network(rng, dom, n0)
        r_eq = float(np.sqrt(dom.area / (np.pi * len(st0))))  # radius if area split evenly
        cfg = [(n, *METHODS[n]) for n in METHODS]
        cfg = [(n, f, r_eq if r == "eq" else r) for n, f, r in cfg]
        rec = {"trial": t, "region": reg, "area": dom.area, "n0": len(st0), "r_eq": r_eq}
        for name, f, r in cfg:
            st, added = greedy(dom, st0, f, r, n_add)
            rec[f"{name}_n"] = len(added)
            m = metrics(dom, st, radii=(50,), Ls=(25, 100))
            for k, v in m.items():
                rec[f"{name}_{k}"] = v
        rows.append(rec)
        if t % 20 == 0:
            print(f"  trial {t}/{n_trials}  {reg:12s} area {dom.area:8.0f} n0={len(st0):2d}  ({time.time()-t0:.0f}s)")
    return rows


def summarize_mc(rows, n_add=10):
    """Pairwise LEA vs each buffer radius, on trials where the buffer run placed
    all stations (it stops early when the buffers already cover the domain)."""
    keys = ["max_d", "mean_d", "p95_d", "cov50", "ksd25", "ksd100"]
    better_high = {"cov50"}
    out = {"n_trials": len(rows)}
    for b in [n for n in METHODS if n != "LEA"]:
        sub = [r for r in rows if r[f"{b}_n"] == n_add]
        res = {"stopped_early_%": round(100 * (1 - len(sub) / len(rows)), 1), "n_compared": len(sub)}
        for k in keys:
            lea = np.array([r[f"LEA_{k}"] for r in sub]); buf = np.array([r[f"{b}_{k}"] for r in sub])
            diff = (buf / lea - 1) * 100 if len(sub) else np.array([np.nan])
            gain = diff if k in better_high else -diff          # >0 means buffer did better
            res[k] = {"LEA_mean": round(float(lea.mean()), 3) if len(sub) else None,
                      "BUF_mean": round(float(buf.mean()), 3) if len(sub) else None,
                      "LEA_better_%": round(float((gain < -0.5).mean() * 100), 1),
                      "tie_%": round(float((np.abs(gain) <= 0.5).mean() * 100), 1),
                      "BUF_better_%": round(float((gain > 0.5).mean() * 100), 1),
                      "median_%_BUF_vs_LEA": round(float(np.median(diff)), 1)}
        out[b] = res
    return out


if __name__ == "__main__":
    n_trials = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    ks = kansas_study()
    rows = monte_carlo(n_trials=n_trials)
    summ = summarize_mc(rows)
    json.dump({"kansas": ks, "monte_carlo": summ, "n_trials": len(rows)},
              open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=float)
    with open(os.path.join(OUT, "monte_carlo_trials.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(json.dumps(summ, indent=1))
