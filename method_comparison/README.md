# Largest Empty Area vs. station buffer

This compares two ways to choose where new stations go. Both place one station at a time and then repeat:

- **LEA (Largest Empty Area)** from Patrignani et al. (2020):
  - Candidate circles are centred on Voronoi vertices and on the points where Voronoi edges cross the boundary.
  - Each circle's radius is the distance to the nearest station.
  - Circles are clipped to the boundary and ranked by clipped area.
  - The new station goes at the centroid of the largest clipped circle.
- **Buffer** (what `index.html` does now): buffer each station by r km, remove the buffers from the boundary, and place the new station at the centroid of the largest uncovered piece. Tested at r = 25 and 50 km.

All distances and areas are computed in a local equal-area projection centred on the boundary, so the code works anywhere in the world.

## Scripts

Run both with `D:\anaconda3\python.exe`.

| Script | Scenario | Outputs |
|---|---|---|
| `compare_methods.py 300` | Adding stations to an existing network | `results/` |
| `scratch_start.py 150` | Building a network from zero | `results_scratch/` |

## 1. Adding to an existing network

**Kansas: 8 example stations, 10 added** (`results/kansas_sequences.png`)

| Method | Largest gap (km) | Mean distance (km) | Kriging SD, L = 50 km |
|---|---|---|---|
| **LEA** | **150** | **47** | **0.866** |
| Buffer 25 km | 191 | 67 | 0.909 |
| Buffer 50 km | 191 | 61 | 0.890 |

With both buffer sizes, the uncovered area is one connected blob whose centroid is in central Kansas. All 10 buffer stations cluster there.

**Worldwide: 300 random concave boundaries in six regions from 45°S to 64°N, each with 4–40 existing stations, 10 added**

A win means a result more than 0.5% better.

| Buffer | Buffer gives no answer* | Largest gap: LEA / buffer wins | Mean distance: LEA / buffer wins | Median change in largest gap, buffer vs. LEA |
|---|---|---|---|---|
| 25 km | 0% | 88% / 9% | 93% / 5% | +51% |
| 50 km | 28% | 37% / 59% | 54% / 31% | −6% |

\* The buffers cover the whole domain before all stations are placed, so the method has nothing left to suggest.

The 50 km row only counts trials where the buffer method finished. Those are mostly sparse networks, where 50 km happens to match the station spacing.

## 2. Building from scratch

LEA needs a starting rule because there are no Voronoi vertices yet:

- **First station:** the point farthest from the boundary.
- **Second station:** candidate centres are the boundary vertices, with circles clipped and ranked the same way.
- **From the third station on:** the published method, unchanged.

Every network is scored against an **all-at-once layout** with the same number of stations. That layout is a centroidal Voronoi (k-means) layout: near-hexagonal inside, adapted to the boundary. It is what you would install if all stations went in at once and siting order did not matter.

**Kansas from scratch** (`results_scratch/kansas_scratch_maps.png`). Each cell shows largest gap / mean distance, in km.

| Stations | LEA | Buffer 25 km | Buffer 50 km | All at once |
|---|---|---|---|---|
| 5 | 280 / 109 | 347 / 172 | 335 / 153 | 182 / 83 |
| 10 | 191 / 65 | 335 / 161 | 309 / 133 | 106 / 57 |
| 20 | 119 / 45 | 315 / 143 | 270 / 100 | 77 / 40 |
| 40 | 108 / 31 | 287 / 118 | 146 / 45 | 57 / 28 |

**Worldwide: 150 random boundaries.** Each cell shows the median percentage above the all-at-once layout (lower is better), for largest gap / mean distance.

| Stations | LEA | Buffer 25 km | Buffer 50 km (no answer) |
|---|---|---|---|
| 5 | +45% / +18% | +74% / +84% | +58% / +52% (0%) |
| 10 | +49% / +13% | +132% / +130% | +101% / +68% (1%) |
| 20 | +47% / +12% | +204% / +170% | +129% / +60% (13%) |
| 30 | +49% / +7% | +246% / +176% | +137% / +50% (33%) |

- **LEA vs. buffer:** LEA beats the buffer method in 79–100% of trials at every network size.
- **Buffer 25 km:** it never leaves the central cluster, so it gets worse as more stations are added.
- **LEA vs. the all-at-once layout:**
  - Mean distance ends up within 7% at 30 stations.
  - The largest gap stays about 50% larger, because edges and corners are filled last.

## Verdict

**Adding to an existing network:** use LEA.

- It needs no parameters and always gives an answer.
- The clipped circle protects it against long, thin gaps.
- The buffer method is only competitive at about 50 km, in sparse networks, and it often gives no answer at all.
- At 25 km the buffer method is clearly worse.

**Building from scratch:** the buffer method is not usable, because it grows a single cluster from the centre. LEA works, but it under-serves edges and corners. If all stations will be planned together, start from an all-at-once (hexagonal/centroidal) layout, and use LEA to rank the siting order or to add stations later.

**LEA's weak spot in both cases is corners.** They are never candidate points and their clipped circles are small, so they are filled last. The effect is small in dense existing networks and large when starting from zero.

Buffers are still useful as a coverage display ("% of the state within r km of a station"), but not for choosing where stations go.
