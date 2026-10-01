/*
 * Largest Empty Area (LEA) siting engine.
 *
 * Patrignani, A., Mohankumar, N., Redmond, C., Santos, E. A., & Knapp, M. (2020).
 * Optimizing the spatial configuration of mesoscale environmental monitoring
 * networks using a geometric approach. J. Atmos. Oceanic Technol., 37(5), 943-956.
 *
 * All geometry is computed in km on a spherical Lambert azimuthal equal-area
 * projection centred on the region, so areas and distances hold anywhere on Earth.
 *
 * Requires d3-delaunay (global `d3`) and polygon-clipping (global `polygonClipping`).
 * Geometries use the polygon-clipping layout: MultiPolygon = [Polygon], Polygon = [outer, ...holes],
 * ring = [[x, y], ...] (closed).
 */
(function (global) {
  'use strict';

  const R_EARTH = 6371.0088;            // mean Earth radius, km
  const D2R = Math.PI / 180;
  const CIRCLE_VERTICES = 128;          // empty-circle polygon resolution
  const MAX_BOUNDARY_VERTICES = 5000;   // boundaries are simplified above this (as in the paper)

  // --------------------------------------------------------------------------
  // Projection
  // --------------------------------------------------------------------------
  function makeProjection(lon0, lat0) {
    const l0 = lon0 * D2R, p0 = lat0 * D2R;
    const sp0 = Math.sin(p0), cp0 = Math.cos(p0);
    return {
      lon0, lat0,
      fwd(lon, lat) {
        const l = lon * D2R - l0, p = lat * D2R;
        const sp = Math.sin(p), cp = Math.cos(p), cl = Math.cos(l);
        const k = Math.sqrt(2 / (1 + sp0 * sp + cp0 * cp * cl));
        return [R_EARTH * k * cp * Math.sin(l), R_EARTH * k * (cp0 * sp - sp0 * cp * cl)];
      },
      inv(x, y) {
        const rho = Math.hypot(x, y);
        if (rho < 1e-12) return [lon0, lat0];
        const c = 2 * Math.asin(Math.min(1, rho / (2 * R_EARTH)));
        const sc = Math.sin(c), cc = Math.cos(c);
        const lat = Math.asin(cc * sp0 + (y * sc * cp0) / rho);
        const lon = l0 + Math.atan2(x * sc, rho * cp0 * cc - y * sp0 * sc);
        return [((lon / D2R + 540) % 360) - 180, lat / D2R];
      }
    };
  }

  // --------------------------------------------------------------------------
  // Planar geometry helpers
  // --------------------------------------------------------------------------
  function closeRing(r) {
    const a = r[0], b = r[r.length - 1];
    return a[0] === b[0] && a[1] === b[1] ? r : r.concat([[a[0], a[1]]]);
  }

  function ringCentroidArea(r) {
    let a = 0, cx = 0, cy = 0;
    for (let i = 0; i < r.length - 1; i++) {
      const f = r[i][0] * r[i + 1][1] - r[i + 1][0] * r[i][1];
      a += f; cx += (r[i][0] + r[i + 1][0]) * f; cy += (r[i][1] + r[i + 1][1]) * f;
    }
    a /= 2;
    return a ? { a, x: cx / (6 * a), y: cy / (6 * a) } : { a: 0, x: r[0][0], y: r[0][1] };
  }

  /** Area (km2) and centroid of a MultiPolygon (outer rings add, holes subtract). */
  function areaCentroid(mp) {
    let A = 0, X = 0, Y = 0;
    for (const poly of mp) {
      poly.forEach((ring, j) => {
        const c = ringCentroidArea(ring);
        const w = (j === 0 ? 1 : -1) * Math.abs(c.a);
        A += w; X += w * c.x; Y += w * c.y;
      });
    }
    return A > 0 ? { area: A, x: X / A, y: Y / A } : { area: 0, x: NaN, y: NaN };
  }

  function bboxOf(mp) {
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    for (const poly of mp) for (const ring of poly) for (const [x, y] of ring) {
      if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y;
    }
    return [x0, y0, x1, y1];
  }

  function segmentsOf(mp) {
    let n = 0;
    for (const poly of mp) for (const ring of poly) n += ring.length - 1;
    const s = new Float64Array(4 * n);
    let k = 0;
    for (const poly of mp) for (const ring of poly) for (let i = 0; i < ring.length - 1; i++) {
      s[k++] = ring[i][0]; s[k++] = ring[i][1]; s[k++] = ring[i + 1][0]; s[k++] = ring[i + 1][1];
    }
    return s;
  }

  function segDist2(px, py, ax, ay, bx, by) {
    const dx = bx - ax, dy = by - ay, L = dx * dx + dy * dy;
    let t = L ? ((px - ax) * dx + (py - ay) * dy) / L : 0;
    t = t < 0 ? 0 : t > 1 ? 1 : t;
    const ex = ax + t * dx - px, ey = ay + t * dy - py;
    return ex * ex + ey * ey;
  }

  function minDistSegs(px, py, s) {
    let m = Infinity;
    for (let i = 0; i < s.length; i += 4) {
      const d = segDist2(px, py, s[i], s[i + 1], s[i + 2], s[i + 3]);
      if (d < m) m = d;
    }
    return Math.sqrt(m);
  }

  function pointInSegs(px, py, s) {
    let inside = false;
    for (let i = 0; i < s.length; i += 4) {
      const ay = s[i + 1], by = s[i + 3];
      if ((ay > py) !== (by > py)) {
        const xi = s[i] + ((py - ay) * (s[i + 2] - s[i])) / (by - ay);
        if (xi > px) inside = !inside;
      }
    }
    return inside;
  }

  function segIntersect(ax, ay, bx, by, cx, cy, dx, dy) {
    const rx = bx - ax, ry = by - ay, sx = dx - cx, sy = dy - cy;
    const den = rx * sy - ry * sx;
    if (Math.abs(den) < 1e-15) return null;
    const t = ((cx - ax) * sy - (cy - ay) * sx) / den;
    const u = ((cx - ax) * ry - (cy - ay) * rx) / den;
    if (t < 0 || t > 1 || u < 0 || u > 1) return null;
    return [ax + t * rx, ay + t * ry];
  }

  function circleRing(cx, cy, r, n = CIRCLE_VERTICES) {
    const ring = [];
    for (let i = 0; i < n; i++) {
      const t = (2 * Math.PI * i) / n;
      ring.push([cx + r * Math.cos(t), cy + r * Math.sin(t)]);
    }
    ring.push([ring[0][0], ring[0][1]]);
    return ring;
  }

  /**
   * Exact area and centroid of (polygon ∩ disk), in one pass over the boundary edges.
   * Each edge a->b forms a triangle with the disk centre; the part of that triangle inside the disk
   * is a mix of triangles (edge pieces inside) and circular sectors (edge pieces outside).
   * Signed pieces add up to the intersection, provided outer rings run counter-clockwise and holes clockwise.
   */
  function diskIntersection(cx, cy, r, segs, ranges = [[0, segs.length / 4]]) {
    const r2 = r * r;
    let A = 0, Mx = 0, My = 0;
    const piece = (px, py, qx, qy) => {
      const mx = (px + qx) / 2, my = (py + qy) / 2;
      if (mx * mx + my * my <= r2) {                               // inside: triangle (0, p, q)
        const w = (px * qy - qx * py) / 2;
        A += w; Mx += (w * (px + qx)) / 3; My += (w * (py + qy)) / 3;
      } else {                                                     // outside: sector from p to q
        const th = Math.atan2(px * qy - py * qx, px * qx + py * qy);
        if (th === 0) return;
        const w = (r2 * th) / 2, a = Math.abs(th);
        const d = (4 * r * Math.sin(a / 2)) / (3 * a), m = Math.atan2(py, px) + th / 2;
        A += w; Mx += w * d * Math.cos(m); My += w * d * Math.sin(m);
      }
    };
    for (const [k0, k1] of ranges) for (let i = 4 * k0; i < 4 * k1; i += 4) {
      const ax = segs[i] - cx, ay = segs[i + 1] - cy, bx = segs[i + 2] - cx, by = segs[i + 3] - cy;
      const dx = bx - ax, dy = by - ay;
      const qa = dx * dx + dy * dy, qb = 2 * (ax * dx + ay * dy), qc = ax * ax + ay * ay - r2;
      const disc = qb * qb - 4 * qa * qc;
      let t0 = 0;
      if (qa > 0 && disc > 0) {
        const sq = Math.sqrt(disc);
        for (const t of [(-qb - sq) / (2 * qa), (-qb + sq) / (2 * qa)]) {
          if (t > 0 && t < 1) {
            piece(ax + dx * t0, ay + dy * t0, ax + dx * t, ay + dy * t);
            t0 = t;
          }
        }
      }
      piece(ax + dx * t0, ay + dy * t0, bx, by);
    }
    return A > 0 ? { area: A, x: cx + Mx / A, y: cy + My / A } : { area: 0, x: cx, y: cy };
  }

  /** Outer rings counter-clockwise, holes clockwise (needed by diskIntersection). */
  function orient(mp) {
    return mp.map(poly => poly.map((ring, j) => {
      const ccw = ringCentroidArea(ring).a > 0;
      return (j === 0) === ccw ? ring : ring.slice().reverse();
    }));
  }

  /** Uniform-grid index over boundary segments: fast point-in-polygon and proximity tests. */
  class SegIndex {
    constructor(segs, bbox) {
      this.segs = segs;
      const n = segs.length / 4;
      const [x0, y0, x1, y1] = bbox;
      this.x0 = x0; this.y0 = y0; this.x1 = x1; this.y1 = y1;
      this.G = Math.max(1, Math.min(512, Math.ceil(Math.sqrt(n))));
      this.cw = (x1 - x0) / this.G || 1;
      this.ch = (y1 - y0) / this.G || 1;
      this.cells = Array.from({ length: this.G * this.G }, () => []);
      for (let s = 0; s < n; s++) {
        const ax = segs[4 * s], ay = segs[4 * s + 1], bx = segs[4 * s + 2], by = segs[4 * s + 3];
        const [i0, i1, j0, j1] = this._range(Math.min(ax, bx), Math.max(ax, bx), Math.min(ay, by), Math.max(ay, by));
        for (let j = j0; j <= j1; j++) for (let i = i0; i <= i1; i++) this.cells[j * this.G + i].push(s);
      }
      this.stamp = new Uint32Array(n);
      this.tick = 0;
    }
    _range(xa, xb, ya, yb) {
      const c = (v, lo) => Math.max(0, Math.min(this.G - 1, Math.floor(v / lo)));
      return [c(xa - this.x0, this.cw), c(xb - this.x0, this.cw), c(ya - this.y0, this.ch), c(yb - this.y0, this.ch)];
    }
    /** Calls fn(segmentIndex) once per segment near the box; stops early if fn returns true. */
    each(xa, ya, xb, yb, fn) {
      if (xb < this.x0 || xa > this.x1 || yb < this.y0 || ya > this.y1) return false;
      if (++this.tick === 0xffffffff) { this.stamp.fill(0); this.tick = 1; }
      const [i0, i1, j0, j1] = this._range(xa, xb, ya, yb);
      for (let j = j0; j <= j1; j++) for (let i = i0; i <= i1; i++) {
        for (const s of this.cells[j * this.G + i]) {
          if (this.stamp[s] === this.tick) continue;
          this.stamp[s] = this.tick;
          if (fn(s) === true) return true;
        }
      }
      return false;
    }
    pointIn(px, py) {
      if (px < this.x0 || px > this.x1 || py < this.y0 || py > this.y1) return false;
      const s = this.segs;
      let inside = false;
      this.each(px, py, this.x1, py, k => {
        const ay = s[4 * k + 1], by = s[4 * k + 3];
        if ((ay > py) !== (by > py)) {
          const xi = s[4 * k] + ((py - ay) * (s[4 * k + 2] - s[4 * k])) / (by - ay);
          if (xi > px) inside = !inside;
        }
      });
      return inside;
    }
    /** True if any boundary segment lies closer than r to (px, py). */
    anyWithin(px, py, r) {
      const s = this.segs, r2 = r * r;
      return this.each(px - r, py - r, px + r, py + r,
        k => segDist2(px, py, s[4 * k], s[4 * k + 1], s[4 * k + 2], s[4 * k + 3]) < r2);
    }
  }

  // Max-heap for the pole-of-inaccessibility search
  class Heap {
    constructor() { this.a = []; }
    get size() { return this.a.length; }
    push(v) {
      const a = this.a; a.push(v);
      let i = a.length - 1;
      while (i > 0) { const p = (i - 1) >> 1; if (a[p].max >= a[i].max) break; [a[p], a[i]] = [a[i], a[p]]; i = p; }
    }
    pop() {
      const a = this.a, top = a[0], last = a.pop();
      if (a.length) {
        a[0] = last;
        let i = 0;
        for (;;) {
          const l = 2 * i + 1, r = l + 1;
          let m = i;
          if (l < a.length && a[l].max > a[m].max) m = l;
          if (r < a.length && a[r].max > a[m].max) m = r;
          if (m === i) break;
          [a[m], a[i]] = [a[i], a[m]]; i = m;
        }
      }
      return top;
    }
  }

  /** Point inside the polygon farthest from its boundary (polylabel algorithm). */
  function poleOfInaccessibility(mp, segs, inside) {
    const [x0, y0, x1, y1] = bboxOf(mp);
    const w = x1 - x0, h = y1 - y0, size = Math.min(w, h);
    if (!(size > 0)) return [x0, y0];
    const precision = size / 2000;
    const pin = inside || ((x, y) => pointInSegs(x, y, segs));
    const sd = (x, y) => (pin(x, y) ? 1 : -1) * minDistSegs(x, y, segs);
    const cell = (x, y, hh) => { const d = sd(x, y); return { x, y, h: hh, d, max: d + hh * Math.SQRT2 }; };
    const heap = new Heap();
    const hh = size / 2;
    for (let x = x0; x < x1; x += size) for (let y = y0; y < y1; y += size) heap.push(cell(x + hh, y + hh, hh));
    const c = areaCentroid(mp);
    let best = cell(c.x, c.y, 0);
    const bc = cell(x0 + w / 2, y0 + h / 2, 0);
    if (bc.d > best.d) best = bc;
    let iter = 0;
    while (heap.size && iter++ < 50000) {
      const q = heap.pop();
      if (q.d > best.d) best = q;
      if (q.max - best.d <= precision) continue;
      const h2 = q.h / 2;
      heap.push(cell(q.x - h2, q.y - h2, h2)); heap.push(cell(q.x + h2, q.y - h2, h2));
      heap.push(cell(q.x - h2, q.y + h2, h2)); heap.push(cell(q.x + h2, q.y + h2, h2));
    }
    return [best.x, best.y];
  }

  function densifyLonLat(ring, maxDeg) {
    const out = [];
    for (let i = 0; i < ring.length - 1; i++) {
      const [a, b] = [ring[i], ring[i + 1]];
      const k = Math.max(1, Math.ceil(Math.max(Math.abs(b[0] - a[0]), Math.abs(b[1] - a[1])) / maxDeg));
      for (let j = 0; j < k; j++) out.push([a[0] + ((b[0] - a[0]) * j) / k, a[1] + ((b[1] - a[1]) * j) / k]);
    }
    out.push(ring[ring.length - 1]);
    return out;
  }

  function simplifyRing(ring, tol) {
    const n = ring.length;
    if (n <= 4) return ring;
    const keep = new Uint8Array(n);
    keep[0] = keep[n - 1] = 1;
    const tol2 = tol * tol, stack = [[0, n - 1]];
    while (stack.length) {
      const [a, b] = stack.pop();
      let m = -1, md = tol2;
      for (let i = a + 1; i < b; i++) {
        const d = segDist2(ring[i][0], ring[i][1], ring[a][0], ring[a][1], ring[b][0], ring[b][1]);
        if (d > md) { md = d; m = i; }
      }
      if (m > 0) { keep[m] = 1; stack.push([a, m], [m, b]); }
    }
    const out = ring.filter((_, i) => keep[i]);
    return out.length >= 4 ? out : ring;
  }

  // --------------------------------------------------------------------------
  // Domain (boundary) set-up
  // --------------------------------------------------------------------------
  /**
   * boundaryLonLat: MultiPolygon coordinates in [lon, lat].
   * Returns the projected domain used by every other function.
   */
  function buildDomain(boundaryLonLat, opts = {}) {
    let lo0 = Infinity, lo1 = -Infinity, la0 = Infinity, la1 = -Infinity;
    let wo0 = Infinity, wo1 = -Infinity;                 // longitudes wrapped to 0..360
    for (const poly of boundaryLonLat) for (const [lon, lat] of poly[0]) {
      if (lon < lo0) lo0 = lon; if (lon > lo1) lo1 = lon; if (lat < la0) la0 = lat; if (lat > la1) la1 = lat;
      const w = lon < 0 ? lon + 360 : lon;
      if (w < wo0) wo0 = w; if (w > wo1) wo1 = w;
    }
    // regions that straddle the 180° meridian (e.g., Alaska) are centred using wrapped longitudes
    let lonC = (lo0 + lo1) / 2;
    if (wo1 - wo0 < lo1 - lo0) lonC = (((wo0 + wo1) / 2 + 180) % 360) - 180;
    const proj = makeProjection(lonC, (la0 + la1) / 2);

    let mp = boundaryLonLat.map(poly => poly.map(ring =>
      densifyLonLat(closeRing(ring), 0.05).map(([lon, lat]) => proj.fwd(lon, lat))));

    let nv = 0;
    for (const poly of mp) for (const ring of poly) nv += ring.length;
    if (nv > MAX_BOUNDARY_VERTICES) {
      const [x0, y0, x1, y1] = bboxOf(mp);
      let tol = Math.max(x1 - x0, y1 - y0) / 20000;
      let simplified = mp;
      while (nv > MAX_BOUNDARY_VERTICES) {
        simplified = mp.map(poly => poly.map(ring => simplifyRing(ring, tol)));
        nv = 0;
        for (const poly of simplified) for (const ring of poly) nv += ring.length;
        tol *= 1.6;
      }
      mp = simplified;
    }
    mp = orient(polygonClipping.union(mp));

    const bbox = bboxOf(mp);
    const segs = segmentsOf(mp);
    // segment range and bbox of each polygon part (islands, exclaves)
    const parts = [];
    let k0 = 0;
    for (const poly of mp) {
      let n = 0;
      for (const ring of poly) n += ring.length - 1;
      parts.push({ s0: k0, s1: k0 + n, bbox: bboxOf([poly]), area: areaCentroid([poly]).area });
      k0 += n;
    }
    const index = new SegIndex(segs, bbox);
    const { area } = areaCentroid(mp);

    // Regular sample of points inside the region, used for network statistics
    const target = opts.samples || 15000;
    const step = Math.sqrt(area / target);
    const pts = [];
    for (let y = bbox[1] + step / 2; y < bbox[3]; y += step)
      for (let x = bbox[0] + step / 2; x < bbox[2]; x += step)
        if (index.pointIn(x, y)) pts.push(x, y);

    return { proj, mp, bbox, segs, parts, index, area, samples: Float64Array.from(pts), sampleStep: step };
  }

  function contains(dom, x, y) { return dom.index.pointIn(x, y); }

  function nearestFn(pts) {
    if (pts.length >= 3) {
      const del = d3.Delaunay.from(pts);
      let hint = 0;
      return { del, nearest(x, y) { hint = del.find(x, y, hint); return Math.hypot(pts[hint][0] - x, pts[hint][1] - y); } };
    }
    return { del: pts.length === 2 ? d3.Delaunay.from(pts) : null,
      nearest: (x, y) => Math.min(...pts.map(p => Math.hypot(p[0] - x, p[1] - y))) };
  }

  function voronoiBounds(pts, dom) {
    let [x0, y0, x1, y1] = dom.bbox;
    for (const [x, y] of pts) { if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y; }
    const pad = 10 * Math.max(x1 - x0, y1 - y0, 1);
    return [x0 - pad, y0 - pad, x1 + pad, y1 + pad];
  }

  // --------------------------------------------------------------------------
  // One LEA step
  // --------------------------------------------------------------------------
  /**
   * Finds the largest empty area for the current network `pts` ([[x, y], ...] in km).
   * Returns { x, y }         suggested station (centroid of the largest empty area)
   *         { cx, cy, r }    centre and radius of the empty circle
   *         { area, shape }  clipped empty area (km2) and its geometry
   * Starting rules when Voronoi vertices do not exist yet:
   *   0 stations -> point farthest from the boundary
   *   1 station  -> candidate centres are the boundary vertices
   */
  function leaStep(pts, dom) {
    const n = pts.length;
    if (n === 0) {
      const p = poleOfInaccessibility(dom.mp, dom.segs, (x, y) => contains(dom, x, y));
      return { x: p[0], y: p[1], cx: p[0], cy: p[1], r: minDistSegs(p[0], p[1], dom.segs), area: dom.area, shape: dom.mp };
    }

    // (i-ii) candidate centres: Voronoi vertices inside + Voronoi edges crossing the boundary
    const cand = [];
    const seen = new Set();
    const add = (x, y) => {
      const key = Math.round(x * 1e4) + ',' + Math.round(y * 1e4);
      if (!seen.has(key)) { seen.add(key); cand.push(x, y); }
    };
    const [bx0, by0, bx1, by1] = dom.bbox;
    const s = dom.segs;
    const nf = nearestFn(pts);

    // parts (islands, exclaves) touched by at least one Voronoi vertex or edge
    const touched = new Uint8Array(dom.parts.length);
    const partOfSeg = k => dom.parts.findIndex(p => k >= p.s0 && k < p.s1);
    const partOfPoint = (x, y) => dom.parts.findIndex(p => x >= p.bbox[0] && x <= p.bbox[2] && y >= p.bbox[1] && y <= p.bbox[3] &&
      pointInSegs(x, y, s.subarray(4 * p.s0, 4 * p.s1)));

    if (n > 1) {
      const vor = nf.del.voronoi(voronoiBounds(pts, dom));
      const cc = vor.circumcenters;
      for (let i = 0; i < cc.length; i += 2) {
        const x = cc[i], y = cc[i + 1];
        if (x >= bx0 && x <= bx1 && y >= by0 && y <= by1 && contains(dom, x, y)) {
          add(x, y);
          const p = partOfPoint(x, y);
          if (p >= 0) touched[p] = 1;
        }
      }
      for (let i = 0; i < n; i++) {
        const cell = vor.cellPolygon(i);
        if (!cell) continue;
        for (let j = 0; j < cell.length - 1; j++) {
          const [ax, ay] = cell[j], [bx, by] = cell[j + 1];
          const ex0 = Math.min(ax, bx), ex1 = Math.max(ax, bx), ey0 = Math.min(ay, by), ey1 = Math.max(ay, by);
          if (ex1 < bx0 || ex0 > bx1 || ey1 < by0 || ey0 > by1) continue;
          dom.index.each(ex0, ey0, ex1, ey1, k => {
            const p = segIntersect(ax, ay, bx, by, s[4 * k], s[4 * k + 1], s[4 * k + 2], s[4 * k + 3]);
            if (p) { add(p[0], p[1]); touched[partOfSeg(k)] = 1; }
          });
        }
      }
    }
    // A part that lies inside a single Voronoi cell (always the case with one station, and for small
    // islands) has no vertex or edge candidates. Its point farthest from the cell's station is a vertex
    // of its convex hull, so those vertices become its candidates.
    dom.parts.forEach((p, j) => {
      if (touched[j]) return;
      const vx = [];
      for (let k = p.s0; k < p.s1; k++) vx.push([s[4 * k], s[4 * k + 1]]);
      const hull = vx.length > 3 ? d3.Delaunay.from(vx).hull : vx.map((_, i) => i);
      for (const i of hull) add(vx[i][0], vx[i][1]);
    });

    // (iii) empty-circle radius = distance to nearest station
    const m = cand.length / 2;
    const R = new Float64Array(m);
    for (let i = 0; i < m; i++) R[i] = nf.nearest(cand[2 * i], cand[2 * i + 1]);
    const order = Array.from({ length: m }, (_, i) => i).sort((a, b) => R[b] - R[a]);

    // (iv-vi) clip circles by the boundary, rank by clipped area
    // (a clipped circle can never exceed the full circle, so stop once no circle can win)
    let best = null;
    for (const i of order) {
      const r = R[i];
      if (!(r > 0)) continue;
      if (best && Math.PI * r * r <= best.area) break;
      const x = cand[2 * i], y = cand[2 * i + 1];
      let area, gx = x, gy = y, clipped = false;
      if (contains(dom, x, y) && !dom.index.anyWithin(x, y, r)) {
        area = Math.PI * r * r;                                      // circle fully inside
      } else {
        // only parts whose bounding box reaches the circle can overlap it
        const near = dom.parts.filter(p => p.bbox[0] <= x + r && p.bbox[2] >= x - r && p.bbox[1] <= y + r && p.bbox[3] >= y - r);
        if (best && near.reduce((a, p) => a + p.area, 0) <= best.area) continue;
        const d = diskIntersection(x, y, r, s, near.map(p => [p.s0, p.s1]));
        area = d.area; gx = d.x; gy = d.y; clipped = true;
        if (!(area > 0)) continue;
      }
      if (!best || area > best.area) best = { cx: x, cy: y, r, area, x: gx, y: gy, clipped };
    }
    if (!best) return null;

    // geometry of the winning empty area (for display)
    best.shape = best.clipped
      ? polygonClipping.intersection(dom.mp, [circleRing(best.cx, best.cy, best.r)])
      : [[circleRing(best.cx, best.cy, best.r)]];
    delete best.clipped;
    if (!contains(dom, best.x, best.y)) {
      // centroid of an oddly shaped clipped area fell outside the region: use its most interior point
      const part = best.shape.reduce((a, b) => (areaCentroid([b]).area > areaCentroid([a]).area ? b : a));
      const p = poleOfInaccessibility([part], segmentsOf([part]));
      best.x = p[0]; best.y = p[1];
    }
    return best;
  }

  // --------------------------------------------------------------------------
  // Network diagnostics
  // --------------------------------------------------------------------------
  /** Distance statistics over the region for network `pts`, footprint radius r (km). */
  function networkStats(pts, dom, r) {
    if (!pts.length) return null;
    const nf = nearestFn(pts);
    const S = dom.samples;
    const n = S.length / 2;
    const dist = new Float64Array(n);
    let within = 0, sum = 0, max = 0;
    for (let i = 0; i < n; i++) {
      const d = nf.nearest(S[2 * i], S[2 * i + 1]);
      dist[i] = d; sum += d; if (d > max) max = d; if (d <= r) within++;
    }
    dist.sort();
    const median = n % 2 ? dist[(n - 1) / 2] : (dist[n / 2 - 1] + dist[n / 2]) / 2;
    return { coveredPct: (100 * within) / n, meanDist: sum / n, medianDist: median, maxDist: max };
  }

  /** Voronoi cells (station support areas) clipped to the region. */
  function voronoiCells(pts, dom) {
    if (pts.length < 2) return pts.length ? [dom.mp] : [];
    const vor = d3.Delaunay.from(pts).voronoi(voronoiBounds(pts, dom));
    const out = [];
    for (let i = 0; i < pts.length; i++) {
      const cell = vor.cellPolygon(i);
      if (!cell) continue;
      const clipped = polygonClipping.intersection(dom.mp, [cell]);
      if (clipped.length) out.push(clipped);
    }
    return out;
  }

  /**
   * Support area of each station (km2): its Voronoi cell inside the region (0 if none).
   * Only cells that cross the boundary are clipped; cells fully inside use their polygon area.
   */
  function supportAreas(pts, dom) {
    const n = pts.length;
    if (n === 0) return [];
    if (n === 1) return [dom.area];
    const vor = d3.Delaunay.from(pts).voronoi(voronoiBounds(pts, dom));
    const s = dom.segs;
    const out = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      const cell = vor.cellPolygon(i);
      if (!cell) continue;
      let crosses = false;
      for (let j = 0; j < cell.length - 1 && !crosses; j++) {
        const [ax, ay] = cell[j], [bx, by] = cell[j + 1];
        crosses = dom.index.each(Math.min(ax, bx), Math.min(ay, by), Math.max(ax, bx), Math.max(ay, by),
          k => !!segIntersect(ax, ay, bx, by, s[4 * k], s[4 * k + 1], s[4 * k + 2], s[4 * k + 3]));
      }
      out[i] = !crosses && contains(dom, cell[0][0], cell[0][1])
        ? Math.abs(ringCentroidArea(cell).a)
        : areaCentroid(polygonClipping.intersection(dom.mp, [cell])).area;
    }
    return out;
  }

  /**
   * Gini coefficient of station support areas (Patrignani et al., 2020): G = (A - B) / A, with A = 0.5 the
   * area under the line of perfect equality and B the area under the Lorenz curve (trapezoidal rule).
   * 0 = all stations support the same area. Stations with no area inside the region are ignored.
   */
  function gini(areas) {
    const a = areas.filter(v => v > 0).sort((x, y) => x - y);
    if (a.length < 2) return 0;
    const total = a.reduce((p, q) => p + q, 0);
    let cum = 0, prev = 0, B = 0;
    for (const v of a) { cum += v / total; B += (prev + cum) / 2 / a.length; prev = cum; }
    return (0.5 - B) / 0.5;
  }

  /** Voronoi cell of station i clipped to the region (MultiPolygon, km). */
  function supportCell(pts, dom, i) {
    if (pts.length === 1) return dom.mp;
    const cell = d3.Delaunay.from(pts).voronoi(voronoiBounds(pts, dom)).cellPolygon(i);
    return cell ? polygonClipping.intersection(dom.mp, [cell]) : [];
  }

  global.LEA = {
    supportAreas, supportCell, gini,
    makeProjection, buildDomain, leaStep, networkStats, voronoiCells, areaCentroid, contains, diskIntersection
  };
})(typeof window !== 'undefined' ? window : globalThis);
