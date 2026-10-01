/*
 * EASE-Grid 2.0 global grids (EPSG:6933, Lambert cylindrical equal-area, standard parallel 30°, WGS84).
 * Brodzik, M. J., Billingsley, B., Haran, T., Raup, B., & Savoie, M. H. (2012). EASE-Grid 2.0:
 * Incremental but significant improvements for Earth-gridded data sets. ISPRS Int. J. Geo-Inf., 1(1), 32-45.
 *
 * Rows and columns are 0-based from the upper-left corner, as in NSIDC and SMAP products.
 * Cells are rectangles in latitude/longitude because the projection is cylindrical.
 */
(function (global) {
  'use strict';

  const A = 6378137.0;                    // WGS84 semi-major axis (m)
  const E2 = 0.0066943799901413165;       // WGS84 first eccentricity squared
  const E = Math.sqrt(E2);
  const D2R = Math.PI / 180;
  const PHI1 = 30 * D2R;                  // standard parallel
  const K0 = Math.cos(PHI1) / Math.sqrt(1 - E2 * Math.sin(PHI1) ** 2);

  const q = phi => {
    const s = Math.sin(phi);
    return (1 - E2) * (s / (1 - E2 * s * s) - (1 / (2 * E)) * Math.log((1 - E * s) / (1 + E * s)));
  };
  const QP = q(Math.PI / 2);

  // Authalic-latitude series for the inverse
  const E4 = E2 * E2, E6 = E4 * E2;
  const C2 = E2 / 3 + (31 * E4) / 180 + (517 * E6) / 5040;
  const C4 = (23 * E4) / 360 + (251 * E6) / 3780;
  const C6 = (761 * E6) / 45360;

  function fwd(lon, lat) {
    return [A * K0 * lon * D2R, (A * q(lat * D2R)) / (2 * K0)];
  }

  function inv(x, y) {
    const beta = Math.asin(Math.max(-1, Math.min(1, (2 * K0 * y) / (A * QP))));
    const phi = beta + C2 * Math.sin(2 * beta) + C4 * Math.sin(4 * beta) + C6 * Math.sin(6 * beta);
    return [x / (A * K0) / D2R, phi / D2R];
  }

  const GRIDS = {
    M36: { name: 'EASE-Grid 2.0 Global, 36 km', cell: 36032.220840584, cols: 964, rows: 406 }
  };

  function makeGrid(id) {
    const g = GRIDS[id];
    const x0 = (-g.cols * g.cell) / 2, y0 = (g.rows * g.cell) / 2;   // upper-left corner

    /** { row, col } of the cell containing lon/lat, or null outside the grid. */
    function cellOf(lon, lat) {
      const [x, y] = fwd(lon, lat);
      const col = Math.floor((x - x0) / g.cell), row = Math.floor((y0 - y) / g.cell);
      return row >= 0 && row < g.rows && col >= 0 && col < g.cols ? { row, col } : null;
    }

    /** [[south, west], [north, east]] of a cell. */
    function bounds(row, col) {
      const [w, n] = inv(x0 + col * g.cell, y0 - row * g.cell);
      const [e, s] = inv(x0 + (col + 1) * g.cell, y0 - (row + 1) * g.cell);
      return [[s, w], [n, e]];
    }

    function center(row, col) {
      return inv(x0 + (col + 0.5) * g.cell, y0 - (row + 0.5) * g.cell);
    }

    /**
     * Cells that touch a region. `boundary` is a MultiPolygon in [lon, lat];
     * `contains(lon, lat)` tests whether a point is inside it.
     * A cell counts if its center is inside or the boundary line passes through it.
     */
    function cellsInRegion(boundary, contains) {
      let lo0 = Infinity, lo1 = -Infinity, la0 = Infinity, la1 = -Infinity;
      for (const poly of boundary) for (const [lon, lat] of poly[0]) {
        lo0 = Math.min(lo0, lon); lo1 = Math.max(lo1, lon); la0 = Math.min(la0, lat); la1 = Math.max(la1, lat);
      }
      const a = cellOf(lo0, la1) || { row: 0, col: 0 };
      const b = cellOf(lo1, la0) || { row: g.rows - 1, col: g.cols - 1 };
      const keys = new Set();
      for (let r = a.row; r <= b.row; r++) for (let c = a.col; c <= b.col; c++) {
        const [lon, lat] = center(r, c);
        if (contains(lon, lat)) keys.add(r * g.cols + c);
      }
      // walk the boundary in steps well under a cell to catch partly covered cells
      const step = (g.cell / 111320) / 4;
      for (const poly of boundary) for (const ring of poly) for (let i = 0; i < ring.length - 1; i++) {
        const [x1, y1] = ring[i], [x2, y2] = ring[i + 1];
        const k = Math.max(1, Math.ceil(Math.max(Math.abs(x2 - x1), Math.abs(y2 - y1)) / step));
        for (let j = 0; j <= k; j++) {
          const c = cellOf(x1 + ((x2 - x1) * j) / k, y1 + ((y2 - y1) * j) / k);
          if (c) keys.add(c.row * g.cols + c.col);
        }
      }
      return [...keys].sort((p, q2) => p - q2).map(k => ({ row: Math.floor(k / g.cols), col: k % g.cols }));
    }

    return { ...g, id, cellOf, bounds, center, cellsInRegion };
  }

  global.EASE = { fwd, inv, makeGrid };
})(typeof window !== 'undefined' ? window : globalThis);
