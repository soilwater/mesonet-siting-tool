# Mesonet Siting Tool

Find where your next monitoring stations should go. Each new station goes in the largest empty area left by the
stations you already have, following the Largest Empty Area method of Patrignani et al. (2020).

It works for any region on Earth: a US state, a country, an uploaded GeoJSON, or a region drawn on the map.

## What it does

1. **Region:** pick a US state or country, upload a GeoJSON, or draw the region.
2. **Existing stations:** upload a CSV (`lat`, `lon`, optional name) or GeoJSON, or click to add stations. You can also start with no stations.
3. **New stations:** get the location and build order of 1 to 250 new stations.
4. **Overlays (optional):** show the 36-km EASE-Grid 2.0, US counties, or your own boundaries, and count the stations in each cell or area.

The results include the largest empty area, the median, mean and longest distance to a station, the largest area
represented by one station, and the Gini coefficient of the areas represented by each station. New stations,
station counts per area, and stations with their area can be downloaded as CSV or GeoJSON.

The example loads the Kansas Mesonet as of 2016 (56 stations), the network analyzed in Patrignani et al. (2020).

## Run it

Open `index.html` in a web browser. No installation is needed; the page loads Leaflet, d3-delaunay and
polygon-clipping from public CDNs.

To serve it locally instead:

```bash
python -m http.server 8765
```

and open http://localhost:8765.

## Files

| File | Contents |
|---|---|
| `index.html` | The app |
| `lea.js` | Largest Empty Area method, support areas and Gini coefficient |
| `ease.js` | EASE-Grid 2.0 global 36-km grid (EPSG:6933) |
| `data/` | US states and counties, world countries, and the Kansas example, loaded by the app |
| `method_comparison/` | Comparison of the Largest Empty Area method with fixed-radius station buffers |

## References

- Patrignani, A., Mohankumar, N., Redmond, C., Santos, E. A., & Knapp, M. (2020). Optimizing the spatial configuration
  of mesoscale environmental monitoring networks using a geometric approach. *Journal of Atmospheric and Oceanic
  Technology*, 37(5), 943–956. https://doi.org/10.1175/JTECH-D-19-0167.1
- Baker, C. B., Cosh, M., Bolten, J., et al. (2022). Working toward a National Coordinated Soil Moisture Monitoring
  Network: Vision, progress, and future directions. *Bulletin of the American Meteorological Society*, 103(12),
  E2719–E2732. https://doi.org/10.1175/BAMS-D-21-0178.1
- Brodzik, M. J., Billingsley, B., Haran, T., Raup, B., & Savoie, M. H. (2012). EASE-Grid 2.0: Incremental but
  significant improvements for Earth-gridded data sets. *ISPRS International Journal of Geo-Information*, 1(1), 32–45.

## Data sources

- US states and counties: US Census Bureau 2010 cartographic boundaries (1:20M), GeoJSON by Eric Celeste
  (https://eric.clst.org/tech/usgeojson/).
- World countries: `countries.geojson`.
- Basemap: Esri World Dark Gray Canvas.
