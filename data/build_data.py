"""
Builds the data files used by the web app. Run with D:\\anaconda3\\python.exe data\\build_data.py

Inputs
  data/gz_2010_us_040_00_20m.json   US Census 2010 cartographic boundaries, states (1:20M),
  data/gz_2010_us_050_00_20m.json   and counties, as GeoJSON from https://eric.clst.org/tech/usgeojson/
                                    (50 states, DC and Puerto Rico; Latin-1 encoded)
  ks_mesonet_stations.xlsx          Kansas Mesonet stations today (optional; data/ks_mesonet_stations.csv is used if absent)
  KS_mesonet_geoinfo.csv            Kansas Mesonet stations in 2016 (Patrignani et al., 2020)
  countries.geojson                 World country boundaries

Outputs
  data/us_states.js                 window.US_STATES (FeatureCollection, UTF-8, 5-decimal coordinates)
  data/us_states.geojson            same, as GeoJSON for loading through the app or a GIS
  data/us_counties.geojson          counties re-encoded as UTF-8
  data/us_counties.js               window.US_COUNTIES (id = state+county FIPS, name = "Riley County, Kansas")
  data/kansas_example.js            window.KANSAS_EXAMPLE = { boundary, stations }
  data/ks_mesonet_stations.csv      Kansas Mesonet stations as CSV
  data/countries.js                 window.COUNTRIES (FeatureCollection, 4-decimal coordinates)
"""
import csv, json, os, zipfile
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def rounded(coords, nd=5):
    if isinstance(coords[0], (int, float)):
        return [round(coords[0], nd), round(coords[1], nd)]
    return [rounded(c, nd) for c in coords]


def load_census(name):
    with open(os.path.join(HERE, name), encoding="latin-1") as fh:
        return json.load(fh)


def read_xlsx(path):
    """Minimal .xlsx reader (the file's autofilter trips openpyxl)."""
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", ns):
            shared.append("".join(t.text or "" for t in si.iter("{%s}t" % ns["m"])))
    def col_index(ref):                      # "C12" -> 2 (empty cells are omitted from the XML)
        n = 0
        for ch in ref:
            if not ch.isalpha():
                break
            n = n * 26 + ord(ch.upper()) - 64
        return n - 1

    rows = []
    for r in ET.fromstring(z.read("xl/worksheets/sheet1.xml")).iter("{%s}row" % ns["m"]):
        vals = {}
        for c in r.findall("m:c", ns):
            v = c.find("m:v", ns)
            vals[col_index(c.get("r"))] = "" if v is None else shared[int(v.text)] if c.get("t") == "s" else v.text
        rows.append([vals.get(i, "") for i in range(max(vals) + 1)] if vals else [])
    head, body = rows[0], rows[1:]
    return [dict(zip(head, r)) for r in body if any(r)]


# ---- states -----------------------------------------------------------------
states = load_census("gz_2010_us_040_00_20m.json")
features = sorted(({
    "type": "Feature",
    "properties": {"name": f["properties"]["NAME"], "fips": f["properties"]["STATE"]},
    "geometry": {"type": f["geometry"]["type"], "coordinates": rounded(f["geometry"]["coordinates"])}
} for f in states["features"]), key=lambda f: f["properties"]["name"])
fc = {"type": "FeatureCollection", "features": features}
with open(os.path.join(HERE, "us_states.js"), "w", encoding="utf-8") as fh:
    fh.write("// US states, DC and Puerto Rico. US Census 2010 cartographic boundaries (1:20M), via eric.clst.org.\n")
    fh.write("window.US_STATES = " + json.dumps(fc, separators=(",", ":"), ensure_ascii=False) + ";\n")
with open(os.path.join(HERE, "us_states.geojson"), "w", encoding="utf-8") as fh:
    json.dump(fc, fh, separators=(",", ":"), ensure_ascii=False)

# ---- counties --------------------------------------------------------------------
counties = load_census("gz_2010_us_050_00_20m.json")
with open(os.path.join(HERE, "us_counties.geojson"), "w", encoding="utf-8") as fh:
    json.dump(counties, fh, separators=(",", ":"), ensure_ascii=False)
state_name = {f["properties"]["STATE"]: f["properties"]["NAME"] for f in states["features"]}
county_features = [{
    "type": "Feature",
    "properties": {
        "id": f["properties"]["STATE"] + f["properties"]["COUNTY"],
        "name": f"{f['properties']['NAME']} {f['properties']['LSAD']}".strip() + f", {state_name[f['properties']['STATE']]}"
    },
    "geometry": {"type": f["geometry"]["type"], "coordinates": rounded(f["geometry"]["coordinates"], 4)}
} for f in counties["features"]]
with open(os.path.join(HERE, "us_counties.js"), "w", encoding="utf-8") as fh:
    fh.write("// US counties and equivalents, incl. Alaska, Hawaii and Puerto Rico. US Census 2010 (1:20M), via eric.clst.org.\n")
    fh.write("window.US_COUNTIES = " + json.dumps({"type": "FeatureCollection", "features": county_features},
                                                 separators=(",", ":"), ensure_ascii=False) + ";\n")

# ---- Kansas example ------------------------------------------------------------
kansas = next(f for f in features if f["properties"]["name"] == "Kansas")

# current network: from the spreadsheet if present, otherwise from the CSV copy made earlier
xlsx = os.path.join(ROOT, "ks_mesonet_stations.xlsx")
current_csv = os.path.join(HERE, "ks_mesonet_stations.csv")
if os.path.exists(xlsx):
    stations = read_xlsx(xlsx)
    with open(current_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "county", "lat", "lon", "elevation_m", "network", "abbr"])
        for s in stations:
            w.writerow([s["NAME"].strip(), s["COUNTY"], round(float(s["LATITUDE"]), 6), round(float(s["LONGITUDE"]), 6),
                        round(float(s["ELEVATION"]), 1), s["NETWORK"], s["ABBR"]])
with open(current_csv, encoding="utf-8") as fh:
    st = [[r["name"], round(float(r["lat"]), 5), round(float(r["lon"]), 5)] for r in csv.DictReader(fh)]

# 2016 network used in Patrignani et al. (2020): KS_mesonet_geoinfo.csv (56 stations)
with open(os.path.join(ROOT, "KS_mesonet_geoinfo.csv"), encoding="utf-8-sig") as fh:
    st2016 = [[r["station_name"].strip(), round(float(r["lat"]), 5), round(float(r["lon"]), 5)] for r in csv.DictReader(fh)]

with open(os.path.join(HERE, "kansas_example.js"), "w", encoding="utf-8") as fh:
    fh.write("// Kansas examples: state boundary (US Census 2010, 1:20M); Kansas Mesonet today (ks_mesonet_stations.csv)\n")
    fh.write("// and in 2016 (KS_mesonet_geoinfo.csv), the network analyzed in Patrignani et al. (2020).\n")
    fh.write("window.KANSAS_EXAMPLE = " + json.dumps({"boundary": kansas["geometry"], "stations": st, "stations2016": st2016},
                                                    separators=(",", ":"), ensure_ascii=False) + ";\n")

# ---- world countries -------------------------------------------------------------
with open(os.path.join(ROOT, "countries.geojson"), encoding="utf-8") as fh:
    world = json.load(fh)
countries = sorted(({
    "type": "Feature",
    "properties": {"name": f["properties"]["name"], "continent": f["properties"].get("continent", "")},
    "geometry": {"type": f["geometry"]["type"], "coordinates": rounded(f["geometry"]["coordinates"], 4)}
} for f in world["features"] if f.get("geometry")), key=lambda f: f["properties"]["name"].lower())
with open(os.path.join(HERE, "countries.js"), "w", encoding="utf-8") as fh:
    fh.write("// World country boundaries (countries.geojson).\n")
    fh.write("window.COUNTRIES = " + json.dumps({"type": "FeatureCollection", "features": countries},
                                               separators=(",", ":"), ensure_ascii=False) + ";\n")

print(f"{len(features)} states, {len(counties['features'])} counties, {len(st)} + {len(st2016)} (2016) Kansas stations, {len(countries)} countries")
for name in ["us_states.js", "us_states.geojson", "us_counties.geojson", "us_counties.js", "kansas_example.js", "ks_mesonet_stations.csv", "countries.js"]:
    print(f"  {name:24s} {os.path.getsize(os.path.join(HERE, name)) / 1e6:6.2f} MB")
