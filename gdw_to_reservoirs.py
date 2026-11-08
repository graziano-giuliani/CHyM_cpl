#!/usr/bin/env python3
"""
Build the CHyM reservoir table (reservoirs.txt) from the Global Dam Watch
barrier database (GDW_barriers_v1_0, CSV text export).

For every GDW dam with at least MIN_CAP_MCM of storage:
  * the dam is placed on a CHyM cell: the routed land cell within 25 km
    minimising |log(CHyM upstream area / GDW catchment area)| + DIST_W * km,
    i.e. the best area match, nearer cells preferred on near-ties. The
    upa_km2 column then holds that cell's CHyM upstream area (the GDW value
    is kept in the comment), so mod_reservoir picks exactly the same cell;
  * the cell is followed downstream; only dams draining to the coupled sea
    (msk = 1: Mediterranean, Marmara, Black Sea, Azov) are kept;
  * if two dams fall on the same CHyM cell, the larger one is kept;
  * dams in the basins listed in EXCLUDE are written commented out, and so
    are dams for which no cell within 25 km has an upstream area within 30 %
    of the GDW catchment (to be placed by hand) and dams without a
    construction year in GDW (to be dated by hand).

Columns written (see mod_reservoir.F90):
  name lon lat upa_km2 cap_km3 surf_km2 evap_mm_yr use_m3s year0 opmon qmean_m3s
from GDW:  lon/lat = LONG_RIV/LAT_RIV (position on the river network,
           LONG_DAM/LAT_DAM if missing), upa = CATCH_SKM, cap = CAP_MCM/1000,
           surf = AREA_SKM, year0 = YEAR_DAM, qmean = DIS_AVG_LS/1000
           (long-term naturalized mean discharge, HydroRIVERS).
not in GDW: evap_mm_yr (net open-water evaporation; EVAP_MED / EVAP_BS
           below), use_m3s (consumptive withdrawal, 0), opmon (first month of
           the operational year: OPMON_BS for snowmelt-fed Black Sea rivers,
           OPMON_MED for rain-fed Mediterranean ones).

Usage:
  python3 gdw_to_reservoirs.py GDW_barriers_v1_0.txt chym_stk.nc reservoirs.txt
          [--min-cap 1000] [--only NAME1,NAME2,...]
"""

import argparse
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from netCDF4 import Dataset

OCEAN = 15
SEARCH_KM = 25.0
DIST_W = 0.02                      # area-mismatch units per km
EVAP_BS, EVAP_MED = 400, 800       # mm/yr, net evaporation (rough defaults)
OPMON_BS, OPMON_MED = 3, 10        # before snowmelt / before the wet season
# Basins left in the table but commented out: (label, lon, lat of the mouth)
EXCLUDE = [("Danube", 29.7, 45.2)]

IR = np.array([-1, 0, 1, 1, 1, 0, -1, -1, 0])
JR = np.array([1, 1, 1, 0, -1, -1, -1, 0, 0])


def gc_km(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lo2 - lo1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(np.minimum(a, 1.0)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gdw")
    ap.add_argument("stk")
    ap.add_argument("out")
    ap.add_argument("--min-cap", type=float, default=1000.0, help="minimum capacity, MCM")
    ap.add_argument("--only", default="", help="comma-separated GDW DAM_NAMEs to keep")
    a = ap.parse_args()

    d = pd.read_csv(a.gdw, low_memory=False, encoding="latin-1").replace(-99, np.nan)
    d["lon"] = d["LONG_RIV"].fillna(d["LONG_DAM"])
    d["lat"] = d["LAT_RIV"].fillna(d["LAT_DAM"])
    d = d[(d["CAP_MCM"] >= a.min_cap) & d["lon"].notna()]
    if a.only:
        names = [n.strip().lower() for n in a.only.split(",")]
        d = d[d["DAM_NAME"].fillna("").str.lower().isin(names)]

    with Dataset(a.stk) as s:
        lon = np.asarray(s["lon"][:], float)
        lat = np.asarray(s["lat"][:], float)
        lus = np.asarray(s["lus"][:])
        msk = np.asarray(s["msk"][:])
        fdm = np.asarray(s["fdm"][:])
        up = np.asarray(s["dra"][:], float) + np.asarray(s["aer"][:], float)
    ny, nx = fdm.shape
    routed = (lus != OCEAN) & (fdm >= 1) & (fdm <= 8)
    inner = np.zeros_like(routed)
    inner[1:-1, 1:-1] = True
    routed &= inner
    d = d[(d["lon"] >= lon.min()) & (d["lon"] <= lon.max()) &
          (d["lat"] >= lat.min()) & (d["lat"] <= lat.max())]

    def outlet(j, i):
        """Follow the flow to the last land cell; return it and the sea mask."""
        for _ in range(ny * nx):
            f = fdm[j, i]
            if lus[j, i] == OCEAN or not (1 <= f <= 8):
                return (j, i), None           # sink
            jj, ii = j + JR[f - 1], i + IR[f - 1]
            if lus[jj, ii] == OCEAN:
                return (j, i), msk[jj, ii]
            j, i = jj, ii
        return (j, i), None

    excl = []
    for label, elo, ela in EXCLUDE:
        dd = gc_km(ela, elo, lat, lon)
        dd = np.where(routed, dd, np.inf)
        cand = dd <= 50
        j, i = np.unravel_index(np.argmax(np.where(cand, up, -1)), up.shape)
        excl.append((label, outlet(j, i)[0]))

    rows, seen = [], {}
    for _, r in d.iterrows():
        dist = gc_km(r["lat"], r["lon"], lat, lon)
        cand = (dist <= SEARCH_KM) & routed
        if not cand.any():
            continue
        if np.isfinite(r["CATCH_SKM"]) and r["CATCH_SKM"] > 0:
            score = np.where(cand, np.abs(np.log(up / r["CATCH_SKM"])) + DIST_W * dist, np.inf)
        else:
            score = np.where(cand, -up + 1e3 * dist, np.inf)
        j, i = np.unravel_index(np.argmin(score), score.shape)
        (oj, oi), sea = outlet(j, i)
        if sea != 1:
            continue                           # not draining to the coupled sea
        blacksea = lon[oj, oi] > 27.3 and lat[oj, oi] > 40.95 and not (
            lon[oj, oi] < 29.3 and lat[oj, oi] < 41.2)
        basin = next((lab for lab, o in excl if o == (oj, oi)), "")
        label = next((str(x) for x in (r["DAM_NAME"], r["RES_NAME"], r["RIVER"])
                      if isinstance(x, str) and x.strip()),
                     f"unnamed_{r['COUNTRY']}" if isinstance(r["COUNTRY"], str) else "unnamed")
        name = re.sub(r"[^A-Za-z0-9_]+", "_", label).strip("_")
        row = dict(
            name=f"{name}_{int(r['GDW_ID'])}", lon=r["lon"], lat=r["lat"],
            upa=up[j, i], gupa=r["CATCH_SKM"],
            cap=r["CAP_MCM"] / 1000.0,
            surf=r["AREA_SKM"] if np.isfinite(r["AREA_SKM"]) else 0.0,
            evap=EVAP_BS if blacksea else EVAP_MED, use=0.0,
            year=int(r["YEAR_DAM"]) if np.isfinite(r["YEAR_DAM"]) else 1900,
            opmon=OPMON_BS if blacksea else OPMON_MED,
            qmean=r["DIS_AVG_LS"] / 1000.0 if np.isfinite(r["DIS_AVG_LS"]) else 0.0,
            river=str(r["RIVER"]) if isinstance(r["RIVER"], str) else "",
            use_main=str(r["MAIN_USE"]) if isinstance(r["MAIN_USE"], str) else "",
            dor=r["DOR_PC"], sea="Black Sea" if blacksea else "Mediterranean",
            cell=(j, i), cupa=up[j, i], excluded=basin, noyear=not np.isfinite(r["YEAR_DAM"]))
        row["dist"] = dist[j, i]
        if (j, i) in seen:
            k = seen[(j, i)]
            if row["cap"] > rows[k]["cap"]:
                rows[k] = row
            continue
        seen[(j, i)] = len(rows)
        rows.append(row)

    rows.sort(key=lambda r: (r["sea"], r["river"], -r["upa"]))
    for r in rows:
        if np.isfinite(r["gupa"]) and abs(r["upa"] / r["gupa"] - 1) > 0.3 and not r["excluded"]:
            # no CHyM cell within 25 km has a comparable upstream area:
            # written commented out, to be placed by hand
            r["excluded"] = (f"CHECK CELL: CHyM upa {r['upa']:.0f} vs GDW {r['gupa']:.0f}")
            print(f"commented out {r['name']}: CHyM upa {r['upa']:.0f} vs GDW {r['gupa']:.0f}")
        elif r["noyear"] and not r["excluded"]:
            # no construction year (mostly satellite-detected GOODD/GROD
            # entries): year0 = 1900 would make it active all the time
            r["excluded"] = "CHECK: NO YEAR IN GDW, set year0 by hand"
            print(f"commented out {r['name']}: no construction year in GDW")
    with open(a.out, "w") as f:
        f.write("# CHyM reservoir table generated from Global Dam Watch GDW_barriers_v1_0\n"
                f"# dams with capacity >= {a.min_cap:.0f} MCM draining to the Mediterranean or Black Sea\n"
                "# from GDW: lon lat (on the river network) cap_km3 surf_km2 year0 qmean_m3s\n"
                "# upa_km2: CHyM upstream area of the chosen cell (GDW catchment in the comment)\n"
                f"# defaults, NOT from GDW: evap_mm_yr ({EVAP_BS} Black Sea, {EVAP_MED} Mediterranean),\n"
                f"#   use_m3s 0, opmon ({OPMON_BS} Black Sea, {OPMON_MED} Mediterranean)\n"
                "# year0 = 1900 where GDW gives no year; lines starting with '#' are ignored\n"
                "#\n# name lon lat upa_km2 cap_km3 surf_km2 evap_mm_yr use_m3s year0 opmon qmean_m3s\n")
        for r in rows:
            note = (f"   # {r['river']}, {r['sea']}, main use {r['use_main'] or '?'}, "
                    f"DOR {r['dor']:.0f}%, GDW catchment {r['gupa']:.0f} km2, "
                    f"cell {r['dist']:.1f} km away")
            line = (f"{r['name']:32s} {r['lon']:8.3f} {r['lat']:7.3f} {r['upa']:8.0f} "
                    f"{r['cap']:7.3f} {r['surf']:8.1f} {r['evap']:5.0f} {r['use']:4.0f} "
                    f"{r['year']:5d} {r['opmon']:3d} {r['qmean']:8.1f}")
            if r["excluded"]:
                tag = r["excluded"] if r["excluded"].startswith("CHECK") else f"{r['excluded']} basin, left out"
                line = "# " + line + f"   [{tag}]"
            f.write(line + note + "\n")
    n_act = sum(1 for r in rows if not r["excluded"])
    print(f"{len(rows)} dams written to {a.out} ({n_act} active, "
          f"{len(rows) - n_act} commented out)")


if __name__ == "__main__":
    main()
