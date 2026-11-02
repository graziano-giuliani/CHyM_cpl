#!/usr/bin/env python3
"""
Fix routing errors in the CHyM static file (chym_stk.nc) for the MED-12
coupled configuration.

What it does
  1. Endorheic basins: the terminal-lake cell of each closed basin listed in
     SINKS becomes a sink (fdm = 9). Water routed into it leaves the model
     (it evaporates in reality) instead of being sent to the sea through a
     DEM breach. Needs the patched CHyM code (fdm = 9 handled as a sink).
  2. Mis-routed outlets: every catchment outlet inside a REROUTE box that
     drains into msk = 0 (outside the Med/Black Sea catchments) is re-routed
     along the valley floor (least-cost path on the DEM) to the nearest river
     channel of the given target river. Used for the upper Rhone: Valais,
     Lake Geneva and the Arve currently flow north out of the catchment.
  3. Recomputes the drainage area 'dra' from the edited flow directions
     (convention of the original file: upstream area, excluding the cell).
  4. Fixes the units label of 'alf' (Manning velocity, m/s, not km/h).
  5. Checks: no loops, every land cell ends in the sea, outside the
     catchment or in a sink; prints the large basins before and after.

Fields not used by the coupled model ('run', 'acc') are left untouched.

Usage:  python3 fix_chym_stk.py chym_stk.nc chym_stk_fixed.nc
"""

import heapq
import shutil
import sys

import numpy as np
from netCDF4 import Dataset

IFILE = sys.argv[1] if len(sys.argv) > 1 else "chym_stk.nc"
OFILE = sys.argv[2] if len(sys.argv) > 2 else "chym_stk_fixed.nc"

OCEAN = 15          # land-use code CHyM treats as water
SINK = 9            # flow direction meaning "terminal sink" (patched code)
ALFAMIN = 0.10      # minimum velocity used by the CHyM preprocessor (m/s)

# Closed basins: (lon, lat) of the terminal lake, search radius (deg) and the
# largest upstream area (km2) the sink cell may have. The sink is the cell
# with the largest drainage area within the radius below that limit, i.e. the
# point the routed water passes through.
SINKS = {
    "Tuz Golu + Konya closed basin": dict(lon=33.40, lat=38.80, radius=0.30, max_area=60000.0),
    "Akarcay (Eber/Aksehir lakes)":  dict(lon=31.15, lat=38.55, radius=0.25, max_area=12000.0),
    "Burdur lake":                   dict(lon=30.20, lat=37.72, radius=0.15, max_area=8000.0),
    "Beysehir lake":                 dict(lon=31.50, lat=37.75, radius=0.20, max_area=8000.0),
    "Acigol":                        dict(lon=29.85, lat=37.82, radius=0.12, max_area=3000.0),
}

# Outlets into msk = 0 inside 'box' are re-routed to the target river: any
# cell draining to the target mouth with an upstream area >= min_dra.
REROUTE = {
    "Upper Rhone (Valais, Lake Geneva, Arve)": dict(
        box=(5.6, 7.9, 45.6, 46.7), mouth=(4.80, 43.40), min_dra=1000.0),
}

IR = np.array([-1, 0, 1, 1, 1, 0, -1, -1, 0])   # as in mod_chym_param.F90
JR = np.array([1, 1, 1, 0, -1, -1, -1, 0, 0])
DIRCODE = {(int(JR[k]), int(IR[k])): k + 1 for k in range(8)}   # (dj, di) -> fdm


def load(path):
    with Dataset(path) as d:
        return {v: np.asarray(d[v][:]) for v in
                ("lon", "lat", "aer", "msk", "lus", "fdm", "dra", "dem", "alf")}


def downstream(fdm, lus):
    """Flat index of the downstream cell; -1 for non-routed cells."""
    ny, nx = fdm.shape
    J, I = np.mgrid[0:ny, 0:nx]
    ok = (lus != OCEAN) & (fdm >= 1) & (fdm <= 8)
    jn = np.where(ok, J + JR[np.clip(fdm, 1, 8) - 1], 0)
    inn = np.where(ok, I + IR[np.clip(fdm, 1, 8) - 1], 0)
    ok &= (jn >= 0) & (jn < ny) & (inn >= 0) & (inn < nx)
    return np.where(ok, jn * nx + inn, -1).ravel()


def terminal(nxt, land):
    """For each cell: its last land cell (outlet) and where it ends."""
    n = nxt.size
    step = np.where((nxt >= 0) & land[np.maximum(nxt, 0)], nxt, np.arange(n))
    out = step.copy()
    for _ in range(30):
        out = out[out]
    if np.any(step[out] != out):
        raise RuntimeError("routing loop detected")
    return out


def drainage(nxt, area):
    """Upstream area excluding the cell itself, by topological accumulation."""
    n = nxt.size
    indeg = np.bincount(nxt[nxt >= 0], minlength=n)
    acc = area.astype(np.float64).copy()
    queue = list(np.flatnonzero(indeg == 0))
    while queue:
        k = queue.pop()
        t = nxt[k]
        if t >= 0:
            acc[t] += acc[k]
            indeg[t] -= 1
            if indeg[t] == 0:
                queue.append(t)
    if np.any(indeg > 0):
        raise RuntimeError("routing loop detected")
    return acc - area


def nearest(lon, lat, lo, la):
    d = ((lon - lo) * np.cos(np.radians(la))) ** 2 + (lat - la) ** 2
    return np.unravel_index(np.argmin(d), lon.shape)


def report(title, g, nxt):
    land = (g["lus"] != OCEAN).ravel()
    out = terminal(nxt, land)
    nn = nxt[out]
    ends = np.where(nn >= 0, g["msk"].ravel()[np.maximum(nn, 0)], -1)  # -1: sink
    a = g["aer"].ravel()
    print(f"\n{title}")
    for code, name in ((1, "to Med/Black Sea"), (0, "out of the catchments"), (-1, "into sinks")):
        sel = land & (ends == code)
        print(f"   land {name:22s} {a[sel].sum():12.0f} km2")
    return out, ends


g = load(IFILE)
lon, lat = g["lon"], g["lat"]
ny, nx = g["fdm"].shape
fdm = g["fdm"].copy()
alf = g["alf"].copy()
lus = g["lus"]
land2 = lus != OCEAN
land = land2.ravel()
area = g["aer"].ravel().astype(np.float64)

nxt0 = downstream(fdm, lus)
dra0 = drainage(nxt0, area)
out0, ends0 = report("BEFORE", g, nxt0)
print(f"   max |dra recomputed - dra file| at routed cells: "
      f"{np.abs(dra0 - g['dra'].ravel())[land].max():.1f} km2")

# ----------------------------------------------------------- 1. sinks
print("\nEndorheic sinks:")
for name, s in SINKS.items():
    r = np.hypot((lon - s["lon"]) * np.cos(np.radians(s["lat"])), lat - s["lat"])
    cand = land2 & (r <= s["radius"]) & (dra0.reshape(ny, nx) <= s["max_area"])
    if not cand.any():
        print(f"   {name}: no candidate cell found, skipped")
        continue
    k = np.flatnonzero(cand.ravel())[np.argmax(dra0[cand.ravel()])]
    j, i = divmod(k, nx)
    old_out = out0[k]
    fdm[j, i] = SINK
    print(f"   {name:32s} sink at {lon[j, i]:.2f}E {lat[j, i]:.2f}N, "
          f"{dra0[k] + area[k]:8.0f} km2 cut from the basin draining at "
          f"{lon.ravel()[old_out]:.2f}E {lat.ravel()[old_out]:.2f}N")

# ----------------------------------------------------------- 2. re-routing
print("\nRe-routing outlets that leave the catchment:")
for name, rr in REROUTE.items():
    lo0, lo1, la0, la1 = rr["box"]
    inbox = (lon >= lo0) & (lon <= lo1) & (lat >= la0) & (lat <= la1)
    nxt = downstream(fdm, lus)
    out = terminal(nxt, land)
    dra = drainage(nxt, area)
    mj, mi = nearest(lon, lat, *rr["mouth"])
    mouth = mj * nx + mi
    target = (land & (out == mouth) & (dra >= rr["min_dra"])).reshape(ny, nx) & inbox
    msk_flat = g["msk"].ravel()
    outlets = [k for k in np.flatnonzero(inbox.ravel() & land)
               if nxt[k] >= 0 and not land[nxt[k]] and msk_flat[nxt[k]] == 0]
    outlets.sort(key=lambda k: -dra[k])
    dem = g["dem"].astype(np.float64)
    for k0 in outlets:
        if dra[k0] + area[k0] < 50.0:          # tiny edge cells: leave them
            continue
        # Dijkstra over land cells in the box; cost = distance x (1 + dem/200)
        start = divmod(k0, nx)
        best = {start: 0.0}
        prev = {}
        heap = [(0.0, start)]
        goal = None
        while heap:
            c, (j, i) = heapq.heappop(heap)
            if target[j, i] and (j, i) != start:
                goal = (j, i)
                break
            if c > best[(j, i)]:
                continue
            for dj in (-1, 0, 1):
                for di in (-1, 0, 1):
                    jj, ii = j + dj, i + di
                    if (dj or di) and inbox[jj, ii] and land2[jj, ii] and fdm[jj, ii] != SINK:
                        step = np.hypot(dj, di) * (1.0 + dem[jj, ii] / 200.0)
                        nc = c + step
                        if nc < best.get((jj, ii), np.inf):
                            best[(jj, ii)] = nc
                            prev[(jj, ii)] = (j, i)
                            heapq.heappush(heap, (nc, (jj, ii)))
        if goal is None:
            print(f"   {name}: no path from {lon.ravel()[k0]:.2f}E {lat.ravel()[k0]:.2f}N")
            continue
        path = [goal]
        while path[-1] != start:
            path.append(prev[path[-1]])
        path = path[::-1]
        for (ja, ia), (jb, ib) in zip(path[:-1], path[1:]):
            fdm[ja, ia] = DIRCODE[(jb - ja, ib - ia)]
            alf[ja, ia] = max(alf[ja, ia], ALFAMIN)
            target[ja, ia] = True                 # later outlets may join this path
        print(f"   outlet {lon[start]:.2f}E {lat[start]:.2f}N ({dra[k0] + area[k0]:6.0f} km2) "
              f"-> {len(path) - 1} cells -> joins at {lon[goal]:.2f}E {lat[goal]:.2f}N")

# ----------------------------------------------------------- 3. checks, dra
nxt1 = downstream(fdm, lus)
dra1 = drainage(nxt1, area)                      # raises on loops
g1 = dict(g, fdm=fdm)
out1, ends1 = report("AFTER", g1, nxt1)

print("\nLargest basins draining to the Med/Black Sea (km2):")
mouths = np.unique(out1[land & (ends1 == 1)])
for k in mouths[np.argsort(dra1[mouths])[::-1]][:15]:
    print(f"   {lon.ravel()[k]:6.2f}E {lat.ravel()[k]:5.2f}N  before {g['dra'].ravel()[k]:9.0f}"
          f"  after {dra1[k]:9.0f}")

# ----------------------------------------------------------- 4. write
shutil.copyfile(IFILE, OFILE)
with Dataset(OFILE, "r+") as d:
    d["fdm"][:] = fdm
    d["fdm"].setncattr("comment", "1-8: flow direction; 9: terminal sink (endorheic)")
    d["dra"][:] = dra1.reshape(ny, nx).astype(np.float32)
    d["alf"][:] = alf.astype(np.float32)
    d["alf"].setncattr("units", "m s-1")
    d.history = (f"{' '.join(sys.argv)}: endorheic sinks, upper Rhone re-routed, "
                 f"dra recomputed\n" + d.history)
print(f"\nwrote {OFILE}")
