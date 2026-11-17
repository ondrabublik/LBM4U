"""Builds a cLBMnet (LBM4U) case from a JSON description, without LBM Studio.

usage: python make_case.py <case.json>            write the case and print the derived quantities
       python make_case.py <case.json> --check    only print the derived quantities (tau, Mach, dt, ...)

The JSON (all lengths in metres, SI units):
{
  "geometry_dir": "simulations/damBreak",        # gets mesh/blocks.txt and mesh/boundary/
  "simulation": "sim1",                          # simulations/damBreak/sim1/{parameters.txt, init/}
  "blocks": [4, 2, 4],                           # level-0 blocks in x, y, z
  "cells": 16,                                   # ni = nj = nk (even if refined)
  "h": 0.0005,                                   # cell size of level 0 [m]; H = cells * h is written for you
  "refine": [{"box": [0, 0, 0, 0.016, 0.016, 0.008], "level": 1}],   # optional, 2:1 balance is enforced
  "walls": "z < 0.002",                          # optional: python/numpy expression in x, y, z -> wall cells (type 1)
  "inlet":  {"where": "(z > 0.03) & ((x-0.01)**2 + (y-0.004)**2 < 0.002**2)"},  # optional, type 2 cells
  "outlet": {"where": "x > 0.0315"},             # optional, type 3 cells
  "fluid": "(x < 0.01) & (z < 0.02)",            # liquid region (flag 1), the rest is gas (flag 2); "true" = full
  "velocity": [0, 0, 0],                         # optional initial velocity of the liquid [m/s]
  "params": {"problemID": 1, "uMax": 0.6, ...},  # anything for parameters.txt (H, ni, nj, nk are added)
  "particles": {"sphere": [x, y, z, R], "d": 1e-4, "rho": 2500, "phi": 0.005, "seed": 1},   # problemID 2
  # problemID 22 (ellipsoids): "shape": [a1, a2, a3] (semi-axes [m]) instead of "d", optional
  # "orientation": "random" (default) or [qw, qx, qy, qz] and "drag_factors": [f1, f2, f3]
  "concentration": {"sphere": [x, y, z, R], "d": 1e-4, "phi": 0.005}                      # problemID 3, 4, 33
  # problemID 33 (spheroids, params particleAspectRatio): optional "orientation": "random" (default, A = I/3),
  # [px, py, pz] (all particles along this axis) or [A11, A22, A33, A12, A13, A23] -> init/orientation*.txt
}
"particles" may also be {"list": [[x, y, z, d, rho], ...]}; "sphere" may be replaced by "box": [x0,y0,z0,x1,y1,z1]
or "where": "<expression>". The domain border needs no wall cells: a missing neighbour block is a bounce-back wall.
"""
import json
import math
import os
import sys

import numpy as np


def die(msg):
    sys.exit("make_case: " + msg)


def leaves(nb, cells, h0, refine):
    """octree leaves (bi, bj, bk, level), 2:1 balanced over faces, edges and corners"""
    lmax = max([r["level"] for r in refine], default=0)
    target = {}  # requested level of each block (l, i, j, k)

    def wanted(l, i, j, k):
        size = cells * h0 / (1 << l)
        lo = np.array([i, j, k]) * size
        hi = lo + size
        need = 0
        for r in refine:
            b = r["box"]
            if all(lo[d] < b[d + 3] and hi[d] > b[d] for d in range(3)):
                need = max(need, r["level"])
        return need

    out = set()

    def add(l, i, j, k):
        if l < lmax and wanted(l, i, j, k) > l:
            for o in range(8):
                add(l + 1, 2 * i + (o & 1), 2 * j + ((o >> 1) & 1), 2 * k + (o >> 2))
        else:
            out.add((i, j, k, l))

    for k in range(nb[2]):
        for j in range(nb[1]):
            for i in range(nb[0]):
                add(0, i, j, k)
    if lmax == 0:
        return sorted(out, key=lambda b: (b[2], b[1], b[0]))
    # balance: level map on the finest block grid, split blocks touching a block two levels finer
    n = [x << lmax for x in nb]
    while True:
        lev = np.zeros(n, dtype=np.int8)
        for (i, j, k, l) in out:
            s = 1 << (lmax - l)
            lev[i * s:(i + 1) * s, j * s:(j + 1) * s, k * s:(k + 1) * s] = l
        split = []
        for (i, j, k, l) in out:
            s = 1 << (lmax - l)
            sl = tuple(slice(max(c * s - 1, 0), min((c + 1) * s + 1, n[d])) for d, c in enumerate((i, j, k)))
            if lev[sl].max() > l + 1:
                split.append((i, j, k, l))
        if not split:
            break
        for (i, j, k, l) in split:
            out.remove((i, j, k, l))
            for o in range(8):
                out.add((2 * i + (o & 1), 2 * j + ((o >> 1) & 1), 2 * k + (o >> 2), l + 1))
    return sorted(out, key=lambda b: (b[3], b[2], b[1], b[0]))


def evaluate(expr, x, y, z):
    if expr is None:
        return np.zeros(x.shape, dtype=bool)
    if isinstance(expr, bool) or str(expr).strip().lower() == "true":
        return np.ones(x.shape, dtype=bool) if expr in (True, "true", "True") else np.zeros(x.shape, dtype=bool)
    env = {"x": x, "y": y, "z": z, "np": np, "pi": math.pi, "sqrt": np.sqrt, "abs": np.abs}
    r = eval(expr, {"__builtins__": {}}, env)
    return np.broadcast_to(np.asarray(r, dtype=bool), x.shape)


def orientation_tensor(o):
    """orientation tensor A11 A22 A33 A12 A13 A23 of the Eulerian spheroids, None for random orientation"""
    if o in (None, "random"):
        return None
    o = [float(v) for v in o]
    if len(o) == 3:
        p = np.array(o) / np.linalg.norm(o)
        return np.array([p[0] * p[0], p[1] * p[1], p[2] * p[2], p[0] * p[1], p[0] * p[2], p[1] * p[2]])
    if len(o) == 6:
        A = np.array(o)
        return A / (A[0] + A[1] + A[2])
    die("concentration orientation: 'random', [px, py, pz] or [A11, A22, A33, A12, A13, A23]")


def region(spec, x, y, z):
    if "sphere" in spec:
        cx, cy, cz, R = spec["sphere"]
        return (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 < R * R
    if "box" in spec:
        b = spec["box"]
        return (x >= b[0]) & (x < b[3]) & (y >= b[1]) & (y < b[4]) & (z >= b[2]) & (z < b[5])
    return evaluate(spec.get("where"), x, y, z)


def derived(c):
    """physical <-> lattice scales as in Block::initCommon, warnings for typical problems"""
    p = c["params"]
    h0 = c["h"]
    uMax = float(p["uMax"])
    Cu = uMax / 0.1
    Ct = h0 / Cu
    nu = float(p["vKin"])
    tau = 0.5 + 3 * nu * Ct / h0 ** 2
    g = [float(p.get(a, 0)) for a in ("gx", "gy", "gz")]
    gl = math.sqrt(sum(v * v for v in g)) * Ct * Ct / h0
    lmax = max([r["level"] for r in c.get("refine", [])], default=0)
    steps = int(float(p.get("maxIter", 0)))
    print("cell size h0 = %.4g m, H = %.4g m, Cu = %.4g m/s, dt (level 0) = %.4g s" % (h0, c["cells"] * h0, Cu, Ct))
    print("maxIter = %d -> physical time %.4g s, saveRate %s" % (steps, steps * Ct, p.get("saveRate", "-")))
    for l in range(lmax + 1):
        print("level %d: h = %.4g m, tau = %.4f" % (l, h0 / (1 << l), 0.5 + (tau - 0.5) * (1 << l)))
    print("lattice gravity |g| dt^2/h = %.3g (keep below ~1e-3 .. 5e-3)" % gl)
    warn = []
    if tau < 0.51 and float(p.get("smagorinsky", 0)) == 0:
        warn.append("tau = %.4f close to 0.5: use smagorinsky = 0.1-0.17 and tauMinusMax = 1, or larger uMax/h" % tau)
    if tau > 2:
        warn.append("tau = %.2f is large (inaccurate walls, slow): reduce uMax or refine" % tau)
    if gl > 5e-3:
        warn.append("lattice gravity %.3g is large: hydrostatic density change 3 g L / Cu^2 may exceed a few %%" % gl)
    height = c["blocks"][2] * c["cells"] * h0
    drho = 3 * math.sqrt(sum(v * v for v in g)) * height / Cu ** 2 * float(p.get("fluidGravity", 1))
    print("hydrostatic density variation over the domain height 3 g Lz / Cu^2 = %.3g" % drho)
    if drho > 0.1:
        warn.append("density variation %.2f > 0.1: increase uMax (Cu) or reduce the height" % drho)
    vm = int(p.get("viscosityModel", 0))
    if vm:
        tmin, tmax = float(p.get("tauMin", 0.5001)), float(p.get("tauMax", 1.5))
        Cvis = h0 * h0 / Ct
        print("non-Newtonian: representable kinematic viscosity on level 0: %.3g .. %.3g m^2/s (tauMin %.4g, tauMax %.4g)"
              % ((tmin - 0.5) / 3 * Cvis, (tmax - 0.5) / 3 * Cvis, tmin, tmax))
        for l in range(1, lmax + 1):
            Cl = Cvis / (1 << l)
            print("   level %d: %.3g .. %.3g m^2/s" % (l, (tmin - 0.5) / 3 * Cl, (tmax - 0.5) / 3 * Cl))
    for w in warn:
        print("warning: " + w)


def main():
    if len(sys.argv) < 2:
        die(__doc__)
    with open(sys.argv[1], encoding="utf-8") as f:
        c = json.load(f)
    check_only = "--check" in sys.argv
    for key in ("geometry_dir", "simulation", "blocks", "cells", "h", "params"):
        if key not in c:
            die("missing '%s' in %s" % (key, sys.argv[1]))
    n = int(c["cells"])
    h0 = float(c["h"])
    p = dict(c["params"])
    p.update({"ni": n, "nj": n, "nk": n, "H": repr(n * h0)})
    c["params"] = p
    for req in ("problemID", "model", "uMax", "density", "vKin", "maxIter"):
        if req not in p:
            die("params needs '%s'" % req)
    refine = c.get("refine", [])
    if refine and n % 2:
        die("refinement needs an even number of cells per block")
    derived(c)
    if check_only:
        return

    geom = c["geometry_dir"]
    sim = os.path.join(geom, c["simulation"])
    os.makedirs(os.path.join(geom, "mesh", "boundary"), exist_ok=True)
    os.makedirs(os.path.join(sim, "init"), exist_ok=True)
    blocks = leaves(c["blocks"], n, h0, refine)
    nb = c["blocks"]
    with open(os.path.join(geom, "mesh", "blocks.txt"), "w") as f:
        for (i, j, k, l) in blocks:
            f.write("%d %d %d %d\n" % (i, j, k, l) if refine else "%d %d %d\n" % (i, j, k))
    with open(os.path.join(sim, "parameters.txt"), "w") as f:
        for kv in p.items():
            f.write("%s = %s\n" % kv)

    # cell centres of one block in the file order (i outer, k fastest)
    ii, jj, kk = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij")
    ii, jj, kk = ii.ravel(), jj.ravel(), kk.ravel()
    vel = c.get("velocity", [0, 0, 0])
    conc = c.get("concentration")
    nWall = nIn = nOut = nLiq = 0
    total_conc = 0.0
    centres = []
    for (bi, bj, bk, l) in blocks:
        h = h0 / (1 << l)
        x = (bi * n + ii + 0.5) * h
        y = (bj * n + jj + 0.5) * h
        z = (bk * n + kk + 0.5) * h
        ct = np.zeros(x.shape, dtype=int)
        ct[evaluate(c.get("walls"), x, y, z)] = 1
        if "inlet" in c:
            ct[region(c["inlet"], x, y, z) & (ct == 0)] = 2
        if "outlet" in c:
            ct[region(c["outlet"], x, y, z) & (ct == 0)] = 3
        tag = ("L%d_%d_%d_%d" % (l, bi, bj, bk)) if l > 0 else str(nb[0] * nb[1] * bk + nb[0] * bj + bi)
        bpath = os.path.join(geom, "mesh", "boundary", "blockBoundary%s.txt" % tag)
        if ct.any():
            np.savetxt(bpath, ct, fmt="%d")
        elif os.path.exists(bpath):
            os.remove(bpath)
        liquid = evaluate(c.get("fluid", "true"), x, y, z) & (ct == 0)
        flag = np.where(liquid, 1, 2)
        u = np.where(liquid, 1.0, 0.0)
        data = np.column_stack([np.ones(x.shape), u * vel[0], u * vel[1], u * vel[2], flag])
        np.savetxt(os.path.join(sim, "init", "blockInit%s.txt" % tag), data, fmt=["%.0f", "%.10g", "%.10g", "%.10g", "%d"])
        if conc:
            vp = math.pi / 6 * conc["d"] ** 3
            cc = np.where(region(conc, x, y, z) & liquid, conc["phi"] * h ** 3 / vp, 0.0)
            total_conc += cc.sum()
            np.savetxt(os.path.join(sim, "init", "concentration%s.txt" % tag), cc, fmt="%.10g")
            A = orientation_tensor(conc.get("orientation", "random"))
            opath = os.path.join(sim, "init", "orientation%s.txt" % tag)
            if A is not None:
                np.savetxt(opath, np.tile(A, (x.size, 1)), fmt="%.10g")
            elif os.path.exists(opath):
                os.remove(opath)
        nWall += (ct == 1).sum()
        nIn += (ct == 2).sum()
        nOut += (ct == 3).sum()
        nLiq += liquid.sum()
        centres.append((l, h))
    print("blocks: %d (%s), liquid cells %d, wall %d, inlet %d, outlet %d"
          % (len(blocks), ", ".join("level %d: %d" % (L, sum(1 for b in blocks if b[3] == L)) for L in sorted({b[3] for b in blocks})),
             nLiq, nWall, nIn, nOut))
    if conc:
        print("Eulerian particles: %.6g in total (particleDiameter/particleDensity in params must match d = %g)" % (total_conc, conc["d"]))

    part = c.get("particles")
    if part:
        if "list" in part:
            pts = [list(map(float, q)) for q in part["list"]]
        else:
            rng = np.random.default_rng(part.get("seed", 1))
            L = [nb[d] * n * h0 for d in range(3)]
            shape = part.get("shape")
            vp = 4.0 / 3.0 * math.pi * shape[0] * shape[1] * shape[2] if shape else math.pi / 6 * part["d"] ** 3
            # sample candidates in the bounding box of the domain, keep those in the region
            cand = rng.random((200000, 3)) * L
            inside = region(part, cand[:, 0], cand[:, 1], cand[:, 2])
            vol = inside.mean() * L[0] * L[1] * L[2]
            count = int(round(part["phi"] * vol / vp)) if "count" not in part else int(part["count"])
            sel = cand[inside]
            while len(sel) < count:
                more = rng.random((200000, 3)) * L
                sel = np.vstack([sel, more[region(part, more[:, 0], more[:, 1], more[:, 2])]])
            if shape:
                ori = part.get("orientation", "random")
                fac = list(part.get("drag_factors", [1, 1, 1]))
                pts = []
                for q in sel[:count]:
                    if ori == "random":  # uniform random rotation (Shoemake)
                        u1, u2, u3 = rng.random(3)
                        qq = [math.sqrt(u1) * math.cos(2 * math.pi * u3), math.sqrt(1 - u1) * math.sin(2 * math.pi * u2),
                              math.sqrt(1 - u1) * math.cos(2 * math.pi * u2), math.sqrt(u1) * math.sin(2 * math.pi * u3)]
                    else:
                        qq = list(ori)
                    pts.append([q[0], q[1], q[2], shape[0], shape[1], shape[2], part["rho"]] + qq + fac)
            else:
                pts = [[q[0], q[1], q[2], part["d"], part["rho"]] for q in sel[:count]]
        with open(os.path.join(sim, "init", "particles.txt"), "w") as f:
            f.write("%d\n" % len(pts))
            for q in pts:
                f.write(" ".join("%.10g" % v for v in q) + "\n")
        print("Lagrangian particles: %d -> init/particles.txt" % len(pts))
    print("written: %s" % sim)


if __name__ == "__main__":
    main()
