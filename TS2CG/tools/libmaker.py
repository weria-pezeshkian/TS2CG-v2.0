"""
Helper tool to generate a lib file section
"""

import argparse
from pathlib import Path
import numpy as np
import logging
from dataclasses import dataclass
from typing import List, Optional, Sequence, Dict, Tuple, Set, Protocol
import MDAnalysis as mda
import networkx as nx
import random

logger = logging.getLogger(__name__)

class ProtocolVirtualSite3(Protocol):
    vid: int         
    funct: int
    a1: int
    a2: int
    a3: int

def write_file(universe: mda.Universe, filename: str = 'lib.txt',name:str="lipid") -> str:
    """
    Write a simple DAPC-style text file from an MDAnalysis Universe.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        The MDAnalysis Universe containing your atoms.
    filename : str
        Path where the output text file will be written.

    Returns
    -------
    str
        The filename of the written file.
    """
    lines = [f"[{name}]"]
    for i, atom in enumerate(universe.atoms, start=1):
        x, y, z = atom.position
        lines.append(f"{i} {atom.name} {x:.1f} {y:.1f} {z:.1f}")

    text = "\n".join(lines) + "\n"
    with open(filename, 'w') as f:
        f.write(text)

def stretch_positions(pos: np.ndarray, z_factor: float, xy_factor: float) -> np.ndarray:
    """
    Stretch a set of 3D positions.

    Parameters
    ----------
    pos : np.ndarray, shape (N,3)
        Original coordinates.
    z_factor : float
        Multiplicative factor to apply to all z-coordinates.
    xy_factor : float
        Multiplicative factor to apply to all x- and y-coordinates.

    Returns
    -------
    np.ndarray, shape (N,3)
        New coordinates with:
         - x -> x * xy_factor
         - y -> y * xy_factor
         - z -> z * z_factor
    """
    # Make a copy so we don’t overwrite the input
    stretched = pos.copy()
    # Stretch z
    stretched[:, 2] *= z_factor
    # Stretch x and y
    stretched[:, :2] *= xy_factor
    return stretched

def rotate_to_z(v):
    """
    Build a rotation matrix that rotates vector v to align with the positive z-axis.
    """
    v = v / np.linalg.norm(v)
    z = np.array([0.0, 0.0, 1.0])

    axis = np.cross(v, z)
    axis_len = np.linalg.norm(axis)
    if axis_len < 1e-8:
        return np.eye(3)
    axis = axis / axis_len

    angle = np.arccos(np.clip(np.dot(v, z), -1.0, 1.0))

    K = np.array([[    0,     -axis[2],  axis[1]],
                  [ axis[2],     0,     -axis[0]],
                  [-axis[1],  axis[0],     0   ]])
    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)
    return R

def maker_gro(file: str, flipper: bool=False, sz: float=.5, sxy: float=.5, base: str="PO4",output: str="lib.txt") -> None:
    u = mda.Universe(file)
    resname=u.residues[0].resname
    pos = u.atoms.positions.copy()  # (N, 3)

    am2 = u.select_atoms(f"name {base}")
    if len(am2) != 1:
        raise ValueError("Did not find exactly one atom named AM2.")
    am2_pos = am2.positions[0]
    pos -= am2_pos

    cov = np.cov(pos.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    # sort eigenvectors by descending eigenvalue
    order = np.argsort(eigvals)[::-1]
    pc1 = eigvecs[:, order[0]]  # first principal component

    R = rotate_to_z(pc1)

    pos = pos.dot(R.T)

    pos=stretch_positions(pos,sz,sxy)
    if flipper:
        pos[:,2]*=-1
    u.atoms.positions = pos
    write_file(u,output,resname)
    print(f"Wrote lib file entry into {output}")

def mol_to_graph(mol, constraints: bool = False):
    G = nx.Graph()

    # keep only non-hydrogen atoms as graph nodes
    heavy = {a.nr for a in mol.atoms if not a.atom.upper().startswith("H")}

    # nodes
    for a in mol.atoms:
        if a.nr not in heavy:
            continue
        G.add_node(a.nr, name=a.atom, atype=a.atype, charge=a.charge)

    # bonds as edges
    for b in getattr(mol, "bonds", []):
        if b.ai in heavy and b.aj in heavy:
            G.add_edge(b.ai, b.aj)

    # optionally treat constraints as edges too
    if constraints:
        for c in getattr(mol, "constraints", []):
            if c.ai in heavy and c.aj in heavy:
                G.add_edge(c.ai, c.aj)

    return G

def straight_roads_from_source(G: nx.Graph, source: int) -> List[List[int]]:
    """
    For each neighbor of source, walk forward while the current node has exactly one
    'forward' neighbor (degree-2 interior along that direction).
    """
    roads: List[List[int]] = []
    for nb in G.neighbors(source):
        road: List[int] = []
        prev, cur = source, nb
        while True:
            road.append(cur)
            nxts = [x for x in G.neighbors(cur) if x != prev]
            if len(nxts) != 1:
                break
            prev, cur = cur, nxts[0]
        roads.append(road)
    return roads

def decide_road(mol, road: List[int], *, seed: int | None = None) -> Tuple[str, int, bool]:
    """
    Decide how to place one road starting from the base bead.

    Rules:
    1) If there are >=3 continuous carbon-prefix atom types (atype startswith 'C')
       at the beginning of the road => direction = "neg"
    2) If a charge exists:
         - if not neg by rule 1 => direction = "pos"
         - if neg by rule 1 => keep neg and shift so first charged atom is at z=0
    3) Fallback: randomize ("neg" or "pos")

    Returns:
      (direction, z_shift_steps, radial)
        - direction: "neg" or "pos"
        - z_shift_steps: integer dz steps to shift the whole road upward
        - radial: True if the road should be placed off-axis
    """
    atoms_by_nr = {a.nr: a for a in mol.atoms}
    rng = random.Random(seed)
    eps = 1e-12

    # --- rule 1: continuous carbon prefix ---
    c_prefix = 0
    for n in road:
        if str(atoms_by_nr[n].atype).startswith("C"):
            c_prefix += 1
        else:
            break

    carbon_neg = c_prefix >= 3

    # --- find first charged atom ---
    first_charged_step = None
    for step, n in enumerate(road, start=1):
        if abs(float(atoms_by_nr[n].charge)) > eps:
            first_charged_step = step
            break

    direction = "neg" if carbon_neg else None
    z_shift_steps = 0
    radial = False

    # --- rule 2: charge handling ---
    if first_charged_step is not None:
        radial = True
        if direction is None:
            direction = "pos"
        else:
            # negative due to carbon prefix: shift upward so charge lands at z=0
            z_shift_steps = first_charged_step

    # --- rule 4: fallback ---
    if direction is None:
        direction = rng.choice(["neg", "pos"])

    return direction, z_shift_steps, radial

def walk_straight_chain(G: nx.Graph, fork: int, first: int, blocked: Set[int]) -> List[int]:
    """
    Walk a 'straight' branch chain starting at `first` coming from `fork`,
    stopping if:
      - dead end (no forward candidates)
      - junction (more than one forward candidate)
      - next step would enter `blocked`
    """
    chain: List[int] = []
    prev, cur = fork, first

    while True:
        if cur in blocked:
            break
        chain.append(cur)
        blocked.add(cur)  # mark as used immediately to avoid collisions across branches

        candidates = [x for x in G.neighbors(cur) if x != prev and x not in blocked]
        if len(candidates) != 1:
            break
        prev, cur = cur, candidates[0]

    return chain

def apply_virtual_sites3(coords: np.ndarray, vs3: list[ProtocolVirtualSite3]) -> np.ndarray:
    # coords indexed by atomnr-1
    out = coords.copy()
    for v in vs3:
        out[v.vid - 1] = (out[v.a1 - 1] + out[v.a2 - 1] + out[v.a3 - 1]) / 3.0
    return out

def layout_xyz(G, mol, source: int, dz: float = 1.0, seed: int | None = None) -> List[Tuple[float, float, float]]:
    """
    Heuristic layout with global carbon constraint:

    - source is immutable at (0,0,0)
    - never overwrite coords for already-placed nodes (negative roads remain fixed once placed)
    - carbon-run detection: ANY run of >=3 consecutive atype starting with 'C' anywhere in the chain
    - rule:
        * if carbon-run exists -> chain wants NEG; if charge exists, allow upward shift so charge can reach z~0
        * if no carbon-run but charge exists -> POS
        * else fallback prefers inherited sign
      additionally:
        * GLOBAL constraint: for any chain with carbon-run, ensure all C-type atoms in that chain have z <= 0
          (shift the newly placed chain down if needed)
    - break rule:
        * if chain would start POS but has carbon-run: break at start of first carbon run:
            - beads before anchor stay on parent x,y above anchor
            - anchor+below goes down (negative) on a side lane
          then apply global carbon constraint (so C never ends up above 0 even if fork is on +z trail)
    - collision:
        * distance-based
        * radii tried: 2.0, 1.5, 1.0
        * cutoff starts > 1*dz and decays each try
    """
    rng = np.random.default_rng(seed)
    atoms_by_nr = {a.nr: a for a in mol.atoms}
    eps = 1e-12

    # distance collision params
    max_tries = 40
    radii = (2.0, 1.5, 1.0)
    cutoff0 = 1.35 * dz
    cutoff_decay = 0.92

    coords: Dict[int, Tuple[float, float, float]] = {source: (0.0, 0.0, 0.0)}
    placed: Set[int] = {source}
    placed_pts: List[Tuple[float, float, float]] = [(0.0, 0.0, 0.0)]  # for distance checks

    def is_C(nr: int) -> bool:
        return str(atoms_by_nr[nr].atype).strip().upper().startswith("C")

    # -------------------------
    # SOURCE ROADS
    # -------------------------
    roads = straight_roads_from_source(G, source)
    source_deg = G.degree(source)

    # Build road meta
    road_meta: List[Dict[str, object]] = []
    for road in roads:
        # carbon-run anywhere + first run start index
        run = 0
        max_run = 0
        carbon_run_start = None  # 1-based
        for i, n in enumerate(road, start=1):
            if is_C(n):
                run += 1
                if run > max_run:
                    max_run = run
                if run == 3 and carbon_run_start is None:
                    carbon_run_start = i - 2
            else:
                run = 0
        carbon_run_exists = (max_run >= 3)

        # charge + first charge step
        first_charge_step = None
        for step, n in enumerate(road, start=1):
            if abs(float(atoms_by_nr[n].charge)) > eps:
                first_charge_step = step
                break
        has_charge = (first_charge_step is not None)

        road_meta.append(
            {
                "road": road,
                "carbon_run_exists": carbon_run_exists,
                "carbon_run_start": carbon_run_start,
                "has_charge": has_charge,
                "first_charge_step": first_charge_step,
            }
        )

    # Partition: negatives first, positives second
    # would-start at source: charge => pos else neg
    neg_indices = []
    pos_indices = []
    for i, m in enumerate(road_meta):
        would_sgn = 1.0 if m["has_charge"] else -1.0
        if (m["carbon_run_exists"] and would_sgn > 0.0) or (would_sgn > 0.0 and not m["carbon_run_exists"]):
            pos_indices.append(i)
        else:
            neg_indices.append(i)

    # degree-1 failsafe
    if source_deg == 1 and road_meta:
        neg_indices = [0]
        pos_indices = []

    # choose one main negative on-axis
    main_neg_idx = neg_indices[0] if neg_indices else None

    # helper: propose chain coords then apply global-carbon constraint
    def apply_global_carbon_constraint(proposed: List[Tuple[int, float, float, float]]) -> List[Tuple[int, float, float, float]]:
        # if no carbon in proposed, nothing to do
        max_z_carbon = None
        for nr, x, y, z in proposed:
            if is_C(nr):
                if max_z_carbon is None or z > max_z_carbon:
                    max_z_carbon = z
        if max_z_carbon is None:
            return proposed
        if max_z_carbon <= 0.0:
            return proposed
        # shift down so max carbon ends at z=0
        shift = max_z_carbon
        return [(nr, x, y, z - shift) for (nr, x, y, z) in proposed]

    # helper: collision test for a list of candidate points
    def collides(candidate: List[Tuple[int, float, float, float]], cutoff: float) -> bool:
        c2 = cutoff * cutoff
        for _nr, x, y, z in candidate:
            for px, py, pz in placed_pts:
                dx = x - px
                dy = y - py
                dz_ = z - pz
                if (dx * dx + dy * dy + dz_ * dz_) < c2:
                    return True
        return False

    # ---- place negative roads (fixed once placed) ----
    for idx in neg_indices:
        road = road_meta[idx]["road"]
        carbon_run_exists = bool(road_meta[idx]["carbon_run_exists"])
        has_charge = bool(road_meta[idx]["has_charge"])
        first_charge_step = road_meta[idx]["first_charge_step"]

        # decide sgn/shift
        # carbon => neg; if charged, shift up so first charge at z=0 (then carbon constraint may shift back down)
        sgn = -1.0 if carbon_run_exists else -1.0
        z_shift_steps = int(first_charge_step) if (carbon_run_exists and has_charge and first_charge_step is not None) else 0

        # choose lane
        if idx == main_neg_idx:
            bx, by = 0.0, 0.0
        else:
            bx, by = 0.0, 0.0
            found = False
            for rad in radii:
                cutoff = cutoff0
                for _ in range(max_tries):
                    th = float(rng.uniform(0.0, 2.0 * np.pi))
                    bx_try = rad * float(np.cos(th))
                    by_try = rad * float(np.sin(th))

                    proposed = []
                    for i_step, nr in enumerate(road, start=1):
                        if nr in placed:
                            continue
                        z = (z_shift_steps * dz) + (sgn * dz * float(i_step))
                        proposed.append((nr, bx_try, by_try, z))
                    proposed = apply_global_carbon_constraint(proposed)

                    if not collides(proposed, cutoff):
                        bx, by = bx_try, by_try
                        found = True
                        break
                    cutoff *= cutoff_decay
                if found:
                    break

        proposed = []
        for i_step, nr in enumerate(road, start=1):
            if nr in placed:
                continue
            z = (z_shift_steps * dz) + (sgn * dz * float(i_step))
            proposed.append((nr, bx, by, z))
        proposed = apply_global_carbon_constraint(proposed)

        for nr, x, y, z in proposed:
            coords[nr] = (float(x), float(y), float(z))
            placed.add(nr)
            placed_pts.append((float(x), float(y), float(z)))

    # ---- place positive roads (allowed to kink to stay connected to neg) ----
    for idx in pos_indices:
        road = road_meta[idx]["road"]
        carbon_run_exists = bool(road_meta[idx]["carbon_run_exists"])
        carbon_run_start = road_meta[idx]["carbon_run_start"]
        has_charge = bool(road_meta[idx]["has_charge"])
        first_charge_step = road_meta[idx]["first_charge_step"]

        would_sgn = 1.0 if has_charge else -1.0  # source fallback neg
        break_mode = bool(carbon_run_exists and would_sgn > 0.0 and carbon_run_start is not None)

        # non-break positive: sgn=+1, no shift
        # (if it had carbon_run_exists and would_sgn>0 we'd be in break_mode)
        if not break_mode:
            sgn = 1.0
            z_shift_steps = 0
        else:
            sgn = 1.0
            z_shift_steps = 0

        # choose side lane for tail (bead 1 stays on-axis for neat connection)
        bx, by = 0.0, 0.0
        found = False
        for rad in radii:
            cutoff = cutoff0
            for _ in range(max_tries):
                th = float(rng.uniform(0.0, 2.0 * np.pi))
                bx_try = rad * float(np.cos(th))
                by_try = rad * float(np.sin(th))

                proposed = []
                if break_mode:
                    k = int(carbon_run_start)
                    for i_step, nr in enumerate(road, start=1):
                        if nr in placed:
                            continue
                        if i_step < k:
                            # pre-anchor on-axis above anchor plane (anchor plane at 0)
                            z = dz * float(k - i_step)
                            proposed.append((nr, 0.0, 0.0, z))
                        else:
                            # anchor+below on lane, anchor at z=0
                            z = -dz * float(i_step - k)
                            proposed.append((nr, bx_try, by_try, z))
                else:
                    for i_step, nr in enumerate(road, start=1):
                        if nr in placed:
                            continue
                        z = (z_shift_steps * dz) + (sgn * dz * float(i_step))
                        if i_step == 1:
                            proposed.append((nr, 0.0, 0.0, z))
                        else:
                            proposed.append((nr, bx_try, by_try, z))

                proposed = apply_global_carbon_constraint(proposed)

                if not collides(proposed, cutoff):
                    bx, by = bx_try, by_try
                    found = True
                    break
                cutoff *= cutoff_decay
            if found:
                break

        # commit
        proposed = []
        if break_mode:
            k = int(carbon_run_start)
            for i_step, nr in enumerate(road, start=1):
                if nr in placed:
                    continue
                if i_step < k:
                    z = dz * float(k - i_step)
                    proposed.append((nr, 0.0, 0.0, z))
                else:
                    z = -dz * float(i_step - k)
                    proposed.append((nr, bx, by, z))
        else:
            for i_step, nr in enumerate(road, start=1):
                if nr in placed:
                    continue
                z = (z_shift_steps * dz) + (sgn * dz * float(i_step))
                if i_step == 1:
                    proposed.append((nr, 0.0, 0.0, z))
                else:
                    proposed.append((nr, bx, by, z))

        proposed = apply_global_carbon_constraint(proposed)

        for nr, x, y, z in proposed:
            coords[nr] = (float(x), float(y), float(z))
            placed.add(nr)
            placed_pts.append((float(x), float(y), float(z)))

    # -------------------------
    # PROPAGATION: forks anywhere
    # -------------------------
    queue: List[int] = list(placed)

    while queue:
        u = queue.pop(0)
        ux, uy, uz = coords[u]
        inherited_sgn = -1.0 if uz <= 0.0 else 1.0

        unplaced_neighbors = [v for v in G.neighbors(u) if v not in placed]
        if not unplaced_neighbors:
            continue
        unplaced_neighbors = sorted(unplaced_neighbors)

        # build straight chain for each neighbor
        chains: Dict[int, List[int]] = {}
        lengths: Dict[int, int] = {}
        for v in unplaced_neighbors:
            prev = u
            cur = v
            chain: List[int] = []
            seen_local = {u}
            while True:
                if cur in placed:
                    break
                chain.append(cur)
                seen_local.add(cur)
                nxts = [x for x in G.neighbors(cur) if x != prev and x not in placed and x not in seen_local]
                if len(nxts) != 1:
                    break
                prev, cur = cur, nxts[0]
            chains[v] = chain
            lengths[v] = len(chain)

        main_v = max(unplaced_neighbors, key=lambda vv: lengths.get(vv, 0))
        radial_vs = [v for v in unplaced_neighbors if v != main_v]

        # place main straight then radials
        for kind, vs in (("straight", [main_v]), ("radial", radial_vs)):
            if not vs:
                continue

            if kind == "radial":
                k_out = len(vs)
                base_angles = (2.0 * np.pi * np.arange(k_out)) / k_out
                chosen_rot = 0.0
                for attempt in range(max_tries):
                    rot = 0.0 if attempt == 0 else float(rng.uniform(0.0, 2.0 * np.pi))
                    cutoff = cutoff0 * (cutoff_decay ** attempt)
                    ok = True
                    for j, v in enumerate(vs):
                        chain = chains[v]
                        if not chain:
                            continue

                        # carbon run + start
                        run = 0
                        max_run = 0
                        carbon_run_start = None
                        for i_step, nr in enumerate(chain, start=1):
                            if is_C(nr):
                                run += 1
                                if run > max_run:
                                    max_run = run
                                if run == 3 and carbon_run_start is None:
                                    carbon_run_start = i_step - 2
                            else:
                                run = 0
                        carbon_run_exists = (max_run >= 3)

                        # charge
                        first_charge_step = None
                        for step, nr in enumerate(chain, start=1):
                            if abs(float(atoms_by_nr[nr].charge)) > eps:
                                first_charge_step = step
                                break
                        has_charge = (first_charge_step is not None)

                        would_sgn = 1.0 if has_charge else inherited_sgn
                        break_mode = bool(carbon_run_exists and would_sgn > 0.0 and carbon_run_start is not None)

                        theta = float(base_angles[j] + rot)
                        x = ux + float(np.cos(theta))
                        y = uy + float(np.sin(theta))

                        # probe point
                        z_probe = uz if break_mode else (uz + (would_sgn * dz))
                        # distance check
                        c2 = cutoff * cutoff
                        for px, py, pz in placed_pts:
                            dx = x - px
                            dy = y - py
                            dz_ = z_probe - pz
                            if (dx * dx + dy * dy + dz_ * dz_) < c2:
                                ok = False
                                break
                        if not ok:
                            break
                    if ok:
                        chosen_rot = rot
                        break

            for j, v in enumerate(vs):
                chain = chains[v]
                if not chain:
                    continue

                # carbon run + start
                run = 0
                max_run = 0
                carbon_run_start = None
                for i_step, nr in enumerate(chain, start=1):
                    if is_C(nr):
                        run += 1
                        if run > max_run:
                            max_run = run
                        if run == 3 and carbon_run_start is None:
                            carbon_run_start = i_step - 2
                    else:
                        run = 0
                carbon_run_exists = (max_run >= 3)

                # charge + first charge step
                first_charge_step = None
                for step, nr in enumerate(chain, start=1):
                    if abs(float(atoms_by_nr[nr].charge)) > eps:
                        first_charge_step = step
                        break
                has_charge = (first_charge_step is not None)

                would_sgn = 1.0 if has_charge else inherited_sgn
                break_mode = bool(carbon_run_exists and would_sgn > 0.0 and carbon_run_start is not None)

                # decide non-break sign/shift
                if not break_mode:
                    if carbon_run_exists:
                        sgn = -1.0
                        z_shift_steps = int(first_charge_step) if (has_charge and first_charge_step is not None) else 0
                    elif has_charge:
                        sgn = 1.0
                        z_shift_steps = 0
                    else:
                        sgn = inherited_sgn
                        z_shift_steps = 0
                else:
                    sgn = 1.0
                    z_shift_steps = 0

                # choose lane coords
                # straight: prefer keep (ux,uy); for positive / break, allow kink
                if kind == "radial":
                    theta0 = float(base_angles[j] + chosen_rot)
                    base_bx = ux + float(np.cos(theta0))
                    base_by = uy + float(np.sin(theta0))
                else:
                    base_bx, base_by = ux, uy

                allow_kink = True
                keep_first_on_parent_xy = (kind == "straight")

                bx, by = base_bx, base_by

                # choose a side lane for tails if break or positive (kink)
                want_side_lane = bool(break_mode or (allow_kink and sgn > 0.0))
                if want_side_lane:
                    found = False
                    for rad in radii:
                        cutoff = cutoff0
                        for _ in range(max_tries):
                            th = float(rng.uniform(0.0, 2.0 * np.pi))
                            bx_try = ux + rad * float(np.cos(th))
                            by_try = uy + rad * float(np.sin(th))

                            proposed = []
                            if break_mode:
                                k = int(carbon_run_start)
                                for i_step, nr in enumerate(chain, start=1):
                                    if nr in placed:
                                        continue
                                    if i_step < k:
                                        z = uz + dz * float(k - i_step)
                                        proposed.append((nr, ux, uy, z))
                                    else:
                                        z = uz - dz * float(i_step - k)
                                        proposed.append((nr, bx_try, by_try, z))
                            else:
                                for i_step, nr in enumerate(chain, start=1):
                                    if nr in placed:
                                        continue
                                    z = uz + (z_shift_steps * dz) + (sgn * dz * float(i_step))
                                    if keep_first_on_parent_xy and sgn > 0.0 and i_step == 1:
                                        proposed.append((nr, ux, uy, z))
                                    else:
                                        proposed.append((nr, bx_try, by_try, z))

                            proposed = apply_global_carbon_constraint(proposed)

                            if not collides(proposed, cutoff):
                                bx, by = bx_try, by_try
                                found = True
                                break
                            cutoff *= cutoff_decay
                        if found:
                            break

                # commit placement
                proposed = []
                if break_mode:
                    k = int(carbon_run_start)
                    for i_step, nr in enumerate(chain, start=1):
                        if nr in placed:
                            continue
                        if i_step < k:
                            z = uz + dz * float(k - i_step)
                            proposed.append((nr, ux, uy, z))
                        else:
                            z = uz - dz * float(i_step - k)
                            proposed.append((nr, bx, by, z))
                else:
                    for i_step, nr in enumerate(chain, start=1):
                        if nr in placed:
                            continue
                        z = uz + (z_shift_steps * dz) + (sgn * dz * float(i_step))
                        if keep_first_on_parent_xy and sgn > 0.0 and allow_kink and i_step == 1:
                            proposed.append((nr, ux, uy, z))
                        else:
                            proposed.append((nr, bx, by, z))

                proposed = apply_global_carbon_constraint(proposed)

                for nr, x, y, z in proposed:
                    coords[nr] = (float(x), float(y), float(z))
                    placed.add(nr)
                    placed_pts.append((float(x), float(y), float(z)))
                    queue.append(nr)

    # hard enforce base at origin (and never moved)
    coords[source] = (0.0, 0.0, 0.0)

    # final failsafe
    for a in mol.atoms:
        if a.nr not in coords:
            coords[a.nr] = (0.0, 0.0, 0.0)

    return np.asarray([coords[a.nr] for a in mol.atoms])

def universe_from_mol_and_coords(mol, coords, *, resname: str = "MOL", resid: int = 1):
    """
    Build an MDAnalysis Universe from a MoleculeType and coordinates.

    - Atom names come from AtomRow.atom
    - coords must be list/array of shape (N, 3) in the same order as mol.atoms
    - All atoms are put into a single residue (required for clean GRO writing)
    """
    n_atoms = len(mol.atoms)

    # One residue, all atoms assigned to residue index 0
    u = mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=1,
        atom_resindex=np.zeros(n_atoms, dtype=int),
        trajectory=True,
    )

    # Atom-level attributes
    u.add_TopologyAttr("name", [a.atom for a in mol.atoms])
    u.add_TopologyAttr("id", np.arange(1, n_atoms + 1, dtype=int))

    # Residue-level attributes (length must equal n_residues == 1)
    u.add_TopologyAttr("resname", [resname])
    u.add_TopologyAttr("resid", [resid])

    # Positions
    u.atoms.positions = np.asarray(coords, dtype=np.float32)

    return u

def maker_itp(file: str, flipper: bool=False, sz: float=1, sxy: float=1, base: str="PO4",output: str="lib.txt", constraints: bool=False) -> None:
    from ..core.itp_parser import read_itp_molecules as itp

    mols=itp(file)
    all_failed=True
    for mol in mols:
        G=mol_to_graph(mol,constraints)
        try:
            source = next(n for n, d in G.nodes(data=True) if d["name"] == base)
        except StopIteration:
            print(f"base was not in {mol.name} and thus {mol.name} was skipped")
            continue
        all_failed=False
        G.nodes[source]["source"] = True

        coords = layout_xyz(G,mol, source, dz=1.0)
        coords=apply_virtual_sites3(coords,mol.virtual_sites)

        scale = np.array([sxy, sxy, -sz if flipper else sz], dtype=float)
        coords=coords*scale

        u=universe_from_mol_and_coords(mol,coords)
        write_file(u,mol.name+"_"+output,mol.name)
        
        u.dimensions = np.array([10.0, 10.0, 10.0, 90.0, 90.0, 90.0], dtype=np.float32)

        with mda.coordinates.GRO.GROWriter(mol.name+"_"+"layout.gro") as W:
            W.write(u.atoms)
        print(f"Wrote lib file entry into {output}")
    if all_failed:
        print("The base bead was not present in any of the molecules, no file was written.")
    else:
        print("""File(s) created!
    Please note that the libmaker is a crude helper and the correctness of the generated lipid cannot be guaranteed. A visualization aid has been written as a .gro. Check before you run it.
    To use the new lipid, copy the contents of the generated file into the LIB file.
            """)

def library_file_preparer(args: List[str]) -> None:
    """Experimental tool to generate LIB entry from gro or pdb"""
    parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-p','--path',default="structure.gro",help="Single lipid structure file.")
    parser.add_argument('-f','--flip',default=False,action='store_true',help="Flip the results along the xy plane, so +z -> -z")
    parser.add_argument('-sz','--scalingz',default=1,type=float,help="(0 to 1), default 1,Move beads closer together in z direction, easier placement, more minimization")
    parser.add_argument('-sxy','--scalingxy',default=1,type=float,help="(0 to 1), default 1, Move beads closer together in x-y direction, easier placement, more minimization")
    parser.add_argument('-b','--base',default="PO4",type=str,help="Name of the bead that serves as a reference with coordinates 0.0 0.0 0.0")
    parser.add_argument('-o','--output',default="lib.txt",type=str,help="Output text file to be copied into the LIB file.")
    parser.add_argument('--consider_constraints',default=False,action='store_true',help="Uses constraints in the itp as bonds")
    
    args = parser.parse_args(args)
    logging.basicConfig(level=logging.INFO)
    file=Path(args.path)

    try:
        if file.suffix == ".gro" or file.suffix == ".pdb":
            maker_gro(file=file, flipper=args.flip, sz=args.scalingz, sxy=args.scalingxy, base=args.base,output=args.output)
            print("""File created!
    Please note that the libmaker is a crude helper and the correctness of the generated lipid cannot be guaranteed.
    To use the new lipid, copy the contents of the generated file into the LIB file.
            """)
        elif file.suffix == ".itp":
            maker_itp(file=file, flipper=args.flip, sz=args.scalingz, sxy=args.scalingxy, base=args.base,output=args.output, constraints=args.consider_constraints)
            
        else:
            raise ValueError(f"Unsupported file type: {file.suffix}. Only .gro and .itp are supported.")

    except Exception as e:
        logger.error(f"Error: {e}")
        raise
