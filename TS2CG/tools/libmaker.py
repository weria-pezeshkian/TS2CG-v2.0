"""
Helper tool to generate a lib file section
"""

import argparse
from pathlib import Path
import numpy as np
import logging
from dataclasses import dataclass
from typing import List, Optional, Sequence, Dict
import MDAnalysis as mda
import networkx as nx
import random

logger = logging.getLogger(__name__)


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

def mol_to_graph(mol) -> nx.Graph:
    G = nx.Graph()

    # nodes: use atom index (nr) as node id, store atom name as attribute
    for a in mol.atoms:
        G.add_node(a.nr, name=a.atom)

    # edges: bonds
    for b in mol.bonds:
        G.add_edge(b.ai, b.aj)

    return G

def straight_roads_from_source(G: nx.Graph, source: int):
    """
    Returns a list of roads. Each road is a list of node ids starting AFTER source:
    [v1, v2, v3, ...] where v1 is a neighbor of source.
    """
    roads = []
    for nb in G.neighbors(source):
        road = []
        prev = source
        cur = nb
        while True:
            road.append(cur)
            # forward candidates excluding where we came from
            nxts = [x for x in G.neighbors(cur) if x != prev]
            if len(nxts) != 1:
                # stop if dead end (0) or branch/cycle junction (>1)
                break
            prev, cur = cur, nxts[0]
        roads.append(road)
    return roads

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


def choose_neg_pos_roads(mol, roads: List[List[int]], seed: int | None = None) -> Tuple[List[int], List[int], List[List[int]]]:
    """
    Returns (neg_road, pos_road, rest_roads)

    Cascade rules (only proceed if unresolved):
      1) more *type* starting with 'C' => negative
      2) if charges appear, chain with charge => positive
      3) longer => negative
      4) randomize
    """
    roads_sorted = sorted(roads, key=len, reverse=True)
    main = roads_sorted[:2]
    rest = roads_sorted[2:]

    if not main:
        return [], [], rest
    if len(main) == 1:
        return main[0], [], rest

    r1, r2 = main[0], main[1]
    atoms_by_nr = {a.nr: a for a in mol.atoms}

    def count_type_C(road: List[int]) -> int:
        return sum(1 for n in road if str(atoms_by_nr[n].atype).startswith("C"))

    def has_charge(road: List[int], eps: float = 1e-12) -> bool:
        return any(abs(float(atoms_by_nr[n].charge)) > eps for n in road)

    # 1) more type starting with C => negative
    c1, c2 = count_type_C(r1), count_type_C(r2)
    if c1 != c2:
        neg, pos = (r1, r2) if c1 > c2 else (r2, r1)
        return neg, pos, rest

    # 2) chain with charge => positive (if only one has charge)
    q1, q2 = has_charge(r1), has_charge(r2)
    if q1 != q2:
        pos, neg = (r1, r2) if q1 else (r2, r1)
        return neg, pos, rest

    # 3) longer => negative
    if len(r1) != len(r2):
        neg, pos = (r1, r2) if len(r1) > len(r2) else (r2, r1)
        return neg, pos, rest

    # 4) randomize
    rng = random.Random(seed)
    return (r1, r2, rest) if rng.random() < 0.5 else (r2, r1, rest)


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


# -------- only function you call --------
def layout_xyz(G, mol, source: int, dz: float = 1.0, seed: int | None = None) -> List[Tuple[float, float, float]]:
    rng = np.random.default_rng(seed)

    roads = straight_roads_from_source(G, source)

    # --- failsafe: source degree 1 => only negative branch, no positive branch ---
    # (Explicitly override the neg/pos decision)
    if G.degree(source) == 1:
        if roads:
            neg_road, pos_road, _rest = roads[0], [], roads[1:]
        else:
            neg_road, pos_road, _rest = [], [], []
    else:
        neg_road, pos_road, _rest = choose_neg_pos_roads(mol, roads, seed=seed)

    coords: Dict[int, Tuple[float, float, float]] = {source: (0.0, 0.0, 0.0)}
    backbone: Set[int] = {source}

    # place backbone on z-axis
    for k, n in enumerate(neg_road, start=1):
        coords[n] = (0.0, 0.0, -dz * k)
        backbone.add(n)
    for k, n in enumerate(pos_road, start=1):
        coords[n] = (0.0, 0.0, dz * k)
        backbone.add(n)

    blocked: Set[int] = set(backbone)

    tol = 1e-6
    occupied: Set[Tuple[int, int, int]] = set()
    for x, y, z in coords.values():
        occupied.add((round(x / tol), round(y / tol), round(z / tol)))

    max_tries = 20
    r = 1.0

    for fork in sorted(backbone):
        fork_x, fork_y, fork_z = coords[fork]

        outs = [nb for nb in G.neighbors(fork) if nb not in blocked]
        if not outs:
            continue
        outs = sorted(outs)
        k_out = len(outs)

        base_angles = (2.0 * np.pi * np.arange(k_out)) / k_out

        chosen_rot = 0.0
        for attempt in range(max_tries):
            rot = 0.0 if attempt == 0 else float(rng.uniform(0.0, 2.0 * np.pi))

            ok = True
            for idx in range(k_out):
                if fork == source and idx in (0, 1):
                    x, y = fork_x, fork_y
                    sgn = -1.0 if idx == 0 else 1.0
                else:
                    theta = float(base_angles[idx] + rot)
                    x = fork_x + r * float(np.cos(theta))
                    y = fork_y + r * float(np.sin(theta))
                    if fork == source:
                        sgn = -1.0 if (idx % 2 == 0) else 1.0
                    else:
                        sgn = -1.0 if fork_z <= 0 else 1.0

                z1 = fork_z + sgn * dz * 1.0
                key = (round(x / tol), round(y / tol), round(z1 / tol))
                if key in occupied:
                    ok = False
                    break

            if ok:
                chosen_rot = rot
                break

        for idx, nb in enumerate(outs):
            if nb in blocked:
                continue

            if fork == source and idx in (0, 1):
                bx, by = fork_x, fork_y
                branch_sgn = -1.0 if idx == 0 else 1.0
            else:
                theta = float(base_angles[idx] + chosen_rot)
                bx = fork_x + r * float(np.cos(theta))
                by = fork_y + r * float(np.sin(theta))
                if fork == source:
                    branch_sgn = -1.0 if (idx % 2 == 0) else 1.0
                else:
                    branch_sgn = -1.0 if fork_z <= 0 else 1.0

            chain = walk_straight_chain(G, fork, nb, blocked)

            bx_try, by_try = bx, by
            for _ in range(max_tries):
                collided = False
                for i, node in enumerate(chain):
                    z = fork_z + branch_sgn * dz * float(i + 1)
                    key = (round(bx_try / tol), round(by_try / tol), round(z / tol))
                    if key in occupied:
                        collided = True
                        break
                if not collided:
                    break

                theta = float(rng.uniform(0.0, 2.0 * np.pi))
                rr = 0.0 if (fork == source and idx in (0, 1)) else r
                bx_try = fork_x + rr * float(np.cos(theta))
                by_try = fork_y + rr * float(np.sin(theta))

            for i, node in enumerate(chain):
                z = fork_z + branch_sgn * dz * float(i + 1)
                x = float(bx_try)
                y = float(by_try)
                coords[node] = (x, y, float(z))
                occupied.add((round(x / tol), round(y / tol), round(z / tol)))

    return [coords[a.nr] for a in mol.atoms]

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

def maker_itp(file: str, flipper: bool=False, sz: float=.5, sxy: float=.5, base: str="PO4",output: str="lib.txt") -> None:
    from ..core.itp_parser import read_itp_molecules as itp

    mols=itp(file)
    for mol in mols:
        G=mol_to_graph(mol)
        try:
            source = next(n for n, d in G.nodes(data=True) if d["name"] == base)
        except StopIteration:
            print(f"base was not in {mol.name} and thus {mol.name} was skipped")
            continue
        G.nodes[source]["source"] = True
        coords = layout_xyz(G,mol, source, dz=1.0)
        u=universe_from_mol_and_coords(mol,coords)
        write_file(u,mol.name+"_"+output,mol.name)
        print(f"Wrote lib file entry into {output}")
        u.dimensions = np.array([10.0, 10.0, 10.0, 90.0, 90.0, 90.0], dtype=np.float32)

        with mda.coordinates.GRO.GROWriter(mol.name+"_"+"layout.gro") as W:
            W.write(u.atoms)
        


    """    
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
    """


def library_file_preparer(args: List[str]) -> None:
    """Experimental tool to generate LIB entry from gro or pdb"""
    parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-p','--path',default="structure.gro",help="Single lipid structure file.")
    parser.add_argument('-f','--flip',default=False,action='store_true',help="Flip the results along the xy plane, so +z -> -z")
    parser.add_argument('-sz','--scalingz',default=.5,type=float,help="Move beads closer together in z direction, easier placement, more minimization")
    parser.add_argument('-sxy','--scalingxy',default=.5,type=float,help="Move beads closer together in x-y direction, easier placement, more minimization")
    parser.add_argument('-b','--base',default="PO4",type=str,help="Name of the bead that serves as a reference with coordinates 0.0 0.0 0.0")
    parser.add_argument('-o','--output',default="lib.txt",type=str,help="Output text file to be copied into the LIB file.")
    
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
            maker_itp(file=file, flipper=args.flip, sz=args.scalingz, sxy=args.scalingxy, base=args.base,output=args.output)
            print("""File created!
    Please note that the libmaker is a crude helper and the correctness of the generated lipid cannot be guaranteed.
    To use the new lipid, copy the contents of the generated file into the LIB file.
            """)
        else:
            raise ValueError(f"Unsupported file type: {file.suffix}. Only .gro and .itp are supported.")

    except Exception as e:
        logger.error(f"Error: {e}")
        raise
