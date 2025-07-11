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

def maker(file: str, flipper: bool=False, sz: float=.5, sxy: float=.5, base: str="PO4",output: str="lib.txt") -> None:
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

    try:
        maker(file=args.path, flipper=args.flip, sz=args.scalingz, sxy=args.scalingxy, base=args.base,output=args.output)
        print("""File created!
Please note that the libmaker is a crude helper and the correctness of the generated lipid cannot be guaranteed.
To use the new lipid, copy the contents of the generated file into the LIB file.
        """)

    except Exception as e:
        logger.error(f"Error: {e}")
        raise
