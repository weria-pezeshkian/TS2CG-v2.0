"""
CLI tool to place lipids to assign circular domains around inclusions or points.
"""

import argparse
from pathlib import Path
import numpy as np
import logging
from dataclasses import dataclass
from typing import List, Optional, Sequence, Dict
from scipy.special import logsumexp
import networkx as nx
from scipy.spatial.distance import cdist

from ..core.point import Point

logger = logging.getLogger(__name__)

@dataclass
class LipidSpec:
    """Specification for a lipid type and its properties"""
    domain_id: int
    name: str
    percentage: float
    curvature: float
    density: float

def parse_lipid_file(file_path: Path) -> List[LipidSpec]:
    """Parse lipid specification file into structured data"""
    lipids = []
    with open(file_path) as f:
        for line in f:
            if line.strip() and not line.startswith(';'):
                domain_id, name, percentage, curvature, density = line.split()
                lipids.append(LipidSpec(
                    domain_id=int(domain_id),
                    name=name,
                    percentage=float(percentage),
                    curvature=float(curvature),
                    density=float(density)
                ))

    total = sum(lipid.percentage for lipid in lipids)
    if not np.isclose(total, 1.0, atol=0.01):
        raise ValueError(f"Lipid percentages must sum to 1.0 (got {total:.2f})")

    return lipids

def write_input_str(lipids: Sequence[LipidSpec], output_file: Path, old_input: Optional[Path] = None) -> None:
    """
    Write input.str file for TS2CG, preserving all comments and sections except [Lipids List].
    Maintains exact formatting of the original file.
    """
    # If no old input exists, just write the lipids section
    if not old_input or not old_input.exists():
        with open(output_file, 'w') as f:
            f.write("[Lipids List]\n")
            for lipid in lipids:
                f.write(f"Domain {lipid.domain_id}\n")
                f.write(f"{lipid.name} 1 1 {lipid.density}\n")
                f.write("End\n")
        return

    # Read the old file and process it section by section
    with open(old_input) as f:
        content = f.read()

    # Split the content into sections
    sections = []
    current_section = []
    in_lipids_section = False

    for line in content.split('\n'):
        # Detect section headers
        if line.strip().startswith('['):
            # If we were building a section, save it
            if current_section:
                sections.append('\n'.join(current_section) + '\n')
                current_section = []

            # Check if we're entering the lipids section
            if '[Lipids List]' in line:
                in_lipids_section = True
                # Create new lipids section
                new_section = ['[Lipids List]']
                for lipid in lipids:
                    new_section.extend([
                        f"Domain {lipid.domain_id}",
                        f"{lipid.name} 1 1 {lipid.density}",
                        "End"
                    ])
                sections.append('\n'.join(new_section) + '\n')
            else:
                in_lipids_section = False
                current_section.append(line)
        # For non-header lines
        elif not in_lipids_section:
            current_section.append(line)

    # Add the last section if it exists
    if current_section and not in_lipids_section:
        sections.append('\n'.join(current_section))

    # Write everything back to the file
    with open(output_file, 'w') as f:
        f.write('\n'.join(sections))

def _get_centers(membrane,Type,dummys):
    from_type=[inc["point_id"] for inc in membrane.inclusions.get_by_type(Type)]
    from_dummy=dummys.strip().replace(" ","").split(",")

    to_set=list(set(from_type+from_dummy))

    return list(filter(None,to_set))

def _dijkstra_within_radius(distance_matrix, start_index, radius,percent=100):
    if percent != 100:
        distance_matrix = distance_matrix[distance_matrix < np.percentile(distance_matrix,percent)] = 0
    G = nx.from_numpy_array(distance_matrix)
    shortest_paths = nx.single_source_dijkstra_path_length(G, start_index)
    reachable_nodes = [node for node, dist in shortest_paths.items() if dist < radius]
    return reachable_nodes

def _within_radius(distance_matrix,radius):
    reachable_nodes=[]
    for i,value in enumerate(distance_matrix):
        if value <= radius:
            reachable_nodes.append(i)

    return reachable_nodes


def circular_domains(membrane: Point, radius: float, pointid: list, domain: int, path_dist: bool=False, percent: float=100.0, layer: str = "both") -> None:
    """Assign lipids to domains based on curvature preferences"""
    layers = [membrane.outer]

    if membrane.monolayer or layer.lower() == "outer":
        layers = [membrane.outer]
    elif layer.lower() == "inner":
        layers = [membrane.inner]
    else:
        layers = [membrane.outer, membrane.inner]

    for layer in layers:
        if path_dist:
            dist_matrix = cdist(layer.coordinates, layer.coordinates)
        for point in pointid:
            if path_dist:                
                nodes=_dijkstra_within_radius(dist_matrix,int(point),radius,percent)
            else:
                dist_matrix=cdist(layer.coordinates,[layer.coordinates[int(point)]])
                nodes=_within_radius(dist_matrix,radius)
            for index in nodes:
                layer.domain_ids[index]=domain


def DAI(args: List[str]) -> None:
    """Main entry point for Domain Placer tool"""
    parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-p','--point_dir',default="point/",help="Specify the path to the point folder")
    parser.add_argument('-r','--radius',default=1,type=float,help="The radius around a protein in which domains should be changed.")
    parser.add_argument('-t','--type',default=0,type=int,help="The protein type around which domains should be changed.")
    parser.add_argument('-d','--Domain',default=1,type=int,help="The domain number that should be set around the protein.")
    parser.add_argument('-l','--leaflet',default="both",help="Choose which membrane leaflet to alter. Default is both")
    parser.add_argument('-dummy','--dummy',default="",help="Create a dummy protein to place a circular domain around it. Excepts pointids like 3,7,22")
    parser.add_argument('-pd','--path_distance',default=False,type=bool,help="Slower execution, but needed for higher curvature membranes to assign the domain to only one membrane part")
    parser.add_argument('-pdP','--path_distance_percentile',default=100.0,type=float,help="Manipulates neighbors in path distance, tests have shown that 2 works well")

    args = parser.parse_args(args)
    logging.basicConfig(level=logging.INFO)

    try:
        membrane = Point(args.point_dir)

        centers=_get_centers(membrane, args.type,args.dummy)
        circular_domains(membrane=membrane,radius=args.radius,pointid=centers,domain=args.Domain, layer=args.leaflet,path_dist=args.path_distance,percent=args.path_distance_percentile)

        output_dir = args.point_dir
        membrane.save(output_dir)
        logger.info(f"Updated membrane domains in {output_dir}")

    except Exception as e:
        logger.error(f"Error: {e}")
        raise
