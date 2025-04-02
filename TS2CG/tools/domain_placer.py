"""
CLI tool to place lipids in membrane domains based on local curvature preferences.
Uses domain_input.txt format for lipid specifications and generates input.str for next steps.

Example domain_input.txt:
; domain lipid percentage c0 density
0 POPC .5 0.179 0.64
2 POPG .5 0.629 0.64
"""

import argparse
from pathlib import Path
import numpy as np
import logging
from dataclasses import dataclass
from typing import List, Optional, Sequence, Dict
from scipy.special import logsumexp

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

def calculate_curvature_weights(local_curvature: float, lipids: Sequence[LipidSpec],
                              k_factor: float, area:float=1.0, max_delta: float = 5.0) -> np.ndarray:
    """Calculate Boltzmann weights for each lipid type at given curvature"""
    # Calculate curvature differences
    delta_curvatures = np.array([abs(2 * local_curvature - lipid.curvature) for lipid in lipids])

    # Calculate log weights using the log-sum-exp trick for numerical stability
    log_weights = -k_factor * (delta_curvatures**2)*area
    max_log_weight = np.max(log_weights)
    exp_weights = np.exp(log_weights - max_log_weight)
    weights = exp_weights / np.sum(exp_weights)

    return weights

def assign_domains(membrane: Point, lipids: Sequence[LipidSpec], layer: str = "both",
                  k_factor: float = 1.0,area:bool=False) -> None:
    """Assign lipids to domains based on curvature preferences"""
    layers = [membrane.outer]

    if membrane.monolayer or layer.lower() == "outer":
        layers = [membrane.outer]
    elif layer.lower() == "inner":
        layers = [membrane.inner]
    else:
        layers = [membrane.outer, membrane.inner]

    for membrane_layer in layers:
        layer_name = "outer" if membrane_layer is membrane.outer else "inner"
        logger.info(f"Processing {layer_name} layer")

        n_points = len(membrane_layer.ids)
        curvatures = membrane_layer.mean_curvature

        # Flip curvature sign for inner membrane
        if layer_name == "inner":
            curvatures = -curvatures

        # Initialize domain assignments and tracking
        new_domains = np.full(n_points, -1)
        remaining_counts = {i: int(lipid.percentage * n_points)
                          for i, lipid in enumerate(lipids)}
        remaining_counts[len(lipids)-1] += n_points - sum(remaining_counts.values())
        available_lipids = set(range(len(lipids)))

        # Randomly process points
        perm=rng.permuted(range(len(membrane_layer.ids)))
        permuted_ids=membrane_layer.ids[perm]
        permuted_areas=membrane_layer.area[perm]
        for c,idx in enumerate(permuted_ids):
            local_curv = curvatures[idx]

            # Calculate weights only for available lipid types
            valid_lipids = [i for i in available_lipids
                          if remaining_counts[i] > 0]

            if not valid_lipids:
                logger.warning(f"No lipids available for point {idx}")
                valid_lipids = list(range(len(lipids)))

            valid_specs = [lipids[i] for i in valid_lipids]
            weights = calculate_curvature_weights(local_curv, valid_specs, k_factor,permuted_areas[c])

            # Choose lipid type and update bookkeeping
            chosen_idx = rng.choice(valid_lipids, p=weights)
            new_domains[idx] = lipids[chosen_idx].domain_id
            remaining_counts[chosen_idx] -= 1

            if remaining_counts[chosen_idx] == 0:
                available_lipids.remove(chosen_idx)

        # Update membrane
        membrane_layer.domain_ids = new_domains

        # Log results
        for lipid in lipids:
            actual_count = np.sum(new_domains == lipid.domain_id)
            logger.info(f"{lipid.name}: {actual_count/n_points*100:.1f}% "
                       f"(target: {lipid.percentage*100:.1f}%)")

def DOP(args: List[str]) -> None:
    """Main entry point for Domain Placer tool"""
    parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-p', '--point-dir', type=Path, default="point",
                       help='Path to membrane point directory (default: point/)')
    parser.add_argument('-i', '--input', type=Path, default="domain_input.txt",
                       help='Path to lipid specification file (default: domain_input.txt)')
    parser.add_argument('-l', '--leaflet', choices=['both', 'inner', 'outer'],
                       default='both', help='Which membrane leaflet(s) to modify (default: both)')
    parser.add_argument('-k', '--k-factor', type=float, default=1.0,
                       help='Scaling factor for curvature preference strength (default: 1.0)')
    parser.add_argument('-o', '--output', type=Path,
                       help='Output directory (defaults to input directory)')
    parser.add_argument('-ni', '--new-input', type=Path, default="input_DOP.str",
                       help='Path for output input.str file (default: input.str)')
    parser.add_argument('-oi', '--old-input', type=Path,
                       help='Path to existing input.str to preserve additional sections')
    parser.add_argument('--seed', type=int, help='Random seed for reproducibility')
    parser.add_argument('--area',default=False,action='store_true',help="Considers the area of each point for weight calculation")

    args = parser.parse_args(args)
    logging.basicConfig(level=logging.INFO)

    # setup numpy random number generator
    global rng
    rng = np.random.default_rng(args.seed)

    try:
        membrane = Point(args.point_dir)
        lipids = parse_lipid_file(args.input)

        assign_domains(membrane, lipids, args.leaflet, args.k_factor,args.area)

        write_input_str(lipids, args.new_input, args.old_input)
        logger.info(f"Created input file: {args.new_input}")

        output_dir = args.output or args.point_dir
        membrane.save(output_dir)
        logger.info(f"Updated membrane domains in {output_dir}")

    except Exception as e:
        logger.error(f"Error: {e}")
        raise
