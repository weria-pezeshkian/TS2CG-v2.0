"""
Tool to place lipids in membrane domains based on local curvature preferences.
"""

import argparse
from pathlib import Path
import numpy as np
import logging
from dataclasses import dataclass
from typing import List, Optional, Sequence

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
                              k_factor: float, area: float = 1.0) -> np.ndarray:
    """Calculate Boltzmann weights for each lipid type at given curvature"""
    # Calculate curvature differences
    delta_curvatures = np.array([abs(2 * local_curvature - lipid.curvature) for lipid in lipids])

    # Calculate log weights using the log-sum-exp trick for numerical stability
    log_weights = -k_factor * (delta_curvatures**2) * area
    max_log_weight = np.max(log_weights)
    exp_weights = np.exp(log_weights - max_log_weight)
    weights = exp_weights / np.sum(exp_weights)

    return weights

def assign_domains(membrane: Point, lipids: Sequence[LipidSpec], leaflet: str = "both",
                  k_factor: float = 1.0, area_weighted: bool = False, seed: Optional[int] = None) -> None:
    """Assign lipids to domains based on curvature preferences"""

    # Set random seed
    rng = np.random.default_rng(seed)

    # Determine leaflets to process
    leaflets = []
    if membrane.monolayer or leaflet.lower() == "outer":
        leaflets = [("outer", membrane.outer)]
    elif leaflet.lower() == "inner":
        leaflets = [("inner", membrane.inner)]
    elif leaflet.lower() == "both":
        leaflets = [("outer", membrane.outer)]
        if not membrane.monolayer:
            leaflets.append(("inner", membrane.inner))

    for leaflet_name, membrane_leaflet in leaflets:
        logger.info(f"Processing {leaflet_name} leaflet")

        n_points = len(membrane_leaflet.ids)
        curvatures = membrane_leaflet.mean_curvature

        # Flip curvature sign for inner membrane
        if leaflet_name == "inner":
            curvatures = -curvatures

        # Initialize domain assignments and tracking
        new_domains = np.full(n_points, -1)
        remaining_counts = {i: int(lipid.percentage * n_points)
                          for i, lipid in enumerate(lipids)}
        remaining_counts[len(lipids)-1] += n_points - sum(remaining_counts.values())
        available_lipids = set(range(len(lipids)))

        # Randomly process points
        perm = rng.permuted(range(len(membrane_leaflet.ids)))
        permuted_ids = membrane_leaflet.ids[perm]
        permuted_areas = membrane_leaflet.area[perm]

        for c, idx in enumerate(permuted_ids):
            local_curv = curvatures[idx]

            # Calculate weights only for available lipid types
            valid_lipids = [i for i in available_lipids
                          if remaining_counts[i] > 0]

            if not valid_lipids:
                logger.warning(f"No lipids available for point {idx}")
                valid_lipids = list(range(len(lipids)))

            valid_specs = [lipids[i] for i in valid_lipids]
            if area_weighted:
                weights = calculate_curvature_weights(local_curv, valid_specs, k_factor, permuted_areas[c])
            else:
                weights = calculate_curvature_weights(local_curv, valid_specs, k_factor)

            # Choose lipid type and update bookkeeping
            chosen_idx = rng.choice(valid_lipids, p=weights)
            new_domains[idx] = lipids[chosen_idx].domain_id
            remaining_counts[chosen_idx] -= 1

            if remaining_counts[chosen_idx] == 0:
                available_lipids.remove(chosen_idx)

        # Update membrane
        membrane_leaflet.domain_ids = new_domains

        # Log results
        for lipid in lipids:
            actual_count = np.sum(new_domains == lipid.domain_id)
            logger.info(f"{lipid.name}: {actual_count/n_points*100:.1f}% "
                       f"(target: {lipid.percentage*100:.1f}%)")

def DOP(args: List[str]) -> None:
    """Main entry point for Domain Placer tool"""
    parser = argparse.ArgumentParser(description="Place lipids in membrane domains based on curvature preferences")

    parser.add_argument('--point-dir', '-p', default="point/",
                       help="Path to point folder (default: point/)")
    parser.add_argument('--lipid-specs', '-s', default="domain_input.txt",
                       help="Path to lipid specification file (default: domain_input.txt)")
    parser.add_argument('--leaflet', '-l', choices=['both', 'inner', 'outer'], default='both',
                       help="Which membrane leaflet to modify (default: both)")
    parser.add_argument('--k-factor', '-k', type=float, default=1.0,
                       help="Curvature preference strength (default: 1.0)")
    parser.add_argument('--output-dir', '-o',
                       help="Output directory (default: overwrite input with backup)")
    parser.add_argument('--no-backup', action='store_true',
                       help="Skip creating backup when overwriting input")
    parser.add_argument('--new-input', '-ni', default="input_DOP.str",
                       help="Path for output input.str file (default: input_DOP.str)")
    parser.add_argument('--old-input', '-no',
                       help="Path to existing input.str to preserve additional sections")
    parser.add_argument('--area', action='store_true',
                       help="Consider area of each point for weight calculation")
    parser.add_argument('--seed', type=int,
                       help="Random seed for reproducibility")

    args = parser.parse_args(args)

    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    # Load membrane and lipid specifications
    logger.info(f"Loading membrane from {args.point_dir}")
    membrane = Point(args.point_dir)

    logger.info(f"Loading lipid specifications from {args.lipid_specs}")
    lipids = parse_lipid_file(Path(args.lipid_specs))

    # Assign domains
    assign_domains(membrane, lipids, args.leaflet, args.k_factor, args.area, args.seed)

    # Write input.str file
    old_input_path = Path(args.old_input) if args.old_input else None
    write_input_str(lipids, Path(args.new_input), old_input_path)
    logger.info(f"Created input file: {args.new_input}")

    # Save results
    output_dir = args.output_dir if args.output_dir else args.point_dir
    backup = not args.no_backup
    membrane.save(output_dir, backup=backup)

    logger.info(f"Finished updating membrane domains!")
