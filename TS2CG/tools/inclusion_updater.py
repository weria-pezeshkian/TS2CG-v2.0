"""
CLI tool to place protein inclusions in membrane.
Reads input.str for protein definitions and places new inclusions in the membrane.

Usage:
    TS2CG INU -p point -t 0 -r 2 -N 10 -c 0.1 -o point_new -l both
"""

import argparse
from pathlib import Path
import numpy as np
import logging
from dataclasses import dataclass
from typing import List, Optional, Set, Dict, Sequence, Literal
from scipy.special import logsumexp
from scipy.spatial.distance import cdist

from ..core.point import Point

logger = logging.getLogger(__name__)

def pbc_wrap(membrane):
    """Wrap membrane coordinates into the primary box"""
    box = membrane.box
    membrane.outer.coordinates = membrane.outer.coordinates - box * np.round(membrane.outer.coordinates / box)
    if not membrane.monolayer:
        membrane.inner.coordinates = membrane.inner.coordinates - box * np.round(membrane.inner.coordinates / box)
    return membrane

def get_nearby_points_both_leaflets(membrane: Point, leaflet: str, point_idx: int,
                                  radius: float) -> Dict[str, np.ndarray]:
    """
    Find points within radius in both leaflets from a point in the specified leaflet.

    Args:
        membrane: Membrane Point object
        leaflet: Which leaflet the center point is in ('inner' or 'outer')
        point_idx: Index of the center point
        radius: Exclusion radius

    Returns:
        Dict with excluded points for each leaflet
    """
    # Get coordinates of the center point
    source_layer = membrane.outer if leaflet == 'outer' else membrane.inner
    center_coords = source_layer.coordinates[point_idx]

    excluded = {}

    # Check both leaflets
    for target_leaflet in ['outer', 'inner']:
        if target_leaflet == 'inner' and membrane.monolayer:
            continue

        target_layer = membrane.outer if target_leaflet == 'outer' else membrane.inner

        # Calculate displacement vectors
        displacement = target_layer.coordinates - center_coords

        # Apply minimum image convention
        displacement = displacement - membrane.box * np.round(displacement / membrane.box)

        # Calculate distances
        distances = np.sqrt(np.sum(displacement**2, axis=1))

        # Get points within radius
        excluded[target_leaflet] = np.where(distances <= radius)[0]

    return excluded

def calculate_curvature_weights(curvatures: np.ndarray, target_curvature: Optional[float],
                              k_factor: float) -> np.ndarray:
    """Calculate Boltzmann weights based on curvature preference"""
    if target_curvature is None:
        # If no curvature preference, return uniform weights
        return np.ones_like(curvatures) / len(curvatures)

    # Calculate weights using the log-sum-exp trick for numerical stability
    deltas = curvatures - target_curvature
    log_weights = -k_factor * deltas**2
    log_weights -= logsumexp(log_weights)  # normalize
    return np.exp(log_weights)

def get_points_near_existing_proteins(membrane: Point, radius: float) -> Dict[str, Set[int]]:
    """Get points that are too close to existing proteins in any leaflet"""
    excluded_points = {'outer': set(), 'inner': set()}

    # Skip if no existing proteins
    if not membrane.inclusions:
        return excluded_points

    # Check each existing protein
    for inclusion in membrane.inclusions:
        point_id = inclusion['point_id']
        # Determine which leaflet the protein is in (assuming outer if point exists in both)
        if point_id in membrane.outer.ids:
            source_leaflet = 'outer'
        else:
            source_leaflet = 'inner'

        # Get points too close in both leaflets
        nearby = get_nearby_points_both_leaflets(
            membrane,
            source_leaflet,
            point_id,
            radius
        )

        # Update excluded points for each leaflet
        for leaflet in nearby:
            excluded_points[leaflet].update(nearby[leaflet])

    return excluded_points

def place_proteins(membrane: Point, type_id: int, radius: float,
                  num_proteins: Optional[int] = None, target_curvature: Optional[float] = None,
                  k_factor: float = 1.0, leaflet: str = 'both') -> Dict[str, int]:
    """Place proteins in membrane with given constraints"""

    # Get next available type_id
    existing_type_ids = [inc['type_id'] for inc in membrane.inclusions]
    type_id = type_id or max(existing_type_ids, default=0) + 1

    # Get points excluded by existing proteins
    excluded_by_existing = get_points_near_existing_proteins(membrane, radius)
    logger.info("Excluding points near existing proteins:")
    for l, points in excluded_by_existing.items():
        if points:
            logger.info(f"  {l} leaflet: {len(points)} points excluded")

    # Initialize available points for each leaflet
    available_points = {
        'outer': set(membrane.outer.ids) - excluded_by_existing['outer']
                if leaflet in ['both', 'outer'] else set(),
        'inner': (set(membrane.inner.ids) - excluded_by_existing['inner']
                if not membrane.monolayer and leaflet in ['both', 'inner']
                else set())
    }

    total_proteins = num_proteins if num_proteins else len(membrane.outer.ids)
    results = {'outer': 0, 'inner': 0}

    logger.info(f"Attempting to place {total_proteins} of protein type {type_id}")
    logger.info(f"Available points for placement:")
    for l, points in available_points.items():
        if points:
            logger.info(f"  {l} leaflet: {len(points)} points")
    logger.info(f"Radius: {radius:.1f} nm" +
               (f", Target curvature: {target_curvature:.3f}" if target_curvature is not None else ""))

    placed = 0
    while placed < total_proteins:
        # Determine valid leaflets (those with available points)
        valid_leaflets = []
        if leaflet in ['both', 'outer'] and available_points['outer']:
            valid_leaflets.append('outer')
        if leaflet in ['both', 'inner'] and available_points['inner']:
            valid_leaflets.append('inner')

        if not valid_leaflets:
            logger.warning("No more valid points available for placement")
            break

        # Randomly choose leaflet
        current_leaflet = rng.choice(valid_leaflets)
        membrane_layer = membrane.outer if current_leaflet == 'outer' else membrane.inner

        # Get valid indices for current leaflet
        valid_indices = np.array(list(available_points[current_leaflet]))

        # Calculate weights based on curvature preference
        curvatures = membrane_layer.mean_curvature[valid_indices]
        if current_leaflet == 'inner':
            curvatures = -curvatures  # Flip curvature for inner leaflet

        weights = calculate_curvature_weights(
            curvatures,
            target_curvature,
            k_factor
        )

        # Choose placement point
        chosen_idx = rng.choice(valid_indices, p=weights)

        # Add protein inclusion
        membrane.inclusions.add_protein(
            type_id=type_id,
            point_id=chosen_idx
        )

        # Update available points accounting for exclusion radius in both leaflets
        nearby = get_nearby_points_both_leaflets(
            membrane,
            current_leaflet,
            chosen_idx,
            radius
        )

        # Update available points for both leaflets
        for leaflet_name, excluded in nearby.items():
            available_points[leaflet_name] -= set(excluded)

        placed += 1
        results[current_leaflet] += 1
        logger.debug(f"Placed protein at point {chosen_idx} in {current_leaflet} leaflet")

    return results

def INU(args: List[str]) -> None:
    """Main entry point for protein inclusion tool"""
    parser = argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-p', '--point-dir', type=Path, default="point",
                       help='Path to membrane point directory (default: point/)')
    parser.add_argument('-t', '--type-id', type=int, required=False,
                       help='Specify which protein type to add into the membrane ')
    parser.add_argument('-r', '--radius', type=float, required=True,
                       help='Exclusion radius for protein placement (point to point distance)')
    parser.add_argument('-c', '--curvature', type=float,
                       help='Target curvature for placement (optional)')
    parser.add_argument('-n', '--num-proteins', type=int,
                       help='Number of proteins to place (optional)')
    parser.add_argument('-k', '--k-factor', type=float, default=1.0,
                       help='Scaling factor for curvature preference strength (default: 1.0)')
    parser.add_argument('-l', '--leaflet', choices=['both', 'inner', 'outer'],
                       default='both', help='Which membrane leaflet(s) to modify (default: both)')
    parser.add_argument('-o', '--output', type=Path,
                       help='Output directory (defaults to input directory)')
    parser.add_argument('--seed', type=int,
                       help='Random seed for reproducibility')

    args = parser.parse_args(args)
    logging.basicConfig(level=logging.INFO)

    if args.num_proteins is None and args.curvature is not None:
        logger.warning(
            "Curvature preference specified without number of proteins to place. "
            "This will attempt to place proteins at all valid points, which may not be desired. "
            "Consider specifying -N/--num-proteins to limit the number of proteins placed."
        )

    # setup numpy random number generator
    global rng
    rng = np.random.default_rng(args.seed)

    try:
        # Load membrane
        membrane = Point(args.point_dir)

        # wrap membrane inside box
        membrane = pbc_wrap(membrane)

        # Place proteins
        results = place_proteins(
            membrane=membrane,
            type_id=args.type_id,
            radius=args.radius,
            num_proteins=args.num_proteins,
            target_curvature=args.curvature,
            k_factor=args.k_factor,
            leaflet=args.leaflet
        )

        # Log results
        total_placed = sum(results.values())
        logger.info(f"Successfully placed {total_placed} of protein type {args.type_id}:")
        for leaflet, count in results.items():
            if count > 0:
                logger.info(f"  {leaflet} leaflet: {count} proteins")

        # Save membrane
        output_dir = args.output or args.point_dir
        membrane.save(output_dir)
        logger.info(f"Updated membrane in {output_dir}")

    except Exception as e:
        logger.error(f"Error: {e}")
        raise
