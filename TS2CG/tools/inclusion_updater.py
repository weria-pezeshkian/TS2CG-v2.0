"""
Tool to place protein inclusions in membrane.
"""

import argparse
import numpy as np
import logging
from typing import List, Optional, Set
from scipy.special import logsumexp
from scipy.spatial import KDTree

from ..core.point import Point

logger = logging.getLogger(__name__)


def get_excluded_points(membrane: Point, radius: float) -> Set[int]:
    """Get points too close to existing proteins."""
    excluded = set()

    if not membrane.inclusions:
        return excluded

    box_dims = np.array(membrane.box)

    # Shift and wrap coordinates for KDTree
    coordinates = membrane.outer.coordinates + box_dims / 2
    coordinates = coordinates % box_dims

    # Build KDTree with periodic boundary conditions
    tree = KDTree(coordinates, boxsize=box_dims)

    # Check exclusions around each existing protein
    for inclusion in membrane.inclusions:
        point_id = inclusion['point_id']
        if point_id < len(coordinates):
            center_coord = coordinates[point_id]
            nearby_indices = tree.query_ball_point(center_coord, radius)
            excluded.update(nearby_indices)

    return excluded


def calculate_placement_weights(curvatures: np.ndarray, target_curvature: Optional[float],
                              k_factor: float) -> np.ndarray:
    """Calculate placement weights based on curvature preference."""
    if target_curvature is None:
        return np.ones_like(curvatures) / len(curvatures)

    deltas = curvatures - target_curvature
    log_weights = -k_factor * deltas**2
    log_weights -= logsumexp(log_weights)
    return np.exp(log_weights)


def place_proteins(membrane: Point, protein_type: int, radius: float, num_proteins: int,
                  target_curvature: Optional[float] = None, k_factor: float = 1.0,
                  seed: Optional[int] = None) -> None:
    """Place proteins in outer membrane leaflet."""

    # Set random seed
    rng = np.random.default_rng(seed)

    # Get box dimensions
    box_dims = np.array(membrane.box)

    # Shift and wrap coordinates for KDTree
    coordinates = membrane.outer.coordinates + box_dims / 2
    coordinates = coordinates % box_dims

    # Get excluded points
    excluded = get_excluded_points(membrane, radius)
    logger.info(f"Placing {num_proteins} proteins of type {protein_type} with radius {radius}")

    placed_total = 0

    while placed_total < num_proteins:
        # Find available placement points
        all_indices = set(range(len(coordinates)))
        available_indices = list(all_indices - excluded)

        if not available_indices:
            logger.warning(f"No more valid placement points. Placed {placed_total}/{num_proteins}")
            break

        # Calculate placement weights
        available_indices = np.array(available_indices)
        curvatures = membrane.outer.mean_curvature[available_indices]
        weights = calculate_placement_weights(curvatures, target_curvature, k_factor)

        # Choose placement point
        chosen_idx = rng.choice(available_indices, p=weights)

        # Add protein
        membrane.inclusions.add_protein(
            type_id=protein_type,
            point_id=chosen_idx
        )

        # Update excluded points using KDTree
        tree = KDTree(coordinates, boxsize=box_dims)
        center_coord = coordinates[chosen_idx]
        nearby_indices = tree.query_ball_point(center_coord, radius)
        excluded.update(nearby_indices)

        placed_total += 1
        logger.info(f"Placed protein {placed_total} at point {chosen_idx}")

    logger.info(f"Finished placing {placed_total} proteins")


def INU(args: List[str]) -> None:
    """Main entry point for protein inclusion tool"""
    parser = argparse.ArgumentParser(description="Place protein inclusions in outer membrane")

    parser.add_argument('--point-dir', '-p', default="point/",
                       help="Path to point folder (default: point/)")
    parser.add_argument('--protein-type', '-t', type=int, required=True,
                       help="Protein type ID to place")
    parser.add_argument('--radius', '-r', type=float, required=True,
                       help="Exclusion radius for protein placement")
    parser.add_argument('--num-proteins', '-n', type=int, required=True,
                       help="Number of proteins to place")
    parser.add_argument('--curvature', '-c', type=float,
                       help="Target curvature for placement (optional)")
    parser.add_argument('--k-factor', '-k', type=float, default=1.0,
                       help="Curvature preference strength (default: 1.0)")
    parser.add_argument('--output-dir', '-o',
                       help="Output directory (default: overwrite input with backup)")
    parser.add_argument('--no-backup', action='store_true',
                       help="Skip creating backup when overwriting input")
    parser.add_argument('--seed', type=int,
                       help="Random seed for reproducibility")

    args = parser.parse_args(args)

    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    # Load membrane
    logger.info(f"Loading membrane from {args.point_dir}")
    membrane = Point(args.point_dir)

    # Place proteins
    place_proteins(
        membrane=membrane,
        protein_type=args.protein_type,
        radius=args.radius,
        num_proteins=args.num_proteins,
        target_curvature=args.curvature,
        k_factor=args.k_factor,
        seed=args.seed
    )

    # Save results
    output_dir = args.output_dir if args.output_dir else args.point_dir
    backup = not args.no_backup
    membrane.save(output_dir, backup=backup)

    logger.info(f"Finished updating membrane inclusions!")
