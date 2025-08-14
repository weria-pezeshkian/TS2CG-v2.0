"""
Tool to assign circular domains around inclusions or specific points.
"""

import argparse
import numpy as np
import logging
from typing import List, Optional
from scipy.spatial import KDTree

from ..core.point import Point

logger = logging.getLogger(__name__)


def get_domain_centers(membrane: Point, protein_type: Optional[int] = None,
                      manual_points: Optional[List[int]] = None) -> List[int]:
    """Get domain center points from protein inclusions and/or manually specified points."""
    domain_centers = []

    if protein_type is not None:
        protein_points = [inc["point_id"] for inc in membrane.inclusions.get_by_type(protein_type)]
        domain_centers.extend(protein_points)
        logger.info(f"Added {len(domain_centers)} domain centers on proteins of type {protein_type}")

    if manual_points:
        domain_centers.extend(manual_points)
        logger.info(f"Added {len(manual_points)} manually specified domain centers")

    unique_centers = list(set(domain_centers))
    logger.info(f"Total domain centers: {len(unique_centers)}")

    return unique_centers

def assign_circular_domains(membrane: Point, radius: float, domain_centers: List[int],
                          domain_id: int, leaflet: str = "both") -> None:
    """Assign domain IDs to all points within radius of domain centers."""

    # Get box dimensions
    box_dims = np.array(membrane.box)

    # Determine which leaflets to process
    leaflets_to_process = []
    if membrane.monoleaflet or leaflet.lower() == "outer":
        leaflets_to_process = [("outer", membrane.outer)]
    elif leaflet.lower() == "inner":
        leaflets_to_process = [("inner", membrane.inner)]
    elif leaflet.lower() == "both":
        leaflets_to_process = [("outer", membrane.outer), ("inner", membrane.inner)]

    for leaflet_name, membrane_leaflet in leaflets_to_process:
        # Shift and wrap coordinates for KDTree
        coordinates = membrane_leaflet.coordinates + box_dims / 2
        coordinates = coordinates % box_dims

        # Build KDTree with periodic boundary conditions
        tree = KDTree(coordinates, boxsize=box_dims)

        for center_idx in domain_centers:
            if center_idx >= len(membrane_leaflet.coordinates):
                logger.warning(f"Domain center {center_idx} does not exist in {leaflet_name} leaflet")
                continue

            # Use wrapped coordinates for center query
            center_coord = coordinates[center_idx]

            # Query ball points around center with periodic boundaries
            within_radius_indices = tree.query_ball_point(center_coord, radius)

            # Assign domain ID
            if len(within_radius_indices) > 0:
                within_radius_indices = np.array(within_radius_indices)
                membrane_leaflet.domain_ids[within_radius_indices] = domain_id

        logger.info(f"Finished processing {leaflet_name} leaflet")

def DAI(args: List[str]) -> None:
    """Main entry point for circular domain assignment tool"""
    parser = argparse.ArgumentParser(description="Assign circular domains around inclusions or points")

    parser.add_argument('--point-dir', '-p', default="point/",
                       help="Path to point folder (default: point/)")
    parser.add_argument('--radius', '-r', type=float, required=True,
                       help="Radius for circular domain assignment")
    parser.add_argument('--domain-id', '-d', type=int,
                       help="Domain ID to assign to points within radius (required for domain assignment)")
    parser.add_argument('--protein-type', '-t', type=int,
                       help="Protein type ID to use as domain centers")
    parser.add_argument('--manual-points', '-m',
                       help="Comma-separated point IDs to use as domain centers (e.g., '3,7,22')")
    parser.add_argument('--leaflet', '-l', choices=['outer', 'inner', 'both'], default='both',
                       help="Which membrane leaflet to modify (default: both)")
    parser.add_argument('--output-dir', '-o',
                       help="Output directory (default: overwrite input with backup)")
    parser.add_argument('--no-backup', action='store_true',
                       help="Skip creating backup when overwriting input")

    args = parser.parse_args(args)

    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    # Load membrane
    logger.info(f"Loading membrane from {args.point_dir}")
    membrane = Point(args.point_dir)

    # Validate inputs for domain assignment
    if args.protein_type is None and args.manual_points is None:
        raise ValueError("Must specify either --protein-type or --manual-points for domain assignment")

    if args.domain_id is None:
        raise ValueError("Must specify --domain-id for domain assignment")

    # Parse manual points
    manual_points = None
    if args.manual_points:
        try:
            manual_points = [int(p.strip()) for p in args.manual_points.split(',')]
        except ValueError:
            raise ValueError("Manual points must be comma-separated integers")

    # Get domain centers
    domain_centers = get_domain_centers(membrane, args.protein_type, manual_points)

    if not domain_centers:
        raise ValueError("No domain centers found")

    # Assign domains
    assign_circular_domains(membrane, args.radius, domain_centers, args.domain_id, args.leaflet)

    # Save results
    output_dir = args.output_dir if args.output_dir else args.point_dir
    backup = not args.no_backup
    membrane.save(output_dir, backup=backup)

    logger.info(f"Finished with updating membrane domains!")
