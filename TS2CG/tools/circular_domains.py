"""
Tool to assign circular domains around inclusions or specific points.
"""

import argparse
import numpy as np
import logging
from typing import List, Optional
from scipy.spatial import KDTree
import networkx as nx

from ..core.point import Point

logger = logging.getLogger(__name__)


def get_domain_centers(membrane: Point, protein_type: Optional[int] = None,
                      manual_points: Optional[List[int]] = None) -> List[int]:
    """Get domain center points from protein inclusions and/or manually specified points."""
    domain_centers = []

    if protein_type is not None:
        protein_points = [inc["point_id"] for inc in membrane.inclusions.get_by_type(protein_type)]
        domain_centers.extend(protein_points)
        logger.info(f"Added {len(protein_points)} domain centers on proteins of type {protein_type}")

    if manual_points:
        domain_centers.extend(manual_points)
        logger.info(f"Added {len(manual_points)} manually specified domain centers")

    unique_centers = list(set(domain_centers))
    return unique_centers


def assign_circular_domains(membrane: Point, radius: float, domain_centers: List[int],
                          domain_id: int, leaflet: str = "both", use_path_distance: bool = False,
                          edge_cutoff: float = 5) -> None:
    """Assign domain IDs to all points within radius of domain centers."""

    # Determine which leaflets to process
    leaflets_to_process = []
    if membrane.monolayer or leaflet == "outer":
        leaflets_to_process = [("outer", membrane.outer)]
    elif leaflet == "inner":
        leaflets_to_process = [("inner", membrane.inner)]
    elif leaflet == "both":
        leaflets_to_process = [("outer", membrane.outer)]
        if not membrane.monolayer:
            leaflets_to_process.append(("inner", membrane.inner))

    for leaflet_name, membrane_leaflet in leaflets_to_process:
        logger.info(f"Processing {leaflet_name} leaflet")

        # Build KDTree from point coordinates
        box_dims = np.array(membrane.box)
        coordinates = membrane_leaflet.coordinates
        tree = KDTree(coordinates, boxsize=box_dims)

        for center_idx in domain_centers:
            if center_idx >= len(coordinates):
                logger.warning(f"Domain center {center_idx} does not exist in {leaflet_name} leaflet")
                continue

            if use_path_distance:
                logger.info(f"Using path distance (Dijkstra) for domain assignment")

                # Create sparse distance matrix
                dist_matrix = KDTree.sparse_distance_matrix(tree, tree, max_distance=edge_cutoff)

                # Convert to networkx graph
                G = nx.from_scipy_sparse_array(dist_matrix)

                # Find shortest paths within radius
                shortest_paths = nx.single_source_dijkstra_path_length(
                    G, center_idx, cutoff=radius
                )
                # Convert back to original indices
                points_within_radius = list(shortest_paths.keys())
            else:
                logger.info("Using euclidean distance (KDTree) for domain assignment")
                # Simple KDTree query for euclidean distance
                points_within_radius = tree.query_ball_point(coordinates[center_idx], radius)

            # Assign domain ID
            if len(points_within_radius) > 0:
                membrane_leaflet.domain_ids[points_within_radius] = domain_id
                logger.info(f"Assigned domain {domain_id} to {len(points_within_radius)} points")

        logger.info(f"Finished processing {leaflet_name} leaflet")


def DAI(args: List[str]) -> None:
    """Main entry point for circular domain assignment tool"""
    parser = argparse.ArgumentParser(description="Assign circular domains around inclusions or points")

    parser.add_argument('--point-dir', '-p', default="point/",
                       help="Path to point folder (default: point/)")
    parser.add_argument('--radius', '-r', type=float, required=True,
                       help="Radius for circular domain assignment")
    parser.add_argument('--domain-id', '-d', type=int, required=True,
                       help="Domain ID to assign to points within radius")
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
    parser.add_argument('--path-distance', action='store_true',
                       help="Use path distance (Dijkstra) instead of euclidean distance for curved membranes")
    parser.add_argument('--edge-cutoff', type=float, default=5.0,
                        help="Maximum distance (nm) for graph edges in path distance mode. "
                             "Edges longer than this are excluded, preventing shortcuts across membrane folds")

    args = parser.parse_args(args)

    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    # Load membrane
    logger.info(f"Loading membrane from {args.point_dir}")
    membrane = Point(args.point_dir)

    # Parse manual points
    manual_points = None
    if args.manual_points:
        try:
            manual_points = [int(p.strip()) for p in args.manual_points.split(',')]
        except ValueError:
            parser.error("Manual points must be comma-separated integers")

    # Get domain centers
    domain_centers = get_domain_centers(membrane, args.protein_type, manual_points)

    if not domain_centers:
        parser.error("No domain centers found. You must specify either --protein-type or --manual-points for domain assignment")

    # Assign domains
    assign_circular_domains(membrane, args.radius, domain_centers, args.domain_id,
                           args.leaflet, args.path_distance, args.edge_cutoff)

    # Save results
    output_dir = args.output_dir if args.output_dir else args.point_dir
    backup = not args.no_backup
    membrane.save(output_dir, backup=backup)

    logger.info(f"Finished with updating membrane domains!")
