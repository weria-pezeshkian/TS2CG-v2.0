from typing import Optional, Dict, Union, List

import numpy as np

class Membrane:
    """Represents a membrane layer with associated properties."""

    def __init__(self, data: np.ndarray):
        """Initialize membrane layer from raw data."""
        self._process_data(data)

    def _process_data(self, data: np.ndarray):
        """Convert raw data into accessible properties."""
        # Basic properties
        self.ids = data[0].astype(int)
        self.domain_ids = data[1].astype(int)
        self.area = data[2]

        # Coordinates and vectors
        self.coordinates = data[3:6].T  # X, Y, Z
        self.normals = data[6:9].T      # Nx, Ny, Nz
        self.principal_vectors = {
            'p1': data[9:12].T,   # P1x, P1y, P1z
            'p2': data[12:15].T    # P2x, P2y, P2z
        }
        self.curvature = {
            'c1': data[15],        # First principal curvature
            'c2': data[16]         # Second principal curvature
        }
        self.edges = data[17].astype(bool) # Boolean array if point is edge

    @property
    def mean_curvature(self) -> np.ndarray:
        """Calculate mean curvature for all points."""
        return (self.curvature['c1'] + self.curvature['c2']) / 2

    @property
    def gaussian_curvature(self) -> np.ndarray:
        """Calculate Gaussian curvature for all points."""
        return self.curvature['c1'] * self.curvature['c2']

    def get_points_by_domain(self, domain_id: int) -> np.ndarray:
        """Get coordinates of all points in a specific domain."""
        mask = self.domain_ids == domain_id
        return self.coordinates[mask]

    def get_edge_ids(self) -> np.ndarray:
        """Get ids of all points that are an edge."""
        return np.where(self.edges)[0]
