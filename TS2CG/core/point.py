import logging
from typing import Optional, Dict, Union, List

import shutil
from io import StringIO
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

def loadtxt_fix(filename, skiprows):
    # needed because of bad formatting out of PLM
    with open(filename, 'r') as f:
        lines = [' '.join(line.replace('-', ' -').split()) + '\n' for line in f]
        return np.loadtxt(StringIO(''.join(lines)), skiprows=skiprows)

class Point:
    """
    A class representing a membrane structure with inclusions and exclusions.
    Can be initialized from a point folder or built from scratch.

    Examples:
        # Create bilayer
        >>> membrane = Point.create_empty()
        >>> membrane.add_membrane_points(coordinates, normals)

        # Create monolayer
        >>> monolayer = Point.create_empty(monolayer=True)
        >>> monolayer.add_membrane_points(coordinates, normals)
    """

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
        @classmethod
        def create_empty(cls):
            """Create an empty membrane layer."""
            instance = cls.__new__(cls)
            instance.ids = np.array([], dtype=int)
            instance.domain_ids = np.array([], dtype=int)
            instance.area = np.array([])
            instance.coordinates = np.array([]).reshape(0, 3)
            instance.normals = np.array([]).reshape(0, 3)
            instance.principal_vectors = {
                'p1': np.array([]).reshape(0, 3),
                'p2': np.array([]).reshape(0, 3)
            }
            instance.curvature = {
                'c1': np.array([]),
                'c2': np.array([])
            }
            return instance

        def add_points(self, coordinates: np.ndarray, normals: Optional[np.ndarray] = None,
                      domain_ids: Optional[np.ndarray] = None, areas: Optional[np.ndarray] = None):
            """
            Add points to the membrane layer.

            Args:
               coordinates: Nx3 array of point coordinates
                normals: Nx3 array of normal vectors (optional)
                domain_ids: Array of domain IDs (optional)
                areas: Array of point areas (optional)
            """
            n_points = len(coordinates)

            # Set IDs
            self.ids = np.arange(n_points, dtype=int)

            # Set coordinates
            self.coordinates = np.asarray(coordinates)

            # Set normals or generate default
            if normals:
                self.normals = np.asarray(normals)
            else:
                self.normals = np.zeros_like(coordinates)
                self.normals[:, 2] = 1.0  # Default normal pointing up

            # Set domain IDs or default to 0
            if domain_ids:
                self.domain_ids = np.asarray(domain_ids)
            else:
                self.domain_ids = np.zeros(n_points, dtype=int)

            # Set areas or default to 1.0
            if areas:
                self.area = np.asarray(areas)
            else:
                self.area = np.ones(n_points)

            # Initialize principal vectors and curvatures
            self.principal_vectors['p1'] = np.zeros_like(coordinates)
            self.principal_vectors['p2'] = np.zeros_like(coordinates)
            self.curvature['c1'] = np.zeros(n_points)
            self.curvature['c2'] = np.zeros(n_points)

        @property
        def mean_curvature(self) -> np.ndarray:
            """Calculate mean curvature for all points."""
            return (self.curvature['c1'] + self.curvature['c2']) / 2


        def gaussian_curvature(self) -> np.ndarray:
            """Calculate Gaussian curvature for all points."""
            return self.curvature['c1'] * self.curvature['c2']


        def get_points_by_domain(self, domain_id: int) -> np.ndarray:
            """Get coordinates of all points in a specific domain."""
            mask = self.domain_ids == domain_id
            return self.coordinates[mask]

    class Inclusion:
        """Manages protein inclusions in the membrane."""

        def __init__(self, data: Optional[np.ndarray] = None):
            self.points = []
            self._process_data(data)

        def _process_data(self, data: np.ndarray):
            """Process inclusion data."""
            if data is None or len(data) == 0:
                return
            else:
                if data.shape[0] != 6:
                    msg = f"Inclusion definition in IncData.dat has wrong size. Expected (6,N), got {data.shape}"
                    raise ValueError(msg)

            for i in range(data.shape[1]):
                point = {
                    'id': int(data[0, i]),
                    'type_id': int(data[1, i]),
                    'point_id': int(data[2, i]),
                    'orientation': data[3:6, i]
                }
                self.points.append(point)

        @property
        def count(self) -> int:
            """Number of protein inclusions."""
            return len(self.points)

        def get_all(self) -> List[dict]:
            """Get all points with protein inclusions."""
            return [p['point_id'] for p in self.points]

        def get_by_type(self, type_id: int) -> List[dict]:
            """Get all inclusions of a specific type."""
            return [p for p in self.points if p['type_id'] == type_id]

        def add_protein(self, type_id: int, point_id: int,
                       orientation: Optional[np.ndarray] = np.array([1, 0, 0])):
            """
            Add a protein inclusion.

            Args:
                type_id: Type identifier for the protein
                point_id: Point ID where protein should be placed
                orientation: Vector specifying protein orientation
            """
            if orientation is None:
                orientation = np.array([0, 0, 1])

            orientation=orientation/np.linalg.norm(orientation)

            point = {
                'id': len(self.points),
                'type_id': type_id,
                'point_id': point_id,
                'orientation': orientation
            }
            self.points.append(point)

    class Exclusion:
        """Manages lipid exclusions in the membrane."""

        def __init__(self, data: Optional[np.ndarray] = None):
            self.points = []
            if data is not None:
                self._process_data(data)

        def _process_data(self, data: np.ndarray):
            """Process exclusion data."""
            if data is None or len(data) == 0:
                return
            if len(data.shape)==1:
                data.shape = (3,1)
            for i in range(data.shape[1]):
                point = {
                    'id': int(data[0,i]),
                    'point_id': int(data[1,i]),
                    'radius': float(data[2,i])
                }
                self.points.append(point)

        @property
        def count(self) -> int:
            """Number of exclusion points."""
            return len(self.points)

        def get_all(self) -> List[int]:
            """Get all excluded points."""
            return [p['point_id'] for p in self.points]

        def add_pore(self, point_id: int, radius: float = 1.0):
            """
            Add a pore in the lipid membrane.

            Args:
                point_id: Point ID where lipids should be excluded
                radius: Radius of exclusion zone
            """
            point = {
                'id': len(self.points),
                'point_id': point_id,
                'radius': radius
            }
            self.points.append(point)

    def __init__(self, path: Union[str, Path]):
        """
        Initialize point class from a point folder.

        Args:
            path: Path to the point folder
        """
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"Point folder not found: {self.path}")

        self._load_data()

    def _load_data(self):
        """Load all data from the point folder."""
        try:
            # Try to load outer membrane (required)
            outer_data = self._load_membrane_file(self.path / "OuterBM.dat")
            if outer_data is None:
                raise FileNotFoundError("OuterBM.dat can not be parsed!")

            # Set monolayer flag based on presence of InnerBM.dat
            inner_data = self._load_membrane_file(self.path / "InnerBM.dat")
            self.monolayer = inner_data is None

            # Create membrane instances
            self.outer = self.Membrane(outer_data)
            self.inner = None if self.monolayer else self.Membrane(inner_data)

            # Load modifications data
            inc_data = self._load_modification_file("IncData.dat")
            exc_data = self._load_modification_file("ExcData.dat")
            self.inclusions = self.Inclusion(inc_data)
            self.exclusions = self.Exclusion(exc_data)

        except Exception as e:
            logger.error("Failed to load membrane data", exc_info=True)
            raise

    def _load_membrane_file(self, file_path: Path) -> Optional[np.ndarray]:
        """
        Load membrane definition file.
        Returns None if file doesn't exist.
        """

        if not file_path.exists():
            logger.info(f"Membrane file {file_path.name} not found")
            return None

        try:
            # Read the first few lines to check for Box information
            with open(file_path) as f:
                first_lines = [next(f) for _ in range(4)]

            # Store box dimensions if this is OuterBM.dat
            if "OuterBM" in file_path.name:
                self.box = self._parse_box_line(first_lines[0])
                skiprows = 4
            else:
                skiprows = 3

            return loadtxt_fix(file_path, skiprows).T

        except Exception as e:
            logger.warning(f"Error loading {file_path.name}: {e}")
            return None

    def _parse_box_line(self, line: str) -> tuple:
        """Parse box dimensions from header line."""
        parts = line.split()
        return np.array([float(x) for x in parts[1:4]])

    def _load_modification_file(self, filename: str) -> Optional[np.ndarray]:
        """Load modification (inclusion/exclusion) file."""
        try:
            return loadtxt_fix(self.path / filename, skiprows=2).T
        except (ValueError, FileNotFoundError):
            return None

    @staticmethod
    def _ensure_path(path: Optional[Union[str, Path]]) -> Optional[Path]:
        """Convert string path to Path object if needed."""
        if path is None:
            return None
        return Path(path) if isinstance(path, str) else path

    def save(self, output_path: Optional[Union[str, Path]] = None):
        """
        Save membrane structure to files.

        Args:
            output_path: Path where to save the point folder. If None, saves to original location.
                        Backup is only created if saving to the original location.
        """
        output_path = self._ensure_path(output_path) if output_path else self.path

        # Create backup only if we're overwriting the original folder
        if output_path == self.path:
            self._create_backup()

        # Create output directory if it doesn't exist
        output_path.mkdir(parents=True, exist_ok=True)

        self._save_membranes(output_path)
        self._save_modifications(output_path)

        logger.info(f"Saved point data to: {output_path}")

    def _create_backup(self):
        """
        Create a backup of the point folder.
        If #folder# exists, creates ##folder##, etc.
        """
        def get_backup_path(n_hashes: int) -> Path:
            """Generate backup path with specified number of hashes."""
            hashes = '#' * n_hashes
            return self.path.parent / f"{hashes}{self.path.name}{hashes}"

        # Start with one hash on each side
        n_hashes = 1
        backup_path = get_backup_path(n_hashes)

        # Keep incrementing hashes until we find a non-existing path
        while backup_path.exists():
            n_hashes += 1
            backup_path = get_backup_path(n_hashes)

        shutil.copytree(self.path, backup_path)
        logger.info(f"Created backup at: {backup_path}")

        return backup_path

    def _save_membranes(self, output_path: Path):
        """Save membrane data to files."""
        # Always save outer membrane
        if len(self.outer.ids) > 0:
            self._save_single_membrane(output_path / "OuterBM.dat", self.outer)

        # Save inner membrane only for bilayers
        if not self.monolayer and self.inner is not None and len(self.inner.ids) > 0:
            self._save_single_membrane(output_path / "InnerBM.dat", self.inner)

    def _save_single_membrane(self, output_path: Path, membrane):
        """Helper method to save a single membrane layer."""
        data = np.zeros((17, len(membrane.ids)))
        data[0] = membrane.ids
        data[1] = membrane.domain_ids
        data[2] = membrane.area
        data[3:6] = membrane.coordinates.T
        data[6:9] = membrane.normals.T
        data[9:12] = membrane.principal_vectors['p1'].T
        data[12:15] = membrane.principal_vectors['p2'].T
        data[15] = membrane.curvature['c1']
        data[16] = membrane.curvature['c2']

        # Create header
        headers = []

        # Add box dimensions for OuterBM.dat
        if "Outer" in output_path.name:
            headers.append(f"Box     {self.box[0]:.3f}     {self.box[1]:.3f}     {self.box[2]:.3f}")

        headers.extend([
            f"< Point NoPoints     {len(membrane.ids)}>",
            "< id domain_id area X Y Z Nx Ny Nz P1x P1y P1z P2x P2y P2z C1 C2  >",
            f"< {'Outer' if 'Outer' in output_path.name else 'Inner'} >"
        ])

        header = '\n'.join(headers)

        np.savetxt(
            output_path,
            data.T,
            header=header,
            comments='',
            fmt = ['%10d', '%4d', '%9.3f'] + ['%9.3f']*3 + ['%7.3f']*11
        )
        logger.info(f"Saved {len(membrane.ids)} points to {output_path.name}")

    def _save_modifications(self, output_path: Path):
        """Save inclusion and exclusion data to files."""
        # Save inclusions
        if self.inclusions.points:
            data = np.zeros((6, len(self.inclusions.points)))
            for i, point in enumerate(self.inclusions.points):
                data[0:3, i] = [point['id'], point['type_id'], point['point_id']]
                data[3:6, i] = point['orientation']

            header = f"< Inclusion NoInc {len(self.inclusions.points)} >\n"
            header += "< id typeid pointid lx ly lz >"

            np.savetxt(
                output_path / "IncData.dat",
                data.T,
                header=header,
                comments='',
                fmt=['%12d', '%12d', '%12d', '%8.3f', '%8.3f', '%8.3f']
            )

        # Save exclusions
        if self.exclusions.points:
            data = np.zeros((3, len(self.exclusions.points)))
            for i, point in enumerate(self.exclusions.points):
                data[0:3, i] = [point['id'], point['point_id'], point['radius']]

            header = f"< Exclusion NoExc {len(self.exclusions.points)} >\n"
            header += "< id typeid radius >"

            np.savetxt(
                output_path / "ExcData.dat",
                data.T,
                header=header,
                comments='',
                fmt=['%12d', '%12d', '%12d']
            )

    def update_domains(self, domain_ids: Optional[np.ndarray] = None):
        """
        Update domain assignments for membrane layer(s).
        For bilayers, updates both leaflets. For monolayers, updates only the outer leaflet.

        Args:
            domain_ids: New domain assignments as numpy array
        """
        if domain_ids is None:
            logger.warning("No domain IDs provided for update")
            return

        # Update outer membrane (always present)
        if len(domain_ids) != len(self.outer.ids):
            raise ValueError(
                f"Domain IDs length ({len(domain_ids)}) does not match "
                f"number of membrane points ({len(self.outer.ids)})"
            )

        self.outer.domain_ids = np.asarray(domain_ids, dtype=int)

        # Update inner membrane if this is a bilayer
        if not self.monolayer and self.inner is not None:
            self.inner.domain_ids = np.asarray(domain_ids, dtype=int)
            logger.debug("Updated domains for both leaflets")
        else:
            logger.debug("Updated domains for outer leaflet only")

    @classmethod
    def create_empty(cls, box, monolayer=False):
        """
        Create an empty Point instance.

        Args:
            box: Tuple of (x, y, z) box dimensions
            monolayer: If True, only creates outer membrane layer
        """
        instance = cls.__new__(cls)
        instance.path = None
        instance.box = box
        instance.monolayer = monolayer

        # Initialize outer membrane (always present)
        instance.outer = cls.Membrane.create_empty()

        # Initialize inner membrane only for bilayers
        instance.inner = None if monolayer else cls.Membrane.create_empty()

        # Initialize modifications
        instance.inclusions = cls.Inclusion()
        instance.exclusions = cls.Exclusion()

        return instance

    def add_membrane_points(self, coordinates: np.ndarray, normals: Optional[np.ndarray] = None,
                          domain_ids: Optional[np.ndarray] = None, areas: Optional[np.ndarray] = None,
                          bilayer_spacing: float = 4.0):
        """
        Add points to membrane layer(s) with proper bilayer spacing.

        Args:
            coordinates: Nx3 array of point coordinates (midplane coordinates)
            normals: Nx3 array of normal vectors (optional)
            domain_ids: Array of domain IDs (optional)
            areas: Array of point areas (optional)
            bilayer_spacing: Distance between leaflets in nm (default=4.0, only used for bilayers)
        """
        # Generate default normals if not provided (pointing up)
        if normals is None:
            normals = np.zeros_like(coordinates)
            normals[:, 2] = 1.0

        # Normalize the normal vectors
        normals = normals / np.linalg.norm(normals, axis=1)[:, np.newaxis]

        if self.monolayer:
            # For monolayer, use coordinates directly for outer membrane
            self.outer.add_points(coordinates, normals, domain_ids, areas)
        else:
            # For bilayer, offset both leaflets from the midplane
            half_spacing = bilayer_spacing / 2

            # Add outer leaflet
            outer_coordinates = coordinates + normals * half_spacing
            self.outer.add_points(outer_coordinates, normals, domain_ids, areas)

            # Add inner leaflet
            inner_coordinates = coordinates - normals * half_spacing
            self.inner.add_points(inner_coordinates, -normals, domain_ids, areas)

    def add_lipids(self, coordinates: np.ndarray, normals: np.ndarray,
                domain_ids: Optional[np.ndarray] = None, areas: Optional[np.ndarray] = None,
                bilayer_spacing: float = 4.0):
        """
        convenience method to add lipids to the membrane.
        adds points to both leaflets offset by bilayer_spacing along the normal vectors.

        args:
            coordinates: nx3 array of midplane lipid positions
            domain_ids: array of domain ids (optional)
            layer: which layer to add lipids to ("both", "inner", or "outer")
            bilayer_spacing: distance between inner and outer leaflets in nm (default=4.0)
        """
        self.add_membrane_points(
            coordinates,
            normals=normals,
            domain_ids=domain_ids,
            areas=areas,
            bilayer_spacing=bilayer_spacing
        )
