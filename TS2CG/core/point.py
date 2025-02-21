import logging
from typing import Optional, Dict, Union, List

import shutil
from io import StringIO
from pathlib import Path

import numpy as np

from .inclusion import Inclusion
from .exclusion import Exclusion
from .membrane import Membrane

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
    """

    @staticmethod
    def _ensure_path(path: Optional[Union[str, Path]]) -> Optional[Path]:
        """Convert string path to Path object if needed."""
        if path is None:
            return None
        else:
            return Path(path) if isinstance(path, str) else path

    def __init__(self, path: Union[str, Path]):
        """
        Initialize point class from a point folder.

        Args:
            path: pathlib
                Path to the point folder
        """
        self.path = self._ensure_path(path)
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
            self.outer = Membrane(outer_data)
            self.inner = None if self.monolayer else Membrane(inner_data)

            # Load modifications data
            inc_data = self._load_modification_file("IncData.dat")
            exc_data = self._load_modification_file("ExcData.dat")
            self.inclusions = Inclusion(inc_data)
            self.exclusions = Exclusion(exc_data)

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
            return
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
        data = np.zeros((18, len(membrane.ids)))
        data[0] = membrane.ids
        data[1] = membrane.domain_ids
        data[2] = membrane.area
        data[3:6] = membrane.coordinates.T
        data[6:9] = membrane.normals.T
        data[9:12] = membrane.principal_vectors['p1'].T
        data[12:15] = membrane.principal_vectors['p2'].T
        data[15] = membrane.curvature['c1']
        data[16] = membrane.curvature['c2']
        data[17] = membrane.edges.astype(int)

        # Create header
        headers = []

        # Add box dimensions for OuterBM.dat
        if "Outer" in output_path.name:
            headers.append(f"Box     {self.box[0]:.3f}     {self.box[1]:.3f}     {self.box[2]:.3f}")

        headers.extend([
            f"< Point NoPoints     {len(membrane.ids)}>",
            "< id domain_id area X Y Z Nx Ny Nz P1x P1y P1z P2x P2y P2z C1 C2 vtype >",
            f"< {'Outer' if 'Outer' in output_path.name else 'Inner'} >"
        ])

        header = '\n'.join(headers)

        np.savetxt(
            output_path,
            data.T,
            header=header,
            comments='',
            fmt = ['%10d', '%4d', '%9.3f'] + ['%9.3f']*3 + ['%7.3f']*11 + ['%9d']
        )
        logger.info(f"Saved {len(membrane.ids)} points to {output_path.name}")

    def _save_modifications(self, output_path: Path):
        """Save inclusion and exclusion data to files."""
        # Save inclusions
        if self.inclusions:
            data = np.zeros((6, len(self.inclusions)))
            for i, point in enumerate(self.inclusions):
                data[0:3, i] = [point['id'], point['type_id'], point['point_id']]
                data[3:6, i] = point['orientation']

            header = f"< Inclusion NoInc {len(self.inclusions)} >\n"
            header += "< id typeid pointid lx ly lz >"

            np.savetxt(
                output_path / "IncData.dat",
                data.T,
                header=header,
                comments='',
                fmt=['%12d', '%12d', '%12d', '%8.3f', '%8.3f', '%8.3f']
            )

        # Save exclusions
        if self.exclusions:
            data = np.zeros((3, len(self.exclusions)))
            for i, point in enumerate(self.exclusions):
                data[0:3, i] = [point['id'], point['point_id'], point['radius']]

            header = f"< Exclusion NoExc {len(self.exclusions)} >\n"
            header += "< id typeid radius >"

            np.savetxt(
                output_path / "ExcData.dat",
                data.T,
                header=header,
                comments='',
                fmt=['%12d', '%12d', '%12d']
            )

    def save(self, output_path: Optional[Union[str, Path]] = None, backup: Optional[bool] = True):
        """
        Save membrane structure to files.

        Args:
            output_path: pathlib
                Path where to save the point folder. If None, saves to original location.
                Backup is only created if saving to the original location.
            backup: bool
                wheter or not to write output files
        """
        output_path = self._ensure_path(output_path) if output_path else self.path

        # Create backup only if we're overwriting the original folder
        if output_path == self.path and backup:
            self._create_backup()

        # Create output directory if it doesn't exist
        output_path.mkdir(parents=True, exist_ok=True)

        self._save_membranes(output_path)
        self._save_modifications(output_path)

        logger.info(f"Saved point data to: {output_path}")
