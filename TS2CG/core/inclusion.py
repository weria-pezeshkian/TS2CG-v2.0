from typing import Optional, Dict, Union, List

import numpy as np

class Inclusion:
    """Manages protein inclusions in the membrane."""

    def __init__(self, data: Optional[np.ndarray] = None):
        self.inclusions = []
        self._process_data(data)

    def __getitem__(self, ndx):
        return self.inclusions[ndx]

    def __len__(self):
        return len(self.inclusions)

    def _process_data(self, data: np.ndarray):
        """Process inclusion data."""
        if data is None or len(data) == 0:
            return
        else:
            if data.shape[0] != 6:
                msg = f"Inclusion definition in IncData.dat has wrong size. Expected (6,N), got {data.shape}."
                raise ValueError(msg)

        for i in range(data.shape[1]):
            inclusion = {
                'id': int(data[0, i]),
                'type_id': int(data[1, i]),
                'point_id': int(data[2, i]),
                'orientation': data[3:6, i]
            }
            self.inclusions.append(inclusion)

    def get_all(self) -> List[int]:
        """Get all inclusion with protein inclusions."""
        return [i['point_id'] for i in self.inclusions]

    def get_by_type(self, type_id: int) -> List[dict]:
        """Get all inclusions of a specific type."""
        return [i for i in self.inclusions if i['type_id'] == type_id]

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

        orientation = orientation / np.linalg.norm(orientation)

        inclusion = {
            'id': len(self.inclusions),
            'type_id': type_id,
            'point_id': point_id,
            'orientation': orientation
        }
        self.inclusions.append(inclusion)

    def remove_protein(self, index: int):
        """
        remove a protein inclusion.

        Args:
            index: int
                index of inclusion to remove
        """
        del self.inclusions[index]

    def remove_proteins(self, indices: List[int]):
        """
        remove a list of protein inclusions.

        Args:
            indices: List[int]
                list of indices of inclusions to remove
        """
        for index in indices:
            self.remove_protein[index]
