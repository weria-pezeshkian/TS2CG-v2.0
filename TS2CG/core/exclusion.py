from typing import Optional, Dict, Union, List

import numpy as np

class Exclusion:
    """Manages lipid exclusions in the membrane."""

    def __init__(self, data: Optional[np.ndarray] = None):
        self.exclusions = []
        self._process_data(data)

    def __getitem__(self, ndx):
        return self.exclusions[ndx]

    def __len__(self):
        return len(self.exclusions)

    def _process_data(self, data: np.ndarray):
        """Process exclusion data."""

        if data is None or len(data) == 0:
            return
        else:
            if data.shape[0] != 3:
                msg = f"Exclusion definition in ExcData.dat has wrong size. Expected (3,N), got {data.shape}."
                raise ValueError(msg)
        try:
            data.shape[1]
        except IndexError:
            data=data.reshape(3,1)
        for i in range(data.shape[1]):
            point = {
                'id': int(data[0,i]),
                'point_id': int(data[1,i]),
                'radius': float(data[2,i])
            }
            self.exclusions.append(point)

    def get_all(self) -> List[int]:
        """Get all excluded exclusions."""
        return [e['point_id'] for e in self.exclusions]

    def add_pore(self, point_id: int, radius: float = 1.0):
        """
        Add a pore in the lipid membrane.

        Args:
            point_id: Point ID where lipids should be excluded
            radius: Radius of exclusion zone
        """
        exclusion = {
            'id': len(self.exclusions),
            'point_id': point_id,
            'radius': radius
        }
        self.exclusions.append(exclusion)

    def remove_pore(self, index: int):
        """
        remove a pore in the lipid membrane.

        Args:
            index: int
                index of exclusion to remove
        """
        del self.exclusions[index]

    def remove_pores(self, indices: List[int]):
        """
        remove a list of pores from the lipid membrane.

        Args:
            indices: List[int]
                list of indices of inclusions to remove
        """
        for index in indices:
            self.remove_pore[index]
