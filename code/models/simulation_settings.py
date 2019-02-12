from typing import List, Tuple

import numpy as np


class SimulationSettings:
    # Tuple[species_name, label]
    def __init__(self, start_time: float, end_time: float, precision: int,
                 plotted_species: List[str]):
        self.plotted_species = plotted_species
        self.start_time = start_time
        self.end_time = end_time
        self.precision = precision

    def generate_time_space(self):
        return np.linspace(self.start_time, self.end_time, self.precision)
