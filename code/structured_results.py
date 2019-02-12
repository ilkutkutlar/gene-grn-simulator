from typing import Dict, List

import numpy as np


class StructuredResults:
    species: Dict[str, float]

    def __init__(self, unstructured_results: np.ndarray,
                 species: List[str], time_space: list):
        self.species = StructuredResults.label_results(unstructured_results, species)
        self.time_space = time_space

    def results_between_times(self, species: str, t1: float, t2: float):
        results: list

        i = 0
        for x in self.time_space:
            if t1 <= x <= t2:
                results.append(self.species[species])
            i += 1

        return results

    @staticmethod
    def label_results(results: np.ndarray, labels: List[str]):
        labeled_results = dict()

        # Attach species names to results
        i = 0
        for s in labels:
            # Syntax meaning: a list consisting of the ith element of each list in results
            labeled_results[s] = results[:, i]
            i += 1

        return labeled_results
