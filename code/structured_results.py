from typing import Dict, List

import numpy as np


class StructuredResults:
    species: Dict[str, np.ndarray]

    def __init__(self, unstructured_results: np.ndarray,
                 names_of_species: List[str], time_space):

        self.species = StructuredResults.label_results(unstructured_results, names_of_species)
        self.time_space = time_space

    """
    :param str species: The name of the species for which results will be returned.
    :param float t1: The lower bound for time period
    :param float t2: The uppor bound for time period
    """
    def results_between_times(self, species: str, t1: float, t2: float):
        results = list()

        i = 0
        for x in self.time_space:
            if t1 <= x <= t2:

                results.append(self.species[species][i])
            i += 1

        return results

    """
    :param np.ndarray results: Raw, unstructured results output by np.odeint.
    :param List[str] labels: The labels to use for labeling the unstructured results
    """
    @staticmethod
    def label_results(results: np.ndarray, labels: List[str]) -> Dict[str, np.ndarray]:
        labeled_results = dict()

        # Attach species names to results
        i = 0
        for s in labels:
            # Syntax meaning: a list consisting of the ith element of each list in results
            labeled_results[s] = results[:, i]
            i += 1

        return labeled_results
