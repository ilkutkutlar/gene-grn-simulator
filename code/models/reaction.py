from typing import List, Dict, Callable

from models.models import NamedVector
from models.network import Network


class Reaction:
    def __init__(self, left: List[str], right: List[str], rate_fn: Callable[[Network], float]):
        self.left = left
        self.right = right
        self.rate_fn = rate_fn

    def rate_function(self, n: Network) -> float:
        return self.rate_fn(n)

    def change_vector(self, n: Network) -> NamedVector:
        # fpm + MmyR -> [k1] fpm:MmyR
        # ---------------------------
        # MmyR -= k1[MmyR][fpm]     | General: if MmyR in left, then MmyR -= (k) * (left1) * (left2)
        # fpm -= k1[MmyR][fpm]
        # fpm:MmyR += k1[MmyR][fpm] | General: if fpm:MmyR in right, then fpm:MmyR += (k_) * (left1) * (left2)

        # fpm:MmyR -> [k_1] fpm + MmyR
        # ---------------------------
        # MmyR += k_1 [fpm:MmyR]
        # fpm += k_1 [fpm:MmyR]
        # fpm:MmyR -= k_1[fpm:MmyR]

        change: Dict[str, float] = dict()

        for x in n.species:

            if x not in change:
                change[x] = 0

            if x in self.left:
                change[x] -= self.rate_function(n)

            if x in self.right:
                change[x] += self.rate_function(n)

        return change

    def __str__(self) -> str:
        return "Reaction"
