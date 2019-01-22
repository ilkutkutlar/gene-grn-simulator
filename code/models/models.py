from typing import List, Dict

Vector = List[float]
NamedVector = Dict[str, float]


def __str__(self) -> str:
    left = self.left[0] if self.left else "∅"
    right = self.right[0] if self.right else "∅"
    return "Reaction (" + left + " -> " + right + "): " + self.rate_function_ast
