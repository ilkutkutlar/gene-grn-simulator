from typing import NamedTuple

from models.reg_type import RegType


class Regulation(NamedTuple):
    from_gene: str
    to_gene: str
    reg_type: RegType

    def __str__(self) -> str:
        sign = " -> " if self.reg_type == RegType.ACTIVATION else " -| "
        return "Reg: " + self.from_gene + sign + self.to_gene
