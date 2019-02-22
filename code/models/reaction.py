from typing import List, Dict

from models.formulae import Formula, TranscriptionFormula, TranslationFormula, DegradationFormula, CustomFormula
from models.reg_type import RegType


class Reaction:
    def __init__(self, name: str, left: List[str], right: List[str],
                 rate_fn: Formula):
        self.name = name
        self.left = left
        self.right = right
        self.rate_function = rate_fn

    def rate(self, n: Dict[str, float]) -> float:
        return self.rate_function.compute(n)

    def change_vector(self, n: Dict[str, float]) -> Dict[str, float]:
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

        for x in n:

            if x not in change:
                change[x] = 0

            if x in self.left:
                change[x] -= self.rate(n)

            if x in self.right:
                change[x] += self.rate(n)

        return change

    def __str__(self) -> str:
        left = self.left[0] if self.left else "∅"
        right = self.right[0] if self.right else "∅"

        if isinstance(self.rate_function, TranscriptionFormula):
            # left = "∅"
            # right = self.rate_fn.transcribed_species if self.rate_fn.transcribed_species else "∅"

            trans_rate = str(self.rate_function.rate)
            kd = str(self.rate_function.kd)
            hill_coeff = str(self.rate_function.hill_coeff)

            string = "Type: Transcription" + "\n"
            string += "Reaction: " + left + " -> " + right + "\n\n"

            string += "Rate: " + trans_rate + "\n"
            string += "Kd: " + kd + "\n"
            string += "n: " + hill_coeff + "\n\n"

            string += "== Regulation == \n"

            if self.rate_function.regulators:
                from_gene = self.rate_function.regulators[0].from_gene
                sign = " -> " if self.rate_function.regulators[0].reg_type == RegType.ACTIVATION else " -| "
                to_gene = self.rate_function.regulators[0].to_gene
                string += from_gene + sign + to_gene
            else:
                string += "Not regulated"

            return string

        elif isinstance(self.rate_function, TranslationFormula):
            rate = str(self.rate_function.rate)

            string = "Type: Translation" + "\n"
            string += "Reaction: " + left + " -> " + right + "\n\n"
            string += "Rate: " + rate + "\n"

            return string
        elif isinstance(self.rate_function, DegradationFormula):
            rate = str(self.rate_function.rate)

            string = "Type: Degradation" + "\n"
            string += "Reaction: " + left + " -> " + right + "\n\n"
            string += "Rate: " + rate + "\n"

            return string
        elif isinstance(self.rate_function, CustomFormula):
            rate_function_ast = str(self.rate_function.rate_function)

            params = ""
            for p in self.rate_function.parameters:
                params += "\n       • " + p + ": " + str(self.rate_function.parameters[p])

            string = "Type: Custom Reaction" + "\n"
            string += "Reaction: " + left + " -> " + right + "\n"
            string += "Rate function: " + rate_function_ast + "\n\n"
            string += "== Parameters == \n"
            string += params

            return string
