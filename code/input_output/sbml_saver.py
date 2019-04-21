import libsbml
from libsbml._libsbml import parseL3Formula

from models.formulae.custom_formula import CustomFormula
from models.formulae.degradation_formula import DegradationFormula
from models.formulae.transcription_formula import TranscriptionFormula
from models.network import Network
from models.reaction import Reaction
from models.reg_type import RegType
from models.regulation import Regulation


class SbmlSaver:
    @staticmethod
    def network_to_sbml(net):
        sbml_document = libsbml.SBMLDocument(2, 3)  # SBML level 2, version 3
        model = sbml_document.createModel()

        for name, initial_amount in net.species.items():
            spec = model.createSpecies()
            spec.setId(name)
            spec.setInitialAmount(initial_amount)

        # parameters = {p.getId(): p.getValue() for p in x.getKineticLaw().getListOfParameters()}
        # sbo = x.getSBOTerm()

        symbols = []

        for r in net.reactions:
            reaction = model.createReaction()
            reaction.setName(r.name)

            for x in r.left:
                z = reaction.createReactant()
                z.setSpecies(x)

            for x in r.right:
                z = reaction.createProduct()
                z.setSpecies(x)

            # TODO: http://sbml.org/Software/libSBML/5.17.0/docs//python-api/create_simple_model_8py-example.html
            law = reaction.createKineticLaw()
            law.setMath(parseL3Formula(r.rate_function.get_formula_string()))

            if isinstance(r.rate_function, CustomFormula):
                params = r.rate_function.parameters
                for p in params:
                    param = law.createParameter()
                    param.setId(p)
                    param.setValue(params[p])

                # TODO
                # if not symbols:
                #     syms = r.rate_function.symbols
                #     for s in syms:
                #         p = model.createParameter()
                #         p.setId(s)
                #         p.setValue(syms[s])
                #         symbols.append(p)

        return sbml_document

    @staticmethod
    def save_network_to_file(net, filename):
        sbml_document = SbmlSaver.network_to_sbml(net)
        s = libsbml.SBMLWriter()
        s.writeSBMLToFile(sbml_document, filename)


def test():
    species = {"x": 0, "y": 20}

    reg = Regulation("y", "x", RegType.REPRESSION, 40)
    x_trans = TranscriptionFormula(5, "x")
    x_trans.set_regulation(2, [reg])

    reactions = [Reaction("", [], ["x"], x_trans),
                 Reaction("", ["y"], [], DegradationFormula(0.3, "y")),
                 Reaction("", ["y"], [], CustomFormula("10*2", {'r': 4}, {'t': 5}))]

    net: Network = Network()
    net.species = species
    net.reactions = reactions

    m = SbmlSaver.network_to_sbml(net)
    print(m.getListOfReactions()[2].getKineticLaw().getListOfParameters()[0].getValue())
    # print(helper.ast_to_string(m.getListOfReactions()[0].getKineticLaw().getMath()))


if __name__ == '__main__':
    test()
