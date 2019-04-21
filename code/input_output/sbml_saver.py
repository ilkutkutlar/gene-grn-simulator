import libsbml
from libsbml._libsbml import parseL3Formula

from models.formulae.custom_formula import CustomFormula
from models.formulae.degradation_formula import DegradationFormula
from models.formulae.transcription_formula import TranscriptionFormula
from models.formulae.translation_formula import TranslationFormula
from models.network import Network
from models.reaction import Reaction
from models.reg_type import RegType
from models.regulation import Regulation


class SbmlSaver:
    @staticmethod
    def network_to_sbml(net):
        sbml_document = libsbml.SBMLDocument(2, 3)  # SBML level 2, version 3
        model = sbml_document.createModel()
        compartment = model.createCompartment()
        compartment.setId("cell")
        compartment.setSize(1)

        substance_def = model.createUnitDefinition()
        substance_def.setId("substance")
        substance = substance_def.createUnit()
        substance.setKind(libsbml.UNIT_KIND_ITEM)
        substance.setMultiplier(1)

        for name, initial_amount in net.species.items():
            spec = model.createSpecies()
            spec.setId(name)
            spec.setInitialAmount(initial_amount)
            spec.setCompartment("cell")
            spec.setHasOnlySubstanceUnits(True)

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

                # TODO -> Problem is, won't properly save symbols of an
                #  imported SBML file back to the exported SBML file.
                # if not symbols:
                #     syms = r.rate_function.symbols
                #     for s in syms:
                #         p = model.createParameter()
                #         p.setId(s)
                #         p.setValue(syms[s])
                #         symbols.append(p)
            elif isinstance(r.rate_function, TranscriptionFormula):
                if r.rate_function.regulators:
                    for reg in r.rate_function.regulators:
                        m = model.createModifier()
                        m.setSpecies(reg.from_gene)
            elif isinstance(r.rate_function, TranslationFormula):
                m = model.createModifier()
                m.setSpecies(r.rate_function.mrna_species)

        return sbml_document

    @staticmethod
    def save_network_to_file(net, filename):
        sbml_document = SbmlSaver.network_to_sbml(net)
        s = libsbml.SBMLWriter()
        s.writeSBMLToFile(sbml_document, filename)
