import libsbml

# 1. Core objects of libsbml:
# http://sbml.org/Software/libSBML/5.17.0/docs//python-api/group__core.html
# 2. Classes:
# http://sbml.org/Software/libSBML/5.17.0/docs//python-api/annotated.html
from models import Network


class SbmlParser:
    def __init__(self, filename):
        self.filename = filename

    def parse(self):
        net: Network = Network()

        sbml = libsbml.SBMLReader()
        parsed = sbml.readSBML(self.filename)
        model: libsbml.Model = parsed.getModel()

        for species in model.getListOfSpecies():
            net.species[species.getId()] = species.getInitialAmount()

        # for x in model.getListOfReactions():
        #     print(x.getKineticLaw().getMath())
        print(net.species)


def parsing():
    p = SbmlParser("other_files/BIOMD0000000012.xml")
    p.parse()


parsing()
