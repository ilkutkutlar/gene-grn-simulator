# We want to extract these values from a given SBML (for each cassette):
# transcription rate (ps)       -> SBO:0000183
# translation rate              -> SBO:0000184
# mRNA degradation rate         -> SBO:0000179
# protein degradation rate      -> SBO:0000179
# basal promoter strength       -> SBO:0000485
# initial concentrations

# SBOs:
# Protein   -> SBO:0000252
# mRNA      -> SBO:0000250

# In the SBML files, the interesting parts are:
# listOfSpecies     -> Contains proteins, mRNAs and initial concentrations
# listOfParameters  -> Rates
# listOfReactions   -> Some other rates

import xml.etree.ElementTree as ET
import re


def remove_namespaces(text):
    print(text.replace('xmlns=".*"', ""))


def parse_sbml(file_name):
    data: ET.ElementTree = ET.parse(file_name)

    print(ET.QName(data.getroot().tag).text)

    # for top_level in data.getroot():
    # print(top_level.tag)

    # list_of_species =


parse_sbml('BIOMD0000000012.xml')
remove_namespaces('<sbml xmlns="http://www.sbml.org/sbml/level2/version3" level="2" metaid="_153818" version="3">')
