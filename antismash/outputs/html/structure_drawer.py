# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.


import shutil
import os
import logging

from helperlibs.wrappers.io import TemporaryDirectory

import antismash.common.deprecated as utils
import antismash.common.path as path

from .indigo import Indigo
from .indigo_renderer import IndigoRenderer

def _update_sec_met_entry(clusterfeature, smiles_string):
    clusterfeature.qualifiers['structure'] = [smiles_string]


def generate_chemical_structure_preds(pksnrpsvars, seq_record, options):
    #Create directory to store structures
    options.structuresfolder = os.path.abspath(os.path.join(options.output_dir, "structures"))
    if not os.path.exists(options.structuresfolder):
        os.mkdir(options.structuresfolder)

    #Combine predictions into a prediction of the final chemical structure and generate images
    for genecluster in seq_record.get_clusters():
        geneclusternr = genecluster.get_cluster_number()
        smiles_string = ""
        if pksnrpsvars.compound_pred_dict.has_key(geneclusternr):

            #print "output_modules/html/pksnrpsvars.compound_pred_dict:"
            #print pksnrpsvars.compound_pred_dict

            residues = pksnrpsvars.compound_pred_dict[geneclusternr].replace("(", "").replace(")", "").replace(" + ", " ").replace("-", " ")

            #Now generates SMILES of predicted secondary metabolites without NP.searcher
            residuesList = residues.split(" ")

            #Counts the number of malonate and its derivatives in polyketides
            mal_count = 0
            for i in residuesList:
                if "mal" in i:
                    mal_count += 1

            nrresidues = len(residuesList)

            #Reflecting reduction states of ketide groups starting at beta carbon of type 1 polyketide
            if "pk" in residuesList and "mal" in residuesList[-1]:
                residuesList.pop(residuesList.index('pk')+1)
                residuesList.append('pks-end1')
            elif mal_count == len(residuesList):
                if residuesList[0] == "mal":
                    residuesList[0] = "pks-start1"
                if residuesList[-1] == "ccmal":
                    residuesList.append('pks-end2')

            if nrresidues > 1:
                #Conventionally used aaSMILES was used;
                #chirality expressed with "@@" causes indigo error
                aa_smiles = load_smiles()

                for monomer in residuesList:
                    if monomer in aa_smiles:
                        smiles_string += aa_smiles[monomer]
                    elif '|' in monomer:
                        logging.debug("Substituting 'nrp' for combined monomer %r", monomer)
                        smiles_string += aa_smiles['nrp']
                    else:
                        logging.debug("No SMILES mapping for unknown monomer %r", monomer)
                logging.debug("Cluster %s: smiles_string: %s", geneclusternr, smiles_string)
                with TemporaryDirectory(change=True):
                    smilesfile = open("genecluster" + str(geneclusternr) + ".smi", "w")
                    smilesfile.write(smiles_string)
                    smilesfile.close()
                    depictstatus = depict_smile(geneclusternr, options.structuresfolder)
                if depictstatus == "failed":
                    pksnrpsvars.failedstructures.append(geneclusternr)
        elif utils.get_cluster_type(genecluster) == "ectoine":
            smiles_string = "CC1=NCCC(N1)C(=O)O"
            with TemporaryDirectory(change=True):
                smilesfile = open("genecluster" + str(geneclusternr) + ".smi", "w")
                smilesfile.write(smiles_string)
                smilesfile.close()
                depictstatus = depict_smile(geneclusternr, options.structuresfolder)
            if depictstatus == "failed":
                pksnrpsvars.failedstructures.append(geneclusternr)
            elif genecluster in pksnrpsvars.failedstructures:
                del pksnrpsvars.failedstructures[pksnrpsvars.failedstructures.index(geneclusternr)]
            pksnrpsvars.compound_pred_dict[geneclusternr] = "ectoine"
        _update_sec_met_entry(genecluster, smiles_string)


def load_smiles():
    """Load smiles from a dictionary mapping residues to SMILES string"""
    aa_smiles = {}

    smiles_monomer = open(path.get_full_path(__file__, 'aaSMILES.txt'), 'r')

    for line in smiles_monomer.readlines():
        line = line.strip()
        if not line or line.startswith('#') or line == "END":
            continue
        smiles = line.split()
        assert len(smiles) == 2, "Invalid smiles line {!r}".format(line)

        aa_smiles[smiles[0]] = smiles[1]

    smiles_monomer.close()
    return aa_smiles


def depict_smile(genecluster, structuresfolder):
    indigo = Indigo()
    renderer = IndigoRenderer(indigo)
    query = indigo.loadMoleculeFromFile("genecluster" + str(genecluster) + ".smi")
    indigo.setOption("render-coloring", True)
    renderer.renderToFile(query, "genecluster" + str(genecluster) + ".png")

    indigo.setOption("render-image-size", 200, 150)
    renderer.renderToFile(query, "genecluster" + str(genecluster) + "_icon.png")
    dircontents = os.listdir(os.getcwd())
    geneclusterstring = "genecluster" + str(genecluster) + ".png"
    if geneclusterstring in dircontents:
        shutil.copy("genecluster" + str(genecluster) + ".png", structuresfolder)
        shutil.copy("genecluster" + str(genecluster) + "_icon.png", structuresfolder)
        shutil.copy("genecluster" + str(genecluster) + ".smi", structuresfolder)
        os.remove("genecluster" + str(genecluster) + ".png")
        os.remove("genecluster" + str(genecluster) + "_icon.png")
        os.remove("genecluster" + str(genecluster) + ".smi")
        smiles_input = os.path.join('SMILES', 'input')
        if os.path.exists(smiles_input):
            os.remove(smiles_input)
        return "success"
    else:
        return "failed"
