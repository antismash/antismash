# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""
Identify conserved active site residues in PFAM_Doman / aSDomain features
"""

from io import StringIO
import logging
import os
import re
from tempfile import NamedTemporaryFile

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

# Ignore Biopython experimental warning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO

from antismash.common import path, subprocessing, module_results, secmet
from antismash.config import get_config
from antismash.config.args import ModuleArgs

from .analysis import run_all

NAME = "ActiveSiteFinder"
SHORT_DESCRIPTION = "ActiveSiteFinder identifies conserved active sites in PFAM_Domain/aSDomain features"


class ASFResults(module_results.ModuleResults):
    def to_json(self, *args):
        logging.critical("skipping ASFResults.to_json()")
        return {}

    def from_json(self, *args):
        logging.critical("skipping ASFResults.from_json()")
        return None

    def add_to_record(self, *args):
        logging.critical("skipping ASFResults.add_to_record()")


def check_prereqs():
    """Check the prerequisites"""
    # create a throwaway active_site_finder object without a record
    my_asf = active_site_finder(None, get_config())
    return my_asf.check_prereqs()


def check_options(options):
    return []  # TODO: maybe bail if no full_hmmer?


def get_arguments():
    args = ModuleArgs('Additional analysis', 'asf')
    args.add_analysis_toggle('--asf',
                             dest='asf',
                             action='store_true',
                             default=False,
                             help="Run active site finder analysis.")
    return args


def is_enabled(options):
    return options.asf


def run_on_record(record, results, options):
    run_all(record)
    asf = active_site_finder(record, options)
    asf.execute()
    return ASFResults(record.id)


def get_prediction_annotation(result, predictionChoicesXML):
    "gererate annotation from choices/prediciton information"

    choices = []
    descriptions = []
    query_seq = result.hsps[0].aln[0].seq
    hmm_seq = result.hsps[0].aln[1].seq

    for choice in predictionChoicesXML:
        extracted_aa_List = []
        emissions = []
        matches = []
        skip = False
        predictionOffsetList = choice.find('./offset').text.split(',')
        predictionValueList = choice.find('./value').text.split(',')
        predictionResult = choice.attrib['result']
        predicitionComment = choice.find('./comment').text
        if choice.find('./valueEmission'):
            predictionValueEmissionList = choice.find('./valueEmission').text.split(',')
        else:
            predictionValueEmissionList = []

        choiceOverallMatch = True

        logging.debug("testing %s (%s):", predictionResult, predicitionComment)

        for i in range(0, len(predictionOffsetList)):
            predictionOffset = int(predictionOffsetList[i]) - 1
            predictionValue = predictionValueList[i]
            if predictionValueEmissionList:
                predictionValueEmission = int(predictionValueEmissionList[i])
            else:
                predictionValueEmission = "n.d."
            emissions.append(predictionValueEmission)

            # Check whether coordinates are within HSP match of the HMM profile
            if (result.hsps[0].hit_start > predictionOffset) or (result.hsps[0].hit_end < predictionOffset):
                logging.warning("choice/prediction coordinate %s outside hsp!", predictionOffset)
                choiceOverallMatch = False
                skip = True
                break

            try:
                fixed_predictionOffset = fix_coordinates(predictionOffset - result.hsps[0].hit_start, hmm_seq)
            except ValueError:
                logging.warning("gap-fixed choice/prediction coordinate %s outside hsp! for result: %s, which has length: %s",
                                predictionOffset, result.id, len(query_seq))
                choiceOverallMatch = False
                skip = True
                break

            # Check whether gap-fixed coordinate still is within returned query sequence
            if len(hmm_seq) < fixed_predictionOffset + 1:
                logging.warning("gap-fixed choice/prediction coordinate %s outside hsp! for result: %s, which has length: %s",
                                fixed_predictionOffset, result.id, len(query_seq))
                choiceOverallMatch = False
                skip = True
                break

            try:
                extracted_aa = query_seq[fixed_predictionOffset]
            except IndexError:
                # This error now actually should never be raised
                logging.error("gap-fixed choice/prediction coordinate %s outside hsp! for result: %s, which has length: %s; this should never happen",
                              fixed_predictionOffset, result.id, len(query_seq))
                choiceOverallMatch = False
                skip = True
                break

            extracted_aa_List.append(extracted_aa)

            choiceMatch = bool(re.match("(?i)" + predictionValue, extracted_aa))
            choiceOverallMatch = choiceOverallMatch and choiceMatch

            matches.append(str(choiceMatch))
            logging.debug("Offset %s; fixed offset %s  Expected: %s; observed in query %s: observed in hmm %s; Emission %s; match %s",
                          predictionOffset,
                          fixed_predictionOffset,
                          predictionValue,
                          extracted_aa,
                          hmm_seq[fixed_predictionOffset],
                          predictionValueEmission,
                          choiceMatch)

        logging.debug("Overall Match for prediction %s: %s", predictionResult, str(choiceOverallMatch).upper())
        logging.debug("================================")

        if not skip:
            ASF_string = "Description: %s, choice result: %s, choice coordinates: (%s); residues: (%s); " \
                         "expected for choice: (%s); matchArray: (%s); emission probability array (%s); overall match: %s" % \
                          (predicitionComment,
                           predictionResult,
                           ",".join(predictionOffsetList),
                           ",".join(extracted_aa_List),
                           ",".join(predictionValueList),
                           ",".join(matches),
                           ",".join(emissions),
                           str(choiceOverallMatch).upper())

            descriptions.append(ASF_string)

        if choiceOverallMatch:
            choice_string = "Full match for prediction: %s" % predictionResult
            choices.append(choice_string)
    return descriptions, choices


def fix_coordinates(coordinate, seq):
    "Fix coordinates in hmm reference seq if gaps are present"

    # Check if coordinate is within string range
    if len(seq) - seq.count('.') - 1 < coordinate:
        logging.error('tried to fix coordinate failed; coordinate exceeds sequence length. This should never happen!')
        raise ValueError('Coordinate not within sequence')

    numberOfGaps = seq[:coordinate].count('.')

    while coordinate < len(seq) and seq[coordinate] == ".":
        coordinate += 1
        logging.debug("increase coordinate by 1")
    temp_coordinate = coordinate
    new_coordinate = coordinate + numberOfGaps

    # if seq[coordinate:new_coordinate] also contains gaps, we also have to add these;
    while seq[temp_coordinate:new_coordinate].count('.') > 0:
        t = new_coordinate
        new_coordinate += seq[temp_coordinate:new_coordinate].count('.')
        temp_coordinate = t

    return new_coordinate


def filter_features_by_qualifier(features, query_tag, query_value):
    """ Return all features which contain a 'query_tag' with value 'query_value'
        Note: query has to be exact!
    """
    logging.critical("attempting to filter features by: %s in qualifier %s", query_value, query_tag)
    assert query_tag == "domain"
    filtered_features = []
    for feature_to_test in features:
        if query_value == feature_to_test.domain:
            filtered_features.append(feature_to_test)
    return filtered_features


def get_scaffold_annotation(result, scaffold_xml):
    "generate annotation from scaffold information"

    query_seq = result.hsps[0].aln[0].seq
    hmm_seq = result.hsps[0].aln[1].seq
    scaffoldPos = scaffold_xml.find('./scaffoldOffset').text
    scaffoldValue = scaffold_xml.find('./scaffoldValue').text
    if scaffold_xml.find('./scaffoldEmission'):
        scaffoldEmissionList = scaffold_xml.find('./scaffoldEmission').text.split(',')
    else:
        scaffoldEmissionList = []
    scaffoldPosList = scaffoldPos.split(',')
    scaffoldValueList = scaffoldValue.split(',')

    overallMatch = True

    # Calculate and print match overview line and print alignment and match line
    matchLineStr = []
    for i in range(0, len(hmm_seq)):
        offset = i + result.hsps[0].hit_start + 1
        if hmm_seq[i] == ".":
            matchLineStr.append(" ")
        if str(offset) in scaffoldPosList:
            matchLineStr.append("*")
        else:
            matchLineStr.append(" ")
    logging.debug("%s %s..%s", query_seq, result.hsps[0].query_start, result.hsps[0].query_end)
    logging.debug("%s %s..%s", hmm_seq, result.hsps[0].hit_start, result.hsps[0].hit_end)
    logging.debug("".join(matchLineStr))

    # Check scaffold matches
    extracted_aa_List = []
    matches = []
    emissions = []

    for i in range(0, len(scaffoldPosList)):
        scafPos = int(scaffoldPosList[i]) - 1
        scafValue = scaffoldValueList[i]

        # Check whether scafPos is within HSP coordinates
        if (result.hsps[0].hit_start > scafPos) or (result.hsps[0].hit_end < scafPos):
            logging.warning("scaffold coordinate %s outside hsp!", scafPos)
            logging.debug("Overall Scaffold Match: FALSE\n")
            return None
        # fix position for gap characters in hmm_hit
        try:
            fixedScafPos = fix_coordinates(scafPos - result.hsps[0].hit_start, hmm_seq)
        except ValueError:
            logging.error("gap-fixed scaffold coordinate %s outside hsp for original position: ", scafPos)
            logging.debug("Overall Scaffold Match: FALSE\n")
            return None

        # Check whether amino acid in hmm_seq[fixedScafPos] equals predifined aa from XML
        # (should always match) and thus may be removed when the module is thoroughly tested
        # This statement fails if the fixedScafPos is outside alignment due to gaps..., so I catch the exception...
        try:
            if hmm_seq[fixedScafPos].lower() != scafValue.lower():
                logging.warning("ASF: aa extracted from hmm profile does not match predifined aa in XML config file!")
        except IndexError:
            logging.warning("gap-fixed scaffold coordinate %s outside hsp!", fixedScafPos)
            logging.debug("Overall Scaffold Match: FALSE\n")
            return None

        extracted_aa = query_seq[fixedScafPos]
        extracted_aa_List.append(extracted_aa)
        if scaffoldEmissionList:
            extracted_aa_Emission = int(scaffoldEmissionList[i])
        else:
            extracted_aa_Emission = "n.d."
        emissions.append(extracted_aa_Emission)

        # We have to use a RegEx here to allow negations and more complex queries; ignore case (?i)
        match = bool(re.match("(?i)" + scafValue, query_seq[fixedScafPos]))
        overallMatch = overallMatch and match
        matches.append(str(match))
        logging.debug("Scaffold coordinate %s; fixed scaffold coordinate %s, query aa %s; hmm aa %s; "
                      "scaffold value %s; emission probability %s; match %s",
                      scafPos, fixedScafPos, extracted_aa, hmm_seq[fixedScafPos], scafValue, extracted_aa_Emission, match)

    logging.debug("Overall Scaffold Match: %s\n", str(overallMatch).upper())

    # Generate Feature qualifiers
    ASF_string = ("Scaffold coordinates: (%s); scaffold residues: (%s); expected: (%s); matchArray: (%s); "
                  "emission probability array (%s); overall match: %s") % (
                 ",".join(scaffoldPosList), ",".join(extracted_aa_List), ",".join(scaffoldValueList),
                 ",".join(matches), ",".join(emissions), str(overallMatch).upper())

    return ASF_string


def execute_tool(analysisResource, fileName=None, stdin_data=None):
    "Perform the external program execution"

    cmdlineList = []

    # Assemble commad line list
    # extract program name from XML
    executeObj = analysisResource.find('./Execute')
    cmdlineList.append(executeObj.attrib['program'])

    # Cycle through parameters in XML
    for parameter in list(analysisResource.findall('./Execute/parameters/parameter')):

        if 'prefix' in parameter.attrib:
            cmdlineList.append(parameter.attrib['prefix'])
        cmdlineList.append(parameter.text)

    # Get database name
    database = analysisResource.find('./Execute/database')
    if 'prefix' in database.attrib:
        cmdlineList.append(database.attrib['prefix'])
    # Add searchpath
    cmdlineList.append(path.locate_file(path.get_full_path(__file__, "data", database.text)))

    if fileName:
        # Get (optional) input file prefix (e.g. -query in blast)
        if 'inputfile_prefix' in executeObj.attrib:
            cmdlineList.append(executeObj.attrib['inputfile_prefix'])
        cmdlineList.append(fileName)

    if stdin_data:
        # Get (optional) prefix for stdin (e.g. "-" for hmmpfam / hmmscan
        if 'STDINprefix' in executeObj.attrib:
            cmdlineList.append(executeObj.attrib['STDINprefix'])

    logging.debug("ASF: %s; external program call:\n%s", analysisResource.attrib['name'], " ".join(cmdlineList))

    try:
        if fileName:
            logging.debug("Executing tool with file input")
            result = subprocessing.execute(cmdlineList)
        else:
            logging.debug("Executing tools with STDIN input: %s", stdin_data)
            result = subprocessing.execute(cmdlineList, stdin=stdin_data)
    except OSError:
        logging.warning('OS error on execution of: %s', " ".join(cmdlineList))
        return []
    if not result.successful():
        logging.warning('%s returned %s', cmdlineList[0], result.exit_code)
        return []
    res_stream = StringIO(result.stdout)

    # Get Biopython parser information from XML
    biopython_parser = analysisResource.find('./Execute/BioPythonParser')
    try:
        results = list(SearchIO.parse(res_stream, biopython_parser.text))
    except Exception as e:
        logging.warning('Error parsing results for active site finder analysis: %s ; no hits will be reported', e)
        results = []

    return results


def run_external_tool(analysisResource, features):
    "Generate tempfile containing the extracted Feature sequences and run tool defined in XML file"

    # write fasta file for the features to tempfile
    # FastA headers will be a simple numerical index to avoid names that
    # are too long for tools

    tempfile = ""
    fastafile = []
    for n, feature in enumerate(features):
        fastaHeader = "feature%d" % n
        fastaSeq = feature.translation
        # Never write empty fasta entries
        if not fastaSeq:
            logging.warning("No translation for %s, skipping", fastaHeader)
            continue
        fastafile.append(">%s\n" % fastaHeader)
        fastafile.append("%s\n" % fastaSeq)
    querydata = "".join(fastafile)

    UseSTDIN = "False"
    executeObj = analysisResource.find('./Execute')
    if 'UseSTDIN' in executeObj.attrib:
        UseSTDIN = executeObj.attrib['UseSTDIN']

    results = []
    if not fastafile:
        logging.warning("ASP: No features found containing feature/tag/value %s / %s / %s",
                        analysisResource.find('./Prerequisite/primary_tag_type').text,
                        analysisResource.find('./Prerequisite/tag').text,
                        analysisResource.find('./Prerequisite/tag_value').text)
    elif UseSTDIN == "True":
        results = execute_tool(analysisResource, stdin_data=querydata)
    else:
        with NamedTemporaryFile(prefix='antiSMASH_ASP') as tempfile:
            out_file = open(tempfile.name, "w")
            out_file.write(querydata)
            out_file.close()
            results = execute_tool(analysisResource, fileName=tempfile.name)

    return results


class active_site_finder(object):
    """Active site finder class; perfoms analysis of active sites
    identified from reference aa positions of HMM profiles"""

    def __init__(self, record, options):
        "Initialize ASF object"

        # Set options
        conf = path.get_full_path(__file__, "config", "SignatureResources.xml")

        # Assign variables
        try:
            XMLtree = ET.parse(conf)
        except ET.ParseError:
            logging.exception("Could not load/parse ActiveSiteFinder configuration file %s.", conf)
            raise

        XMLroot = XMLtree.getroot()

        HmmProfilesFilenameObj = XMLroot.findall(".//Execute/database")

        self.record = record
        self.options = options
        self.XMLtree = XMLtree
        self.XMLroot = XMLroot
        self.HmmProfilesFilenameObj = HmmProfilesFilenameObj

    def execute(self):
        """Execute complete analysis; includes XML parsing, export of feature aa sequences,
        external tool calling, evaluation of the results and addition of new annotation to the SeqFeature"""

        XMLroot = self.XMLroot
        record = self.record

        # get dictionary of all features matching primary_tag_type, tag and tag value obtained from XML file

        analysisResources = XMLroot.findall('./analysis')

        for analysisResource in analysisResources:
            # set variables
            analysisResourceName = analysisResource.attrib['name']
            analysisResourceType = analysisResource.attrib['type']
            logging.debug("**********************************************\nresource name %s\n*********************************", analysisResourceName)

            # Get information on features to search for from XML file
            primaryTagType = analysisResource.find('./Prerequisite/primary_tag_type').text
            tag = analysisResource.find('./Prerequisite/tag').text
            tagValue = analysisResource.find('./Prerequisite/tag_value').text

            # get scaffold / choice subtrees form XML file
            scaffold_xml = analysisResource.find('./Alignment/scaffold')
            predictionChoicesXML = analysisResource.findall('./Alignment/choice')

            # Get pre-annoated feature where we have rules defined in XML file
            if primaryTagType == "PFAM_domain":
                SeqFeatureList = record.get_pfam_domains()
                SeqFeatureList = filter_features_by_qualifier(SeqFeatureList, tag, tagValue)
            elif primaryTagType == "aSDomain":
                SeqFeatureList = [domain for domain in record.get_antismash_domains() if tagValue == domain.domain]
            else:
                raise ValueError("Unknown tag type: %s", primaryTagType)

            SeqFeature_byID = {}
            for n, feature in enumerate(SeqFeatureList):
                SeqFeature_byID["feature%d" % n] = feature

            # write multi-fasta tempfile with domain aa sequences and execute tool as defined in XML file
            # and parse with Biopython SearchIO
            results = run_external_tool(analysisResource, SeqFeatureList)
            logging.debug("found %s hsps in hmmer results", len(results))

            # Now cycle through all results for this analysis and associate the results
            # with the corresponding SeqFeature
            for result in results:
                logging.debug("found hit with %s", result.id)
                if not result.hsps:
                    continue

                SeqFeature = SeqFeature_byID[result.id]

                if result.hsps[0].aln[0].id != result.id:
                    # This actually should never happen...
                    logging.exception("Result ID: %s", result.id)
                    logging.exception("Was looking for hit %s but got hit for %s instead",
                                      result.hsps[0].aln[0].id, result.id)
                    break

                # identify scaffolds and annotate
                ASF_string = get_scaffold_annotation(result, scaffold_xml)

                note = None
                scaffold = []
                if ASF_string:
                    scaffold = [ASF_string]
                    logging.debug(ASF_string)
                    note = ["ASF analyisis with definition %s (type %s)" %
                            (analysisResourceName, analysisResourceType)]

                # identify predictions / active sites and annotate
                descriptions, choiceResultList = get_prediction_annotation(result, predictionChoicesXML)
                choices = list(descriptions)

                for description in descriptions:
                    logging.debug("adding ASF choice info to %s %s..%s:", SeqFeature.type, SeqFeature.location.start, SeqFeature.location.end)
                    logging.debug(description)
                    note = ["ASF analyisis with definition %s (type %s)" %
                            (analysisResourceName, analysisResourceType)]

                prediction = list(choiceResultList)
                for choiceResult in choiceResultList:
                    logging.debug("adding ASF choiceResult info to %s %s..%s:", SeqFeature.type, SeqFeature.location.start, SeqFeature.location.end)
                    logging.debug(choiceResult)

                    # Also annotate choiceResult as sec_met qualifier to corresponding CDS feature

                    correspondingCDSFeature = record.get_cds_name_mapping()[SeqFeature.locus_tag]

                    logging.debug("adding ASF-prediction data to sec_met qualifier of %s", correspondingCDSFeature.get_name())
                    logging.critical("skipping addition of ASF prediction to feature.sec_met")
                    if False:
                        sec_met_string = "ASF-prediction: "

                        # Calculate relative locations of domains

                        if SeqFeature.strand == 1:
                            start = ((SeqFeature.location.start - correspondingCDSFeature.location.start + 3) / 3) - 1
                            end = ((SeqFeature.location.end - correspondingCDSFeature.location.start + 3) / 3) - 1
                        else:
                            start = ((correspondingCDSFeature.location.end - SeqFeature.location.end + 3) / 3) - 1
                            end = ((correspondingCDSFeature.location.end - SeqFeature.location.start + 3) / 3) - 1
                        sec_met_string += SeqFeature.qualifiers['domain'][0] + " (" + str(start) + ".." + str(end) + "): "
                        sec_met_string += choiceResult
                        correspondingCDSFeature.qualifiers['sec_met'].append(sec_met_string)
                if choices or scaffold or prediction:
                    cds = record.get_cds_name_mapping()[SeqFeature.locus_tag]
                    cds.asf = secmet.feature.ActiveSiteFinderQualifier(choice=choices,
                                 scaffold=scaffold, note=note, prediction=prediction)
                    if analysisResourceName == "ASP_Thioesterase":
                        print("cds:", cds.get_name(), choices, scaffold, note, prediction)
        return True

    def check_prereqs(self):
        "Check if all required files and applications are around"

        # Tuple is ( binary_name, optional)
        _required_binaries = [
            ('blastp', False),
            ('hmmpfam2', False),
            ('hmmscan', False),
            ('hmmpress', False),
        ]

        _binary_extensions = [
            '.h3f',
            '.h3i',
            '.h3m',
            '.h3p'
        ]

        failure_messages = []

        for binary_name, optional in _required_binaries:
            if path.locate_executable(binary_name) is None and not optional:
                failure_messages.append("Failed to locate file: %r" % binary_name)

        # Get all HMM profile names from XML file
        paths_checked = set()
        for profile in self.HmmProfilesFilenameObj:
            full_hmm_path = path.get_full_path(__file__, "data", profile.text)
            if full_hmm_path in paths_checked:
                continue
            paths_checked.add(full_hmm_path)

            if path.locate_file(full_hmm_path) is None:
                failure_messages.append("Failed to locate file: %s" % profile.text)
                continue

            if profile.text.endswith(".hmm2"):
                continue

            for ext in _binary_extensions:
                binary = "{hmm}{ext}".format(hmm=full_hmm_path, ext=ext)
                if not path.locate_file(binary):
                    result = subprocessing.run_hmmpress(full_hmm_path)
                    if not result.successful():
                        failure_messages.append("Failed to hmmpress {!r}: {!r}".format(profile.text, result.stderr))

                    # hmmpress generates _all_ binary files in one go, so stop the loop
                    break

                binary_mtime = os.path.getmtime(binary)
                hmm_mtime = os.path.getmtime(full_hmm_path)
                if hmm_mtime < binary_mtime:
                    # generated file younger than hmm profile, do nothing
                    continue
                try:
                    from glob import glob
                    for filename in glob("{}.h3?".format(full_hmm_path)):
                        logging.debug("removing outdated file %r", filename)
                        os.remove(filename)
                except OSError as err:
                    failure_messages.append("Failed to remove outdated binary file for %s: %s" %
                                            (profile.text, err))
                    break
                result = subprocessing.run_hmmpress(full_hmm_path)
                if not result.successful():
                    failure_messages.append("Failed to hmmpress %r: %r" % (profile.text, result.stderr))
                    import datetime
                    failure_messages.append("HMM binary files outdated. %s (changed: %s) vs %s (changed: %s)" %
                                            (profile.text, datetime.datetime.fromtimestamp(hmm_mtime),
                                             binary, datetime.datetime.fromtimestamp(binary_mtime)))
                # hmmpress generates _all_ binary files in one go, so stop the loop
                break

        return failure_messages
