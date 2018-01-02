# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging

import Bio

from antismash.common import fasta, path, secmet, subprocessing

from .common import get_scaffold_annotation, get_prediction_annotation, get_signature

def run_all(record):
    for func in [pksi_kr_stereo, asp_thioesterase]:
        func(record)


def pksi_kr_stereo(record):
    """
        Database source: CLUSEAN
        Reference: Reid, R., M. Piagentini, E. Rodriguez, G. Ashley, N.,
                   Viswanathan, J. Carney, D. V. Santi, C. R. Hutchinson,
                   and R. McDaniel. 2003.

                   A model of structure and catalysis for ketoreductase domains
                   in modular polyketide synthases. Biochemistry 42:72-79.
    """
    logging.debug("Pediction of PKSI KR specificities according to Reid et al., Biochemistry 2003, 42, 72-79")

    method_type = "prediction"

    targets = [target for target in record.get_antismash_domains() if target.domain == "PKS_KR"]
    if not targets:
        return
    database = path.get_full_path(__file__, "data", "PKSI-KR.hmm2")
    extra_args = ["-T", "0", # min score
               "-E", "0.1"] # max evalue
    data = fasta.get_fasta_from_features(targets, numeric_names=True)
    results = subprocessing.run_hmmpfam2(database, data, extra_args=extra_args)
    logging.debug("found %s hsps in hmmer results", len(results))  # TODO: temp debug line from old version
    if not results:
        return

#    alignment = ...
    positions = [149, 162, 166]
    expected_values = ["S", "Y", "N"]
    comment = "KR domain putatively catalyzing {}-configuration product formation"
    if get_signature(results[0].aln[0].seq, results[0].aln[1].seq, [102]) == "D":
        comment = comment.format("D")
    else:
        comment = comment.format("L")
    raise RuntimeError("ran pksi_kr_stereo")


def asp_thioesterase(record):
    """ database: Thioesterase.hmm2
        database_source: CLUSEAN
        reference: unknown
    """
    method_type = "active_site"
    targets = [target for target in record.get_antismash_domains() if target.domain == "Thioesterase"]
    database = path.get_full_path(__file__, "data", "Thioesterase.hmm2")
    extra_args = ["-T", "0", # min score
               "-E", "0.1"] # max evalue
    data = fasta.get_fasta_from_features(targets, numeric_names=True)
    results = subprocessing.run_hmmpfam2(database, data, extra_args=extra_args)
    logging.debug("found %s hsps in hmmer results", len(results))  # TODO: temp debug line from old version
    if not results:
        return

    positions = [73, 79, 83, 87, 108]
    expected_values = ["G", "G", "G", "A", "D"]
    emissions = [0.93, 1., 0.99, 0.9, 0.9]

    active_position = 81
    expected_active_value = "S"
    active_emission = 0.93

    hit_string = "active site serine present"

    for result in results:
        antismash_domain = targets[int(result.id)]
        logging.debug("found hit with %s", antismash_domain.locus_tag)
        if not result.hsps:
            continue


        assert result.hsps[0].aln[0].id == result.id

        # identify scaffolds and annotate
        ASF_string = get_scaffold_annotation(result, positions, expected_values, emissions)

        note = None
        scaffold = []
        if ASF_string:
            scaffold = [ASF_string]
            logging.debug(ASF_string)
            note = ["ASF analysis with definition ASP_thioesterase (type active_site)"]

        # identify predictions / active sites and annotate
        description, choice_result = get_prediction_annotation(result, [active_position],
                                             [expected_active_value],
                                             result_label="active site serine present",
                                             comment="active site serine",
                                             expected_emissions=[active_emission])
        choices = [description]

        if description:
            logging.debug("2: adding ASF choice info to %s %s..%s:", antismash_domain.type, antismash_domain.location.start, antismash_domain.location.end)
            logging.debug(description)
            note = ["ASF analyisis with definition ASP_thioesterase (type active_site)"]

        prediction = [choice_result] if choice_result else []
        for choiceResult in prediction:
            logging.debug("adding ASF choiceResult info to %s %s..%s:", antismash_domain.type, antismash_domain.location.start, antismash_domain.location.end)
            logging.debug(choiceResult)

            # Also annotate choiceResult as sec_met qualifier to corresponding CDS feature

            correspondingCDSFeature = record.get_cds_name_mapping()[antismash_domain.locus_tag]

            logging.debug("2: adding ASF-prediction data to sec_met qualifier of %s", correspondingCDSFeature.get_name())
            logging.critical("2: skipping addition of ASF prediction to feature.sec_met")
            if False:
                sec_met_string = "ASF-prediction: "

                # Calculate relative locations of domains

                if antismash_domain.strand == 1:
                    start = ((antismash_domain.location.start - correspondingCDSFeature.location.start + 3) / 3) - 1
                    end = ((antismash_domain.location.end - correspondingCDSFeature.location.start + 3) / 3) - 1
                else:
                    start = ((correspondingCDSFeature.location.end - antismash_domain.location.end + 3) / 3) - 1
                    end = ((correspondingCDSFeature.location.end - antismash_domain.location.start + 3) / 3) - 1
                sec_met_string += antismash_domain.qualifiers['domain'][0] + " (" + str(start) + ".." + str(end) + "): "
                sec_met_string += choiceResult
                correspondingCDSFeature.qualifiers['sec_met'].append(sec_met_string)
        if choices or scaffold or prediction:
            cds = record.get_cds_name_mapping()[antismash_domain.locus_tag]
            cds.asf = secmet.feature.ActiveSiteFinderQualifier(choice=choices,
                         scaffold=scaffold, note=note, prediction=prediction)  # TODO: stop this clobbering
            print("cds:", cds.get_name(), choices, scaffold, note, prediction)
