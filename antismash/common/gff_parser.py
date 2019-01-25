# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Reads GFF files and updates records with the contained information.
"""


import logging
from typing import Dict, IO, List, Set, Union

from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from BCBio import GFF

from antismash.common.errors import AntismashInputError
from antismash.common.secmet import Record
from antismash.common.secmet.errors import SecmetInvalidInputError
from antismash.config import ConfigType


def check_gff_suitability(options: ConfigType, sequences: List[Record]) -> bool:
    """
        Checks that the provided GFF3 file is acceptable

        If only a single record is contained in both sequences and GFF, they
        are assumed to be the same.

        Returns:
            True if only a single entry is contained by both inputs and
                    their sequence coordinates match
    """
    try:
        examiner = GFF.GFFExaminer()
        # file handle is automatically closed by GFF lib
        gff_data = examiner.available_limits(open(options.genefinding_gff3))
        # Check if at least one GFF locus appears in sequence
        gff_ids = set([n[0] for n in gff_data['gff_id']])

        single_entries = False

        if len(gff_ids) == 1 and len(sequences) == 1:
            # If both inputs only have one record, assume is the same,
            # but first check coordinate compatibility
            logging.info("GFF3 and sequence have only one record. Assuming is "
                         "the same as long as coordinates are compatible.")
            limit_info = dict(gff_type=['CDS'])

            record_iter = GFF.parse(open(options.genefinding_gff3), limit_info=limit_info)
            try:
                record = next(record_iter)
            except StopIteration:
                raise AntismashInputError("could not parse records from GFF3 file")

            if not record.features:
                raise AntismashInputError('GFF3 record %s contains no features' % record.id)

            coord_max = max([n.location.end.real for n in record.features])
            if coord_max > len(sequences[0]):
                logging.error('GFF3 record and sequence coordinates are not compatible.')
                raise AntismashInputError('incompatible GFF record and sequence coordinates')

            single_entries = True

        elif not gff_ids.intersection({seq.id for seq in sequences}):
            logging.error('No GFF3 record IDs match any sequence record IDs.')
            raise AntismashInputError("GFF3 record IDs don't match sequence file record IDs.")

        # Check GFF contains CDSs
        if not ('CDS',) in gff_data['gff_type']:
            logging.error('GFF3 does not contain any CDS.')
            raise AntismashInputError("no CDS features in GFF3 file.")

        # Check CDS are childless but not parentless
        if 'CDS' in set([n for key in examiner.parent_child_map(open(options.genefinding_gff3)) for n in key]):
            logging.error('GFF3 structure is not suitable. CDS features must be childless but not parentless.')
            raise AntismashInputError('GFF3 structure is not suitable.')

    except AssertionError as err:
        logging.error('Parsing %r failed: %s', options.genefinding_gff3, err)
        raise AntismashInputError(str(err)) from err
    return single_entries


def get_features_from_file(record: Record, handle: IO,
                           limit_to_seq_id: Union[bool, Dict[str, List[str]]] = False
                           ) -> List[SeqFeature]:
    """ Generates new SeqFeatures from a GFF file.

        Arguments:
            record: the Record that features belong to
            handle: a file handle/stream with the GFF contents
            limit_to_seq_id: False or a dictionary of GFF.parse options

        Returns:
            a list of SeqFeatures parsed from the GFF file
    """
    features = []
    try:
        gff_records = list(GFF.parse(handle, limit_info=limit_to_seq_id))
    except Exception as err:
        raise AntismashInputError("could not parse records from GFF3 file") from err

    for gff_record in gff_records:
        for feature in gff_record.features:
            if feature.type == 'CDS':
                new_features = [feature]
            else:
                new_features = check_sub(feature, record)
                if not new_features:
                    continue

            name = feature.id
            locus_tag = feature.qualifiers.get("locus_tag")

            for qtype in ["gene", "name", "Name"]:
                if qtype in feature.qualifiers:
                    name_tmp = feature.qualifiers[qtype][0]
                    # Assume name/Name to be sane if they don't contain a space
                    if " " in name_tmp:
                        continue
                    name = name_tmp
                    break

            for i, new_feature in enumerate(new_features):
                variant = name
                if len(new_features) > 1:
                    variant = "{0}_{1}".format(name, i)
                new_feature.qualifiers['gene'] = [variant]
                if locus_tag is not None:
                    new_feature.qualifiers["locus_tag"] = locus_tag
                features.append(new_feature)
    return features


def run(record: Record, single_entry: bool, options: ConfigType) -> None:
    """ The entry point of gff_parser.
        Generates new features and adds them to the provided record.

        Arguments:
            record: the secmet.Record instance to alter
            single_entry: if True, record and GFF ids must match
            options: an antismash.Config object, used only for fetching the GFF path

        Returns:
            None
    """
    # If there's only one sequence in both, read all, otherwise, read only appropriate part of GFF3.
    limit_info = False  # type: Union[bool, Dict[str, List[str]]]
    if not single_entry:
        limit_info = {'gff_id': [record.id]}

    with open(options.genefinding_gff3) as handle:
        features = get_features_from_file(record, handle, limit_info)
        for feature in features:
            try:
                record.add_biopython_feature(feature)
            except SecmetInvalidInputError as err:
                raise AntismashInputError from err


def generate_details_from_subfeature(sub_feature: SeqFeature,
                                     existing_qualifiers: Dict[str, List[str]],
                                     locations: List[FeatureLocation],
                                     trans_locations: List[FeatureLocation]) -> Set[str]:
    """ Finds the locations of a subfeature and any mismatching qualifiers

        Arguments:
            sub_feature: the GFF subfeature to work on
            existing_qualifiers: a dict of any existing qualifiers from other
                                 subfeatures
            locations: a list of any existing FeatureLocations from other
                       subfeatures
            trans_locations: a list of any existing FeatureLocations for
                             translations

        Returns:
            a set of qualifiers from the subfeature for which an existing
            qualifier existed but had a different value
    """
    mismatching_qualifiers = set()
    start = sub_feature.location.start.real
    end = sub_feature.location.end.real
    phase = int(sub_feature.qualifiers.get('phase', [0])[0])
    if sub_feature.strand == 1:
        start += phase
    else:
        end -= phase
    locations.append(FeatureLocation(start, end, strand=sub_feature.strand))
    # Make sure CDSs lengths are multiple of three. Otherwise extend to next full codon.
    # This only applies for translation.
    modulus = (end - start) % 3
    if modulus and sub_feature.strand == 1:
        end += 3 - modulus
    elif modulus and sub_feature.strand == -1:
        start -= 3 - modulus
    trans_locations.append(FeatureLocation(start, end, strand=sub_feature.strand))
    # For split features (CDSs), the final feature will have the same qualifiers as the children ONLY if
    # they're the same, i.e.: all children have the same "protein_ID" (key and value).
    for qual in sub_feature.qualifiers:
        if qual not in existing_qualifiers:
            existing_qualifiers[qual] = sub_feature.qualifiers[qual]
        elif existing_qualifiers[qual] != sub_feature.qualifiers[qual]:
            mismatching_qualifiers.add(qual)
    return mismatching_qualifiers


def check_sub(feature: SeqFeature, sequence: Record) -> List[SeqFeature]:
    """ Recursively checks a GFF feature for any subfeatures and generates any
        appropriate SeqFeature instances from them.
    """
    new_features = []
    locations = []  # type: List[FeatureLocation]
    trans_locations = []  # type: List[FeatureLocation]
    qualifiers = {}  # type: Dict[str, List[str]]
    mismatching_qualifiers = set()  # type: Set[str]
    for sub in feature.sub_features:
        if sub.sub_features:  # If there are sub_features, go deeper
            new_features.extend(check_sub(sub, sequence))
        elif sub.type == 'CDS':
            sub_mismatch = generate_details_from_subfeature(sub, qualifiers,
                                                            locations, trans_locations)
            mismatching_qualifiers.update(sub_mismatch)

    for qualifier in mismatching_qualifiers:
        del qualifiers[qualifier]
    if 'Parent' in qualifiers:
        del qualifiers['Parent']

    # if nothing to work on
    if not new_features and not locations:
        return []

    # Only works in tip of the tree, when there's no new_feature built yet. If there is,
    # it means the script just came out of a check_sub and it's ready to return.
    if not new_features:
        new_loc = locations[0]
        # construct a compound location if required
        if len(locations) > 1:
            locations = sorted(locations, key=lambda x: x.start.real)
            trans_locations = sorted(trans_locations, key=lambda x: x.start.real)
            if locations[0].strand == 1:
                new_loc = CompoundLocation(locations)
            else:
                new_loc = CompoundLocation(list(reversed(locations)))
                trans_locations = list(reversed(trans_locations))
        # TODO: use new secmet features
        new_feature = SeqFeature(new_loc)
        new_feature.qualifiers = qualifiers
        new_feature.type = 'CDS'
        new_features.append(new_feature)

    return new_features
