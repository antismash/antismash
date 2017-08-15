import logging
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from BCBio import GFF


def check_gff_suitability(options, sequences):
    """
        Checks that the provided GFF3 file is acceptable

        If only a single record is contained in both sequences and GFF, they
        are assumed to be the same.

        Returns:
            bool - True if only a single entry is contained by both inputs and
                    their sequence coordinates match
    """
    if not options.genefinding_gff3:
        return

    try:
        examiner = GFF.GFFExaminer()
        # file handle is automatically closed by GFF lib
        gff_data = examiner.available_limits(open(options.genefinding_gff3))
        # Check if at least one GFF locus appears in sequence
        gff_ids = set([n[0] for n in gff_data['gff_id']])

        if len(gff_ids) == 1 and len(sequences) == 1:
            # If both inputs only have one record, assume is the same, but first check coordinate compatibility
            logging.info("GFF3 and sequence have only one record. Assuming is "
                         "the same as long as coordinates are compatible.")
            limit_info = dict(gff_type=['CDS'])

            record_iter = GFF.parse(open(options.genefinding_gff3), limit_info=limit_info)
            record = next(record_iter)

            coord_max = max([n.location.end.real for n in record.features])
            if coord_max > len(sequences[0]):
                logging.error('GFF3 record and sequence coordinates are not compatible.')
                raise ValueError('Incompatible GFF record and sequence coordinates')
            else:
                single_entries = True

        elif len(gff_ids.intersection({seq.id for seq in sequences})) == 0:
            logging.error('No GFF3 record IDs match any sequence record IDs.')
            raise ValueError("GFF3 record IDs don't match sequence file record IDs.")

        else:
            single_entries = False

        # Check GFF contains CDSs
        if not ('CDS',) in gff_data['gff_type']:
            logging.error('GFF3 does not contain any CDS.')
            raise ValueError("No CDS features in GFF3 file.")

        # Check CDS are childless but not parentless
        if 'CDS' in set([n for key in examiner.parent_child_map(open(options.genefinding_gff3)) for n in key]):
            logging.error('GFF3 structure is not suitable. CDS features must be childless but not parentless.')
            raise ValueError('GFF3 structure is not suitable.')

    except AssertionError as e:
        logging.error('Parsing %r failed: %s', options.genefinding_gff3, e)
        raise
    return single_entries

def get_features_from_file(seq_record, handle, limit_to_seq_id=False):
    features = []
    for record in GFF.parse(handle, limit_info=limit_to_seq_id):
        for feature in record.features:
            if feature.type == 'CDS':
                new_features = [feature]
            else:
                new_features = check_sub(feature, seq_record)
                if not new_features:
                    continue

            name = feature.id
            locus_tag = None
            if "locus_tag" in feature.qualifiers:
                locus_tag = feature.qualifiers["locus_tag"]

            for qtype in ["gene", "name", "Name"]:
                if qtype in feature.qualifiers:
                    name_tmp = feature.qualifiers[qtype][0]
                    # Assume name/Name to be sane if they don't contain a space
                    if " " in name_tmp:
                        continue
                    name = name_tmp
                    break

            for i, n in enumerate(new_features):
                variant = name
                if len(new_features) > 1:
                    variant = "{0}_{1}".format(name, i)
                n.qualifiers['gene'] = [variant]
                if locus_tag is not None:
                    n.qualifiers["locus_tag"] = locus_tag
                features.append(n)
    return features

def run(sequence, single_entry, options):
    handle = open(options.genefinding_gff3)
    # If there's only one sequence in both, read all, otherwise, read only appropriate part of GFF3.
    if single_entry:
        limit_info = False
    else:
        limit_info = dict(gff_id=[sequence.id])

    features = get_features_from_file(sequence, handle, limit_info)
    sequence.features.extend(features)


def check_sub(feature, sequence):
    new_features = []
    locations = []
    trans_locations = []
    qualifiers = {}
    topop = []
    for sub in feature.sub_features:
        if sub.sub_features:  # If there are sub_features, go deeper
            new_features.extend(check_sub(sub, sequence))
        elif sub.type == 'CDS':
            loc = [sub.location.start.real, sub.location.end.real]
            if 'phase' in sub.qualifiers:
                phase = int(sub.qualifiers['phase'][0])
                if sub.strand == 1:
                    loc[0] += phase
                else:
                    loc[1] -= phase
            locations.append(FeatureLocation(loc[0], loc[1], strand=sub.strand))
            # Make sure CDSs lengths are multiple of three. Otherwise extend to next full codon.
            # This only applies for translation.
            modulus = (loc[1] - loc[0]) % 3
            if modulus == 0:
                trans_locations.append(FeatureLocation(loc[0], loc[1], strand=sub.strand))
            else:
                if sub.strand == 1:
                    trans_locations.append(FeatureLocation(loc[0], loc[1] + (3 - modulus), strand=sub.strand))
                else:
                    trans_locations.append(FeatureLocation(loc[0] - (3 - modulus), loc[1], strand=sub.strand))
            # For split features (CDSs), the final feature will have the same qualifiers as the children ONLY if
            # they're the same, i.e.: all children have the same "protein_ID" (key and value).
            for qual in sub.qualifiers.keys():
                if not qual in qualifiers:
                    qualifiers[qual] = sub.qualifiers[qual]
                if qual in qualifiers and not qualifiers[qual] == sub.qualifiers[qual]:
                    topop.append(qual)

    for n in topop:  # Pop mismatching qualifers over split features
        qualifiers.pop(n, None)
    qualifiers.pop('Parent', None)  # Pop parent.

    # Only works in tip of the tree, when there's no new_feature built yet. If there is,
    # it means the script just came out of a check_sub and it's ready to return.
    if not new_features:
        if len(locations) > 1:
            locations = sorted(locations, key=lambda x: x.start.real)
            trans_locations = sorted(trans_locations, key=lambda x: x.start.real)
            if locations[0].strand == 1:
                new_loc = CompoundLocation(locations)
            else:
                new_loc = CompoundLocation(list(reversed(locations)))
                trans_locations = reversed(trans_locations)
        elif not locations:
            return new_features
        else:
            new_loc = locations[0]
        # TODO: use new secmet features
        new_feature = SeqFeature(new_loc)
        new_feature.qualifiers = qualifiers
        new_feature.type = 'CDS'
        trans = ''.join([n.extract(sequence.seq).translate(stop_symbol='')._data for n in trans_locations])
        new_feature.qualifiers['translation'] = [trans]
        new_features.append(new_feature)

    return new_features
