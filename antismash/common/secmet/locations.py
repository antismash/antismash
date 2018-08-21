# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Helper functions for location operations """

from typing import List, Tuple

from Bio.SeqFeature import FeatureLocation, CompoundLocation

# extend FeatureLocation and CompoundLocation with bio_start / bio_end until upstream implementation
FeatureLocation.bio_start = property(lambda self: self.end if self.strand == -1 else self.start )
FeatureLocation.bio_end = property(lambda self: self.start if self.strand == -1 else self.end )
CompoundLocation.bio_start = property(lambda self: self.parts[0].bio_start )
CompoundLocation.bio_end = property(lambda self: self.parts[-1].bio_end )


def convert_protein_position_to_dna(start: int, end: int, location: FeatureLocation) -> Tuple[int, int]:
    """ Convert a protein position to a nucleotide sequence position for use in generating
        new FeatureLocations from existing FeatureLocations and/or CompoundLocations.

        Arguments:
            position: the position in question, must be contained by the location
            location: the location of the related feature, for handling introns/split locations

        Returns:
            an int representing the calculated DNA location
    """
    if not 0 <= start < end <= len(location) // 3:
        raise ValueError("Protein positions %d and %d must be contained by %s" % (start, end, location))
    if location.strand == -1:
        dna_start = location.start + len(location) - end * 3
        dna_end = location.start + len(location) - start * 3
    else:
        dna_start = location.start + start * 3
        dna_end = location.start + end * 3

    # only CompoundLocations are complicated
    if not isinstance(location, CompoundLocation):
        if not location.start <= dna_start < dna_end <= location.end:
            raise ValueError(("Converted coordinates %d..%d "
                              "out of bounds for location %s") % (dna_start, dna_end, location))
        return dna_start, dna_end

    parts = sorted(location.parts, key=lambda x: x.start)
    gap = 0
    last_end = parts[0].start
    start_found = False
    end_found = False
    for part in parts:
        if start_found and end_found:
            break
        gap += part.start - last_end
        if not start_found and dna_start + gap in part:
            start_found = True
            dna_start = dna_start + gap
        if not end_found and dna_end + gap - 1 in part:
            end_found = True
            dna_end = dna_end + gap

        last_end = part.end

    assert start_found
    assert end_found

    if not location.start <= dna_start < dna_end <= location.end:
        raise ValueError(("Converted coordinates %d..%d "
                          "out of bounds for location %s") % (dna_start, dna_end, location))
    return dna_start, dna_end


def location_bridges_origin(location: CompoundLocation) -> bool:
    """ Determines if a CompoundLocation would cross the origin of a record.

        Arguments:
            location: the CompoundLocation to check

        Returns:
            False if the location does not bridge the origin or if the location
            is of indeterminate strand, otherwise True
    """
    assert isinstance(location, (FeatureLocation, CompoundLocation)), type(location)

    # if it's not compound, it can't bridge at all
    if not isinstance(location, CompoundLocation):
        return False

    # invalid strands mean direction can't be determined, may need to be an error
    if location.strand not in [1, -1]:
        return False

    for i, part in enumerate(location.parts[1:]):
        if location.strand == 1:
            if part.start <= location.parts[i].end:
                return True
        else:
            if part.start >= location.parts[i].end:
                return True
    return False
