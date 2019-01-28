# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The base class of all secmet features """

from collections import OrderedDict
import logging
from typing import Any, Dict, List, Tuple, Union
from typing import Optional  # comment hints  # pylint: disable=unused-import

from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Seq import Seq

from antismash.common.secmet.locations import (
    convert_protein_position_to_dna,
    location_bridges_origin,
    location_contains_other,
    locations_overlap,
)

from ..errors import SecmetInvalidInputError
from ..locations import AfterPosition, BeforePosition, Location


def _adjust_location_by_offset(location: Location, offset: int) -> Location:
    """ Adjusts the given location to account for an offset (e.g. start_codon)
    """
    assert -2 <= offset <= 2, "invalid offset %d" % offset

    def adjust_single_location(part: FeatureLocation) -> FeatureLocation:
        """ only functions on FeatureLocation """
        assert not isinstance(part, CompoundLocation)
        start = part.start
        end = part.end
        if part.strand == -1:
            end = type(end)(end + offset)
        else:
            start = type(start)(start + offset)
        return FeatureLocation(start, end, part.strand)

    if isinstance(location, CompoundLocation):
        part = location.parts[0]
        if location.strand == -1:
            assert part.end == location.end
        else:
            assert part.start == location.start
        location = CompoundLocation([adjust_single_location(part)] + location.parts[1:])
    else:
        location = adjust_single_location(location)

    return location


class Feature:
    """ The base class of any feature. Contains only a location, the label of the
        subclass, the 'notes' qualifier, and other qualifiers not tracked by any
        subclass.
    """
    __slots__ = ["location", "notes", "type", "_qualifiers", "created_by_antismash",
                 "_original_codon_start"]

    def __init__(self, location: Location, feature_type: str,
                 created_by_antismash: bool = False) -> None:
        assert isinstance(location, (FeatureLocation, CompoundLocation)), type(location)
        if location_bridges_origin(location):
            raise ValueError("Features that bridge the record origin cannot be directly created: %s" % location)
        assert location.start <= location.end, "Feature location invalid: %s" % location
        self.location = location
        self.notes = []  # type: List[str]
        assert feature_type
        self.type = str(feature_type)
        self._qualifiers = OrderedDict()  # type: Dict[str, Optional[List[str]]]
        self.created_by_antismash = bool(created_by_antismash)
        self._original_codon_start = None  # type: Optional[int]

    @property
    def strand(self) -> int:
        """ A simple wrapper to access the location strand directly.
        """
        return self.location.strand

    def extract(self, sequence: Seq) -> Seq:
        """ Extracts the section of the given sequence that this feature covers.

            Return type is always a Seq, unlike location.extract.
        """
        assert isinstance(sequence, Seq)
        return self.location.extract(sequence)

    def get_sub_location_from_protein_coordinates(self, start: int, end: int) -> Location:
        """ Generates a FeatureLocation for a protein sequence based on the start
            and end positions within the features protein sequence.

            The start position is inclusive and the end position is exclusive.

            Arguments:
                start: a position between 0 and len(feature.location) // 3 - 1, inclusive
                end: a position between 1 and len(feature.location) // 3, inclusive
        """
        if not 0 <= start <= len(self.location) // 3 - 1:
            raise ValueError("Protein start coordinate must be contained by the feature")

        if not 1 <= end <= len(self.location) // 3:
            if end > 0 and (self.location.strand != -1 and isinstance(self.location.end, AfterPosition)
                            or self.location.strand == -1 and isinstance(self.location.start, BeforePosition)):
                # since some NCBI records (e.g. CP006567.1 and AWEV00000000.1) infer
                # an amino from a partial but don't increase the location to match
                logging.warning("%s protein coordinate %d truncated to match ambiguous end", str(self), end)
                end = len(self.location) // 3
            else:
                raise ValueError("Protein end coordinate must be contained by the feature")

        if start >= end:
            raise ValueError("Protein start coordinate must be less than the end coordinate")

        dna_start, dna_end = convert_protein_position_to_dna(start, end, self.location)
        if not dna_start < dna_end:
            raise ValueError("Invalid protein coordinate conversion (start %d, end %d)" % (dna_start, dna_end))
        if dna_start not in self.location:
            raise ValueError("Protein coordinate start %d (nucl %d) is outside feature %s" % (start, dna_start, self))
        # end check is more complicated as 'in' is inclusive and end is exclusive
        end_contained = dna_end in self.location or dna_end == self.location.end
        for part in self.location.parts:
            if dna_end in part or dna_end == part.end:
                end_contained = True
                break
        if not end_contained:
            raise ValueError("Protein coordinate end %d (nucl %d) is outside feature %s" % (end, dna_end, self))

        if not isinstance(self.location, CompoundLocation):
            return FeatureLocation(dna_start, dna_end, self.location.strand)

        new_locations = []
        for location in sorted(self.location.parts, key=lambda x: x.start):
            if dna_start in location:
                new = FeatureLocation(dna_start, location.end, self.location.strand)
                # the end could also be in this part
                if dna_end - 1 in location:
                    # can't be a compound location with only one, so return a simple one
                    return FeatureLocation(new.start, dna_end, new.strand)
                new_locations.append(new)
            elif dna_end - 1 in location:  # 'in' uses start <= value < end
                new = FeatureLocation(location.start, dna_end, self.location.strand)
                new_locations.append(new)
                break
            elif new_locations:  # found a start, but haven't yet found an end
                new_locations.append(location)

        if not new_locations:
            raise ValueError(("Could not create compound location from"
                              " %s and internal protein coordinates %d..%d (dna %d..%d)") % (
                                str(self.location), start, end, dna_start, dna_end))
        if self.location.strand == -1:
            new_locations.reverse()

        if len(new_locations) == 1:
            return new_locations[0]

        return CompoundLocation(new_locations)

    def get_qualifier(self, key: str) -> Union[None, Tuple[str, ...], bool]:
        """ Fetches a qualifier by key and returns
            - a tuple of items stored under that key,
            - True if the key is present without a value,
            - or None if the key was not present.
        """
        # not present
        if key not in self._qualifiers:
            return None
        qualifier = self._qualifiers[key]
        # present with value(s) (e.g. /key=value)
        if qualifier:
            return tuple(qualifier)
        # present without value (e.g. /pseudo)
        assert qualifier is None
        return True

    def overlaps_with(self, other: Union["Feature", Location]) -> bool:
        """ Returns True if the given feature overlaps with this feature.
            This operation is commutative, a.overlaps_with(b) is equivalent to
            b.overlaps_with(a).
        """
        if isinstance(other, Feature):
            location = other.location
        elif isinstance(other, (CompoundLocation, FeatureLocation)):
            location = other
        else:
            raise TypeError("Container must be a Feature, CompoundLocation, or FeatureLocation, not %s" % type(other))
        return locations_overlap(self.location, location)

    def is_contained_by(self, other: Union["Feature", Location]) -> bool:
        """ Returns True if the given feature is wholly contained by this
            feature.
        """
        if isinstance(other, Feature):
            return location_contains_other(other.location, self.location)
        if isinstance(other, (CompoundLocation, FeatureLocation)):
            return location_contains_other(other, self.location)
        raise TypeError("Container must be a Feature, CompoundLocation or FeatureLocation, not %s" % type(other))

    def to_biopython(self, qualifiers: Dict[str, Any] = None) -> List[SeqFeature]:
        """ Converts this feature into one or more SeqFeature instances.

            Subclasses must manage their own attributes and potential extra
            features.
        """
        feature = SeqFeature(self.location, type=self.type)
        quals = self._qualifiers.copy()
        notes = self._qualifiers.get("note", [])
        assert notes is not None
        notes.extend(self.notes)
        if qualifiers:
            notes += qualifiers.pop("note", [])
            quals.update(qualifiers)
        if notes:
            # sorting helps with consistency and comparison testing
            quals["note"] = sorted(notes)
        if self.created_by_antismash:
            quals["tool"] = ["antismash"]
        if self._original_codon_start is not None:
            start = int(self._original_codon_start)
            quals["codon_start"] = [str(start + 1)]
            # adjust location back if neccessary
            if self.location.strand == -1:
                start *= -1
            if self._original_codon_start != 0:
                feature.location = _adjust_location_by_offset(feature.location, -start)
        # sorted here to match the behaviour of biopython
        for key, val in sorted(quals.items()):
            feature.qualifiers[key] = val
        assert isinstance(feature.qualifiers, dict)
        return [feature]

    def __lt__(self, other: "Feature") -> bool:
        """ Allows sorting Features by location without key complication """
        assert isinstance(other, Feature)
        if self.location.start < other.location.start:
            return True
        elif self.location.start == other.location.start:
            return self.location.end < other.location.end
        return False

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return "%s(%s)" % (self.type, self.location)

    @staticmethod
    def from_biopython(bio_feature: SeqFeature, feature: "Feature" = None,
                       leftovers: Dict[str, Any] = None) -> "Feature":
        """ Converts a SeqFeature into a single Feature instance.

            Arguments:
                bio_feature: the SeqFeature to convert
                feature: a optional Feature instance to update with the values
                         this class tracks
                leftovers: any qualifiers remaining from the original SeqFeature
                           that have not been used by any subclass

            Returns:
                a Feature instance
        """
        if feature is None:
            feature = Feature(bio_feature.location, bio_feature.type)
            if not leftovers:
                assert isinstance(bio_feature.qualifiers, dict)
                leftovers = bio_feature.qualifiers.copy()
            feature.notes = leftovers.pop("note", [])
        else:
            assert isinstance(feature, Feature)
        if leftovers:
            feature.created_by_antismash = leftovers.get("tool") == ["antismash"]
            if "codon_start" in leftovers:
                codon_start = int(leftovers.pop("codon_start")[0]) - 1
                if not 0 <= codon_start <= 2:
                    raise SecmetInvalidInputError("invalid codon_start qualifier: %d" % (codon_start + 1))
                feature._original_codon_start = codon_start  # very much private, so pylint: disable=protected-access
                if feature.location.strand == -1:
                    codon_start *= -1
                if codon_start != 0:
                    # adjust the location for now until converting back to biopython if required
                    feature.location = _adjust_location_by_offset(feature.location, codon_start)

            feature._qualifiers.update(leftovers)  # shouldn't be a public thing, so pylint: disable=protected-access
        else:
            feature.created_by_antismash = False
        return feature

    @staticmethod
    def make_qualifiers_copy(bio_feature: SeqFeature) -> Dict[str, Any]:
        """ Makes a shallow copy of a SeqFeature's qualifiers. Only the 'notes'
            key will have a copy taken at a deeper level.
        """
        qualifiers = bio_feature.qualifiers.copy()
        if "note" in qualifiers:
            qualifiers["note"] = qualifiers["note"].copy()
        return qualifiers
