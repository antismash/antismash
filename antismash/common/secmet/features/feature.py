# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The base class of all secmet features """

from collections import OrderedDict
from typing import Any, Dict, List, Optional, Tuple, Union

from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Seq import Seq

from antismash.common.secmet.locations import (
    convert_protein_position_to_dna,
)


class Feature:
    """ The base class of any feature. Contains only a location, the label of the
        subclass, the 'notes' qualifier, and other qualifiers not tracked by any
        subclass.
    """
    __slots__ = ["location", "notes", "type", "_qualifiers", "created_by_antismash"]

    def __init__(self, location: FeatureLocation, feature_type: str,
                 created_by_antismash: bool = False) -> None:
        assert isinstance(location, (FeatureLocation, CompoundLocation)), type(location)
        assert location.start < location.end, "Feature location invalid"
        self.location = location
        self.notes = []  # type: List[str]
        assert feature_type
        self.type = str(feature_type)
        self._qualifiers = OrderedDict()  # type: Dict[str, List[str]]
        self.created_by_antismash = bool(created_by_antismash)

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

    def get_sub_location_from_protein_coordinates(self, start: int, end: int) -> FeatureLocation:
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
            raise ValueError("Protein end coordinate must be contained by the feature")

        if start >= end:
            raise ValueError("Protein start coordinate must be less than the end coordinate")

        dna_start, dna_end = convert_protein_position_to_dna(start, end, self.location)

        if not 0 <= dna_start - self.location.start < self.location.end - 2:
            raise ValueError("Protein coordinate start %d (nucl %d) is outside feature %s" % (start, dna_start, self))
        if not 2 < dna_end - self.location.start - 1 <= self.location.end:
            raise ValueError("Protein coordinate end %d (nucl %d) is outside feature %s" % (end, dna_end, self))

        if not isinstance(self.location, CompoundLocation):
            return FeatureLocation(dna_start, dna_end, self.location.strand)

        new_locations = []
        # do not sort here either way, as we want biological order. instead, reverse the location parts if needed
        for location in reversed(self.location.parts) if self.location.strand == -1 else self.location.parts:
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

    def get_qualifier(self, key: str) -> Optional[Tuple]:
        """ Fetches a qualifier by key and returns a tuple of items stored under
            that key or None if the key was not present.
        """
        qualifier = self._qualifiers.get(key)
        if qualifier:
            return tuple(qualifier)
        return None

    def overlaps_with(self, other: Union["Feature", FeatureLocation]) -> bool:
        """ Returns True if the given feature overlaps with this feature.
            This operation is commutative, a.overlaps_with(b) is equivalent to
            b.overlaps_with(a).
        """
        if isinstance(other, Feature):
            location = other.location
        elif isinstance(other, FeatureLocation):
            location = other
        else:
            raise TypeError("Container must be a Feature or a FeatureLocation, not %s" % type(other))
        return (self.location.start in location
                or self.location.end - 1 in location
                or location.start in self.location
                or location.end - 1 in self.location)

    def is_contained_by(self, other: Union["Feature", FeatureLocation]) -> bool:
        """ Returns True if the given feature is wholly contained by this
            feature.
        """
        end = self.location.end - 1  # to account for the non-inclusive end
        if isinstance(other, Feature):
            return self.location.start in other.location and end in other.location
        if isinstance(other, FeatureLocation):
            return self.location.start in other and end in other
        raise TypeError("Container must be a Feature or a FeatureLocation, not %s" % type(other))

    def to_biopython(self, qualifiers: Dict[str, Any] = None) -> List[SeqFeature]:
        """ Converts this feature into one or more SeqFeature instances.

            Subclasses must manage their own attributes and potential extra
            features.
        """
        feature = SeqFeature(self.location, type=self.type)
        quals = self._qualifiers.copy()
        notes = self._qualifiers.get("note", []) + self.notes
        if qualifiers:
            notes += qualifiers.pop("note", [])
            quals.update(qualifiers)
        if notes:
            # sorting helps with consistency and comparison testing
            quals["note"] = sorted(notes)
        if self.created_by_antismash:
            quals["tool"] = ["antismash"]
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
            feature._qualifiers.update(leftovers)
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
