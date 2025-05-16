# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature

from antismash.common.secmet.errors import SecmetInvalidInputError
from antismash.common.secmet.features import FeatureLocation, CDSFeature
from antismash.common.secmet.features.cds_feature import (
    _is_valid_translation_length,
    _translation_fits_in_record,
    MAX_TRANSLATION_LENGTH,
)
from antismash.common.secmet.features.feature import CompoundLocation
from antismash.common.secmet.locations import AfterPosition, BeforePosition
from antismash.common.secmet.qualifiers import SecMetQualifier
from antismash.common.secmet.qualifiers.gene_functions import GeneFunction
from antismash.common.test.helpers import DummyRecord


class TestCDSFeature(unittest.TestCase):
    def test_required_identifiers(self):
        with self.assertRaisesRegex(ValueError, "requires at least one of: gene, protein_id, locus_tag"):
            CDSFeature(FeatureLocation(1, 5, 1), translation="A")
        assert CDSFeature(FeatureLocation(1, 5, 1), locus_tag="foo", translation="A")
        assert CDSFeature(FeatureLocation(1, 5, 1), protein_id="foo", translation="A")
        assert CDSFeature(FeatureLocation(1, 5, 1), gene="foo", translation="A")

    def test_bad_strand(self):
        for strand in [0, None]:
            with self.assertRaisesRegex(ValueError, "invalid strand"):
                CDSFeature(FeatureLocation(1, 5, strand), locus_tag="test", translation="A")

    def test_invalid_qualifier(self):
        cds = CDSFeature(FeatureLocation(1, 5, 1), locus_tag="test", translation="A")
        for bad in ["bad", ["stuff"], {}, 1, None]:
            with self.assertRaisesRegex(TypeError, "can only be set to an instance of SecMetQualifier"):
                cds.sec_met = bad

    def test_bad_translation(self):
        loc = FeatureLocation(1, 5, 1)
        for trans in [None, "A?", "A!", ""]:
            with self.assertRaisesRegex(ValueError, "valid translation required|invalid translation characters"):
                CDSFeature(loc, locus_tag="test", translation=trans)

    def test_max_translation_length(self):
        loc = FeatureLocation(0, (MAX_TRANSLATION_LENGTH - 1) * 3, 1)
        CDSFeature(loc, locus_tag="test", translation="A" * (MAX_TRANSLATION_LENGTH - 1))

        loc = FeatureLocation(0, MAX_TRANSLATION_LENGTH * 3, 1)
        with self.assertRaisesRegex(ValueError, "translation too long"):
            CDSFeature(loc, locus_tag="test", translation="A" * MAX_TRANSLATION_LENGTH)

    def test_methionine_starts(self):
        """ All gene translations should start with methionine, even if using
            an alternate start codon
        """
        translation = "LAGIC"
        expected = "MAGIC"
        # via constructor
        cds = CDSFeature(FeatureLocation(0, 15, 1), locus_tag="test", translation=translation)
        assert cds.translation == expected
        # and via setting the property
        cds.translation = translation
        assert cds.translation == expected

    def test_methionine_starts_when_ambiguous(self):
        # ambiguous start coordinates should not be adjusted to methionines
        translation = "LAGIC"
        for location in [
                FeatureLocation(BeforePosition(10), 25, 1),
                FeatureLocation(10, AfterPosition(25), -1),
        ]:
            cds = CDSFeature(location, locus_tag="test", translation=translation)
            assert cds.translation == translation, location


class TestCDSBiopythonConversion(unittest.TestCase):
    def setUp(self):
        self.cds = CDSFeature(FeatureLocation(0, 12, 1),
                              translation="MAAA",
                              locus_tag="loctag",
                              gene="gene",
                              protein_id="prot_id")

    def convert(self):
        bio_features = self.cds.to_biopython()
        assert isinstance(bio_features, list)
        assert len(bio_features) == 1
        return bio_features[0]

    def test_basics(self):
        bio = self.convert()
        assert bio.location == self.cds.location
        assert bio.qualifiers["locus_tag"] == ["loctag"]
        assert bio.qualifiers["gene"] == ["gene"]
        assert bio.qualifiers["protein_id"] == ["prot_id"]
        assert bio.qualifiers["translation"] == ["MAAA"]

        regen = CDSFeature.from_biopython(bio)
        assert regen.location == self.cds.location
        assert regen.locus_tag == self.cds.locus_tag
        assert regen.gene == self.cds.gene
        assert regen.protein_id == self.cds.protein_id

    def test_without_genefunctions(self):
        bio = self.convert()
        assert "gene_functions" not in bio.qualifiers
        assert "gene_kind" not in bio.qualifiers

        regen = CDSFeature.from_biopython(bio)
        assert not regen.gene_functions

    def test_with_genefunctions(self):
        self.cds.gene_functions.add(GeneFunction.ADDITIONAL, "testtool", "dummy")
        bio = self.convert()
        assert "gene_functions" in bio.qualifiers
        assert bio.qualifiers["gene_kind"] == [str(self.cds.gene_function)] == ["biosynthetic-additional"]

        regen = CDSFeature.from_biopython(bio)
        assert regen.gene_function == self.cds.gene_function
        assert regen.gene_functions.get_by_tool("testtool") == self.cds.gene_functions.get_by_tool("testtool")

    def test_without_secmet(self):
        assert not self.cds.sec_met
        bio = self.convert()
        assert "sec_met" not in bio.qualifiers  # for detecting legacy versions
        assert "sec_met_domain" not in bio.qualifiers

        regen = CDSFeature.from_biopython(bio)
        assert not regen.sec_met

    def test_with_secmet(self):
        domains = [SecMetQualifier.Domain("testA", 0.1, 1.1, 3, "test"),
                   SecMetQualifier.Domain("testB", 5.1, 3.9, 5, "dummy")]
        self.cds.sec_met = SecMetQualifier(domains)
        bio = self.convert()
        assert "sec_met" not in bio.qualifiers  # again, detecting leftover legacy versions
        assert len(bio.qualifiers["sec_met_domain"]) == 2
        assert bio.qualifiers["sec_met_domain"] == list(map(str, domains))

        regen = CDSFeature.from_biopython(bio)
        assert regen.sec_met
        assert len(regen.sec_met.domains) == len(domains)
        assert regen.sec_met.domains == domains

    def test_mixed_strand(self):
        bio = self.cds.to_biopython()[0]
        for location in [CompoundLocation([FeatureLocation(1, 5, strand=-1), FeatureLocation(8, 10, strand=1)]),
                         CompoundLocation([FeatureLocation(1, 5, strand=1), FeatureLocation(8, 10, strand=None)])]:
            bio.location = location
            with self.assertRaisesRegex(ValueError, "compound locations with mixed strands"):
                CDSFeature.from_biopython(bio)
        # compound locations starting with an invalid strand will be treated as per a non-compound wtih a bad strand

    def test_translation_outside_record(self):
        rec = DummyRecord(seq="A" * 10)
        for location in [FeatureLocation(0, AfterPosition(6), strand=1),
                         FeatureLocation(BeforePosition(4), 10, strand=-1)]:
            bio = SeqFeature(location, type="CDS")
            bio.qualifiers["translation"] = ["M" * 5]
            with self.assertRaisesRegex(SecmetInvalidInputError, "translation extends out of record"):
                CDSFeature.from_biopython(bio, record=rec)

    def test_invalid_translation_table(self):
        bio = self.cds.to_biopython()[0]
        bio.qualifiers["transl_table"] = ["11a"]
        with self.assertRaisesRegex(SecmetInvalidInputError, "invalid translation table"):
            CDSFeature.from_biopython(bio)

    def test_translation_with_codon_start(self):
        """ Ensures translation extraction takes place *after* location adjustment
            for frame shifts/codon start qualifier
        """
        cds = CDSFeature(FeatureLocation(6, 21, 1), locus_tag="test", translation="MAGIC")
        rec = DummyRecord(seq="xxxxxxATGGCAGGTATTTGTxxxxxx")
        rev = DummyRecord(seq="xxxxxxACAAATACCTGCCATxxxxxx")
        assert cds.location.extract(rec.seq).translate() == cds.translation

        bio = cds.to_biopython()[0]
        bio.qualifiers.pop("translation")

        for record, strand in [(rec, 1), (rev, -1)]:
            for shift in range(0, 3):
                start = cds.location.start - (shift if strand == 1 else 0)
                end = cds.location.end + (shift if strand == -1 else 0)
                bio.qualifiers["codon_start"] = [str(shift + 1)]
                bio.location = FeatureLocation(start, end, strand)
                assert "translation" not in bio.qualifiers
                new = CDSFeature.from_biopython(bio, record=record)
                assert new.translation == cds.translation

    def test_bad_name_generation(self):
        # a CDS with no identifiers but with a pseudo(gene) qualifier
        # used the coordinates, which is very awkward if it's an ambiguous position
        bio = SeqFeature(FeatureLocation(BeforePosition(6), 9, 1), "CDS")
        bio.qualifiers["pseudo"] = True
        bio.qualifiers["translation"] = "M"
        cds = CDSFeature.from_biopython(bio)
        assert cds.get_name() == "pseudo6_9"

        bio.qualifiers.pop("pseudo")
        cds = CDSFeature.from_biopython(bio)
        assert cds.get_name() == "cds6_9"


class TestCDSProteinLocation(unittest.TestCase):
    def setUp(self):
        self.magic_split = Seq("ATGGCAxxxxxxGGTxxxxxxATTTGT")
        self.magic = Seq("ATGGCAGGTATTTGT")
        self.translation = "MAGIC"
        self.sub_locations = [FeatureLocation(0, 6, strand=1),
                              FeatureLocation(12, 15, strand=1),
                              FeatureLocation(21, 27, strand=1)]
        self.location = CompoundLocation(self.sub_locations)
        self.cds = CDSFeature(self.location, locus_tag="compound", translation="A")

    def reverse_strand(self):
        self.magic = self.magic.reverse_complement()
        self.magic_split = self.magic_split.reverse_complement()
        self.sub_locations = [FeatureLocation(loc.start, loc.end, strand=loc.strand*-1) for loc in self.sub_locations]
        self.location = CompoundLocation(self.sub_locations[::self.sub_locations[0].strand])
        self.cds = CDSFeature(self.location, locus_tag="compound", translation="A")

    def test_simple_location_forward_complete(self):
        cds = CDSFeature(FeatureLocation(0, 15, 1), locus_tag="simple", translation="A")
        new = cds.get_sub_location_from_protein_coordinates(0, 5)
        extracted = new.extract(self.magic)
        assert extracted == self.magic
        assert extracted.translate() == self.translation

    def test_simple_location_forward_partial(self):
        cds = CDSFeature(FeatureLocation(0, 15, 1), locus_tag="simple", translation="A")
        for start, end in [(1, 5), (0, 3), (2, 3), (1, 4)]:
            print("testing", start, end)
            new = cds.get_sub_location_from_protein_coordinates(start, end)
            print(new)
            extracted = new.extract(self.magic)
            assert extracted == self.magic[start * 3: end * 3]
            assert extracted.translate() == self.translation[start:end]

    def test_compound_location_forward_full(self):
        new = self.cds.get_sub_location_from_protein_coordinates(0, 5)
        assert isinstance(new, CompoundLocation)
        assert len(new.parts) == 3
        print(list(map(str, self.cds.location.parts)))
        print(list(map(str, new.parts)))
        assert len(new) == len(self.cds.location)
        assert new == self.location, f"{new} != {self.location}"
        extracted = new.extract(self.magic_split)
        assert extracted == self.magic
        assert extracted.translate() == self.translation[0:5]

    def test_compound_forward_within_single(self):
        new = self.cds.get_sub_location_from_protein_coordinates(0, 2)
        assert isinstance(new, FeatureLocation)
        assert len(new) == 6
        assert new.start == 0
        assert new.end == 6
        assert new.extract(self.magic_split).translate() == self.translation[0:2]

        new = self.cds.get_sub_location_from_protein_coordinates(2, 3)
        assert isinstance(new, FeatureLocation)
        assert len(new) == 3
        assert new.start == 12
        assert new.end == 15
        assert new.extract(self.magic_split).translate() == self.translation[2:3]

    def test_compound_forward_over_multiple(self):
        new = self.cds.get_sub_location_from_protein_coordinates(2, 4)
        assert isinstance(new, CompoundLocation)
        print(list(map(str, self.cds.location.parts)))
        print(list(map(str, new.parts)))
        assert len(new.parts) == 2
        assert len(new) == 6
        assert new.parts[0].start == 12
        assert new.parts[0].end == 15
        assert new.parts[1].start == 21
        assert new.parts[1].end == 24
        assert new.extract(self.magic_split).translate() == self.translation[2:4]

    def test_compound_location_reverse_full(self):
        self.reverse_strand()
        cds = CDSFeature(self.location, locus_tag="compound", translation="A")
        new = cds.get_sub_location_from_protein_coordinates(0, 5)
        assert isinstance(new, CompoundLocation)
        assert len(new.parts) == 3
        print(list(map(str, cds.location.parts)))
        print(list(map(str, new.parts)))
        assert len(new) == len(cds.location)
        assert new.extract(self.magic_split).translate() == self.translation[0:5]

    def test_compound_location_reverse_single(self):
        self.reverse_strand()
        cds = CDSFeature(self.location, locus_tag="compound", translation="A")

        new = cds.get_sub_location_from_protein_coordinates(0, 2)
        assert isinstance(new, FeatureLocation)
        assert len(new) == 6
        assert new.start == 21
        assert new.end == 27
        assert new.extract(self.magic_split).translate() == self.translation[0:2]

        new = cds.get_sub_location_from_protein_coordinates(2, 3)
        assert isinstance(new, FeatureLocation)
        assert len(new) == 3
        assert new.start == 12
        assert new.end == 15
        assert new.extract(self.magic_split).translate() == self.translation[2:3]

    def test_compound_location_reverse_multiple(self):
        self.reverse_strand()
        cds = CDSFeature(self.location, locus_tag="compound", translation="A")

        new = cds.get_sub_location_from_protein_coordinates(2, 4)
        assert isinstance(new, CompoundLocation)
        print(list(map(str, cds.location.parts)))
        print(list(map(str, new.parts)))
        assert len(new.parts) == 2
        assert len(new) == 6
        assert new.parts[0].start == 12
        assert new.parts[0].end == 15
        assert new.parts[1].start == 3
        assert new.parts[1].end == 6
        assert new.extract(self.magic_split).translate() == self.translation[2:4]

    def test_frameshifted_location(self):
        location = CompoundLocation([FeatureLocation(3, 9, 1), FeatureLocation(8, 14, 1)])
        assert len(location) == 12
        seq = Seq("ATGATGAGCCCTCGTCTAGACTACAATGA")
        extracted = location.extract(seq)
        assert extracted == "ATGAGCCCCTCG"
        assert len(extracted) == len(location)
        translation = extracted.translate()
        assert translation == "MSPS"

        cds = CDSFeature(location, locus_tag="test", translation=translation)
        new = cds.get_sub_location_from_protein_coordinates(1, 3)
        assert isinstance(new, CompoundLocation)
        assert len(new.parts) == 2
        assert new.start == 6
        assert new.end == 11

    def test_complicated(self):
        parts = [FeatureLocation(121124, 122061, 1), FeatureLocation(122339, 122383, 1),
                 FeatureLocation(122559, 122666, 1), FeatureLocation(122712, 122874, 1),
                 FeatureLocation(123060, 123337, 1), FeatureLocation(123481, 123749, 1),
                 FeatureLocation(123809, 124032, 1), FeatureLocation(124091, 124193, 1),
                 FeatureLocation(124236, 124401, 1), FeatureLocation(124684, 124724, 1)]
        location = CompoundLocation(parts, operator="join")
        cds = CDSFeature(location, locus_tag="complicated", translation="A")
        seq = ("ATGAGCCCTCGTCTAGACTACAATGAAGGATACGATTCCGAAGACGAGGAGATCCCCCGTTACGTACACCAT"
               "TCTAGAGGAAAGAGTCATAGATCCGTGAGGACGTCAGGTCGCTCACGCACGTTGGATTACGACGGGGATGAT"
               "GAAGCTAGTGACCACGCTGCCCCCTCCGGGATTGATCGGGACGCTCGAGCCTGTCCAACATCTCGCAGATAT"
               "ACTGATGACTGCCTTGAGACACATAAATTTCGAGGTGCCCGCTCCTCTCGCTCCCGTGGACGAACCGATGAT"
               "AACAAGGTTTTGTACTACACCAAGTATCGCAGCCCGGCTAAGGACTTGCCTATCGAGCGTGATCCCGAGGGT"
               "ATTAATTTATTCAAGGTCCGACAGCACACACGGCCAAGTGACGCTCATGTGCCCAGTGGATACCGTGAGCCC"
               "TACGAAGTCAAGGTCGACGAGTATGAGGATGATCATCCCCGTACATGCACTAGCCGCCGTGACTCTAGACAG"
               "CCGAAAGTCTACAAGGTCCGGGTTGATGAGTACGAGGATAACCTCCCTGCACGCTCTCACACTGACTTTCGC"
               "GAGTCTCCACGGTCTGAAAGATGCTCTAGCCGCTACACCGAGGACTCGAAGCCTGGGGAGCTTCCTCCCCGC"
               "TCAGGGCCCTGTCGGTCCAGCAGGCCTTCTCCGGTCGATGAGGACGTCGAGTATGAGATCCGTGAGCCCCGA"
               "GGGCATCGCTCCAGTCGACACTCTACAGATGTTGACTTTCAGCCAGTAGAACAACATCCTCGCTTTGGACAA"
               "CGTGGACTCAGCAGACCTTCGCGGGTTGATGAGGAAGTCGATTATGAGATCCGTGAGCCCCGTGGCAATCGT"
               "GTCAGTCACGCTGCTCATGGTGACAGCCCCTGTCAGGACCAAAGCTCCAGGCATATCGGCATTCAATTGTGG"
               "AGTACGCGCGGACCCCGGGCGGCTGGCCGTGGCCGGGGTCCTGATGAGTCTGACGATGTTGAGCCCTAGGCA"
               "GGGAATTGCCGTAATGCTCTTCAAACTGTATAGCAAGCTCAGCATCAATTCTTTAACTGGCAGGCGCTCTGC"
               "TCGCGCGTTTCTCTCTTGGGGTGGTTGGTTTGACTGTAGATTTCCTCTTTCAAGGCTTCTAGATACACCTTT"
               "GGAAGATAGCAACGCTATGCAAGATATTTTTGATAATTCAAATCCTTTTTACACATGGAATAGCTGGTGTTC"
               "CTGTTTTATCTAGGCAATTGACCCACGCCATCTCGGTAGGTACGGTAAAAGCAAGCCGTAATCTCGTATGGC"
               "TTCATCCTTAGCATCGTATAGATCTCCACTCGGGACTCGGCCAGGGATCTTCCATCAATCAACGTGAAGAAG"
               "TCCAGCACCCCGCTGAATCATAATATCCTACCGATTCTGCTCTCTTCACCTCTAGATACCCCTCTAGACTCC"
               "TGTCAACATGTTCCGTACAGTCGAAGACCGCCCGACCCCAAAAGAGGTATATAACTGGCGGCTGTACACCGA"
               "GGCCACCATCATTGCCACTGGTACACTCTTGTGAGTAGGTGCTGTTGTAACGAAAAACATCCAACTGATCCG"
               "CCAGGTTCGGCTATGACTCGGCTTTTGTGGGAACTACCATTGCCCGCCAAAGCTTCGTTGATGCCTTCAACA"
               "TCGTCGAGTCGGAGGCGGCGGATATTTCAAGCAATATCACGTCAACCTTTCAGGCCGGCGCATTTTTCGGCG"
               "CCATCTTCTGCTTCTTGCCTGAGTGAAGCCGTTAGAGACGGTCTCACTGGCTAACCGGACCAAGTGACCGAC"
               "AAAATTGGGCGTAAATGGGCCCTTCAGGCAAACACACTGCTGTTTCTTATTGGCGCGATTGTGATGACGGCT"
               "GCAACACATCACCTTTCCTATATATGTAAGTCATATCCCCGTAGTAGTCAAGGTTGTTAACTAGAGCAGATG"
               "CTGGACGAGCTCTCACCGGCATCGCATGCGGCGCTATCACCGCGACCGTCCCCAGCTATATTGCCGAGCTGT"
               "CAATCGTGTCGATCCGGGGCTTCCTCACCGGGTTCTTCGAAGTCGCATACCAGATTGGTAGCTTGGTTGGAT"
               "TCTGGATCAACTATGGCATTAACGAGAACATGGACAACTCCTCGGCCGCAAGCTGGAGAGTGCCTATGGCAG"
               "TCCAGATCATCCCCGCAGGAGTCCTTTTCATTGGTGGCTTTTCCTCCATGAGAGTCCTCTCTGGCTGATGCG"
               "AAAAGACAGTGAGGATGCCGCGACGGCTGCCCTGGAGGCGTTGAGGAAACTGCCACGGTCTCATCAATGTAA"
               "TCTCCCACCAAGACTCAGGACATAGTCCCATGCTGACTATTTTAGATGTCCAGGAAGACATCGAGATGAACC"
               "GCACCAGGCTGCTGGAGGAAGCTCGGATCGCCGAGAAGTACGGACAAGGTTGGTTGGCATATATCCGAGGCG"
               "CACTCTTCGAGCTCTCGCGCCATGGGATGTGGAATCGTGTTCTGCTCGTCCTCTGTGCCTTTGCACTGCAGA"
               "ATATGTCGGGAGCTGCTGCTATCAACTACTATTCCCCCATACTCTTTGCGTCGTTGGGGATCACTGATGTCG"
               "CTCTGTATACAGGTATTTATGGCCTGGTAAAAGGTAAGTTCTTCTCCTTAAGTATCTCTGGCTGACAATAGG"
               "GATTAACTGATGAGTTTACAGCCGTCGCATCAATTATATTCTACGGCATTCTCATTGATATGTGGGGCCGCC"
               "GACGTCCGACCATTGTTTCGTCACTGGCCTGCCCTCTATGTCTCTGGTTTGTGGGTGCATACGTCAAAGTTG"
               "GGCATCCAGCCGATATCATAGACGCCGGCGGGGAATTGTCCCCCTCCACGGAGGCTGGTGGTAGAGCGGCGA"
               "CTGCGATGATTATGATCTACTCCGTCTTGTAAGTGCCCCTCACTTTTGAATGGGCTTCAGCTTGGAACTCGA"
               "GTAACTGGTATCCAGTTGGTCTTTTGGTCTCAACGGTATCCCCTGGATTGTCTCCGCCGAAATCTTCCCCGG"
               "CGCGCTGCGAAATCTCACGGGGACATGGGCTGCGCTGGTGCAATGGTATGCAATTCCCTTCACCTAGTATCC"
               "ATATCTAAATCAGCAGGTTGATCCAATTCGTTATCACCAAAGCTCTCCCGTACATCTTCAATAGCCTTGGGT"
               "ACGGGACGTGGTTCTTCTTCGCCTCCTGGATGCTGCTCGCTATCATTTGGTCATTCTTTTTTCTCCCGGAAA"
               "CCAAGGGGAAGACTCTCGATGAAATGCATACGATCTTGTACGTTTCTCTCCGTCGAAATGTGGTCTTGGCTA"
               "ATGAATCAGCGGCCATTCTCTCGCCGAAGAGCAGGGTAAGGGTGAGGTTCGAGATAACACTACTAAAAGTGA"
               "TCGGGAGGCTGTCTAGTCCAGTAGTTCTAGAGGACTATTGGCTGGATGATTCCTCTGATGATTTTTGATTGG"
               "TGGTGAAAATGTTGGATGTTTAATGCCAATGTACTGGGAGAGAACATGCCGATAGTACATACCGCTGTGTTG"
               "TATATCGAAGACGGTTGATTTATATATCTTAGTCTTTCAAAAGACGGCACTCACACAATCACACTTCGATGA")
        translation = ("MSPRLDYNEGYDSEDEEIPRYVHHSRGKSHRSVRTSGRSRTLDYDGDDEASDHAAPSGIDRDAR"
                       "ACPTSRRYTDDCLETHKFRGARSSRSRGRTDDNKVLYYTKYRSPAKDLPIERDPEGINLFKVRQ"
                       "HTRPSDAHVPSGYREPYEVKVDEYEDDHPRTCTSRRDSRQPKVYKVRVDEYEDNLPARSHTDFR"
                       "ESPRSERCSSRYTEDSKPGELPPRSGPCRSSRPSPVDEDVEYEIREPRGHRSSRHSTDVDFQPV"
                       "EQHPRFGQRGLSRPSRVDEEVDYEIREPRGNRVSHAAHGDSPCQDQSSRHIGIQLWTGVPVLSR"
                       "QLTHAISTPVNMFRTVEDRPTPKEVYNWRLYTEATIIATGTLLFGYDSAFVGTTIARQSFVDAF"
                       "NIVESEAADISSNITSTFQAGAFFGAIFCFLPEADAGRALTGIACGAITATVPSYIAELSIVSI"
                       "RGFLTGFFEVAYQIGSLVGFWINYGINENMDNSSAASWRVPMAVQIIPAGVLFIGGFSSMREDI"
                       "EMNRTRLLEEARIAEKYGQGWLAYIRGALFELSRHGMWNRVLLVLCAFALQNMSGAAAINYYSP"
                       "ILFASLGITDVALYTGIYGLVKAVASIIFYGILIDMWGRRRPTIVSSLACPLCLWFVGAYVKVG"
                       "HPADIIDAGGELSPSTEAGGRAATAMIMIYSVFWSFGLNGIPWIVSAEIFPGALRNLTGTWAAL"
                       "VQWLIQFVITKALPYIFNSLGYGTWFFFASWMLLAIIWSFFFLPETKGKTLDEMHTIFLSKDGT"
                       "HTITLR")
        new = cds.get_sub_location_from_protein_coordinates(353, 412)
        # pad the beginning to match the location
        assert new.extract(Seq("x" * location.start + seq)).translate() == translation[353:412]

    def test_extends_past_after(self):
        self.sub_locations[-1] = FeatureLocation(21, AfterPosition(29), strand=1)
        self.cds.location = CompoundLocation(self.sub_locations)

        new = self.cds.get_sub_location_from_protein_coordinates(0, 7)
        assert new.end == 27

    def test_extends_past_before(self):
        self.reverse_strand()
        self.sub_locations[0] = FeatureLocation(BeforePosition(2), self.sub_locations[0].end, strand=-1)
        self.cds.location = CompoundLocation(self.sub_locations[::-1])
        new = self.cds.get_sub_location_from_protein_coordinates(0, 7)
        assert new.start == 3


class TestTranslationLength(unittest.TestCase):
    def test_simple(self):
        location = FeatureLocation(0, 9)
        for good in ["AA", "AAA"]:
            assert _is_valid_translation_length(good, location)
        assert not _is_valid_translation_length("AAAA", location)

    def test_compound(self):
        location = CompoundLocation([FeatureLocation(0, 3), FeatureLocation(6, 9)])
        for good in ["A", "AA"]:
            assert _is_valid_translation_length(good, location)
        assert not _is_valid_translation_length("AAA", location)
        # and with an ambiguous end, that becomes ok
        location = CompoundLocation([FeatureLocation(0, 3), FeatureLocation(6, AfterPosition(11))])
        assert _is_valid_translation_length("AAA", location)
        # and reversed ambiguous end
        location = CompoundLocation([FeatureLocation(BeforePosition(0), 3, -1), FeatureLocation(6, 9, -1)])
        for good in ["A", "AA", "AAA"]:
            assert _is_valid_translation_length(good, location)

    def test_ambiguous_coords(self):
        translation = "A" * 10
        for strand in (-1, 1):
            # both ambiguous
            location = FeatureLocation(BeforePosition(10), AfterPosition(20), strand)
            assert _is_valid_translation_length(translation, location)
            # start ambiguous
            location = FeatureLocation(BeforePosition(10), 20, strand)
            assert _is_valid_translation_length(translation, location)
            # end ambiguous
            location = FeatureLocation(10, AfterPosition(20), strand)
            assert _is_valid_translation_length(translation, location)
            # neither
            location = FeatureLocation(10, 20, strand)
            assert not _is_valid_translation_length(translation, location)


class TestTranslationInRecord(unittest.TestCase):
    def setUp(self):
        self.run = _translation_fits_in_record

    def test_simple(self):
        location = FeatureLocation(0, AfterPosition(3), 1)
        size = 9
        assert not self.run(size, location, size - 3)
        assert self.run(size, location, size - 2)  # single ambiguous amino
        assert self.run(size, location, size)

        location = FeatureLocation(BeforePosition(3), 9, -1)
        assert not self.run(size + 3, location, size)
        assert self.run(size + 2, location, size)  # single ambiguous amino
        assert self.run(size, location, size)
