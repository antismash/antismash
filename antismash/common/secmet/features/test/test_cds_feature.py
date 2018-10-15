# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=no-self-use,protected-access,missing-docstring

import unittest

from Bio.Seq import Seq

from antismash.common.secmet.features import FeatureLocation, CDSFeature
from antismash.common.secmet.features.feature import CompoundLocation


class TestCDSFeature(unittest.TestCase):
    def test_required_identifiers(self):
        with self.assertRaises(ValueError):
            CDSFeature(FeatureLocation(1, 5, 1), translation="A")
        assert CDSFeature(FeatureLocation(1, 5, 1), locus_tag="foo", translation="A")
        assert CDSFeature(FeatureLocation(1, 5, 1), protein_id="foo", translation="A")
        assert CDSFeature(FeatureLocation(1, 5, 1), gene="foo", translation="A")

    def test_bad_strand(self):
        with self.assertRaisesRegex(ValueError, "Strand must be"):
            CDSFeature(FeatureLocation(1, 5, 0), locus_tag="test", translation="A")

    def test_invalid_qualifier(self):
        cds = CDSFeature(FeatureLocation(1, 5, 1), locus_tag="test", translation="A")
        for bad in ["bad", ["stuff"], {}, 1]:
            with self.assertRaisesRegex(TypeError, "can only be set to an instance of SecMetQualifier"):
                cds.sec_met = bad

    def test_bad_translation(self):
        loc = FeatureLocation(1, 5, 1)
        for trans in [None, ""]:
            with self.assertRaisesRegex(ValueError, "requires a valid translation"):
                CDSFeature(loc, locus_tag="test", translation=trans)


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
        assert new == self.location, "%s != %s" % (str(new), str(self.location))
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
