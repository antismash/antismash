# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
# pylint: disable=line-too-long
""" A mapping of A-domain substrate names in multiple naming conventions """

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class SubstrateName:
    """ A mapping of substrate names using multiple naming conventions """
    long: str
    short: str
    norine: str

    def to_json(self) -> dict[str, Any]:
        """ Creates a JSON representation of a SubstrateName """
        return dict(vars(self))

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> "SubstrateName":
        """ Creates a SubstrateName from a JSON representation """
        return cls(data["long"], data["short"], data["norine"])

    def __str__(self) -> str:
        """ Use the short name as the canonical string representation """
        return self.short


KNOWN_SUBSTRATES: list[SubstrateName] = [
    SubstrateName("(2S,3R)-2-amino-3-hydroxy-4-(4-nitrophenyl)butanoic acid", "n-oh-nitropheBut", "X"),  # noqa: E501
    SubstrateName("(2S,6R)-diamino-(5R,7)-dihydroxy-heptanoic acid", "dn-doh-Hep", "X"),
    SubstrateName("(4S)-5,5,5-trichloroleucine", "3ClLeu", "X"),
    SubstrateName("(E)-4-methylhex-2-enoic acid", "meHex", "X"),
    SubstrateName("(S,E)-2-amino-4-decenoic acid", "nDece", "X"),
    SubstrateName("1-(1,1-dimethylallyl)-tryptophan", "dmeaTrp", "X"),
    SubstrateName("1-aminocyclopropane-1-carboxylic acid", "ncyc-Prca", "X"),
    SubstrateName("1-pyrroline-5-carboxylic acid", "cPyr", "X"),
    SubstrateName("2,3-diamino-3-methylpropanoic acid", "dnmeProp", "X"),
    SubstrateName("2,3-diaminopropionic acid", "Dpr", "Dpr"),
    SubstrateName("2,3-dihydroxy-para-aminobenzoic acid", "2,3-doh-4-nBza", "X"),
    SubstrateName("2,3-dihydroxybenzoic acid", "2,3-dohBza", "diOH-Bz"),
    SubstrateName("2,3-dihydroxyhexadecanoic acid", "Dhhd", "X"),
    SubstrateName("2,4-diaminobutyric acid", "Dab", "Dab"),
    SubstrateName("2,4-dihydroxypentanoic acid", "DHPea", "X"),
    SubstrateName("2-(1-methylcyclopropyl)-D-glycine", "Cyc-D-Gly", "X"),
    SubstrateName("2-amino-3,5-dimethyl-4-hexenoic Acid", "ndmHex", "X"),
    SubstrateName("2-amino-3-hydroxycyclopent-2-enone", "noh-cycPen", "X"),
    SubstrateName("2-amino-6-hydroxy-4-methyl-8-oxodecanoic acid", "nohme-oxo-Dec", "X"),
    SubstrateName("2-aminoadipic acid", "Aad", "Aad"),
    SubstrateName("2-aminobutyric acid", "Abu", "Abu"),
    SubstrateName("2-aminoisobutyric acid", "Aib", "Aib"),
    SubstrateName("2-carboxy-6-hydroxyoctahydroindole", "cohOhi", "X"),
    SubstrateName("2-chloro-3,5-dihydroxy-4-methylphenylglycine", "Cl-me-Dpg", "X"),
    SubstrateName("2-chlorobenzoic acid", "2ClBza", "X"),
    SubstrateName("2-hydroxy-4-methylpentanoic acid", "ohmePea", "4Me-Hva"),
    SubstrateName("2-hydroxypent-4-enoic acid", "Hpea", "X"),
    SubstrateName("2-ketoglutaric acid", "Kga", "X"),
    SubstrateName("2-ketoisocaproic acid", "Kic", "X"),
    SubstrateName("2-ketoisovaleric acid", "Kiv", "X"),
    SubstrateName("2-methylserine", "meSer", "X"),
    SubstrateName("2-sulfamoylacetic acid", "Sma", "X"),
    SubstrateName("2R-hydroxy-3-methylpentanoic acid", "ohmePen", "Hmp"),
    SubstrateName("2R-hydroxyisovaleric acid", "2R-Hiv", "D-Hiv"),
    SubstrateName("2S,3S-diaminobutyric acid", "3Dab", "L-Dbu"),
    SubstrateName("2S-amino-8-oxodecanoic acid", "n-oxo-Dec", "X"),
    SubstrateName("2S-amino-9,10-epoxy-8-oxodecanoic acid", "n-epox-oxoDec", "C10:0-NH2(2)-Ep(9)-oxo(8)"),  # noqa: E501
    SubstrateName("2S-aminodecanoic acid", "nDec", "X"),
    SubstrateName("2S-aminododecanoic acid", "nDodec", "X"),
    SubstrateName("2S-aminooctanoic acid", "nOct", "X"),
    SubstrateName("2S-hydroxyisocaproic acid", "2S-Hic", "X"),
    SubstrateName("2S-hydroxyisovaleric acid", "2S-Hiv", "Hiv"),
    SubstrateName("2S-methyl-3-oxobutyrine", "me-oxo-But", "X"),
    SubstrateName("3,4-dehydrolysine", "dhLys", "X"),
    SubstrateName("3,4-dihydroxybenzoic acid", "3,4-dohBza", "X"),
    SubstrateName("3,5-dichloro-4-hydroxybenzoylformic acid", "dclohBfa", "X"),
    SubstrateName("3,5-dichloro-4-hydroxyphenylglycine", "DClHpg", "Cl2-Hpg"),
    SubstrateName("3,5-dihydroxyphenylglycine", "dHpg", "Dhpg"),
    SubstrateName("3-(2-nitrocyclopropylalanine)", "niCpa", "X"),
    SubstrateName("3-(3-pyridyl)-alanine", "PyrAla", "X"),
    SubstrateName("3-amino-2,4-dihydroxybenzoic acid", "2,4-doh-3-nBza", "X"),
    SubstrateName("3-amino-4-hydroxybenzoic acid", "nohBza", "X"),
    SubstrateName("3-amino-6-hydroxy-2-piperidone", "nohPid", "Ahp"),
    SubstrateName("3-aminoisobutyric acid", "Ibu", "X"),
    SubstrateName("3-chlorotyrosine", "clTyr", "Cl-Tyr"),
    SubstrateName("3-hydroxy-4-methylproline", "3-oh-4-mePro", "X"),
    SubstrateName("3-hydroxy-O-methyl-5-methyltyrosine", "ohdmTyr", "X"),
    SubstrateName("3-hydroxy-O-methyltyrosine", "ohmeTyr", "X"),
    SubstrateName("3-hydroxy-para-aminobenzoic acid", "ohnBza", "X"),
    SubstrateName("3-hydroxyasparagine", "ohAsn", "OH-Asn"),
    SubstrateName("3-hydroxyaspartic acid", "ohAsp", "OH-Asp"),
    SubstrateName("3-hydroxyglutamine", "ohGln", "bOH-Gln"),
    SubstrateName("3-hydroxykynurenine", "oh-Kyn", "X"),
    SubstrateName("3-hydroxyleucine", "3-oh-Leu", "3OH-Leu"),
    SubstrateName("3-hydroxypicolinic acid", "Hpa", "Hpa"),
    SubstrateName("3-hydroxyquinaldic acid", "Hqa", "X"),
    SubstrateName("3-hydroxytyrosine", "3-ohTyr", "X"),
    SubstrateName("3-hydroxyvaline", "ohVal", "bOH-Val"),
    SubstrateName("3-methoxyanthranilic acid", "mxAnt", "X"),
    SubstrateName("3-methoxyaspartic acid", "mxAsp", "OMe-Asp"),
    SubstrateName("3-methylasparagine", "meAsn", "bMe-Asn"),
    SubstrateName("3-methylaspartic acid", "meAsp", "bMe-Asp"),
    SubstrateName("3-methylglutamic acid", "meGlu", "3Me-Glu"),
    SubstrateName("3-nitrotyrosine", "nitroTyr", "3NO2-Tyr"),
    SubstrateName("3R-chloroproline", "Cl-Pro", "X"),
    SubstrateName("3R-hydroxy-2,4-diaminobutyric acid", "ohDab", "X"),
    SubstrateName("3R-hydroxyasparagine", "3R-ohAsn", "OH-Asn"),
    SubstrateName("3R-hydroxyaspartic acid", "3R-ohAsp", "OH-Asp"),
    SubstrateName("3R-hydroxyhomotyrosine", "oh-hTyr", "X"),
    SubstrateName("3R-hydroxyleucine", "3R-ohLeu", "3OH-Leu"),
    SubstrateName("3R-methyl-D-aspartic acid branched", "D-3R-meAsp", "bMe-Asp"),
    SubstrateName("3R-methylbeta-alanine", "me-bAla", "X"),
    SubstrateName("3R-methylglutamic acid", "3R-meGlu", "3Me-Glu"),
    SubstrateName("3S,4R-dichloroproline", "dClPro", "Cl2-Pro"),
    SubstrateName("3S,4S-dihydroxyhomotyrosine", "doh-Hty", "X"),
    SubstrateName("3S-aminobutyric acid", "3S-Abu", "X"),
    SubstrateName("3S-cyclohex-2-enylalanine", "cycHexAla", "X"),
    SubstrateName("3S-hydroxy-4S-methylproline", "3S-oh-4S-mePro", "X"),
    SubstrateName("3S-hydroxy-6-chlorohistidine", "HClHis", "X"),
    SubstrateName("3S-hydroxyasparagine", "3S-ohAsn", "OH-Asn"),
    SubstrateName("3S-hydroxyleucine", "3S-ohLeu", "3OH-Leu"),
    SubstrateName("3S-hydroxypipecolic acid", "ohPip", "X"),
    SubstrateName("3S-hydroxyproline", "3S-ohPro", "3OH-Pro"),
    SubstrateName("3S-methyl-D-aspartic acid branched", "D-3S-meAsp", "D-bMe-Asp"),
    SubstrateName("3S-methylaspartic acid", "meAsp", "bMe-Asp"),
    SubstrateName("3S-methylleucine", "meLeu", "X"),
    SubstrateName("3S-methylproline", "3S-mePro", "3Me-Pro"),
    SubstrateName("4,5-dehydroarginine", "4,5-dhArg", "X"),
    SubstrateName("4,5-dihydroxyornithine", "Dho", "X"),
    SubstrateName("4-acetamidopyrrole-2-carboxylic acid", "nacPyr", "X"),
    SubstrateName("4-amino-2-hydroxy-3-isopropoxybenzoic acid", "noh-isopox-Bza", "X"),
    SubstrateName("4-aminobutyric acid", "4-Abu", "X"),
    SubstrateName("4-aminophenylalanine", "nPhe", "X"),
    SubstrateName("4-chlorobenzoic acid", "4ClBza", "X"),
    SubstrateName("4-hydroxy-3-nitrobenzoic acid", "oh-nitroBza", "X"),
    SubstrateName("4-hydroxy-D-kynurenine", "oh-D-Kyn", "X"),
    SubstrateName("4-hydroxybenzoic acid", "4-ohBza", "pOH-Bz"),
    SubstrateName("4-hydroxyglutamine", "4-ohGln", "X"),
    SubstrateName("4-hydroxyphenylglycine", "Hpg", "Hpg"),
    SubstrateName("4-hydroxyphenylpyruvic acid", "oh-ph-Pyr", "X"),
    SubstrateName("4-methoxytryptophan", "mxTrp", "X"),
    SubstrateName("4-methylproline", "4mePro", "4Me-Pro"),
    SubstrateName("4-nitrotryptophan", "nitroTrp", "X"),
    SubstrateName("4R-E-butenyl-4R-methylthreonine", "Bmt", "Bmt"),
    SubstrateName("4R-hydroxyproline", "4R-ohPro", "4OH-Pro"),
    SubstrateName("4R-methylproline", "4R-mePro", "4Me-Pro"),
    SubstrateName("4R-propylproline", "proPro", "X"),
    SubstrateName("4S,5-dihydroxy-2S-aminopentanoic acid", "4oh-Orn", "X"),
    SubstrateName("4S-acetyl-5S-methylproline", "ac-mePro", "X"),
    SubstrateName("4S-hydroxylysine", "ohLys", "X"),
    SubstrateName("4S-methylazetidine-2S-carboxylic acid", "meAzca", "X"),
    SubstrateName("4S-methylproline", "4S-mePro", "4Me-Pro"),
    SubstrateName("4S-propenylproline", "pPro", "X"),
    SubstrateName("5,5-dimethylpipecolic acid", "dmePip", "X"),
    SubstrateName("5-aminolevulinic acid", "nLev", "X"),
    SubstrateName("5-chloroanthranilic acid", "clAnt", "X"),
    SubstrateName("5-chlorotryptophan", "ClTrp", "X"),
    SubstrateName("5-methoxytyrosine", "mxTyr", "X"),
    SubstrateName("5-methylorsellinic acid", "meOrs", "X"),
    SubstrateName("5S-methylproline", "5S-mePro", "5Me-Pro"),
    SubstrateName("6,7-dichlorotryptophan", "dclTrp", "X"),
    SubstrateName("6-chloro-4-hydroxy-1-methyl-indole-3-carboxylic acid", "Cl-oh-me-ind-Car", "X"),
    SubstrateName("6-chloro-4-hydroxyindole-3-carboxylic acid", "ClHic", "X"),
    SubstrateName("6-chlorotryptophan", "ClTrp", "X"),
    SubstrateName("6-hydroxy-tetrahydro-isoquinoline-3-carboxylic acid", "Htia", "X"),
    SubstrateName("6-methylsalicylic acid", "meSal", "X"),
    SubstrateName("6S-methyl-pipecolic acid", "mePip", "X"),
    SubstrateName("An acid hydrazine polyene (intermediate 14)", "HyzPoly", "X"),
    SubstrateName("Compound 4 (formed by the decarboxylative condensation of L-Phe and succinyl-CoA)", "PheSuc", "X"),  # noqa: E501
    SubstrateName("D-alanine", "D-Ala", "D-Ala"),
    SubstrateName("D-aspartic acid branched", "D-Asp", "D-Asp"),
    SubstrateName("D-glutamic acid branched", "D-Glu", "D-Glu"),
    SubstrateName("D-leucine", "D-Leu", "D-Leu"),
    SubstrateName("D-lysergic acid", "D-Lya", "X"),
    SubstrateName("D-phenylalanine", "D-phe", "D-Phe"),
    SubstrateName("D-phenyllactic acid", "D-pheLac", "D-Ph-Lac"),
    SubstrateName("D-pipecolic acid", "D-Pip", "D-Hpr"),
    SubstrateName("N-(1-methyl)-tryptophan", "meTrp", "1Me-Trp"),
    SubstrateName("N-(1-propargyl)-tryptophan", "parTrp", "X"),
    SubstrateName("N-formylglycine", "foGly", "X"),
    SubstrateName("N-hydroxyvaline", "ohVal", "X"),
    SubstrateName("N1-methoxytryptophan", "mxTrp", "OMe-Trp"),
    SubstrateName("N5-acetyl-N5-hydroxyornithine", "acohOrn", "Ac-OH-Orn"),
    SubstrateName("N5-cis-anhydromevalonyl-N5-hydroxyornithine", "Z-ahmohOrn", "X"),
    SubstrateName("N5-formyl-N5-hydroxyornithine", "fohOrn", "Fo-OH-Orn"),
    SubstrateName("N5-hydroxyornithine", "ohOrn", "OH-Orn"),
    SubstrateName("N5-trans-anhydromevalonyl-N5-hydroxyornithine", "E-ahmohOrn", "X"),
    SubstrateName("N6-hydroxylysine", "ohLys", "X"),
    SubstrateName("O-methylthreonine", "meTyr", "X"),
    SubstrateName("O-methyltyrosine", "OmeTyr", "X"),
    SubstrateName("R-3-hydroxy-3-methylproline", "3R-oh-3-mePro", "X"),
    SubstrateName("R-aza-beta-tyrosine", "aza-bTyr", "X"),
    SubstrateName("R-beta-hydroxyphenylalanine", "ohPhe", "X"),
    SubstrateName("R-beta-hydroxytyrosine", "R-ohTyr", "bOH-Tyr"),
    SubstrateName("R-beta-methylphenylalanine", "mePhe", "bMe-Phe"),
    SubstrateName("R-beta-methyltryptophan", "meTrp", "X"),
    SubstrateName("R-beta-phenylalanine", "R-bPhe", "bPhe"),
    SubstrateName("R-beta-tyrosine", "R-bTyr", "bTyr"),
    SubstrateName("S-adenosylmethionine", "Sam", "X"),
    SubstrateName("S-beta-hydroxycyclohex-2S-enylalanine", "Hch", "X"),
    SubstrateName("S-beta-hydroxyenduracididine", "ohEnd", "X"),
    SubstrateName("S-beta-hydroxyphenylalanine", "S-ohPhe", "X"),
    SubstrateName("S-beta-tyrosine", "S-bTyr", "bTyr"),
    SubstrateName("Z-dehydroaminobutyric acid", "dhAbu", "dhAbu"),
    SubstrateName("Z-dehydrotyrosine", "dhTyr", "X"),
    SubstrateName("acetic acid", "Ac", "X"),
    SubstrateName("alanine", "Ala", "Ala"),
    SubstrateName("alaninol", "Aol", "X"),
    SubstrateName("allo-isoleucine", "aIle", "aIle"),
    SubstrateName("allo-threonine", "aThr", "aThr"),
    SubstrateName("anthranilic acid", "Ant", "X"),
    SubstrateName("arginine", "Arg", "Arg"),
    SubstrateName("asparagine", "Asn", "Asn"),
    SubstrateName("aspartic acid", "Asp", "Asp"),
    SubstrateName("aspartic acid branched", "Asp", "Asp"),
    SubstrateName("azetidine-2-carboxylic acid", "cAze", "X"),
    SubstrateName("benzoic acid", "Bza", "Bz"),
    SubstrateName("benzoxazolinate", "Boz", "X"),
    SubstrateName("beta-alanine", "bAla", "bAla"),
    SubstrateName("beta-hydroxy-3-hydroxy-O-methyl-5-methyltyrosine", "dohmTyr", "X"),
    SubstrateName("beta-hydroxyarginine", "ohArg", "X"),
    SubstrateName("beta-hydroxyphenylalanine", "ohPhe", "Ph-Ser"),
    SubstrateName("beta-hydroxytyrosine", "ohTyr", "bOH-Tyr"),
    SubstrateName("beta-lysine", "bLys", "bLys"),
    SubstrateName("betaine", "Bet", "X"),
    SubstrateName("butyric acid", "But", "C4:0"),
    SubstrateName("capreomycidine", "Cap", "Cap"),
    SubstrateName("cinnamic acid", "Cin", "X"),
    SubstrateName("citrulline", "Cit", "Cit"),
    SubstrateName("coumaric acid", "Cou", "X"),
    SubstrateName("cysteic acid", "Cya", "CysA"),
    SubstrateName("cysteine", "Cys", "Cys"),
    SubstrateName("cysteine branched", "Cys", "Cys"),
    SubstrateName("dehydroarginine", "dhArg", "X"),
    SubstrateName("dehydrophenylalanine", "dhPhe", "X"),
    SubstrateName("dehydrotryptophan", "dhTrp", "dh-Trp"),
    SubstrateName("dehydrovaline", "dhVal", "X"),
    SubstrateName("dimethylsulfoniopropionic acid", "dmesulf-Pra", "X"),
    SubstrateName("enduracididine", "End", "End"),
    SubstrateName("fatty acid", "fa", "X"),
    SubstrateName("glutamic acid", "Glu", "Glu"),
    SubstrateName("glutamine", "Gln", "Gln"),
    SubstrateName("glycine", "Gly", "Gly"),
    SubstrateName("glycolic acid", "Glyca", "X"),
    SubstrateName("graminine", "Gram", "X"),
    SubstrateName("guanidinoacetic acid", "Gua", "X"),
    SubstrateName("histidine", "His", "His"),
    SubstrateName("homophenylalanine", "hPhe", "Hph"),
    SubstrateName("homoserine", "hSer", "Hse"),
    SubstrateName("homotyrosine", "hTyr", "Hty"),
    SubstrateName("hydroxyproline", "ohPro", "X"),
    SubstrateName("isoleucine", "Ile", "Ile"),
    SubstrateName("isovaline", "Iva", "Ival"),
    SubstrateName("kynurenine", "Kyn", "Kyn"),
    SubstrateName("lactic acid", "Lac", "Lac"),
    SubstrateName("leucine", "Leu", "Leu"),
    SubstrateName("linoleic acid", "Lin", "X"),
    SubstrateName("lysine", "Lys", "Lys"),
    SubstrateName("malic acid", "MalAc", "X"),
    SubstrateName("meta-tyrosine", "mTyr", "X"),
    SubstrateName("methionine", "Met", "Met"),
    SubstrateName("nicotinic acid", "Nic", "X"),
    SubstrateName("norcoronamic acid", "Nca", "norCMA"),
    SubstrateName("ornithine", "Orn", "Orn"),
    SubstrateName("p-hydroxymandelate", "ohMan", "X"),
    SubstrateName("para-aminobenzoic acid", "nBza", "pNH2-Bz"),
    SubstrateName("pentanoic acid", "Pen", "X"),
    SubstrateName("phenazine-1,6-dicarboxylic acid", "dcPhz", "X"),
    SubstrateName("phenylalanine", "Phe", "Phe"),
    SubstrateName("phenylglycine", "phGly", "Ph-Gly"),
    SubstrateName("phenylpyruvic acid", "Ppa", "X"),
    SubstrateName("pipecolic acid", "Pip", "Hpr"),
    SubstrateName("piperazic acid", "Piz", "X"),
    SubstrateName("proline", "Pro", "Pro"),
    SubstrateName("pyrrole-2-carboxylic acid", "Pca", "X"),
    SubstrateName("pyruvic acid", "Pyv", "X"),
    SubstrateName("quinoxaline-2-carboxylic acid", "Qui", "COOH-Qui"),
    SubstrateName("salicylic acid", "Sal", "X"),
    SubstrateName("serine", "Ser", "Ser"),
    SubstrateName("succinic semialdehyde", "Suc", "X"),
    SubstrateName("succinyl-hydrazinoacetic acid", "SucHyzAc", "X"),
    SubstrateName("threonine", "Thr", "Thr"),
    SubstrateName("trans-2-crotylglycine", "croGly", "X"),
    SubstrateName("tryptophan", "Trp", "Trp"),
    SubstrateName("tyrosine", "Tyr", "Tyr"),
    SubstrateName("valine", "Val", "Val"),
    SubstrateName("valinol", "Valol", "Valol"),
]

_LONG_TO_SUBSTRATE: dict[str, SubstrateName] = {
    sub.long: sub for sub in KNOWN_SUBSTRATES
}

_SHORT_TO_SUBSTRATE: dict[str, SubstrateName] = {
    sub.short: sub for sub in KNOWN_SUBSTRATES
}

_LC_NORINE_TO_SUBSTRATE: dict[str, SubstrateName] = {
    sub.norine.lower(): sub for sub in KNOWN_SUBSTRATES
}


def get_substrate_by_name(name: str) -> SubstrateName:
    """ Look up a substrate by long or short name """
    if name in _LONG_TO_SUBSTRATE:
        return _LONG_TO_SUBSTRATE[name]
    if name in _SHORT_TO_SUBSTRATE:
        return _SHORT_TO_SUBSTRATE[name]
    lc_name = name.lower()
    if lc_name in _LC_NORINE_TO_SUBSTRATE:
        return _LC_NORINE_TO_SUBSTRATE[lc_name]
    raise ValueError(f"Substrate {name} not found")


def is_valid_norine_name(name: str) -> bool:
    """ Check if a name is valid for Norine """
    return name.lower() in _LC_NORINE_TO_SUBSTRATE
