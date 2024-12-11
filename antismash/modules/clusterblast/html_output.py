# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML generation for the clusterblast variants """

from collections import defaultdict
from dataclasses import dataclass, field
import math
from typing import Any, List, Self, Set

from antismash.common import path
from antismash.common.html_renderer import (
    docs_link,
    HTMLSections,
    FileTemplate,
    Markup,
    WILDCARD_TEMPLATE,
)
from antismash.common.json import JSONBase
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer
from antismash.common.secmet import CDSFeature, Record, Region
from antismash.config import get_config

from .data_structures import ReferenceCluster, Protein, Query, Subject
from .results import ClusterBlastResults
from .svg_builder import build_colour_groups

ASDB_URL = (
    "https://antismash-db.secondarymetabolites.org/area"
    f"?record={WILDCARD_TEMPLATE.format('accession')}"
    f"&start={WILDCARD_TEMPLATE.format('real_start')}"
    f"&end={WILDCARD_TEMPLATE.format('real_end')}"
)
MIBIG_URL = f"https://mibig.secondarymetabolites.org/go/{WILDCARD_TEMPLATE.format('accession')}"
TITLE_GENERAL = "Similar known clusters"
TITLE_KNOWN = f"{TITLE_GENERAL} from MIBiG %s"
TITLE_SUB = "Similar subclusters"


def will_handle(_products: List[str], _product_categories: Set[str]) -> bool:
    """ Relevant to every region, so return True for every product """
    return True


def generate_html(region_layer: RegionLayer, results: ClusterBlastResults,
                  record_layer: RecordLayer, options_layer: OptionsLayer
                  ) -> HTMLSections:
    """ Generates the HTML sections for all variants of clusterblast

        Arguments:
            region_layer: the region's HTML wrapper
            results: the clusterblast results to display
            record_layer: the record's HTML wrapper for the region and results
            options_layer: the current global options

        Returns:
            the relevant HTML sections
    """

    html = HTMLSections("clusterblast")

    region = region_layer.region_feature
    base_tooltip = ("Shows %s that are similar to the current region. Genes marked with the "
                    "same colour are interrelated. White genes have no relationship.<br>"
                    "Detailed help and explanations are available "
                    f"{docs_link('here', 'modules/clusterblast/')}.<br><br>"
                    "Click on reference genes to show details of similarities to "
                    "genes within the current region.<br>"
                    "Double click on a reference drawing to reverse the display of the genes.<br>"
                    )
    if results.general:
        tooltip = base_tooltip % "regions from the antiSMASH database"
        tooltip += "<br>Click on a reference name to open that entry in the antiSMASH database (if applicable)."
        references = [ref for ref, _ in results.general.region_results[region.get_region_number() - 1].ranking]
        div = generate_div(region_layer, record_layer, options_layer, "clusterblast",
                           tooltip, references, title=TITLE_GENERAL)
        html.add_detail_section("ClusterBlast", div, "clusterblast")

    if results.knowncluster:
        assert results.knowncluster and results.knowncluster.data_version, "missing MIBiG data version"
        tooltip = base_tooltip % "clusters from the MIBiG database"
        tooltip += "<br>Click on an accession to open that entry in the MIBiG database."
        references = []
        region_results = results.knowncluster.region_results[region.get_region_number() - 1]
        references = [ref for ref, _ in region_results.ranking]
        best_match = region_results.get_best_match()
        if best_match:
            region_layer.most_related_area = best_match
        title = TITLE_KNOWN % results.knowncluster.data_version
        div = generate_div(region_layer, record_layer, options_layer, "knownclusterblast",
                           tooltip, references, title=title)
        html.add_detail_section("KnownClusterBlast", div, "knownclusterblast")

    if results.subcluster:
        tooltip = base_tooltip % "sub-cluster units"
        references = [ref for ref, _ in results.subcluster.region_results[region.get_region_number() - 1].ranking]
        div = generate_div(region_layer, record_layer, options_layer, "subclusterblast",
                           tooltip, references, title=TITLE_SUB)
        html.add_detail_section("SubClusterBlast", div, "subclusterblast")

    return html


def generate_div(region_layer: RegionLayer, record_layer: RecordLayer,
                 options_layer: OptionsLayer, search_type: str,
                 tooltip: str, references: list[ReferenceCluster],
                 **kwargs: Any,
                 ) -> Markup:
    """ Generates the specific HTML section of the body for a given variant of
        clusterblast

        Arguments:
            region_layer: the region's HTML wrapper
            record_layer: the record's HTML wrapper for the region and results
            options_layer: the current global options
            search_type: the name of the variant
            tooltip: the help tooltip for the variant
            references: a list of relevant reference clusters to display

        Returns:
            HTML for the given inputs
    """
    references = references[:options_layer.cb_nclusters]
    template = FileTemplate(path.get_full_path(__file__, "templates", "clusterblast.html"))
    return template.render(record=record_layer, region=region_layer, options=options_layer,
                           tooltip=tooltip, references=references, search_type=search_type, **kwargs)


@dataclass
class GeneMatchJSON(JSONBase):
    """ JSON-friendly structure for the match of a reference gene to the given query gene """
    query: str
    pid: float
    coverage: float

    @classmethod
    def from_subject(cls, query_id: str, subject: Subject) -> "GeneMatchJSON":
        """ Create an instance from the given match.

            Arguments:
                query_id: the identifier of the query CDS
                subject: the match information of the reference gene
        """
        return cls(
            query=query_id,
            pid=math.floor(subject.perc_ident),
            coverage=math.floor(subject.perc_coverage),
        )


@dataclass
class _GeneJSON(JSONBase):
    """ JSON-friendly structure for a single gene """
    locus_tag: str
    start: int
    end: int
    strand: int
    colour: str
    real_start: int
    real_end: int

    def __post_init__(self) -> None:
        assert self.start < self.end, f"{self.start=} < {self.end=}, {self.real_start=}, {self.real_end=}"


@dataclass
class QueryGeneJSON(_GeneJSON):
    """ JSON-friendly structure for a single CDS feature """
    @classmethod
    def from_cds(cls, cds: CDSFeature, colour: str = "white", start: int = 0, end: int = 0) -> Self:
        """ Create an instance from the given CDS feature.

            Arguments:
                cds: the query CDS feature
                colour: an override colour to use for this CDS in the visualisation
        """
        return cls(
            locus_tag=cds.get_name(),
            start=start or cds.start,
            end=end or cds.end,
            strand=cds.location.strand,
            colour=colour,
            real_start=cds.start,
            real_end=cds.end,
        )


@dataclass
class ReferenceGeneJSON(_GeneJSON):
    """ JSON-friendly structure for a single reference gene """
    product: str = ""
    matches: list[GeneMatchJSON] = field(default_factory=list)

    @classmethod
    def from_protein(cls, protein: Protein, *, matches: list[GeneMatchJSON] = None,
                     colour: str = "white",
                     ) -> "ReferenceGeneJSON":
        """ Create an instance from the given reference gene.

            Arguments:
                protein: the reference gene
                matches: a list of each match of this reference to query genes
                colour: an override colour to use for this gene in the visualisation
        """
        if matches is None:
            matches = []
        start, end = list(map(int, protein.location.strip("c").split("-")))
        return cls(
            locus_tag=protein.locus_tag,
            start=protein.draw_start,
            end=protein.draw_end,
            real_start=start,
            real_end=end,
            strand=-1 if protein.strand == "-" else 1,
            colour=colour,
            product=protein.annotations,
            matches=matches,
        )


@dataclass
class ReferenceDataJSON(JSONBase):
    """ JSON-friendly structure for a whole reference area """
    accession: str
    unique_name: str
    label: str
    product: str
    start: int
    end: int
    similarity: int
    genes: list[ReferenceGeneJSON]
    reverse: bool
    real_start: int
    real_end: int

    def __post_init__(self) -> None:
        assert self.start < self.end, f"{self.start=} < {self.end=}, {self.real_start=}, {self.real_end=}"

    @classmethod
    def from_reference(cls, reference: ReferenceCluster, similarity: int, start: int, end: int,
                       genes: list[ReferenceGeneJSON] = None, reverse: bool = False) -> "ReferenceDataJSON":
        """ Create an instance from the given reference area.

            Arguments:
                reference: the reference area
                start:
                end:
                similarity: the similarity score for this reference area and the query area
                genes: a list of reference genes contained by this area
                reverse: whether the area should be reversed in the visualisation
        """
        if genes is None:
            genes = []
        # where possible, abbreviate very long text
        label = reference.description.replace("whole genome shotgun sequence", "WGS")
        return cls(
            accession=reference.accession,
            unique_name=reference.get_name(),
            label=label,
            product=reference.cluster_type,
            start=start,
            end=end,
            similarity=similarity,
            genes=genes,
            reverse=reverse,
            real_start=reference.start,
            real_end=reference.end,
        )


@dataclass
class QueryJSON(JSONBase):
    """ JSON-friendly structure for a whole region/query area """
    start: int
    end: int
    genes: list[QueryGeneJSON] = field(default_factory=list)
    gene_by_name: dict[str, QueryGeneJSON] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.genes and self.gene_by_name:
            assert len(self.genes) == len(self.gene_by_name)
            return
        if self.genes:
            self.gene_by_name.update({gene.locus_tag: gene for gene in self.genes})

    def add_cds(self, cds: CDSFeature) -> None:
        """ Adds a JSON-friendly representation of a CDS to the area.

            Arguments:
                cds: the CDS feature to convert and add
        """
        converted = QueryGeneJSON.from_cds(cds)
        self.genes.append(converted)
        self.gene_by_name[cds.get_name()] = converted

    def __len__(self) -> int:
        return self.end - self.start

    @classmethod
    def from_region(cls, region: Region) -> "QueryJSON":
        """ Creates an instance from the given region.

            Arguments:
                region: the query region
        """
        end = region.end
        genes = []
        if region.crosses_origin():
            origin = region.location.end  # always the end of the record
            end += origin
            for cds in region.cds_children.pre_origin:
                genes.append(QueryGeneJSON.from_cds(cds))
            for cds in region.cds_children.cross_origin:
                genes.append(QueryGeneJSON.from_cds(cds, end=cds.end + origin))
            for cds in region.cds_children.post_origin:
                genes.append(QueryGeneJSON.from_cds(cds, start=cds.start + origin, end=cds.end + origin))
        else:
            genes.extend(QueryGeneJSON.from_cds(cds) for cds in region.cds_children)
        return cls(
            start=region.start,
            end=end,
            genes=genes,
        )


@dataclass
class ClusterblastVariantJSON(JSONBase):
    """ JSON-friendly structure for all results of a particular database (e.g. MIBiG/known) """
    variant_name: str
    url: str = ""
    matches: list[ReferenceDataJSON] = field(default_factory=list)
    query_colours: dict[str, str] = field(default_factory=dict)


@dataclass
class ClusterblastJSON(JSONBase):
    """ JSON-friendly structure for all results of all clusterblast variants """
    query: QueryJSON
    references: list[ClusterblastVariantJSON] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.references:
            return
        for name, url in [
            ("clusterblast", ASDB_URL),
            ("knownclusterblast", MIBIG_URL),
            ("subclusterblast", ""),
        ]:
            self.references.append(ClusterblastVariantJSON(variant_name=name, url=url))


def generate_javascript_data(record: Record, region: Region, results: ClusterBlastResults) -> JSONBase:
    """ Generates JSON data for the javascript to draw relevant results in HTML output

        Arguments:
            record: the relevant Record for the results
            region: the specific Region to generate data for
            results: the ClusterBlastResults that need data extracted

        Returns:
            a JSON-friendly dictionary with the relevant data
    """
    data = ClusterblastJSON(QueryJSON.from_region(region))

    variant_results = [results.general, results.knowncluster, results.subcluster]
    for section_json, section in zip(data.references, variant_results):
        if not section:
            continue
        section_json.query_colours.update({cds.get_name(): "white" for cds in region.cds_children})
        output = section_json.matches
        region_results = section.region_results[region.get_region_number() - 1]
        assert region_results.region is region
        if not region_results.total_hits:
            continue
        ranking = region_results.ranking[:get_config().cb_nclusters]
        colours = build_colour_groups(region.cds_children, ranking)

        for ref, score in ranking:
            assert isinstance(ref, ReferenceCluster)
            pairs_per_ref: dict[str, list[tuple[Query, Subject]]] = defaultdict(list)
            for query, reference in score.scored_pairings:
                pairs_per_ref[reference.name].append((query, reference))
            assert pairs_per_ref
            reference_genes = [region_results.displayed_reference_proteins[tag] for tag in ref.tags]

            start = min(gene.draw_start for gene in reference_genes)
            end = max(gene.draw_end for gene in reference_genes)
            similarity = min(100, int(100 * len(pairs_per_ref) / len(ref.proteins)))
            ref_data = ReferenceDataJSON.from_reference(ref, start=start, end=end, similarity=similarity)
            # force the KCB accessions into their descriptions, as their products aren't unique
            if section is results.knowncluster:
                ref_data.label = f"{ref_data.accession}: {ref_data.label}"
            output.append(ref_data)
            average_strand = 0
            for locus_tag in ref.tags:
                protein: Protein = region_results.displayed_reference_proteins[locus_tag]
                colour = colours.get(locus_tag, "white")
                ref_gene = ReferenceGeneJSON.from_protein(protein, colour=colour)
                ref_data.genes.append(ref_gene)
                if colour == "white":
                    continue
                for query, subject in pairs_per_ref[locus_tag]:
                    average_strand += ref_gene.strand * record.get_cds_by_name(query.id).location.strand
                    record.get_cds_by_name(query.id)
                    ref_gene.matches.append(GeneMatchJSON.from_subject(query.id, subject))
                    assert ref_gene.matches[-1].coverage
                    section_json.query_colours[query.id] = colour
                ref_gene.matches.sort(key=lambda x: (x.pid, x.coverage), reverse=True)
            if average_strand < 0:
                ref_data.reverse = True
    return data
