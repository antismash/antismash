# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for converting records, features, and domains into
    JSON for use by the webpage javascript
"""

import string
import os
from typing import Any, Dict, Iterable, List
from typing import Set  # comment hints, # pylint: disable=unused-import

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import CDSFeature, Feature, Record, Region, SuperCluster
from antismash.common.secmet import Cluster  # comment hints, # pylint: disable=unused-import
from antismash.config import ConfigType
from antismash.detection import clusterfinder_probabilistic as clusterfinder
from antismash.modules import clusterblast, tta
from antismash.outputs.html.generate_html_table import generate_html_table

searchgtr_links = {}  # type: Dict[str, str]  # TODO: refactor away from global


def convert_records(records: List[Record], results: List[Dict[str, ModuleResults]],
                    options: ConfigType) -> List[Dict[str, Any]]:
    """ Convert multiple Records to JSON """
    json_records = []
    for record, result in zip(records, results):
        json_records.append(convert_record(record, options, result))
    return json_records


def convert_record(record: Record, options: ConfigType, result: Dict[str, ModuleResults] = None) -> Dict[str, Any]:
    """ Convert a Record to JSON """
    return {'seq_id': record.id,
            'regions': convert_regions(record, options, result)}


def fetch_tta_features(region: Region, result: Dict[str, ModuleResults]) -> List[Feature]:
    """ Returns a list of all TTA features that overlap with the region """
    hits = []  # type: List[Feature]
    tta_results = result.get(tta.__name__)
    if not tta_results:
        return hits

    assert isinstance(tta_results, tta.TTAResults), type(tta_results)
    for feature in tta_results.features:
        if feature.overlaps_with(region):
            hits.append(feature)

    return hits


def convert_regions(record: Record, options: ConfigType, result: Dict[str, ModuleResults]) -> List[Dict[str, Any]]:
    """Convert Region features to JSON"""
    js_regions = []
    mibig_results = {}  # type: Dict[int, Dict[str, List[clusterblast.results.MibigEntry]]]

    clusterblast_results = result.get(clusterblast.__name__)
    if clusterblast_results is not None:
        assert isinstance(clusterblast_results, clusterblast.results.ClusterBlastResults)
        if clusterblast_results.knowncluster:
            mibig_results = clusterblast_results.knowncluster.mibig_entries

    for region in record.get_regions():
        tta_codons = fetch_tta_features(region, result)

        js_region = {}  # type: Dict[str, Any]
        js_region['start'] = int(region.location.start) + 1
        js_region['end'] = int(region.location.end)
        js_region['idx'] = region.get_region_number()
        mibig_entries = mibig_results.get(js_region['idx'], {})
        js_region['orfs'] = convert_cds_features(record, region.cds_children, options, mibig_entries)
        js_region['clusters'] = get_clusters_from_superclusters(region.superclusters)
        js_region['tta_codons'] = convert_tta_codons(tta_codons)
        js_region['type'] = "-".join(region.products)
        js_region['products'] = region.products
        js_region['anchor'] = "r%dc%d" % (record.record_index, region.get_region_number())

        js_regions.append(js_region)

    return js_regions


def convert_cds_features(record: Record, features: Iterable[CDSFeature], options: ConfigType,
                         mibig_entries: Dict[str, List[clusterblast.results.MibigEntry]]
                         ) -> List[Dict[str, Any]]:
    """ Convert CDSFeatures to JSON """
    js_orfs = []
    for feature in features:
        gene_function = str(feature.gene_function)
        js_orfs.append({"start": feature.location.start + 1,
                        "end": feature.location.end,
                        "strand": feature.strand or 1,
                        "locus_tag": feature.get_name(),
                        "type": gene_function,
                        "description": get_description(record, feature, gene_function, options,
                                                       mibig_entries.get(feature.protein_id, []))})
    return js_orfs


def get_clusters_from_superclusters(superclusters: Iterable[SuperCluster]) -> List[Dict[str, Any]]:
    """ Converts all Clusters in a collection of SuperCluster features to JSON """
    unique_clusters = set()  # type: Set[Cluster]
    for supercluster in superclusters:
        unique_clusters.update(supercluster.clusters)
    js_clusters = []
    # clusterfinder's putative regions can never overlap, so if they exist, collapse
    # them into a single row
    putatives = [cluster for cluster in unique_clusters if cluster.product == clusterfinder.PUTATIVE_PRODUCT]
    non_putatives = [cluster for cluster in unique_clusters if cluster.product != clusterfinder.PUTATIVE_PRODUCT]
    clusters = putatives + sorted(non_putatives,
                                  key=lambda x: (x.location.start, -len(x.location), x.product or "unknown"))
    for i, cluster in enumerate(clusters):
        js_cluster = {"start": cluster.core_location.start,
                      "end": cluster.core_location.end,
                      "tool": cluster.tool,
                      "extent": cluster.neighbourhood_range,
                      "product": cluster.product or "unknown"}
        if cluster.product == "cf_putative":
            js_cluster['height'] = 0
        else:
            js_cluster['height'] = i - len(putatives) + 1
        if cluster.tool == "cassis":
            js_cluster['product'] = cluster.get_qualifier("anchor")[0]
        js_clusters.append(js_cluster)
    return js_clusters


def convert_tta_codons(tta_codons: List[Feature]) -> List[Dict[str, Any]]:
    """Convert found TTA codon features to JSON"""
    js_codons = []
    for codon in tta_codons:
        js_codons.append({'start': codon.location.start + 1,
                          'end': codon.location.end,
                          'strand': codon.strand if codon.strand is not None else 1})
    return js_codons


def generate_pfam2go_tooltip(record: Record, feature: CDSFeature) -> List[str]:
    """Create tooltip text for Pfam to Gene Ontologies results."""
    go_notes = []
    unique_pfams_with_gos = {}
    go_url = 'http://amigo.geneontology.org/amigo/term/'
    go_info_line = "{pf_id}: <a href='{url}{go_id}' target='_blank'> {go_id}:</a> {go_desc}"
    for pfam in record.get_pfam_domains_in_cds(feature):
        if pfam.gene_ontologies:
            pfam_ids = ", ".join([xref for xref in pfam.db_xref if xref.startswith('PF')])
            unique_pfams_with_gos[pfam_ids] = pfam.gene_ontologies
    for unique_id, go_qualifier in sorted(unique_pfams_with_gos.items()):
        for go_id, go_description in sorted(go_qualifier.go_entries.items()):
            go_notes.append(go_info_line.format(pf_id=unique_id, url=go_url, go_id=go_id, go_desc=go_description))
    return go_notes


def generate_asf_tooltip_section(record: Record, feature: CDSFeature) -> str:
    """ Construct tooltip text for activesitefinder annotations """
    asf_notes = []
    for domain in feature.nrps_pks.domains:
        for hit in record.get_domain_by_name(domain.feature_name).asf.hits:
            asf_notes.append("%s (%d..%d): %s" % (domain.name, domain.start, domain.end, hit))
    for pfam in record.get_pfam_domains_in_cds(feature):
        for hit in pfam.asf.hits:
            asf_notes.append("%s (%d..%d): %s" % (pfam.domain, pfam.protein_start, pfam.protein_end, hit))
    if not asf_notes:
        return ""
    return '<span class="bold">Active Site Finder results:</span><br>\n%s<br><br>\n' % "<br>".join(asf_notes)


def get_description(record: Record, feature: CDSFeature, type_: str,
                    options: ConfigType, mibig_result: List[clusterblast.results.MibigEntry]) -> str:
    "Get the description text of a CDS feature"

    blastp_url = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&" \
                 "PROGRAM=blastp&BLAST_PROGRAMS=blastp&QUERY=%s&" \
                 "LINK_LOC=protein&PAGE_TYPE=BlastSearch" % feature.translation
    genomic_context_url = "http://www.ncbi.nlm.nih.gov/projects/sviewer/?" \
                          "Db=gene&DbFrom=protein&Cmd=Link&noslider=1&"\
                          "id=%s&from=%s&to=%s"
    template = '<span class="svgene-tooltip-bold">%s</span><br>\n' % feature.product or feature.get_name()
    template += 'Locus-tag: %s; Protein-ID: %s<br>\n' % (feature.locus_tag, feature.protein_id)

    if feature.get_qualifier('EC_number'):
        template += "EC-number(s): %s<br>\n" % ",".join(feature.get_qualifier('EC_number'))

    for gene_function in feature.gene_functions:
        template += "%s<br>\n" % str(gene_function)

    template += "Location: %d - %d<br><br>\n" % (feature.location.start + 1,  # 1-indexed
                                                 feature.location.end)

    if mibig_result:
        region_number = feature.region.get_region_number()
        mibig_homology_file = os.path.join(options.output_dir, "knownclusterblast",
                                           "region%d" % region_number,
                                           feature.get_accession() + '_mibig_hits.html')
        generate_html_table(mibig_homology_file, mibig_result)
        mibig_path = mibig_homology_file[len(options.output_dir) + 1:]
        template += '<br><a href="%s" target="_new">MiBIG Hits</a><br>\n' % mibig_path

    if type_ == 'transport':
        url = "http://blast.jcvi.org/er-blast/index.cgi?project=transporter;" \
              "program=blastp;database=pub/transporter.pep;" \
              "sequence=sequence%%0A%s" % feature.translation
        template += '<a href="%s" target="_new">TransportDB BLAST on this gene<br>' % url

    key = record.id + "_" + feature.get_name()
    if key in searchgtr_links:
        url = searchgtr_links[key]
        template += '<a href="%s" target="_new">SEARCHGTr on this gene<br>\n' % url

    template += '<a href="%s" target="_new">NCBI BlastP on this gene</a><br>\n' % blastp_url

    context = genomic_context_url % (record.id,
                                     max(feature.location.start - 9999, 0),
                                     min(feature.location.end + 10000, len(record)))
    template += """<a href="%s" target="_new">View genomic context</a><br>\n""" % context

    if options.smcog_trees:
        for note in feature.notes:  # TODO find a better way to store image urls
            if note.startswith('smCOG tree PNG image:'):
                url = note.split(':')[-1]
                entry = '<a href="%s" target="_new">View smCOG seed phylogenetic tree with this gene</a><br>\n'
                template += entry % url
                break

    template += generate_asf_tooltip_section(record, feature)

    go_notes = generate_pfam2go_tooltip(record, feature)
    if go_notes:
        template += '<br><span class="bold">Gene Ontology terms for PFAM domains:</span><br>\n' \
                    '%s<br><br>\n' % "<br>".join(go_notes)

    clipboard_fragment = """<a href="javascript:copyToClipboard('%s')">Copy to clipboard</a>"""
    template += "AA sequence: %s<br>\n" % (clipboard_fragment % feature.translation)
    template += "Nucleotide sequence: %s<br>\n" % (clipboard_fragment % feature.extract(record.seq))

    return "".join(char for char in template if char in string.printable)
