# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for converting records, features, and domains into
    JSON for use by the webpage javascript
"""

import os
from typing import Any, Dict, Iterable, List, Optional, Tuple
from typing import Set  # comment hints, # pylint: disable=unused-import

from antismash.common import html_renderer, path
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import CDSFeature, Feature, Record, Region, SuperCluster, SubRegion
from antismash.common.secmet import Cluster  # comment hints, # pylint: disable=unused-import
from antismash.common.secmet.features.cdscollection import CDSCollection
from antismash.config import ConfigType
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


def convert_record(record: Record, options: ConfigType, result: Optional[Dict[str, ModuleResults]] = None
                   ) -> Dict[str, Any]:
    """ Convert a Record to JSON """
    if result is None:
        result = {}
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

    assert record.record_index  # shouldn't get here without ensuring this
    for region in record.get_regions():
        tta_codons = fetch_tta_features(region, result)

        js_region = {}  # type: Dict[str, Any]
        js_region['start'] = int(region.location.start) + 1
        js_region['end'] = int(region.location.end)
        js_region['idx'] = region.get_region_number()
        mibig_entries = mibig_results.get(js_region['idx'], {})
        js_region['orfs'] = convert_cds_features(record, region.cds_children, options, mibig_entries)
        js_region['clusters'] = get_clusters_from_region_parts(region.superclusters, region.subregions)
        js_region['ttaCodons'] = convert_tta_codons(tta_codons, record)
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
        mibig_hits = []  # type: List[clusterblast.results.MibigEntry]
        if feature.protein_id:
            mibig_hits = mibig_entries.get(feature.protein_id, [])
        description = get_description(record, feature, gene_function, options, mibig_hits)
        js_orfs.append({"start": feature.location.start + 1,
                        "end": feature.location.end,
                        "strand": feature.strand or 1,
                        "locus_tag": feature.get_name(),
                        "type": gene_function,
                        "description": description})
    return js_orfs


def _find_non_overlapping_cluster_groups(collections: Iterable[CDSCollection],
                                         padding: int = 100) -> Dict[CDSCollection, int]:
    """ Finds a group number for each given collection for which no collection in one
        group overlaps with any other collection in the same group.

        Group numbers start at 0 and the leftmost cluster will be in group 0.
        Assumes that the collections provided are sorted.

        Args:
            collections: the collections to group
            padding: the number of base pairs to have as a minimum gap between
                     collections in a group

        Returns:
            a dictionary mapping each CDSCollection to its group number
    """
    if padding < 0:
        raise ValueError("padding cannot be negative")
    if not collections:
        return {}
    groups = []  # type: List[List[CDSCollection]]
    for collection in collections:
        found_group = False
        for group in groups:
            if collection.location.start > group[-1].location.end + padding:
                group.append(collection)
                found_group = True
                break
        if not found_group:  # then start a new group
            groups.append([collection])

    results = {}
    for group_number, group in enumerate(groups):
        for collection in group:
            results[collection] = group_number
    return results


def get_clusters_from_region_parts(superclusters: Iterable[SuperCluster],
                                   subregions: Iterable[SubRegion]) -> List[Dict[str, Any]]:
    """ Converts all Clusters in a collection of SuperCluster features to JSON """
    unique_clusters = set()  # type: Set[Cluster]
    for supercluster in superclusters:
        unique_clusters.update(supercluster.clusters)
    js_clusters = []
    superclusters = sorted(superclusters, key=lambda x: (x.location.start, -len(x.location)))
    supercluster_groupings = _find_non_overlapping_cluster_groups(superclusters)
    start_index = 0
    for supercluster in superclusters:
        # if it's the only supercluster in the region and it's single, don't draw it to minimise noise
        parent = supercluster.parent
        assert isinstance(parent, Region), type(parent)
        if len(parent.superclusters) == 1 and not parent.subregions and len(supercluster.clusters) == 1:
            continue
        js_cluster = {"start": supercluster.location.start + 1,
                      "end": supercluster.location.end - 1,
                      "tool": "rule-based-clusters",
                      "neighbouring_start": supercluster.location.start,
                      "neighbouring_end": supercluster.location.end,
                      "product": "SC %d: %s" % (supercluster.get_supercluster_number(), supercluster.kind),
                      "isSuperCluster": True}
        js_cluster['height'] = supercluster_groupings[supercluster]
        js_clusters.append(js_cluster)

    start_index += max(supercluster_groupings.values())
    subregions = sorted(subregions, key=lambda x: (x.location.start, -len(x.location), x.tool))
    for i, subregion in enumerate(subregions):
        js_cluster = {"start": subregion.location.start,
                      "end": subregion.location.end,
                      "tool": subregion.tool,
                      "neighbouring_start": subregion.location.start,
                      "neighbouring_end": subregion.location.end,
                      "product": subregion.anchor,
                      "height": start_index + i}
        js_clusters.append(js_cluster)

    start_index += len(subregions) + 2  # allow for label above
    clusters = sorted(unique_clusters, key=lambda x: (x.location.start, -len(x.location), x.product))
    cluster_groupings = _find_non_overlapping_cluster_groups(clusters)
    for cluster in clusters:
        js_cluster = {"start": cluster.core_location.start,
                      "end": cluster.core_location.end,
                      "tool": cluster.tool,
                      "neighbouring_start": cluster.location.start,
                      "neighbouring_end": cluster.location.end,
                      "product": cluster.product,
                      "height": cluster_groupings[cluster] * 2 + start_index,
                      "isSuperCluster": False}
        js_clusters.append(js_cluster)

    return js_clusters


def convert_tta_codons(tta_codons: List[Feature], record: Record) -> List[Dict[str, Any]]:
    """Convert found TTA codon features to JSON"""
    js_codons = []
    for codon in tta_codons:
        cdses = record.get_cds_features_within_location(codon.location, with_overlapping=True)
        js_codons.append({
            'start': codon.location.start + 1,
            'end': codon.location.end,
            'strand': codon.strand if codon.strand is not None else 1,
            'containedBy': [cds.get_name() for cds in cdses]
        })
    return js_codons


def generate_pfam2go_tooltip(record: Record, feature: CDSFeature) -> List[html_renderer.Markup]:
    """Create tooltip text for Pfam to Gene Ontologies results."""
    go_notes = []
    unique_pfams_with_gos = {}
    go_url = 'http://amigo.geneontology.org/amigo/term/'
    go_info_line = "{pf_id}: <a href='{url}{go_id}' target='_blank'> {go_id}:</a> {go_desc}"
    for pfam in record.get_pfam_domains_in_cds(feature):
        if pfam.gene_ontologies:
            pfam_id = pfam.full_identifier
            unique_pfams_with_gos[pfam_id] = pfam.gene_ontologies
    for unique_id, go_qualifier in sorted(unique_pfams_with_gos.items()):
        for go_id, go_description in sorted(go_qualifier.go_entries.items()):
            go_notes.append(go_info_line.format(pf_id=unique_id, url=go_url, go_id=go_id, go_desc=go_description))
    return list(map(html_renderer.Markup, go_notes))


def generate_asf_tooltip_section(record: Record, feature: CDSFeature) -> Dict[Tuple[str, int, int], List[str]]:
    """ Construct tooltip text for activesitefinder annotations """
    asf_notes = {}
    for domain in feature.nrps_pks.domains:
        hits = record.get_domain_by_name(domain.feature_name).asf.hits
        if hits:
            asf_notes[(domain.name, domain.start, domain.end)] = hits
    for pfam in record.get_pfam_domains_in_cds(feature):
        if not pfam.domain:
            continue
        if pfam.asf.hits:
            asf_notes[(pfam.domain, pfam.protein_start, pfam.protein_end)] = pfam.asf.hits
    return asf_notes


def get_description(record: Record, feature: CDSFeature, type_: str,
                    options: ConfigType, mibig_result: List[clusterblast.results.MibigEntry]) -> str:
    "Get the description text of a CDS feature"

    urls = {
        "blastp": ("http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&"
                   "PROGRAM=blastp&BLAST_PROGRAMS=blastp&QUERY=%s&"
                   "LINK_LOC=protein&PAGE_TYPE=BlastSearch") % feature.translation,
        "mibig": "",
        "transport": "",
        "smcog_tree": ""
    }

    genomic_context_url = "http://www.ncbi.nlm.nih.gov/projects/sviewer/?" \
                          "Db=gene&DbFrom=protein&Cmd=Link&noslider=1&"\
                          "id=%s&from=%s&to=%s"

    if mibig_result:
        assert feature.region
        region_number = feature.region.get_region_number()
        mibig_homology_file = os.path.join(options.output_dir, "knownclusterblast",
                                           "region%d" % region_number,
                                           feature.get_accession() + '_mibig_hits.html')
        generate_html_table(mibig_homology_file, mibig_result)
        urls["mibig"] = mibig_homology_file[len(options.output_dir) + 1:]

    if type_ == 'transport':
        urls["transport"] = ("http://blast.jcvi.org/er-blast/index.cgi?project=transporter;"
                             "program=blastp;database=pub/transporter.pep;"
                             "sequence=sequence%%0A%s") % feature.translation

    urls["context"] = genomic_context_url % (record.id,
                                             max(feature.location.start - 9999, 0),
                                             min(feature.location.end + 10000, len(record)))

    if options.smcog_trees:
        for note in feature.notes:  # TODO find a better way to store image urls
            if note.startswith('smCOG tree PNG image:'):
                urls["smcog_tree"] = note.split(':')[-1]
                break

    asf_notes = generate_asf_tooltip_section(record, feature)
    go_notes = generate_pfam2go_tooltip(record, feature)

    urls["searchgtr"] = searchgtr_links.get("{}_{}".format(record.id, feature.get_name()), "")
    template = html_renderer.FileTemplate(path.get_full_path(__file__, "templates", "cds_detail.html"))
    ec_numbers = ""
    ec_number_qual = feature.get_qualifier("EC_number")
    if isinstance(ec_number_qual, list):
        ec_numbers = ",".join(ec_number_qual)
    return template.render(feature=feature, ec_numbers=ec_numbers, go_notes=go_notes,
                           asf_notes=asf_notes, record=record, urls=urls)
