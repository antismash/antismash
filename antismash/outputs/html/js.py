# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for converting records, features, and domains into
    JSON for use by the webpage javascript
"""

import os
from typing import Any, Dict, Iterable, List, Optional, Tuple

from antismash.common import html_renderer, path
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import CDSFeature, Feature, Record, Region, Source
from antismash.common.secmet.locations import FeatureLocation, location_bridges_origin
from antismash.common.secmet.qualifiers.gene_functions import GeneFunction
from antismash.common.secmet.qualifiers.go import GOQualifier
from antismash.common.secmet.features.cdscollection import CDSCollection
from antismash.common.secmet.features.protocluster import SideloadedProtocluster
from antismash.common.secmet.features.subregion import SideloadedSubRegion
from antismash.common.secmet.locations import locations_overlap, split_origin_bridging_location
from antismash.config import ConfigType
from antismash.detection.tigrfam.tigr_domain import TIGRDomain
from antismash.modules import clusterblast, smcog_trees, tfbs_finder as tfbs, tta
from antismash.outputs.html.generate_html_table import generate_html_table

GO_URL = 'http://amigo.geneontology.org/amigo/term/'


def get_region_css(region: Region) -> str:
    """ Collects the CSS classes for the given region and returns it as a
        single string.

        Arguments:
            region: the region to generate CSS for

        Returns:
            a string of the CSS class(es)
    """
    if len(region.product_categories) > 1:
        return "hybrid"
    if len(region.get_unique_protoclusters()) < 1:
        return "unknown"
    classes = [list(region.product_categories)[0]]
    if len(region.get_unique_protoclusters()) == 1:
        classes.append(region.get_unique_protoclusters()[0].product)
    return " ".join(classes)


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
    return {
        'length': len(record.seq),
        'seq_id': record.id,
        'regions': convert_regions(record, options, result)
    }


def fetch_tta_features(region: Region, result: Dict[str, ModuleResults]) -> List[Feature]:
    """ Returns a list of all TTA features that overlap with the region """
    hits: List[Feature] = []
    tta_results = result.get(tta.__name__)
    if not tta_results:
        return hits

    assert isinstance(tta_results, tta.TTAResults), type(tta_results)
    for feature in tta_results.features:
        if feature.overlaps_with(region):
            hits.append(feature)

    return hits


def convert_source(source: Source, region: Region, name: str = None) -> Dict[str, Any]:
    """ Converts a Source feature to JSON """
    result = {
        "regionStart": max(0, source.location.start - region.location.start),
        "regionEnd": min(region.location.end, source.location.end - region.location.start),
        "recordStart": source.location.start,
        "recordEnd": source.location.end,
    }
    if name:
        result["name"] = name
    return result


def convert_regions(record: Record, options: ConfigType, result: Dict[str, ModuleResults]) -> List[Dict[str, Any]]:
    """Convert Region features to JSON"""
    js_regions = []
    mibig_results: Dict[int, Dict[str, List[clusterblast.results.MibigEntry]]] = {}

    clusterblast_results = result.get(clusterblast.__name__)
    if clusterblast_results is not None:
        assert isinstance(clusterblast_results, clusterblast.results.ClusterBlastResults)
        if clusterblast_results.knowncluster:
            mibig_results = clusterblast_results.knowncluster.mibig_entries

    assert record.record_index  # shouldn't get here without ensuring this
    for region in record.get_regions():
        tta_codons = fetch_tta_features(region, result)

        js_region: Dict[str, Any] = {}
        js_region['start'] = int(region.location.start) + 1
        js_region['end'] = int(region.location.end)
        js_region['idx'] = region.get_region_number()
        mibig_entries = mibig_results.get(js_region['idx'], {})
        if region.crosses_origin():
            last_part = region.location.parts[-1]
            assert last_part.start == 0, region.location.parts
            js_region["start"] = region.start
            js_region["end"] = len(record.seq) + last_part.end
            start = len(record.seq)
            end = len(record.seq) + len(last_part)
            js_region["sources"] = [
                convert_source(Source(region.location.parts[0]), region, record.id),
                convert_source(Source(FeatureLocation(start, end)), region, record.id),
            ]
            js_region["sources"][-1]["regionEnd"] = len(record) + len(last_part)

            upper = []
            lower = []
            for cds in region.cds_children:
                if cds.is_contained_by(region.location.parts[0]):
                    upper.append(cds)
                else:
                    lower.append(cds)
            orfs = convert_cds_features(record, upper, options, mibig_entries)
            orfs.extend(convert_cds_features(record, lower, options, mibig_entries, offset=len(record)))
            assert len(orfs) == len(region.cds_children), f"{len(orfs)} != {len(region.cds_children)}"
            js_region["orfs"] = orfs
        else:
            js_region['orfs'] = convert_cds_features(record, region.cds_children, options, mibig_entries)
        js_region['clusters'] = get_clusters_from_region(region, len(record.seq), circular=record.is_circular())
        sites = {
            "ttaCodons": convert_tta_codons(tta_codons, record),
            "bindingSites": convert_binding_sites(region, result),
        }
        js_region['sites'] = sites
        js_region['type'] = region.get_product_string()
        js_region['products'] = region.products
        js_region['product_categories'] = list(region.product_categories)
        js_region['cssClass'] = get_region_css(region)
        js_region['anchor'] = f"r{record.record_index}c{region.get_region_number()}"
        if record.has_multiple_sources():
            sources = [source for source in record.get_sources() if source.overlaps_with(region)]
            if len(sources) > 1:
                js_region['sources'] = [convert_source(source, region) for source in sources]

        js_regions.append(js_region)

    return js_regions


def convert_cds_features(record: Record, features: Iterable[CDSFeature], options: ConfigType,
                         mibig_entries: Dict[str, List[clusterblast.results.MibigEntry]], offset: int = 0,
                         ) -> List[Dict[str, Any]]:
    """ Convert CDSFeatures to JSON """
    js_orfs = []
    for feature in features:
        gene_function = feature.gene_function
        # resistance genes have special markers, not just a colouring, so revert to OTHER
        if gene_function == GeneFunction.RESISTANCE:
            gene_function = GeneFunction.OTHER
        mibig_hits: List[clusterblast.results.MibigEntry] = []
        mibig_hits = mibig_entries.get(feature.get_name(), [])
        description = get_description(record, feature, str(gene_function), options, mibig_hits)
        start = feature.location.start + 1 + offset
        end = feature.location.end + offset
        if location_bridges_origin(feature.location):
            lower, upper = split_origin_bridging_location(feature.location)
            start = min(part.start for part in upper)
            end = max(part.end for part in lower) + len(record)
        js_orfs.append({
            "start": start,
            "end": end,
            "strand": feature.strand or 1,
            "locus_tag": feature.get_name(),
            "type": str(gene_function),
            "description": description,
            "dna": str(feature.location.extract(record.seq)),
            "translation": feature.translation,
            "product": feature.product or "",
        })
        if feature.gene_functions.get_by_tool("resist"):  # don't add to every gene for size reasons
            js_orfs[-1]["resistance"] = True
    return js_orfs


def _find_non_overlapping_cluster_groups(collections: Iterable[CDSCollection],
                                         padding: int = 100, record_length: int = None) -> Dict[CDSCollection, int]:
    """ Finds a group number for each given collection for which no collection in one
        group overlaps with any other collection in the same group.

        Group numbers start at 0 and the leftmost cluster will be in group 0.
        Assumes that the collections provided are sorted.

        Args:
            collections: the collections to group
            padding: the number of base pairs to have as a minimum gap between
                     collections in a group
            record_length: the length of the record for collection sets that
                           contain an origin-crossing

        Returns:
            a dictionary mapping each CDSCollection to its group number
    """
    collections = sorted(collections)
    if padding < 0:
        raise ValueError("padding cannot be negative")
    if not collections:
        return {}
    any_cross_origin = any(collection.crosses_origin() for collection in collections)
    if any_cross_origin and not record_length:
        raise ValueError("Cannot group cross-origin collections without record length")

    groups: List[List[CDSCollection]] = [[collections[0]]]
    for collection in collections[1:]:
        found_group = False
        for group in groups:
            if group[-1].overlaps_with(collection):
                continue
            end = group[-1].end + padding
            start = collection.start
            if group[-1].crosses_origin() and not collection.crosses_origin():
                assert any_cross_origin and record_length
                start += record_length
            if start > end:
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


def get_clusters_from_region(region: Region, record_length: int, circular: bool = False) -> List[Dict[str, Any]]:
    """ Converts all Protoclusters in a collection of CandidateCluster features to JSON

        Arguments:
            region: the region with protoclusters and/or subregions to convert
            record_length: the length of the record
            circular: whether or not the record is circular

        Returns:
            a list of dictionaries, each representing an area in the region
    """
    js_clusters = []

    def reduction(collection: CDSCollection) -> tuple[int, int]:
        if region.crosses_origin():
            if collection.start < record_length / 2:
                return (collection.start + record_length, -len(collection.location))
        return (collection.start, -len(collection.location))

    candidate_clusters = sorted(region.candidate_clusters, key=reduction)

    candidate_cluster_groupings = _find_non_overlapping_cluster_groups(candidate_clusters, record_length=record_length)
    start_index = 0
    for candidate_cluster in candidate_clusters:
        # if it's the only candidate_cluster in the region and it's single, don't draw it to minimise noise
        parent = candidate_cluster.parent
        assert isinstance(parent, Region), type(parent)
        if len(parent.candidate_clusters) == 1 and not parent.subregions and len(candidate_cluster.protoclusters) == 1:
            continue
        # adjust locations to account for cross-origin segments
        start = candidate_cluster.start
        end = candidate_cluster.end
        if region.crosses_origin():
            offset = record_length
            if locations_overlap(candidate_cluster.location.parts[-1], region.location.parts[-1]):
                end += offset
            if locations_overlap(candidate_cluster.location.parts[0], region.location.parts[-1]):
                start += offset
        js_cluster = {
            "start": start + 1,
            "end": end - 1,
            "tool": "",
            "neighbouring_start": start,
            "neighbouring_end": end,
            "product": f"CC {candidate_cluster.get_candidate_cluster_number()}: {candidate_cluster.kind}",
            "kind": "candidatecluster",
            "prefix": ""
        }
        js_cluster['height'] = candidate_cluster_groupings[candidate_cluster]
        js_clusters.append(js_cluster)

    if candidate_cluster_groupings:
        start_index += max(candidate_cluster_groupings.values())

    for subregion in sorted(region.subregions, key=lambda x: (reduction(x), x.tool)):
        start_index += 1
        prefix = ""
        tool = ""
        if isinstance(subregion, SideloadedSubRegion):
            prefix = subregion.tool + (":" if subregion.label else "")
        else:
            tool = subregion.tool
        start = subregion.start
        end = subregion.end
        if subregion.crosses_origin():
            end += record_length
        js_cluster = {"start": start,
                      "end": end,
                      "tool": tool,
                      "neighbouring_start": start,
                      "neighbouring_end": end,
                      "product": subregion.label,
                      "height": start_index,
                      "prefix": prefix,
                      "kind": "subregion"}
        js_clusters.append(js_cluster)

    start_index += 2  # allow for label above
    clusters = region.get_unique_protoclusters()
    cluster_groupings = _find_non_overlapping_cluster_groups(clusters, record_length=record_length)
    for cluster in clusters:
        prefix = ""
        if isinstance(cluster, SideloadedProtocluster):
            prefix = f"{cluster.tool}:"
        start = cluster.core_location.parts[0].start
        end = cluster.core_location.parts[-1].end
        if region.crosses_origin():
            offset = record_length
            if locations_overlap(cluster.core_location.parts[-1], region.location.parts[-1]):
                end += offset
            if locations_overlap(cluster.core_location.parts[0], region.location.parts[-1]):
                start += offset
        neighbour_start = start - cluster.neighbourhood_range
        neighbour_end = end + cluster.neighbourhood_range
        # if it's not a circular record, the above is too generous, so bound them
        if not circular:
            neighbour_start = max(0, neighbour_start)
            neighbour_end = min(record_length, neighbour_end)
        js_cluster = {
            "start": start,
            "end": end,
            "tool": cluster.tool,
            "neighbouring_start": neighbour_start,
            "neighbouring_end": neighbour_end,
            "product": cluster.product,
            "category": cluster.product_category,
            "height": cluster_groupings[cluster] * 2 + start_index,
            "kind": "protocluster",
            "prefix": prefix
        }
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


def convert_binding_sites(region: Region, results: Dict[str, ModuleResults]) -> List[Dict[str, Any]]:
    """ Constructs required JSON-friendly information for the javascript marking
        binding sites in the gene overview.

        Arguments:
            region: the relevant region to generate data for
            results: a dictionary of all module results for the region's record

        Returns:
            a list of dictionaries, one for each marker in the region
    """
    sites: List[Dict[str, Any]] = []
    tfbs_results = results.get(tfbs.__name__)
    if not tfbs_results:
        return sites
    assert isinstance(tfbs_results, tfbs.TFBSFinderResults)
    confidence = tfbs.tfbs_finder.Confidence.MEDIUM
    for hits in tfbs_results.get_hits_by_region(region.get_region_number(),
                                                confidence=confidence, allow_better=True):
        sites.append({
            "loc": hits.start,
            "len": len(hits.consensus),
        })
    return sites


def build_pfam2go_links(go_qualifier: Optional[GOQualifier], prefix: str = "") -> List[str]:
    """ A helper for generating Pfam2GO HTML fragments with links and descriptions

        Arguments:
            go_qualifier: the GOQualifier to use for building the links
            prefix: an optional string to prefix the link with

        Returns:
            a list of strings, each being an HTML formatted link prefixed by
            the given prefix and followed by the description of the GO term

    """
    if go_qualifier is None:  # a pfam may have no matching GO terms
        return []
    template = "{prefix}<a class='external-link' href='{url}{go_id}' target='_blank'>{go_id}</a>: {desc}"
    return [template.format(prefix=prefix, url=GO_URL, go_id=go_id, desc=desc)
            for go_id, desc in go_qualifier.go_entries.items()]


def generate_pfam2go_tooltip(record: Record, feature: CDSFeature) -> List[html_renderer.Markup]:
    """Create tooltip text for Pfam to Gene Ontologies results."""
    go_notes = []
    unique_pfams_with_gos = {}
    for pfam in record.get_pfam_domains_in_cds(feature):
        if pfam.gene_ontologies:
            pfam_id = pfam.full_identifier
            unique_pfams_with_gos[pfam_id] = pfam.gene_ontologies
    for unique_id, go_qualifier in sorted(unique_pfams_with_gos.items()):
        go_notes.extend(build_pfam2go_links(go_qualifier, prefix=f"{unique_id}: "))
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
            asf_notes[(pfam.domain, pfam.protein_location.start, pfam.protein_location.end)] = pfam.asf.hits
    return asf_notes


def generate_pfam_tooltip(record: Record, feature: CDSFeature) -> List[str]:
    """ Construct tooltip text for PFAM annotations """
    pfam_notes = []
    for pfam in record.get_pfam_domains_in_cds(feature):
        pfam_notes.append(f"{pfam.full_identifier} ({pfam.description}): {pfam.protein_location}"
                          f"(score: {pfam.score}, e-value: {pfam.evalue})")
    return pfam_notes


def generate_tigr_tooltip(record: Record, feature: CDSFeature) -> List[str]:
    """ Construct tooltip text for TIGRFam annotations """
    tigr_notes = []
    for tigr in record.get_antismash_domains_in_cds(feature):
        if not isinstance(tigr, TIGRDomain):
            continue
        tigr_notes.append(f"{tigr.identifier} ({tigr.description}): {tigr.protein_location}"
                          f"(score: {tigr.score}, e-value: {tigr.evalue})")
    return tigr_notes


def get_description(record: Record, feature: CDSFeature, type_: str,
                    options: ConfigType, mibig_result: List[clusterblast.results.MibigEntry]) -> str:
    "Get the description text of a CDS feature"

    urls = {
        "blastp": ("http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&"
                   f"PROGRAM=blastp&BLAST_PROGRAMS=blastp&QUERY={feature.translation}&"
                   "LINK_LOC=protein&PAGE_TYPE=BlastSearch"),
        "mibig": "",
        "transport": "",
        "smcog_tree": "",
        "context": (
            "http://www.ncbi.nlm.nih.gov/projects/sviewer/"
            "?Db=gene&DbFrom=protein&Cmd=Link&noslider=1"
            f"&id={record.id}"
            f"&from={max(feature.location.start - 9999, 0)}"  # TODO handle circularity
            f"&to={min(feature.location.end + 10000, len(record))}"  # TODO handle circularity
        ),
    }

    if mibig_result:
        assert feature.region
        region_number = feature.region.get_region_number()
        mibig_homology_file = os.path.join(options.output_dir, "knownclusterblast",
                                           f"region{region_number}",
                                           f"{feature.get_accession()}_mibig_hits.html")
        generate_html_table(mibig_homology_file, mibig_result)
        urls["mibig"] = mibig_homology_file[len(options.output_dir) + 1:]

    if type_ == 'transport':
        urls["transport"] = ("http://blast.jcvi.org/er-blast/index.cgi?project=transporter;"
                             "program=blastp;database=pub/transporter.pep;"
                             f"sequence=sequence%%0A{feature.translation}")

    if smcog_trees in options.all_enabled_modules:
        for note in feature.notes:  # TODO find a better way to store image urls
            if note.startswith('smCOG tree PNG image:'):
                urls["smcog_tree"] = note.split(':')[-1]
                break

    asf_notes = generate_asf_tooltip_section(record, feature)
    go_notes = generate_pfam2go_tooltip(record, feature)
    pfam_notes = generate_pfam_tooltip(record, feature)
    tigr_notes = generate_tigr_tooltip(record, feature)

    template = html_renderer.FileTemplate(path.get_full_path(__file__, "templates", "cds_detail.html"))
    ec_numbers = ""
    ec_number_qual = feature.get_qualifier("EC_number")
    if isinstance(ec_number_qual, list):
        ec_numbers = ",".join(ec_number_qual)
    return template.render(feature=feature, ec_numbers=ec_numbers, go_notes=go_notes,
                           asf_notes=asf_notes, pfam_notes=pfam_notes, tigr_notes=tigr_notes,
                           record=record, urls=urls, add_ncbi_context=options.html_ncbi_context)
