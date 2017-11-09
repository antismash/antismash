# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import string
import os

from antismash.outputs.html.generate_html_table import generate_html_table

searchgtr_links = {}


def convert_records(seq_records, results, options):
    records = []
    for srec, result in zip(seq_records, results):
        records.append(convert_record(srec, options, result))
    return records


def convert_record(record, options, result=None):
    """Convert a SeqRecord to JSON"""
    js_rec = {'seq_id': record.id,
              'clusters': convert_clusters(record, options, result)}
    # TODO: js_rec["orig_id"] = original id if changed
    return js_rec


def fetch_tta_features(cluster, result):
    """ Returns a list of all TTA features that overlap with the cluster """
    hits = []
    tta_results = result.get("antismash.modules.tta")
    if not tta_results:
        return hits

    for feature in tta_results.features:
        if feature.overlaps_with(cluster):
            hits.append(feature)

    return hits


def convert_clusters(record, options, result):
    """Convert cluster SeqFeatures to JSON"""
    js_clusters = []
    mibig_results = {}

    clusterblast_results = result.get("antismash.modules.clusterblast")
    if clusterblast_results and clusterblast_results.knowncluster:
        mibig_results = clusterblast_results.knowncluster.mibig_entries

    for cluster in record.get_clusters():
        tta_codons = fetch_tta_features(cluster, result)

        js_cluster = {}
        js_cluster['start'] = int(cluster.location.start) + 1
        js_cluster['end'] = int(cluster.location.end)
        js_cluster['idx'] = cluster.get_cluster_number()
        mibig_entries = mibig_results.get(js_cluster['idx'], {})
        js_cluster['orfs'] = convert_cds_features(record, cluster.cds_children, options, mibig_entries)
        js_cluster['borders'] = convert_cluster_border_features(cluster.borders)
        js_cluster['tta_codons'] = convert_tta_codons(tta_codons)
        js_cluster['type'] = cluster.get_product_string()
        js_cluster['products'] = cluster.products
        if cluster.probability is not None:
            js_cluster['probability'] = cluster.probability
        if options.input_type == 'prot':
            js_cluster['unordered'] = True
        js_cluster['knowncluster'] = "-"
        js_cluster['BGCid'] = "-"

        if cluster.knownclusterblast:
            bestcluster = cluster.knownclusterblast[0]
            js_cluster['knowncluster'] = bestcluster[0]
            js_cluster['BGCid'] = bestcluster[1]
        js_clusters.append(js_cluster)

    return js_clusters


def convert_cds_features(record, features, options, mibig_entries):
    """Convert CDS SeqFeatures to JSON"""
    js_orfs = []
    for feature in features:
        js_orf = {}
        js_orf['start'] = int(feature.location.start) + 1
        js_orf['end'] = int(feature.location.end)
        # Fix for files that have their coordinates the wrong way around
        if js_orf['start'] > js_orf['end']:
            js_orf['end'], js_orf['start'] = js_orf['start'], js_orf['end']
        js_orf['strand'] = feature.strand if feature.strand is not None else 1
        js_orf['locus_tag'] = feature.get_name()
        js_orf['type'] = str(feature.gene_function)
        js_orf['description'] = get_description(record, feature, js_orf['type'], options, mibig_entries.get(feature.protein_id, {}))
        js_orfs.append(js_orf)
    return js_orfs


def convert_cluster_border_features(borders):
    js_borders = []
    for border in borders:
        border_note = border.qualifiers['note'][0]
        if not border_note.startswith('best prediction'):
            continue
        js_border = {}
        js_border['start'] = int(border.location.start) + 1
        js_border['end'] = int(border.location.end)
        # Always order the coordinates correctly
        if js_border['start'] > js_border['end']:
            js_border['end'], js_border['start'] = js_border['start'], js_border['end']
        js_border['tool'] = border.qualifiers['tool']
        js_borders.append(js_border)
    return js_borders


def convert_tta_codons(tta_codons):
    """Convert found TTA codon features to JSON"""
    js_codons = []
    for codon in tta_codons:
        js_codon = {}
        js_codon['start'] = int(codon.location.start) + 1
        js_codon['end'] = int(codon.location.end)
        js_codon['strand'] = codon.strand if codon.strand is not None else 1
        js_codons.append(js_codon)

    return js_codons


def get_description(record, feature, type_, options, mibig_result):
    "Get the description text of a CDS feature"

    replacements = {
        'locus_tag': feature.locus_tag,
        'protein_id': feature.protein_id,
        'transport_blast_line': '',
        'smcog_tree_line': '',
        'searchgtr_line': '',
        'model_details': get_model_details(feature),
    }

    smcogs = not options.minimal or options.smcogs_enabled or options.smcogs_trees  # TODO make simpler in args

    blastp_url = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&" \
                 "PROGRAM=blastp&BLAST_PROGRAMS=blastp&QUERY=%s&" \
                 "LINK_LOC=protein&PAGE_TYPE=BlastSearch"
    genomic_context_url = "http://www.ncbi.nlm.nih.gov/projects/sviewer/?" \
                          "Db=gene&DbFrom=protein&Cmd=Link&noslider=1&"\
                          "id=%s&from=%s&to=%s"
    template = '<span class="svgene-tooltip-bold">{product}</span><br>\n'
    template += 'Locus-tag: {locus_tag}; Protein-ID: {protein_id}<br>\n'
    if feature.get_qualifier('EC_number'):
        template += "EC-number(s): {ecnumber}<br>\n"
    for gene_function in feature.gene_functions:
        if gene_function.tool != "cluster_definition":
            template += "%s<br>\n" % str(gene_function)
    if smcogs:
        for note in feature.notes:  # TODO find a better way to store image urls
            if note.startswith('smCOG tree PNG image:'):
                entry = '<a href="%s" target="_new">View smCOG seed phylogenetic tree with this gene</a>'
                url = note.split(':')[-1]
                replacements['smcog_tree_line'] = entry % url
    if options.input_type == 'nucl':
        replacements["start"] = int(feature.location.start) + 1  # 1-indexed
        replacements["end"] = int(feature.location.end)
        template += "Location: {start} - {end}<br><br>\n"

    if feature.sec_met:
        template += '<span class="bold">Signature pHMM hits:</span><br>\n{model_details}<br>\n'

    if mibig_result:
        cluster_number = feature.cluster.get_cluster_number()
        mibig_homology_file = os.path.join(options.output_dir, "knownclusterblast",
                                         "cluster%d" % cluster_number,
                                         feature.get_accession() + '_mibig_hits.html')
        generate_html_table(mibig_homology_file, mibig_result)
        replacements['mibig_homology_path'] = mibig_homology_file[len(options.output_dir) + 1:]
        template += '<br><a href="{mibig_homology_path}" target="_new">MiBIG Hits</a><br>\n'
    template += """
{transport_blast_line}
{searchgtr_line}
<a href="{blastp_url}" target="_new">NCBI BlastP on this gene</a><br>
<a href="{genomic_context_url}" target="_new">View genomic context</a><br>
{smcog_tree_line}<br>"""
    if get_ASF_predictions(feature):
        replacements['asf'] = get_ASF_predictions(feature)
        template += '<span class="bold">Active Site Finder results:</span><br>\n{asf}<br><br>\n'
    template += """AA sequence: <a href="javascript:copyToClipboard('{sequence}')">Copy to clipboard</a><br>"""
    template += """Nucleotide sequence: <a href="javascript:copyToClipboard('{dna_sequence}')">Copy to clipboard</a><br>"""

    replacements['product'] = feature.product
    sequence = feature.translation
    dna_sequence = feature.extract(record.seq)
    replacements['blastp_url'] = blastp_url % sequence
    replacements['sequence'] = sequence
    replacements['dna_sequence'] = dna_sequence
    replacements['genomic_context_url'] = genomic_context_url % (record.id,
                                 max(feature.location.start - 9999, 0),
                                 min(feature.location.end + 10000, len(record)))
    ecnumber = feature.get_qualifier('EC_number')
    if ecnumber:
        replacements['ecnumber'] = ", ".join(ecnumber)

    if type_ == 'transport':
        url = "http://blast.jcvi.org/er-blast/index.cgi?project=transporter;" \
              "program=blastp;database=pub/transporter.pep;" \
              "sequence=sequence%%0A%s" % sequence
        transport_blast_line = '<a href="%s" target="_new">TransportDB BLAST on this gene<br>' % url
        replacements['transport_blast_line'] = transport_blast_line

    key = record.id + "_" + feature.get_name()
    if key in searchgtr_links:
        url = searchgtr_links[key]
        searchgtr_line = '<a href="%s" target="_new">SEARCHGTr on this gene<br>' % url
        replacements['searchgtr_line'] = searchgtr_line

    completed = template.format(**replacements)
    return "".join(char for char in completed if char in string.printable)


def get_model_details(feature):
    if feature.sec_met:
        return "<br>".join(map(str, feature.sec_met.domains))
    return ""


def get_ASF_predictions(feature):
    "check whether predictions from the active site finder module are annotated"
    if not hasattr(get_ASF_predictions, "logged"):
        logging.critical("ASF_predictions being skipped in js.py")
        get_ASF_predictions.logged = True
    return ""  # TODO
#    ASFsec_met_quals = [sec_met_qual[16:] for sec_met_qual in feature.qualifiers.get('sec_met', [""]) if sec_met_qual.startswith("ASF-prediction")]
#    result = "<br>\n".join(ASFsec_met_quals)

#    return result
