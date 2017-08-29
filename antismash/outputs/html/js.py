# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from glob import glob
import logging
import re
import os

import antismash.common.deprecated as utils
import antismash.common.path as path
from antismash.outputs.html.generate_html_table import generate_html_table

searchgtr_links = {}

def convert_records(seq_records, results, options):
    records = []
    annotations = load_cog_annotations()
    for srec, result in zip(seq_records, results):
        records.append(convert_record(srec, annotations, options, result))
    return records

def convert_record(record, annotations, options, result=None):
    """Convert a SeqRecord to JSON"""
    js_rec = {}
    js_rec['seq_id'] = record.id
    if "extrarecord" in options:
        if record.id in options.extrarecord:
            if "extradata" in options.extrarecord[record.id]:
                if "orig_id" in options.extrarecord[record.id].extradata:
                    js_rec['orig_id'] = options.extrarecord[record.id].extradata["orig_id"]
    if 'orig_id' not in js_rec:
        js_rec['orig_id'] = ""
    js_rec['clusters'] = convert_clusters(record, annotations, options, result)

    return js_rec

def convert_clusters(record, annotations, options, result=None):
    """Convert cluster SeqFeatures to JSON"""
    js_clusters = []
    mibig_results = {}
    if result:
        clusterblast_results = result["modules"]["antismash.modules.clusterblast"]
        if clusterblast_results and clusterblast_results.knowncluster:
            mibig_results = clusterblast_results.knowncluster.mibig_entries
    for cluster in record.get_clusters():
        tta_codons = []
        all_misc_features = record.get_generics()
        for feature in all_misc_features:
            if not cluster.overlaps_with(feature):
                continue
            for note in feature.notes:
                if note.startswith('tta leucine codon'):
                    tta_codons.append(feature)
                    break

        js_cluster = {}
        js_cluster['start'] = int(cluster.location.start) + 1
        js_cluster['end'] = int(cluster.location.end)
        js_cluster['idx'] = cluster.get_cluster_number()
        mibig_entries = mibig_results.get(js_cluster['idx'], {})
        js_cluster['orfs'] = convert_cds_features(record, cluster.cds_children, annotations, options, mibig_entries)
        js_cluster['borders'] = convert_cluster_border_features(cluster.borders)
        js_cluster['tta_codons'] = convert_tta_codons(tta_codons)
        js_cluster['type'] = cluster.get_product_string()
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
            logging.debug('Found closest cluster "%s" for cluster no. %s',
                          js_cluster['knowncluster'],
                          cluster.get_cluster_number())
        js_clusters.append(js_cluster)

    return js_clusters

def convert_cds_features(record, features, annotations, options, mibig_entries):
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
        js_orf['type'] = get_biosynthetic_type(feature, annotations)
        js_orf['description'] = utils.ascii_string(get_description(record, feature, js_orf['type'], options, mibig_entries.get(feature.protein_id, {})))
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
    "Get the description text of a feature"

    replacements = {
        'locus_tag': feature.locus_tag,
        'protein_id': feature.protein_id,
        'smcog': '-',
        'ecnumber': '-',
        'transport_blast_line': '',
        'smcog_tree_line': '',
        'searchgtr_line': '',
        'start': int(feature.location.start) + 1,
        'end': int(feature.location.end),
        'model_details': get_model_details(feature),
        'asf': ''
    }

    blastp_url = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&" \
                 "PROGRAM=blastp&BLAST_PROGRAMS=blastp&QUERY=%s&" \
                 "LINK_LOC=protein&PAGE_TYPE=BlastSearch"
    genomic_context_url = "http://www.ncbi.nlm.nih.gov/projects/sviewer/?" \
                          "Db=gene&DbFrom=protein&Cmd=Link&noslider=1&"\
                          "id=%s&from=%s&to=%s"
    template = '<span class="svgene-tooltip-bold">%(product)s</span><br>\n'
    template += 'Locus-tag: %(locus_tag)s; Protein-ID: %(protein_id)s<br>\n'
    if feature.get_qualifier('EC_number'):
        template += "EC-number(s): %(ecnumber)s<br>\n"
    if options.smcogs:
        template += "smCOG: %(smcog)s<br>\n"
    if options.input_type == 'nucl':
        template += "Location: %(start)s - %(end)s<br><br>\n"
    if feature.sec_met:
        template += '<span class="bold">Signature pHMM hits:</span><br>\n%(model_details)s<br>\n'

    if options.cb_knownclusters and mibig_result:
        cluster_number = feature.cluster.get_cluster_number()
        mibig_homology_file = os.path.join(options.output_dir, "knownclusterblast",
                                         "cluster%d" % cluster_number,
                                         feature.get_accession() + '_mibig_hits.html')
        generate_html_table(mibig_homology_file, mibig_result)
        replacements['mibig_homology_path'] = mibig_homology_file[len(options.output_dir) + 1:]
        template += '<br><a href="%(mibig_homology_path)s" target="_new">MiBIG Hits</a><br>\n'
    template += """
%(transport_blast_line)s
%(searchgtr_line)s
<a href="%(blastp_url)s" target="_new">NCBI BlastP on this gene</a><br>
<a href="%(genomic_context_url)s" target="_new">View genomic context</a><br>
%(smcog_tree_line)s<br>"""
    if get_ASF_predictions(feature):
        template += '<span class="bold">Active Site Finder results:</span><br>\n%(asf)s<br><br>\n'
    template += """AA sequence: <a href="javascript:copyToClipboard('%(sequence)s')">Copy to clipboard</a><br>"""
    template += """Nucleotide sequence: <a href="javascript:copyToClipboard('%(dna_sequence)s')">Copy to clipboard</a><br>"""

    if not options.smcogs:
        del replacements['smcog']
    if options.input_type == 'prot':
        del replacements['start']
        del replacements['end']

    replacements['product'] = feature.product
    if feature.translation:
        sequence = feature.translation
    else:
        sequence = str(utils.get_aa_sequence(feature))
    dna_sequence = feature.extract(record.seq)
    replacements['blastp_url'] = blastp_url % sequence
    replacements['sequence'] = sequence
    replacements['dna_sequence'] = dna_sequence
    if len(sequence) > 2000:
        len_seq = 30
    else:
        len_seq = (len(sequence) / 80) + 1
    replacements['len_seq'] = len_seq
    replacements['genomic_context_url'] = genomic_context_url % \
                    (record.id,
                     max(feature.location.start - 9999, 0),
                     min(feature.location.end + 10000, len(record)))
    ecnumber = feature.get_qualifier('EC_number')
    if ecnumber:
        replacements['ecnumber'] = ", ".join(ecnumber)
    else:
        del replacements['ecnumber']

    if options.smcogs:
        for note in feature.qualifiers.get('note', []):
            if note.startswith('smCOG:') and '(' in note:
                text = note[6:].split('(', 1)[0]
                smcog, desc = text.split(':', 1)
                desc = desc.replace('_', ' ')
                replacements['smcog'] = '%s (%s)' % (smcog, desc)
            elif note.startswith('smCOG tree PNG image:'):
                entry = '<a href="%s" target="_new">View smCOG seed phylogenetic tree with this gene</a>'
                url = note.split(':')[-1]
                replacements['smcog_tree_line'] = entry % url

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
    replacements['asf'] = get_ASF_predictions(feature)
    if replacements['asf'] == "":
        del replacements['asf']

    return template % replacements


def get_biosynthetic_type(feature, annotations):
    "Get the biosythetic type of a CDS feature"
    ann = 'other'
    for note in feature.notes:
        if not note.startswith('smCOG:'):
            continue

        smcog = note[7:].split(':')[0]
        ann = annotations.get(smcog, 'other')

    if not feature.sec_met:
        return ann

    return feature.sec_met.kind

def get_model_details(feature):
    if feature.sec_met:
        return "<br>".join(map(str, feature.sec_met.domains))
    return ""
#    result = ""
#    for note in feature.qualifiers.get('sec_met', []):
#        if not note.startswith('Domains detected'):
#            continue
#        note = note[18:]
#        result += note.replace(';', '<br>')

#    return result

def get_ASF_predictions(feature):
    "check whether predictions from the active site finder module are annotated"
    return "" #TODO
#    ASFsec_met_quals = [sec_met_qual[16:] for sec_met_qual in feature.qualifiers.get('sec_met', [""]) if sec_met_qual.startswith("ASF-prediction")]
#    result = "<br>\n".join(ASFsec_met_quals)

#    return result

def load_cog_annotations():
    "Load the smCOG type annotations from a file"
    type_keys = {
        'B': 'biosynthetic-additional',
        'T': 'transport',
        'R': 'regulatory',
        'O': 'other'
    }
    annotations = {}
    for line in open(path.get_full_path(__file__, 'cog_annotations.txt'), 'r'):
        line = line.strip()
        cog, _, type_ = line.split('\t', 3)
        annotations[cog] = type_keys.get(type_, 'other')

    return annotations
