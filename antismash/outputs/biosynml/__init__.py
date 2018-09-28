# Copyright (C) 2010-2018 Srikanth Duddela, Daniel Krug, Nestor Zaburannyi
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""BiosynML output format module

"""
import logging
import os
from xml.etree.ElementTree import Element, ElementTree, SubElement, Comment, tostring, parse
from time import gmtime, strftime
from typing import Dict, List

from antismash.common.module_results import ModuleResults
from antismash.common import path
from antismash.config import ConfigType

NAME = "biosynml"
SHORT_DESCRIPTION = "BiosynML output"

def write(results: List[Dict[str, ModuleResults]], options: ConfigType, destination: str) -> None:
    domains = parse(path.get_full_path(__file__, "domainlist.xml"))
    buildingblocks = parse(path.get_full_path(__file__, "buildingblocks.xml"))
    root = generate_template(results)
    genelist = SubElement(root, "genelist")
    domainlist = SubElement(root, "domainlist")
    motiflist = SubElement(root, "motiflist")
    sequencelist = SubElement(root, "Sequencelist")
    populate_sequencelist(root, results, domains, buildingblocks)
    """ Writes the biosynml file into the output directory """
    ElementTree(root).write(os.path.join(options.output_dir, destination), xml_declaration=True)

def generate_template(results):
    root = Element("root")
    header = SubElement(root, "Header")
    SubElement(header, "system").text = "BiosynML"
    SubElement(header, "version").text = "1.1"
    SubElement(header, "date", attrib={"action": "created"}).text = strftime("%d.%m.%Y", gmtime())
    SubElement(header, "date", attrib={"action": "modified"}).text = strftime("%d.%m.%Y", gmtime())
    SubElement(header, "author").text = "antismash"
    SubElement(header, "description").text = "BiosynML output from antismash analysis of " + results.input_file
    return root

def populate_sequencelist(root, results, domains, buildingblocks):
    sequencelist = root.find("Sequencelist")
    for record in results.records:
        source = SubElement(sequencelist, "source", attrib={"id": str(record.record_index + 1)})
        SubElement(source, "label").text = record.id
        for dbxref in record.dbxrefs:
            SubElement(source, "db_xref").text = str(dbxref)
        SubElement(source, "Sequence").text = str(record.seq)
        SubElement(source, "is_circular").text = "1" if record.is_circular() else "0"
        for cluster in record.get_clusters():
            nodelist = add_model(root, cluster)
        populate_genelist(root, record, nodelist, domains, buildingblocks)

def add_model(root, cluster):
    model = SubElement(root, "model", attrib={"id": str(cluster.parent_record.get_cluster_number(cluster))})
    SubElement(model, "title").text = "Biosynthetic pathway" + str(cluster.parent_record.get_cluster_number(cluster))
    SubElement(model, "confidence").text = "Sequence-based prediction"
    SubElement(model, "generator").text = "antiSMASH"
    SubElement(model, "label").text = "model " + str(cluster.parent_record.get_cluster_number(cluster))
    organism = SubElement(model, "organism")
    SubElement(organism, "name").text = cluster.parent_record.description
    SubElement(organism, "strain").text = cluster.parent_record.annotations.get("organism", "")
    SubElement(organism, "identifier", attrib={"source": ""})
    SubElement(organism, "status")
    compound = SubElement(model, "compound")
    SubElement(compound, "name")
    SubElement(compound, "identifier", attrib={"source": ""})
    SubElement(compound, "status")
    SubElement(compound, "citation", attrib={"type": ""})
    genecluster = SubElement(model, "genecluster")
    SubElement(genecluster, "name").text = "Biosynthetic genecluster"
    SubElement(genecluster, "shortname")
    SubElement(genecluster, "status")
    SubElement(genecluster, "type").text = cluster.product
    SubElement(genecluster, "identifier", attrib={"source": ""})
    SubElement(genecluster, "citation", attrib={"type": ""})
    SubElement(genecluster, "sequence", attrib={"source": "sequencelist"}).text = str(cluster.parent_record.record_index + 1)
    region = SubElement(genecluster, "region")
    SubElement(region, "begin").text = str(int(cluster.location.bio_start) + 1)
    SubElement(region, "end").text = str(int(cluster.location.bio_end))
    nodelist = SubElement(model, "nodelist")
    return nodelist

def populate_genelist(root, record, nodelist, domains, buildingblocks):
    genelist = root.find("genelist")
    for cds_feature in record.get_cds_features():
        cds_biopython_feature = cds_feature.to_biopython()[0]
        gene = SubElement(genelist, "gene", attrib={"id": str(len(list(genelist)) + 1)})
        SubElement(gene, "gene_name").text = cds_feature.get_name()
        SubElement(gene, "sequence", attrib={"source": "sequencelist"}).text = str(record.record_index + 1)
        gene_location = SubElement(gene, "gene_location")
        SubElement(gene_location, "begin").text = str(int(cds_biopython_feature.location.bio_start) + (0 if cds_biopython_feature.location.strand == -1 else 1))
        SubElement(gene_location, "end").text = str(int(cds_biopython_feature.location.bio_end) + (1 if cds_biopython_feature.location.strand == -1 else 0))
        gene_qualifiers = SubElement(gene, "gene_qualifiers")
        for qualifier_type in cds_biopython_feature.qualifiers:
            for qualifier_value in cds_biopython_feature.qualifiers[qualifier_type]:
                SubElement(gene_qualifiers, "qualifier", attrib={"name": qualifier_type, "ori": "auto-annotation", "style": "genbank"}).text = qualifier_value
        SubElement(gene_qualifiers, "qualifier", attrib={"name": "gene_location", "ori": "auto-annotation", "style": "genbank"}).text = str(int(cds_biopython_feature.location.bio_start) + (0 if cds_biopython_feature.location.strand == -1 else 1)) + " - " + str(int(cds_biopython_feature.location.bio_end) + (1 if cds_biopython_feature.location.strand == -1 else 0))
        SubElement(gene, "operon").text = str(1)
        populate_domainlist(root, cds_feature, record, nodelist, domains, buildingblocks)
        populate_motiflist(root, cds_feature, record)

def populate_domainlist(root, cds, record, nodelist, domains, buildingblocks):
    domainlist = root.find("domainlist")
    for domain in cds.nrps_pks.domains:
        domain_feature = record.get_domain_by_name(domain.feature_name)
        domain_data = get_domain_data(domain.name, domains)
        domain_biopython_feature = domain_feature.to_biopython()[0]
        domain = SubElement(domainlist, "domain", attrib={"id": str(len(list(domainlist)) + 1)})
        SubElement(domain, "nodeid").text = str(len(list(domainlist)))
        SubElement(domain, "function").text = getattr(domain_data.find("function"), "text", "unspecified")
        SubElement(domain, "subtype").text = domain_feature.domain_subtype
        dstatus = SubElement(domain, "dstatus")
        dstatus.text = "active"
        SubElement(domain, "label").text = getattr(domain_data.find("function"), "text", "unspecified")
        SubElement(domain, "comment")
        SubElement(domain, "chemistry").text = getattr(domain_data.find("chemistry"), "text", "unspecified")
        SubElement(domain, "substrate").text = getattr(domain_data.find("substrate"), "text", "unspecified")
        SubElement(domain, "evidence").text = "Sequence-based prediction"
        domain_qualifiers = SubElement(domain, "domain_qualifiers")
        for qualifier_type in domain_biopython_feature.qualifiers:
            for qualifier_value in domain_biopython_feature.qualifiers[qualifier_type]:
                SubElement(domain_qualifiers, "qualifier", attrib={"name": qualifier_type, "ori": "auto-annotation", "style": "genbank"}).text = qualifier_value
        location = SubElement(domain, "location")
        gene = SubElement(location, "gene")
        SubElement(gene, "geneid", attrib={"source": "genelist"}).text = str(len(list(root.find("genelist"))))
        position = SubElement(gene, "position")
        SubElement(position, "begin").text = str(abs(int(cds.location.bio_start) - int(domain_feature.location.bio_end if domain_feature.location.strand == -1 else domain_feature.location.bio_start)) - (1 if domain_feature.location.strand == -1 else 0))
        SubElement(position, "end").text =  str(abs(int(cds.location.bio_start) - int(domain_feature.location.bio_start if domain_feature.location.strand == -1 else domain_feature.location.bio_end)) - (0 if domain_feature.location.strand == -1 else 1))
        protein = SubElement(location, "protein")
        SubElement(protein, "name").text = cds.get_name()
        SubElement(protein, "sequence").text = domain_feature.translation
        node = SubElement(nodelist, "node", attrib={"id": str(len(list(domainlist)))})
        SubElement(node, "class").text = getattr(domain_data.find("class"), "text", "unspecified")
        SubElement(node, "context").text = getattr(domain_data.find("default_context"), "text", "unspecified")
        for specificity_index, specificity in enumerate(domain_feature.specificity):
            first, last = specificity.split(":")
            prediction = SubElement(domain, "prediction", attrib={"id": str(specificity_index + 1)})
            SubElement(prediction, "method").text = first.strip()
            SubElement(prediction, "readout").text = "prediction"
            SubElement(prediction, "result").text = last.strip()
            if ( first.strip() == "KR activity" ):
                dstatus.text = last.strip()
            if ( first.strip() == "consensus" ):
                buildingblock = SubElement(node, "buildingblock")
                buildingblock_data = get_buildingblock_data(last.strip(), buildingblocks)
                SubElement(buildingblock, "moiety", attrib={"ratio": "1"})
                SubElement(buildingblock, "name").text = getattr(buildingblock_data.find("fullname"), "text", last.strip())
                SubElement(buildingblock, "category").text = getattr(buildingblock_data.find("context"), "text", "unspecified")
                SubElement(buildingblock, "confidence").text = "Consensus Predictions"
                if ( buildingblock_data.find("parent") ):
                    parent = SubElement(buildingblock, "parent").text
                    SubElement(parent, "name").text = getattr(buildingblock_data.find("parent"), "text", "")
                    SubElement(parent, "transform").text = getattr(buildingblock_data.find("transform"), "text", "")

def populate_motiflist(root, cds, record):
    motiflist = root.find("motiflist")
    for motif_feature in cds.motifs:
        motif_biopython_feature = motif_feature.to_biopython()[0]
        motif = SubElement(motiflist, "motif", attrib={"id": str(len(list(motiflist)) + 1), "domainID": "1", "geneID": str(len(list(root.find("genelist"))))})
        SubElement(motif, "motif_name").text = motif_feature.label
        SubElement(motif, "motif_type").text = "BiosynML motif"
        SubElement(motif, "sequence", attrib={"source": "sequencelist"}).text = str(record.record_index + 1)
        motif_location = SubElement(motif, "motif_location")
        SubElement(motif_location, "begin").text = str(abs(int(cds.location.bio_start) - int(motif_biopython_feature.location.bio_end if motif_biopython_feature.location.strand == -1 else motif_biopython_feature.location.bio_start)) - (1 if motif_biopython_feature.location.strand == -1 else 0))
        SubElement(motif_location, "end").text =  str(abs(int(cds.location.bio_start) - int(motif_biopython_feature.location.bio_start if motif_biopython_feature.location.strand == -1 else motif_biopython_feature.location.bio_end)) - (0 if motif_biopython_feature.location.strand == -1 else 1))
        motif_qualifiers = SubElement(motif, "motif_qualifiers")
        for qualifier_type in motif_biopython_feature.qualifiers:
            for qualifier_value in motif_biopython_feature.qualifiers[qualifier_type]:
                SubElement(motif_qualifiers, "qualifier", attrib={"name": qualifier_type, "ori": "auto-annotation", "style": "genbank"}).text = qualifier_value

def get_domain_data(name, domains):
    domain_data = domains.find("domain[qualifier='"+name+"']")
    return domain_data if domain_data else domains

def get_buildingblock_data(name, buildingblocks):
    buildingblock_data = buildingblocks.find("buildingblock[code='"+name+"']")
    return buildingblock_data if buildingblock_data else buildingblocks
