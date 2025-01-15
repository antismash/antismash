# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Generates HTML output for the terpene module """

from antismash.common import path
from antismash.common.html_renderer import FileTemplate, HTMLSections, Markup
from antismash.common.layers import OptionsLayer, RegionLayer, RecordLayer
from antismash.common.secmet.features.protocluster import Protocluster

from .results import DomainPrediction, ProtoclusterPrediction, TerpeneResults
from .terpene_analysis import load_hmm_properties

_hmm_properties = load_hmm_properties()


def will_handle(_products: list[str], categories: set[str]) -> bool:
    """ Returns true if one or more relevant products or product categories are present """
    return "terpene" in categories


def get_domain_description(prediction: DomainPrediction) -> str:
    """ Returns a description of the main type of a domain

        Arguments:
            prediction: a DomainPrediction object

        Returns:
            a string with the type description
    """
    if prediction.domain_type == "ambiguous":
        return "Domain with ambiguous hits"
    return load_hmm_properties()[prediction.domain_type].description


def format_subtype(prediction: DomainPrediction) -> str:
    """ Creates a description of the subtypes of a domain
        If the main type of the domain doesn't have any subtypes associated to it,
        return "none".
        If there exist subtypes, but for this domain the subtype is missing,
        return "unknown".

        Arguments:
            prediction: a DomainPrediction object

        Returns:
            a string with the subtypes description
    """
    if prediction.domain_type == "ambiguous":
        return "none"
    hmm_properties = load_hmm_properties()
    if not prediction.subtypes:
        if hmm_properties[prediction.domain_type].subtypes:
            return "unknown"
        return "none"
    descriptions = []
    for subtype in prediction.subtypes:
        short_description = _hmm_properties[subtype].description
        if ";" in short_description:
            short_description = "".join(short_description.split("; ")[1:])
        descriptions.append(short_description)
    return " or ".join(descriptions)


def format_reactions(prediction: DomainPrediction) -> Markup:
    """ Creates html description list rows containing the reactions of a domain

        Arguments:
            prediction: a DomainPrediction object

        Returns:
            a Markup object containing description list rows
    """
    reactions = []
    for reaction in prediction.reactions:
        substrates = [substrate.name for substrate in reaction.substrates]
        products = [product.name for product in reaction.products]
        for product in products:
            reactions.append(f"<dd>{', '.join(substrates)} &rarr; {product}</dd>")
    return Markup("".join(reactions))


def format_domain_types(domain_preds: list[DomainPrediction]) -> str:
    """ Returns a summary of the main domain types in a cds

        Arguments:
            domain_preds: a list of DomainPrediction objects

        Returns:
            a string containing a summary of the main types
    """
    types = []
    for domain_pred in domain_preds:
        types.append(domain_pred.domain_type)
    return " + ".join(types)


def get_preds_by_cluster(region_layer: RegionLayer, results: TerpeneResults
                         ) -> dict[Protocluster, ProtoclusterPrediction]:
    """ Creates a dictionary with Protocluster objects as keys and
        ProtoclusterPrediction objects as values

        Arguments:
            region_layer: A RegionLayer object
            results: A TerpeneResults object

        Returns:
            a dictionary of Protocluster to ProtoclusterPrediction
    """
    preds_by_cluster = {}
    for cluster in region_layer.get_unique_protoclusters():
        if cluster.product_category == "terpene":
            prediction = results.cluster_predictions[cluster.get_protocluster_number()]
            preds_by_cluster[cluster] = prediction
    return preds_by_cluster


def get_glossary_data(predictions: list[ProtoclusterPrediction]) -> dict[str, str]:
    """ Creates a dictionary containing compound names and extended names,
        for any compounds associated with the region that have an extended name

        Arguments:
            protoclusters: a list of ProtoclusterPrediction objects

        Returns:
            a dictionary of name to extended name
    """
    name_mappings = {}
    for pred in predictions:
        compounds = pred.get_unique_compounds()
        for compound in compounds:
            if compound.extended_name:
                name_mappings[compound.name] = compound.extended_name.capitalize()
        func_groups = pred.get_functional_groups()
        if "PP" in func_groups:
            name_mappings["PP"] = "Diphosphate"
    return dict(sorted(name_mappings.items()))


def generate_html(region_layer: RegionLayer, results: TerpeneResults,
                  _record_layer: RecordLayer, _options_layer: OptionsLayer
                  ) -> HTMLSections:
    """ Generates HTML output for the module """
    html = HTMLSections("terpene")
    preds_by_cluster = get_preds_by_cluster(region_layer, results)
    glossary_data = get_glossary_data(list(preds_by_cluster.values()))

    details_template = FileTemplate(path.get_full_path(__file__, "templates", "details.html"))
    sidepanel_template = FileTemplate(path.get_full_path(__file__, "templates", "sidepanel.html"))

    details_tooltip = " ".join([
        "Detailed information for terpene domains.",
        "Shows possible products of a terpene cluster and \
        lists the subtype and expected reactions for each terpene-related domain."
        ])
    details = details_template.render(
        preds_by_cluster=preds_by_cluster, tooltip=details_tooltip,
        get_domain_description=get_domain_description,
        format_subtype=format_subtype, format_reactions=format_reactions,
        format_domain_types=format_domain_types)
    html.add_detail_section("Terpene", details, class_name="terpene")

    sidepanel_tooltip = "Glossary showing acronyms and their full names."
    sidepanel = sidepanel_template.render(glossary_data=glossary_data, tooltip=sidepanel_tooltip)
    html.add_sidepanel_section("Terpene", sidepanel, class_name="terpene")

    return html
