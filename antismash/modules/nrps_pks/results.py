# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the nrps_pks module """

import logging
from typing import Any, Dict

from Bio.SeqFeature import FeatureLocation

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import AntismashDomain

from .at_analysis.at_analysis import ATSignatureResults
from .parsers import LONG_TO_SHORT

DOMAIN_TYPE_MAPPING = {'Condensation_DCL': 'Condensation',
                       'Condensation_LCL': 'Condensation',
                       'Condensation_Dual': 'Condensation',
                       'Condensation_Starter': 'Condensation',
                       'CXglyc': 'Condensation',
                       'Cglyc': 'Condensation',
                       'cMT': 'MT',
                       'oMT': 'MT',
                       'nMT': 'MT',
                       'Polyketide_cyc': 'Polyketide_cyc',
                       'Polyketide_cyc2': 'Polyketide_cyc'}


_UNKNOWN = "(unknown)"


class PKSResults:
    """ Results for the PKS section of the nrps_pks module """
    __slots__ = ["method_results"]

    class EnforcedDict(dict):
        """ Enforces the value of the key 'signature' to be a consistent type """
        def __setitem__(self, key, val):
            if key == "signature":
                assert isinstance(val, ATSignatureResults)
            return super().__setitem__(key, val)

    def __init__(self):
        self.method_results = PKSResults.EnforcedDict()

    def to_json(self) -> Dict[str, Any]:
        """ Store results as a JSON object """
        results = {}
        for key, value in self.method_results.items():
            if key == "signature":
                value = value.to_json()
            results[key] = value
        return results

    @staticmethod
    def from_json(json: Dict[str, Any]) -> "PKSResults":
        """ Reconstruct a PKSResults instance from JSON object """
        assert isinstance(json, dict)
        results = PKSResults()
        results.method_results = json
        results.method_results["signature"] = ATSignatureResults.from_json(json.get("signature", {}))
        return results


class NRPS_PKS_Results(ModuleResults):
    """ The combined results of the nrps_pks module """
    _schema_version = 1
    __slots__ = ["_pks", "nrps", "consensus", "consensus_transat", "cluster_predictions"]

    def __init__(self, record_id):
        super().__init__(record_id)
        self._pks = PKSResults()
        self.nrps = {}  # not currently used, but will be
        self.consensus = {}
        self.cluster_predictions = {}  # type: Dict[int, str]
        self.consensus_transat = {}

    @property
    def pks(self) -> PKSResults:
        """ The result of the PKS analyses """
        return self._pks

    @pks.setter
    def pks(self, value: PKSResults):
        assert isinstance(value, PKSResults)
        self._pks = value

    def to_json(self) -> Dict[str, Any]:
        results = {"schema_version": self._schema_version}
        results["record_id"] = self.record_id
        results["pks"] = self.pks.to_json()
        results["nrps"] = self.nrps
        results["consensus"] = self.consensus
        results["cluster_predictions"] = self.cluster_predictions
        return results

    @staticmethod
    def from_json(json: Dict[str, Any], _record) -> "NRPS_PKS_Results":
        assert "record_id" in json
        if json.get("schema_version") != NRPS_PKS_Results._schema_version:
            logging.warning("Mismatching schema version, dropping results")
            return None
        results = NRPS_PKS_Results(json["record_id"])
        pks = json.get("pks")
        if pks:
            results.pks = PKSResults.from_json(pks)
        results.nrps = json.get("nrps", {})
        results.consensus = json["consensus"]
        results.cluster_predictions = json["cluster_predictions"]
        return results

    def add_to_record(self, record) -> None:
        """ Save substrate specificity predictions in NRPS/PKS domain sec_met info of record
        """

        for feature in record.get_nrps_pks_cds_features():
            x_count = 0
            nrps_qualifier = feature.nrps_pks
            new_features = []
            gene_id = feature.get_name()
            for domain in nrps_qualifier.domains:
                domain_type = domain.name
                start_aa = domain.start
                end_aa = domain.end
                evalue = domain.evalue
                score = domain.bitscore

                domain.predictions.clear()

                loc = feature.get_sub_location_from_protein_coordinates(start_aa, end_aa)

                #  set up new CDS_motif feature
                new_feature = AntismashDomain(loc)
                new_feature.domain_subtype = domain_type
                if feature.locus_tag:
                    new_feature.locus_tag = feature.locus_tag
                else:
                    new_feature.locus_tag = gene_id
                new_feature.detection = "hmmscan"
                new_feature.database = "nrpspksdomains.hmm"
                new_feature.evalue = evalue
                new_feature.score = score
                if feature.transl_table:
                    transl_table = feature.transl_table
                else:
                    transl_table = 1
                new_feature.translation = str(new_feature.extract(record.seq).translate(table=transl_table))
                domainname = gene_id + domain.label
                if domain_type == "AMP-binding":
                    new_feature.label = domainname
                    new_feature.domain_id = "nrpspksdomains_" + domainname
                    domain.predictions["consensus"] = "nrp"

                elif domain_type == "PKS_AT":
                    new_feature.label = domainname
                    new_feature.domain_id = "nrpspksdomains_" + domainname

                    # For t1pks, t2pks and t3pks
                    if 'transatpks' not in feature.cluster.products:
                        consensus = self.consensus[domainname]
                    else:  # for transatpks
                        consensus = self.consensus_transat[domainname]
                    pks_sig = self.pks.method_results["signature"][domainname]
                    if pks_sig:
                        domain.predictions["PKS signature"] = pks_sig[0].name.rsplit("_", 1)[1]
                    else:
                        domain.predictions["PKS signature"] = _UNKNOWN
                    minowa = self.pks.method_results["minowa_at"][domainname][0][0]
                    domain.predictions["Minowa"] = LONG_TO_SHORT.get(minowa, minowa)
                    domain.predictions["consensus"] = consensus

                elif domain_type == "CAL_domain":
                    new_feature.label = domainname
                    new_feature.domain_id = "nrpspksdomains_" + domainname
                    minowa = self.pks.method_results["minowa_cal"][domainname][0][0]
                    domain.predictions["Minowa"] = LONG_TO_SHORT.get(minowa, minowa)

                elif domain_type == "PKS_KR":
                    new_feature.label = domainname
                    new_feature.domain_id = "nrpspksdomains_" + domainname

                    domain.predictions["KR activity"] = \
                            "active" if self.pks.method_results["kr_activity"][domainname] else "inactive"
                    domain.predictions["KR stereochemistry"] = \
                            self.pks.method_results["kr_stereochem"].get(domainname, _UNKNOWN)
                else:
                    x_count += 1
                    new_feature.domain_id = "nrpspksdomains_" + gene_id.partition(".")[0] \
                                            + "_Xdom"+'{:02d}'.format(x_count)
#                    updated_nrps_qualifier.append(domain) # TODO weird, but should it be done?
                for method, pred in domain.predictions.items():
                    new_feature.specificity.append("%s: %s" % (method, pred))
                mapping = DOMAIN_TYPE_MAPPING.get(domain_type)
                if mapping:
                    new_feature.domain_subtype = domain_type
                    new_feature.domain = mapping
                new_features.append(new_feature)

            for new_feature in new_features:
                record.add_feature(new_feature)
