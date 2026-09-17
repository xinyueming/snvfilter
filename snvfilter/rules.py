"""Filtering rules module.

Loads rules from config/rules.yaml and provides a ``Rules`` class for
programmatic access to each rule category.
"""

import os

from .config import load_rules as _load_rules_yaml


class Rules:
    """Container for all filtering rules.

    Attributes
    ----------
    common : dict
        Common rules shared by somatic and germline filtering.
    somatic : dict
        Somatic-specific rule parameters.
    germline : dict
        Germline-specific rule parameters.
    mrd : dict
        MRD filtering parameters.
    """

    def __init__(self, rules_file=None):
        raw = _load_rules_yaml(rules_file)
        self.common = raw.get("common", {})
        self.somatic = raw.get("somatic", {})
        self.germline = raw.get("germline", {})
        self.mrd = raw.get("mrd", {})

    # -- convenience accessors ------------------------------------------------

    @property
    def wes_blacklist(self):
        """Gene blacklist for WES panels."""
        return self.common.get("product_gene_filter", {}).get("wes_blacklist", [])

    @property
    def exclude_alts(self):
        """ALT types to exclude (structural variants)."""
        return self.common.get("structural_variant_filter", {}).get("exclude_alts", [])

    @property
    def met_retained_introns(self):
        """MET gene introns to retain."""
        return self.common.get("gene_function_filter", {}).get("met_retained_introns", [])

    @property
    def exclude_clinsig(self):
        """CLNSIG values to filter out as benign."""
        return self.common.get("benign_variant_filter", {}).get("exclude_clinsig", [])

    @property
    def tert_positions(self):
        """TERT hotspot positions."""
        return self.common.get("tert_hotspot", {}).get("positions", [])

    @property
    def large_indel_except_genes(self):
        """Genes excluded from large indel length filter."""
        return self.common.get("large_indel_filter", {}).get("except_genes", [])

    @property
    def large_indel_max_length(self):
        """Max indel length before filtering."""
        return self.common.get("large_indel_filter", {}).get("max_length", 50)

    def somatic_common_pass(self, pass_id, sample_type):
        """Return (af_min, vd_min) for a somatic common_pass_N block."""
        block = self.somatic.get(f"common_pass_{pass_id}", {})
        thresholds = block.get(sample_type, {})
        return thresholds.get("af_min"), thresholds.get("vd_min")

    @property
    def somatic_population_freq(self):
        return self.somatic.get("population_freq", {})

    @property
    def somatic_clonal_hematopoiesis(self):
        return self.somatic.get("clonal_hematopoiesis", {})

    @property
    def somatic_strand_bias(self):
        return self.somatic.get("strand_bias", {})

    @property
    def somatic_local_freq_max(self):
        return self.somatic.get("local_freq_max", 0.15)

    @property
    def germline_af_min(self):
        return self.germline.get("af_min", 0.15)

    @property
    def germline_population_freq(self):
        return self.germline.get("population_freq", {})

    @property
    def germline_required_gene_categories(self):
        return self.germline.get("required_gene_categories", [])
