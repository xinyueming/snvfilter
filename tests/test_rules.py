"""Tests for rules loading."""

import os
import tempfile

import pytest

from snvfilter.rules import Rules


class TestRulesDefaults:
    def test_common_rules(self):
        r = Rules()
        assert "wes_blacklist" in r.wes_blacklist or len(r.wes_blacklist) > 0
        assert "<INV>" in r.exclude_alts
        assert "intron13" in r.met_retained_introns
        assert "Benign" in r.exclude_clinsig
        assert "chr5_1295228" in r.tert_positions

    def test_somatic_rules(self):
        r = Rules()
        assert r.somatic_population_freq.get("onekg_max") == 0.01
        assert r.somatic_local_freq_max == 0.15

    def test_germline_rules(self):
        r = Rules()
        assert r.germline_af_min == 0.15
        pf = r.germline_population_freq
        assert pf.get("onekg_max") == 0.05

    def test_mrd_rules(self):
        r = Rules()
        assert r.mrd.get("af_ratio_max") == 1.0

    def test_somatic_common_pass(self):
        r = Rules()
        af, vd = r.somatic_common_pass(1, "blood")
        assert af == 0.001
        assert vd == 3

        af, vd = r.somatic_common_pass(2, "tissue")
        assert af == 0.01
        assert vd == 10

    def test_large_indel_filter(self):
        r = Rules()
        assert r.large_indel_max_length == 50
        assert "EGFR" in r.large_indel_except_genes
