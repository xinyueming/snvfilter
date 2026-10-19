"""Tests for config loading."""

import os
import tempfile

import pytest
import yaml

from snvfilter.config import load_config, load_rules, _deep_merge, load_yaml


class TestDeepMerge:
    def test_simple_override(self):
        base = {"a": 1, "b": 2}
        override = {"b": 3}
        _deep_merge(base, override)
        assert base == {"a": 1, "b": 3}

    def test_nested_merge(self):
        base = {"input": {"key_site": "default.tsv", "product_gene": "default.list"}}
        override = {"input": {"key_site": "custom.tsv"}}
        _deep_merge(base, override)
        assert base["input"]["key_site"] == "custom.tsv"
        assert base["input"]["product_gene"] == "default.list"

    def test_deeply_nested(self):
        base = {"sample": {"type": "blood"}}
        override = {"sample": {"type": "tissue"}}
        _deep_merge(base, override)
        assert base["sample"]["type"] == "tissue"


class TestLoadYaml:
    def test_load_existing(self, tmp_path):
        p = tmp_path / "cfg.yaml"
        p.write_text("key: value\n")
        assert load_yaml(str(p)) == {"key": "value"}

    def test_load_none(self):
        assert load_yaml(None) == {}

    def test_load_missing(self):
        assert load_yaml("/nonexistent/path.yaml") == {}


class TestLoadConfig:
    def test_defaults_only(self):
        cfg = load_config()
        assert "input" in cfg
        assert "sample" in cfg

    def test_user_config_override(self, tmp_path):
        user = tmp_path / "user.yaml"
        user.write_text("sample:\n  type: tissue\n")
        cfg = load_config(config_file=str(user))
        assert cfg["sample"]["type"] == "tissue"

    def test_cli_overrides_user(self, tmp_path):
        user = tmp_path / "user.yaml"
        user.write_text("sample:\n  type: tissue\n")
        cfg = load_config(config_file=str(user), cli_args={"sample": {"type": "blood"}})
        assert cfg["sample"]["type"] == "blood"

    def test_cli_adds_new_keys(self):
        cfg = load_config(cli_args={"extra_key": "extra_val"})
        assert cfg["extra_key"] == "extra_val"


class TestLoadRules:
    def test_default_rules(self):
        rules = load_rules()
        assert "common" in rules
        assert "somatic" in rules
        assert "germline" in rules
        assert "mrd" in rules

    def test_custom_rules(self, tmp_path):
        r = tmp_path / "rules.yaml"
        r.write_text("common:\n  test_key: value\n")
        rules = load_rules(rules_file=str(r))
        assert rules["common"]["test_key"] == "value"


class TestGenomeConfig:
    def test_hg19_merge(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text(
            "input:\n"
            "  key_site: common.tsv\n"
            "  hg19:\n"
            "    refgene: hg19_refgene.txt\n"
            "    gtf: hg19.gtf\n"
            "  hg38:\n"
            "    refgene: hg38_refgene.txt\n"
            "    gtf: hg38.gtf\n"
        )
        c = load_config(config_file=str(cfg), genome="hg19")
        assert c["input"]["key_site"] == "common.tsv"
        assert c["input"]["refgene"] == "hg19_refgene.txt"
        assert c["input"]["gtf"] == "hg19.gtf"
        # hg38 keys should NOT be merged
        assert "hg38" not in c["input"]

    def test_hg38_merge(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text(
            "input:\n"
            "  key_site: common.tsv\n"
            "  hg19:\n"
            "    refgene: hg19_refgene.txt\n"
            "  hg38:\n"
            "    refgene: hg38_refgene.txt\n"
            "    gtf: hg38.gtf\n"
        )
        c = load_config(config_file=str(cfg), genome="hg38")
        assert c["input"]["refgene"] == "hg38_refgene.txt"
        assert c["input"]["gtf"] == "hg38.gtf"
        assert "hg19" not in c["input"]

    def test_no_genome_no_merge(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text(
            "input:\n"
            "  key_site: common.tsv\n"
            "  hg19:\n"
            "    refgene: hg19_refgene.txt\n"
        )
        c = load_config(config_file=str(cfg), genome=None)
        assert c["input"]["key_site"] == "common.tsv"
        # hg19 sub-dict stays untouched
        assert "hg19" in c["input"]

    def test_cli_overrides_genome(self, tmp_path):
        cfg = tmp_path / "cfg.yaml"
        cfg.write_text(
            "input:\n"
            "  hg19:\n"
            "    refgene: hg19_refgene.txt\n"
        )
        c = load_config(config_file=str(cfg), genome="hg19", cli_args={"input": {"refgene": "cli_refgene.txt"}})
        assert c["input"]["refgene"] == "cli_refgene.txt"
