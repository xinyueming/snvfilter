"""Filtering logic — VariantRow class and rule evaluation."""

import re

from .rules import Rules


class VariantRow:
    """Represents a single variant line from an ANNOVAR-annotated VCF.

    All filtering properties (rule_1 … germline_rule_3, somatic_rule_*,
    common_pass_*, etc.) are evaluated lazily.
    """

    def __init__(self, row: str, rules: Rules, **ctx):
        """
        Parameters
        ----------
        row : str
            One tab-delimited line from the ANNOVAR VCF.
        rules : Rules
            Loaded rule configuration.
        ctx : dict
            Pre-loaded data dicts: product_gene_dict, hotspot_gene_dict,
            key_gene_changes, key_gene_tx, key_gene_tx_supplement,
            gene_categorie_dict, refgene_exon_dict, geneid, varannovar,
            snvdepth, snvdepth_control, gene_alias_dict, local_anno_tuple,
            mrd_filter_dict, plus flags sample_type, force_germline,
            is_brca, is_lynch, is_105, product_gene.
        """
        self.rules = rules
        self._ctx = ctx

        # context
        self.product_gene = ctx.get("product_gene", "")
        self.product_gene_dict = ctx.get("product_gene_dict", {})
        self.sample_type = ctx.get("sample_type", "blood")
        self.hotspot_gene_dict = ctx.get("hotspot_gene_dict", {})
        self.key_gene_changes = ctx.get("key_gene_changes", {})
        self.key_gene_tx = ctx.get("key_gene_tx", {})
        self.key_gene_tx_supplement = ctx.get("key_gene_tx_supplement", {})
        self.gene_categorie_dict = ctx.get("gene_categorie_dict", {})
        self.refgene_exon_dict = ctx.get("refgene_exon_dict", {})
        self.gene_alias_dict = ctx.get("gene_alias_dict", {})
        self.varannovar = ctx.get("varannovar", {})
        self.snvdepth = ctx.get("snvdepth", {})
        self.snvdepth_control = ctx.get("snvdepth_control", {})
        self.mrd_filter_dict = ctx.get("mrd_filter_dict", {})

        # flags
        self.force_germline = ctx.get("force_germline", False)
        self.is_brca = ctx.get("is_brca", False)
        self.is_lynch = ctx.get("is_lynch", False)
        self.is_105 = ctx.get("is_105", False)
        self.hrr_panel = "zy.hrr" in (self.product_gene or "")

        # parse raw line
        self.arr = row.split("\t")
        self.alt = self.arr[4]
        self.filter_snv = self.arr[6]
        self.info = self.arr[7]
        self.info_dict = {}
        self.info_dict2 = {}
        self.parser_info(self.info)
        self.parser_info2(self.info)

        self.key_infos = f"{self.arr[0]}_{self.arr[1]}_{self.arr[3]}_{self.arr[4]}"
        if self.key_infos in self.varannovar:
            v = self.varannovar[self.key_infos]
            self.info_dict["Gene.refGeneWithVer"] = v[0]
            self.info_dict["Func.refGeneWithVer"] = v[1]
            self.info_dict["ExonicFunc.refGeneWithVer"] = v[2]
            self.info_dict["AAChange.refGeneWithVer"] = v[3]
            self.info_dict["gene_details"] = v[4]

        # depth / af / vd
        self.dp = int(self.info_dict["DP"])
        self.vd = int(self.info_dict["VD"])
        self.af = float(self.info_dict["AF"])
        self.somatic_SBF = "."
        self.somatic_ODDRATIO = "."
        try:
            self.pvalue = float(self.info_dict["SSF"])
        except Exception:
            self.pvalue = "NA"
            self.somatic_SBF = self.info_dict.get("SBF", ".")
            self.somatic_ODDRATIO = self.info_dict.get("ODDRATIO", ".")

        self.info_prefix = self.arr[8]
        self.info_sample_1 = self.arr[9]
        self.info_sample_dict_1 = {}
        for k, v in zip(self.info_prefix.split(":"), self.info_sample_1.split(":")):
            self.info_sample_dict_1[k] = v
        self.ald = self.info_sample_dict_1.get("ALD", ".").strip()
        self.ald_f = self.ald.split(",")[0].strip()
        self.ald_r = self.ald.split(",")[-1].strip()
        self.rd = self.info_sample_dict_1.get("RD", ".").strip()
        self.rd_f = self.rd.split(",")[0].strip()
        self.rd_r = self.rd.split(",")[-1].strip()
        self.gt = self.info_sample_dict_1.get("GT", ".").strip().replace("/", "|")
        try:
            if int(self.ald_f) + int(self.ald_r) > 0:
                self.ald = self.ald + ";" + str(
                    round(int(self.ald_f) / (int(self.ald_f) + int(self.ald_r)), 2)
                )
            else:
                self.ald = self.ald + ";."
            if int(self.rd_f) + int(self.rd_r) > 0:
                self.rd = self.rd + ";" + str(
                    round(int(self.rd_f) / (int(self.rd_f) + int(self.rd_r)), 2)
                )
            else:
                self.rd = self.rd + ";."
        except Exception:
            pass

        # paired sample
        self.germline_dp = "."
        self.germline_vd = "."
        self.germline_af = "."
        self.germline_ald = "."
        self.germline_ald_f = "."
        self.germline_ald_r = "."
        self.germline_rd = "."
        self.germline_SBF = "."
        self.germline_ODDRATIO = "."
        self.germline_gt = "."
        if len(self.arr) == 11:  # paired
            self.info_sample_2 = self.arr[10]
            self.info_sample_dict_2 = {}
            for k, v in zip(self.info_prefix.split(":"), self.info_sample_2.split(":")):
                self.info_sample_dict_2[k] = v
            self.germline_ald = self.info_sample_dict_2.get("ALD", ".").strip()
            self.germline_ald_f = self.germline_ald.split(",")[0].strip()
            self.germline_ald_r = self.germline_ald.split(",")[-1].strip()
            self.germline_rd = self.info_sample_dict_2.get("RD", ".").strip()
            self.germline_rd_f = self.germline_rd.split(",")[0].strip()
            self.germline_rd_r = self.germline_rd.split(",")[-1].strip()
            self.germline_gt = self.info_sample_dict_2.get("GT", ".").strip().replace("/", "|")
            self.somatic_SBF = self.info_sample_dict_1.get("SBF", ".").strip()
            self.somatic_ODDRATIO = self.info_sample_dict_1.get("ODDRATIO", ".").strip()
            self.germline_SBF = self.info_sample_dict_2.get("SBF", ".").strip()
            self.germline_ODDRATIO = self.info_sample_dict_2.get("ODDRATIO", ".").strip()
            try:
                self.germline_dp = int(self.info_sample_dict_2.get("DP"))
                self.germline_vd = int(self.info_sample_dict_2.get("VD"))
                self.germline_af = float(self.info_sample_dict_2.get("AF"))
            except Exception:
                pass
            if self.germline_vd == 0:
                self.pvalue = "NA"

        if self.force_germline:
            self.vd, self.germline_vd = self.germline_vd, self.vd
            self.dp, self.germline_dp = self.germline_dp, self.dp
            self.af, self.germline_af = self.germline_af, self.af

        # transcript / gene parsing
        self.parser_common_tx(self.gene, self.key_gene_changes, self.key_gene_tx, self.key_gene_tx_supplement)
        self.parser_gene_detail()
        self.refgene_exon_n = "."
        self.parser_exon(self.new_tx, self.refgene_exon_dict)
        try:
            self.geneid = ctx["geneid"][self.gene]
        except Exception:
            self.geneid = "."

        # local annotation
        local_tuple = ctx.get("local_anno_tuple", ({}, {}, []))
        self.local_anno_dict = local_tuple[0]
        self.local_germline_anno_dict = local_tuple[1]
        self.local_header = local_tuple[2]
        self.comment, self.local_freq = self.anno_file(self.local_anno_dict)
        self.germline_comment, self.germline_local_freq = self.anno_file(self.local_germline_anno_dict)

    # ---- parsing helpers -----------------------------------------------------

    def parser_info(self, info):
        for i in info.split(";"):
            if "=" not in i:
                continue
            k, v = i.split("=")
            if k == "cosmic70":
                k = "cosmic96"
            if k == "gnomad312_AF_eas":
                k = "AF_eas"
            if k == "gnomad312_AF_popmax":
                k = "AF_popmax"
            self.info_dict.setdefault(k, v)

    def parser_info2(self, info):
        for i2 in info.split(";"):
            if "=" not in i2:
                continue
            k, v = i2.split("=")
            if k == "cosmic70":
                k = "cosmic96"
            if k == "gnomad312_AF_eas":
                k = "AF_eas"
            if k == "gnomad312_AF_popmax":
                k = "AF_popmax"
            if k == "gnomad312_AF":
                k = "AF"
            self.info_dict2[k] = v

    def parser_exon(self, new_tx, refgene_exon_dict):
        n_trans = new_tx.split(".")[0]
        try:
            self.refgene_exon_n = refgene_exon_dict[n_trans]
        except Exception:
            pass

    def parser_common_tx(self, gene, key_gene_changes, key_gene_tx, key_gene_tx_supplement):
        self.common_aa_change = "."
        self.common_aa_change = _match_gene_tx(self.gene_detail, key_gene_tx_supplement, default=False)
        if self.common_aa_change == ".":
            if key_gene_changes.get(gene):
                for key_change in key_gene_changes.get(gene):
                    for ichange in self.gene_detail:
                        p_info = ichange.split(":")[-1]
                        if key_change in p_info:
                            self.common_aa_change = ichange
                            break
        if self.common_aa_change == ".":
            self.common_aa_change = _match_gene_tx(self.gene_detail, key_gene_tx, default=True)

        p_list_es = ["ASXL1:p.G646Wfs*10"]
        aas = self.common_aa_change.split(":")
        aaa = aas[-1]
        p_es = f"{gene}:{aaa}"
        if "p." in aaa and "fs" in aaa:
            aaa = aaa.replace("X", "*")
            if "fs*" in aaa:
                aaas = aaa.split("fs*")
                if len(aaas) > 1:
                    if p_es in p_list_es:
                        aaas[-1] = str(int(aaas[-1]) + 2)
                    else:
                        aaas[-1] = str(int(aaas[-1]) + 1)
                    aaa = "fs*".join(aaas)
        aas[-1] = aaa
        self.common_aa_change = ":".join(aas)

    def parser_gene_detail(self):
        self.new_gene, self.new_tx, self.new_exon, self.new_cdna, self.new_amid = ".", ".", ".", ".", "."
        if self.common_aa_change != ".":
            if self.common_aa_change != "UNKNOWN":
                change_arr = self.common_aa_change.split(":")
                try:
                    self.new_gene, self.new_tx, self.new_exon = change_arr[:3]
                except Exception:
                    pass
                if len(change_arr) == 5:
                    self.new_cdna, self.new_amid = change_arr[3:5]
                elif len(change_arr) == 4:
                    self.new_cdna = change_arr[3]
                elif len(change_arr) == 3:
                    self.new_exon = "."
                    self.new_cdna = change_arr[-1]
                if ">" not in self.new_cdna:
                    tem_new_cdna = self.new_cdna.lstrip("c.")
                    a_new_cdna = re.findall("[A-Z]+", tem_new_cdna)
                    n_new_cdna = re.findall("\d+", tem_new_cdna)
                    if len(a_new_cdna) == 2:
                        try:
                            if len(n_new_cdna) == 1 and len(a_new_cdna[0]) == 1 and len(a_new_cdna[1]) == 1:
                                self.new_cdna = f"c.{n_new_cdna[0]}{a_new_cdna[0]}>{a_new_cdna[1]}"
                        except Exception:
                            pass
                if "substitution" in self.info_dict.get("ExonicFunc.refGeneWithVer", "") or "stopgain" in self.info_dict.get("ExonicFunc.refGeneWithVer", ""):
                    if self.arr[3][0] != self.arr[4][0]:
                        try:
                            res = re.match(r"^(c.\d+_\d+)([a-z,A-Z]+)", self.new_cdna)
                            self.new_cdna = f"{res.groups()[0]}delins{res.groups()[1]}"
                        except Exception:
                            pass
                self.new_amid = self.new_amid.replace("X", "*")
            else:
                self.new_gene = self.gene
        else:
            self.new_gene = self.gene
        if self.gene in self.gene_alias_dict:
            self.new_gene = self.gene_alias_dict.get(self.gene)

    # ---- properties ----------------------------------------------------------

    @property
    def gene(self):
        return self.info_dict.get("Gene.refGeneWithVer", ".")

    @property
    def func(self):
        return self.info_dict.get("Func.refGeneWithVer", ".")

    @property
    def gene_detail(self):
        if self.info_dict.get("AAChange.refGeneWithVer", ".") == ".":
            changes = self.info_dict.get("GeneDetail.refGeneWithVer", "").split("\\x3b")
            changes = ["." if i == "." else f"{self.gene}:{i}" for i in changes]
        else:
            changes = self.info_dict["AAChange.refGeneWithVer"].split(",")
        return changes

    @property
    def gene_details_string(self):
        new_changes = []
        p_list_es = ["ASXL1:p.G646Wfs*10"]
        for aa in self.gene_detail:
            aas = aa.split(":")
            aaa = aas[-1]
            p_es = f"{self.gene}:{aaa}"
            if "p." in aaa and "fs" in aaa:
                aaa = aaa.replace("X", "*")
                if "fs*" in aaa:
                    aaas = aaa.split("fs*")
                    if len(aaas) > 1:
                        if p_es in p_list_es:
                            aaas[-1] = str(int(aaas[-1]) + 2)
                        else:
                            aaas[-1] = str(int(aaas[-1]) + 1)
                        aaa = "fs*".join(aaas)
            aas[-1] = aaa
            new_changes.append(":".join(aas))
        return ",".join(new_changes)

    @property
    def exonic_func(self):
        return self.info_dict.get("ExonicFunc.refGeneWithVer", ".")

    @property
    def clinvar_sig(self):
        return self.info_dict.get("CLNSIG", ".")

    @property
    def onco_sig(self):
        return self.info_dict.get("ONCSIG", ".")

    @property
    def onekg_af(self):
        return self.info_dict.get("1000g2015aug_all", ".")

    @property
    def gnomad_af(self):
        return self.info_dict.get("AF_eas", ".")

    @property
    def cosmic(self):
        return self.info_dict.get("cosmic96", ".")

    @property
    def in_pathogenic(self):
        return self.clinvar_sig in ["drug_response"] or self.onco_sig in [
            "Oncogenic/Likely_oncogenic", "Oncogenic", "Likely_oncogenic",
        ]

    @property
    def status(self):
        if self.force_germline:
            return "germline_data"
        return self.info_dict.get("STATUS", ".")

    @property
    def is_somatic(self):
        return self.status in [".", "StrongSomatic", "LikelySomatic"]

    @property
    def is_germline(self):
        return self.status in ["Germline", "AFDiff", "germline_data", "StrongLOH", "LikelyLOH"]

    @property
    def gene_category(self):
        categories = []
        for cat in ["genetic", "drug", "immune_positive", "immune_negative", "immune_hyperprogressive", "chemical", "mmr"]:
            if self.gene_categorie_dict.get(cat, {}).get(self.gene):
                categories.append(cat)
        return ";".join(categories) if categories else "."

    @property
    def repeat_region_anno(self):
        complex_anno = ""
        if self.info_dict.get("simple_repeat", ".") != "." or self.info_dict.get("rmsk", ".") != ".":
            complex_anno = "LCR"
        if self.info_dict.get("genomicSuperDups", ".") != ".":
            complex_anno = complex_anno + ";SDR" if complex_anno else "SDR"
        return complex_anno if complex_anno else "."

    @property
    def mrd_label(self):
        if not self.mrd_filter_dict:
            return "."
        key_info1 = f"{self.arr[0]}:{self.arr[1]}:{self.arr[3]}:{self.arr[4]}"
        key_info2 = f"{self.gene}:{self.new_cdna}" if self.new_cdna != "." else None
        key_info3 = f"{self.gene}:{self.new_amid}" if self.new_amid != "." else None
        for label in ["unique", "general"]:
            for k in ["position", "cHGVS", "pHGVS"]:
                if key_info1 in self.mrd_filter_dict[k][label] or \
                   key_info2 in self.mrd_filter_dict[k][label] or \
                   key_info3 in self.mrd_filter_dict[k][label]:
                    return label
        return "."

    @property
    def end_pos(self):
        pos = int(self.arr[1])
        len_ref = len(self.arr[3])
        len_alt = len(self.arr[4])
        if len_ref > len_alt:
            return pos + (len_ref - len_alt)
        return pos

    @property
    def af_string(self):
        if self.af == ".":
            return "."
        return f"{round(self.af * 100, 2)}%"

    @property
    def germline_af_string(self):
        if self.germline_af == ".":
            return "."
        return f"{round(self.germline_af * 100, 2)}%"

    @property
    def snv_depth(self):
        key_infos = f"{self.arr[0]}_{str(int(self.arr[1]))}"
        if key_infos in self.snvdepth:
            try:
                if int(self.snvdepth[key_infos]) > int(self.dp):
                    return self.snvdepth[key_infos]
            except Exception:
                return self.dp
        try:
            if int(self.snvdepth[f"{self.arr[0]}_{str(int(self.arr[1]) + 1)}"]) > int(self.dp):
                return self.snvdepth[f"{self.arr[0]}_{str(int(self.arr[1]) + 1)}"]
        except Exception:
            return self.dp
        return self.dp

    @property
    def snv_depth_control(self):
        key_infos = f"{self.arr[0]}_{str(int(self.arr[1]))}"
        if key_infos in self.snvdepth_control:
            try:
                if int(self.snvdepth_control[key_infos]) > int(self.germline_dp):
                    return self.snvdepth_control[key_infos]
            except Exception:
                return self.germline_dp
        try:
            if int(self.snvdepth_control[f"{self.arr[0]}_{str(int(self.arr[1]) + 1)}"]) > int(self.germline_dp):
                return self.snvdepth_control[f"{self.arr[0]}_{str(int(self.arr[1]) + 1)}"]
        except Exception:
            return self.germline_dp
        return self.germline_dp

    # ---- annotation ----------------------------------------------------------

    def anno_file(self, anno_db):
        comment = "\t".join(["."] * 9) + "\n"
        local_freq = 0
        try:
            key1 = tuple(self.arr[:2] + self.arr[3:5])
            key2 = tuple(self.common_aa_change)
            if key2 in anno_db:
                comment = "\t".join(anno_db[key2]) + "\n"
                local_freq = float(anno_db[key2][0])
            elif key1 in anno_db:
                comment = "\t".join(anno_db[key1]) + "\n"
                local_freq = float(anno_db[key1][0])
        except Exception:
            pass
        return comment, local_freq

    # ---- common rules --------------------------------------------------------

    @property
    def rule_1(self):
        """Product gene list filter."""
        if "FullRNA.gene.list" in self.product_gene or "zy.wes.gene.list" in self.product_gene:
            if self.gene in self.rules.wes_blacklist:
                return False
            if len(self.gene.split("\\x3b")) > 1:
                return False
            return True
        if self.product_gene_dict.get(self.gene) or self.product_gene_dict.get(self.gene_alias_dict.get(self.gene)):
            return True
        if "\\x3b" in self.gene:
            parts = self.gene.split("\\x3b")
            if self.product_gene_dict.get(parts[0]) or self.product_gene_dict.get(parts[1]):
                return True
        return False

    @property
    def rule_2(self):
        """Structural variant filter."""
        return self.alt not in self.rules.exclude_alts

    @property
    def rule_3(self):
        """Gene function filter."""
        if self.gene == "MET":
            if self.new_exon in self.rules.met_retained_introns:
                return True
        if "intronic" in self.func and "athogenic" in self.info_dict.get("CLNSIGCONF", "."):
            return True
        if ("exonic" in self.func or "splicing" in self.func) and "ncRNA" not in self.func:
            if self.exonic_func != "synonymous_SNV":
                return True
        return False

    @property
    def rule_4(self):
        """Benign variant filter."""
        return self.clinvar_sig not in self.rules.exclude_clinsig

    @property
    def rule_6(self):
        """TERT hotspot positions."""
        tert_var = f"{self.arr[0]}_{self.arr[1]}"
        if self.gene == "TERT" and tert_var in self.rules.tert_positions:
            return True
        return False

    @property
    def rule_7(self):
        """Large indel filter."""
        if (len(self.arr[3]) > self.rules.large_indel_max_length or len(self.arr[4]) > self.rules.large_indel_max_length) \
                and self.gene not in self.rules.large_indel_except_genes:
            return False
        return True

    # ---- somatic common pass conditions -------------------------------------

    @property
    def somatic_common_pass_1(self):
        if "MRD" in self.product_gene:
            return True
        af_min, vd_min = self.rules.somatic_common_pass(1, self.sample_type)
        if af_min is None or vd_min is None:
            return False
        return self.af >= af_min and self.vd >= vd_min

    @property
    def somatic_common_pass_2(self):
        if "MRD" in self.product_gene:
            return True
        af_min, vd_min = self.rules.somatic_common_pass(2, self.sample_type)
        if af_min is None or vd_min is None:
            return False
        return self.af >= af_min and self.vd >= vd_min

    @property
    def somatic_common_pass_3(self):
        return self.filter_snv == "PASS"

    @property
    def somatic_common_pass_4(self):
        if self.is_paired and self.germline_af != "." and float(self.germline_af) > 0 \
                and float(self.af) / float(self.germline_af) < 3:
            return False
        return True

    @property
    def somatic_common_pass_5(self):
        if "LCR" in self.repeat_region_anno and self.gene not in ["BRCA1", "BRCA2", "MLH1", "MSH6", "PMS2", "MSH2"]:
            if len(self.arr[3]) != len(self.arr[4]) or len(self.arr[3]) > 1:
                thresh = self.rules.somatic_common_pass(5, self.sample_type)  # not used, threshold is per-sample-type
                lcr_cfg = self.rules.somatic.get("common_pass_5", {}).get("lcr_af_threshold", {})
                af_thresh = lcr_cfg.get(self.sample_type, 0.01)
                if float(self.af) < af_thresh:
                    return False
        return True

    # ---- somatic rules -------------------------------------------------------

    @property
    def somatic_rule_1(self):
        pf = self.rules.somatic_population_freq
        if self.onekg_af != "." and float(self.onekg_af) > pf.get("onekg_max", 0.01):
            return False
        if self.gnomad_af != "." and float(self.gnomad_af) > pf.get("gnomad_max", 0.01):
            return False
        if self.info_dict2.get("AF", ".") != "." and float(self.info_dict2["AF"]) > pf.get("gnomad_af_max", 0.01):
            return False
        if self.info_dict.get("AF_popmax", ".") != "." and float(self.info_dict["AF_popmax"]) > pf.get("af_popmax_max", 0.01):
            return False
        return True

    @property
    def somatic_rule_2(self):
        ch = self.rules.somatic_clonal_hematopoiesis
        if self.is_paired and self.germline_af != "." and float(self.germline_af) > 0 \
                and float(self.af) / float(self.germline_af) < ch.get("af_ratio_max", 1.0):
            return False
        if self.is_paired and self.germline_af != "." and float(self.germline_af) > 0 \
                and self.pvalue > ch.get("fisher_pvalue_max", 0.05):
            return False
        return True

    @property
    def somatic_rule_3(self):
        sb = self.rules.somatic_strand_bias
        vd = self.vd
        if vd == 0:
            return False
        return sb.get("min_ratio", 0.2) <= float(int(self.ald_f) / vd) <= sb.get("max_ratio", 0.8)

    @property
    def somatic_rule_4_1(self):
        is_hotspot = self.hotspot_gene_dict.get(f"{self.gene}|{self.new_amid}") \
            or self.hotspot_gene_dict.get(f"{self.gene}|{self.new_cdna}")
        is_pathogenic = self.in_pathogenic and self.somatic_common_pass_3
        return (is_hotspot or is_pathogenic) and self.somatic_common_pass_1

    @property
    def somatic_rule_4_2(self):
        return (self.somatic_common_pass_2 and self.somatic_common_pass_3
                and self.local_freq <= self.rules.somatic_local_freq_max
                and self.somatic_common_pass_4 and self.somatic_common_pass_5)

    # ---- germline rules ------------------------------------------------------

    @property
    def germline_rule_1(self):
        if self.germline_af != "." and float(self.germline_af) <= self.rules.germline_af_min:
            return False
        return True

    @property
    def germline_rule_2(self):
        pf = self.rules.germline_population_freq
        if self.onekg_af != "." and float(self.onekg_af) > pf.get("onekg_max", 0.05):
            return False
        if self.gnomad_af != "." and float(self.gnomad_af) > pf.get("gnomad_max", 0.05):
            return False
        if self.info_dict2.get("AF", ".") != "." and float(self.info_dict2["AF"]) > pf.get("gnomad_af_max", 0.05):
            return False
        if self.info_dict.get("AF_popmax", ".") != "." and float(self.info_dict["AF_popmax"]) > pf.get("af_popmax_max", 0.05):
            return False
        return True

    @property
    def germline_rule_3(self):
        genetic_dict = self.gene_categorie_dict.get("genetic", {})
        return bool(genetic_dict.get(self.gene))

    # ---- MRD -----------------------------------------------------------------

    def mrd_filter(self, germline_filter=False):
        key_site1 = f"{self.arr[0]}:{self.arr[1]}:{self.arr[3]}:{self.arr[4]}"
        key_site2 = f"{self.gene}:{self.new_cdna}" if self.new_cdna != "." else None
        if (key_site1 and key_site1 in self.mrd_filter_dict["position"]["all"]) \
                or (key_site2 and key_site2 in self.mrd_filter_dict["cHGVS"]["all"]):
            if not germline_filter and self.is_paired and self.germline_af != "." \
                    and float(self.germline_af) > 0 and float(self.af) / float(self.germline_af) < 1:
                return False
            if germline_filter and self.germline_af != "." and float(self.germline_af) < self.rules.germline.get("germline_af_min", 0.15):
                return False
            return True
        return False

    # ---- output --------------------------------------------------------------

    @property
    def header(self):
        return "\t".join([
            "Chr", "Pos", "End_pos", "Ref", "Alt", "Filter", "GeneID",
            "Gene.refGeneWithVer", "Func.refGeneWithVer", "ExonicFunc.refGeneWithVer",
            "AAChange.refGeneWithVer", "gene_details", "simple_repeat", "rmsk",
            "genomicSuperDups", "Low_complexity_and_SuperDup_regions",
            "1000g2015aug_all", "AF_all", "AF_popmax", "AF_eas", "avsnp150",
            "CLNALLELEID", "CLNDN", "CLNDISDB", "CLNREVSTAT", "CLNSIG",
            "CLNSIGCONF", "InterVar", "ACMG", "ONCDN", "ONCDISDB",
            "ONCREVSTAT", "ONCSIG", "cosmic96", "SIFT_pred",
            "Polyphen2_HDIV_pred", "LRT_pred", "MutationTaster_pred",
            "FATHMM_pred", "PROVEAN_pred", "M-CAP_pred", "PrimateAI_pred",
            "REVEL_score", "REVEL_rankscore", "TAll_depth_raw", "TAll_depth",
            "TAlt_depth", "TAlt_depth_forward", "TAlt_depth_reverse",
            "TAlt_depth_ALD", "TREF_depth_RD", "Somatic_SBF", "Somatic_ODDRATIO",
            "genotype", "Status", "germline_depth_raw", "germline_depth",
            "germline_alt_depth", "germline_alt_depth_forward",
            "germline_alt_depth_reverse", "germline_alt_depth_ALD",
            "germline_alt_depth_RD", "germline_SBF", "germline_ODDRATIO",
            "germline_af", "germline_genotype", "fisher_p_value", "gene_category",
            "mrd_label", "new_transcript", "all_exon_num", "new_gene", "new_exon",
            "new_cdna", "new_amid", "dcs_supporting_reads", "scs_supporting_reads",
            "TAF", "level",
        ] + self.local_header[5:]) + "\n"

    @property
    def row(self):
        new_amid_clean = self.new_amid
        if new_amid_clean.endswith("*fs*1"):
            new_amid_clean = new_amid_clean.replace("*fs*1", "*")
        elif new_amid_clean.endswith("fs*1"):
            new_amid_clean = new_amid_clean.replace("fs*1", "*")

        return "\t".join([str(s) for s in [
            self.arr[0], self.arr[1], self.end_pos, self.arr[3], self.arr[4],
            self.arr[6], self.geneid, self.gene, self.func, self.exonic_func,
            self.common_aa_change, self.gene_details_string,
            self.info_dict.get("simple_repeat", "."),
            self.info_dict.get("rmsk", "."),
            self.info_dict.get("genomicSuperDups", ".").replace("\\x3b", ";").replace("\\x3d", "="),
            self.repeat_region_anno, self.onekg_af, self.info_dict2.get("AF", "."),
            self.info_dict.get("AF_popmax", "."), self.gnomad_af,
            self.info_dict.get("avsnp150", "."), self.info_dict.get("CLNALLELEID", "."),
            self.info_dict.get("CLNDN", "."), self.info_dict.get("CLNDISDB", "."),
            self.info_dict.get("CLNREVSTAT", "."), self.info_dict.get("CLNSIG", "."),
            self.info_dict.get("CLNSIGCONF", "."), self.info_dict.get("InterVar", "."),
            self.info_dict.get("ACMG", "."), self.info_dict.get("ONCDN", "."),
            self.info_dict.get("ONCDISDB", "."), self.info_dict.get("ONCREVSTAT", "."),
            self.info_dict.get("ONCSIG", "."), self.info_dict.get("cosmic96", "."),
            self.info_dict.get("SIFT_pred", "."),
            self.info_dict.get("Polyphen2_HDIV_pred", "."),
            self.info_dict.get("LRT_pred", "."),
            self.info_dict.get("MutationTaster_pred", "."),
            self.info_dict.get("FATHMM_pred", "."),
            self.info_dict.get("PROVEAN_pred", "."),
            self.info_dict.get("M-CAP_pred", "."),
            self.info_dict.get("PrimateAI_pred", "."),
            self.info_dict.get("REVEL_score", "."),
            self.info_dict.get("REVEL_rankscore", "."),
            self.snv_depth, self.dp, self.vd, self.ald_f, self.ald_r, self.ald,
            self.rd, self.somatic_SBF, self.somatic_ODDRATIO, self.gt, self.status,
            self.snv_depth_control, self.germline_dp, self.germline_vd,
            self.germline_ald_f, self.germline_ald_r, self.germline_ald,
            self.germline_rd, self.germline_SBF, self.germline_ODDRATIO,
            self.germline_af_string, self.germline_gt, self.pvalue,
            self.gene_category, self.mrd_label, self.new_tx, self.refgene_exon_n,
            self.new_gene, self.new_exon, self.new_cdna, new_amid_clean,
            self.info_dict.get("DCS", "."), self.info_dict.get("SCS", "."),
            self.af_string, "",
        ]]) + "\n"

    @property
    def is_paired(self):
        return len(self.arr) == 11

    @property
    def pass_common_rule(self):
        return (self.rule_1 and self.rule_2 and self.rule_3 and self.rule_4 and self.rule_7) \
            or (self.rule_1 and self.rule_6 and self.rule_7)

    @property
    def is_somatic_filter(self):
        if "MRD" in self.product_gene:
            if self.is_somatic:
                return self.mrd_filter()
        else:
            if self.pass_common_rule:
                if self.is_somatic:
                    if (self.somatic_rule_1 and self.somatic_rule_2 and self.somatic_rule_3
                            and (self.somatic_rule_4_1 or self.somatic_rule_4_2)) or self.rule_6:
                        return True
        return False

    @property
    def is_germline_filter(self):
        if "MRD" in self.product_gene:
            if self.is_germline:
                return self.mrd_filter(germline_filter=True)
        else:
            if self.pass_common_rule:
                if self.is_germline:
                    if (self.germline_rule_1 and self.germline_rule_2 and self.germline_rule_3) or self.rule_6:
                        return True
        return False

    @property
    def somatic_row(self):
        if self.is_somatic_filter:
            return self.row.strip("\n") + "\t" + self.comment
        return ""

    @property
    def germline_row(self):
        if self.is_germline_filter:
            return self.row.strip("\n") + "\t" + self.germline_comment
        return ""

    @property
    def somatic_row_all(self):
        return self.row.strip("\n") + "\t" + self.comment

    @property
    def germline_row_all(self):
        return self.row.strip("\n") + "\t" + self.germline_comment

    # ---- debug output --------------------------------------------------------

    @property
    def common_debug_header(self):
        return self.header.replace("\n", "\t") + "\t".join(["rule_1", "rule_2", "rule_3", "rule_4"])

    @property
    def debug_somatic_header(self):
        return self.common_debug_header + "\t" + "\t".join([
            "somatic_rule_1", "somatic_rule_2_1", "somatic_rule_2_2",
            "somatic_rule_3", "fanal_rule",
        ]) + "\n"

    @property
    def debug_germline_header(self):
        return self.common_debug_header + "\t" + "\t".join([
            "germline_rule_1", "germline_rule_2", "final_rule",
        ]) + "\n"

    @property
    def common_debug_row(self):
        new_row = self.row.replace("\n", "")
        return f"{new_row}\t{self.rule_1}\t{self.rule_2}\t{self.rule_3}\t{self.rule_4}"

    @property
    def debug_somatic_row(self):
        return f"{self.common_debug_row}\t{self.somatic_rule_1}\t{self.somatic_rule_2}\t{self.somatic_rule_2}\t{self.somatic_rule_3}\t{self.is_somatic_filter}\n"

    @property
    def debug_germline_row(self):
        return f"{self.common_debug_row}\t{self.germline_rule_1}\t{self.germline_rule_2}\t{self.is_germline_filter}\n"


# ---- standalone helper -------------------------------------------------------

def _match_gene_tx(changes, gene_tx_dict, default=False):
    """Match gene transcript from changes list (module-level for tests)."""
    from .utils import match_gene_tx
    return match_gene_tx(changes, gene_tx_dict, default)
