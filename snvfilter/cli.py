"""CLI entry point for snvfilter."""

import argparse
import os
import sys

from . import __version__
from .config import load_config, load_rules
from .rules import Rules
from .utils import (
    detect_panel_flags,
    get_anno_file,
    read_depthtsv,
    read_gene_trans,
    read_gtf,
    read_hotspot_gene_dict,
    read_key_site,
    read_mrd_filter,
    read_product_gene_dict,
    read_refgene,
    read_varannovar,
    load_resource_gene_categories,
)
from .filters import VariantRow


def _build_cli_args_dict(args):
    """Convert argparse Namespace to a flat dict for config merging."""
    d = {}
    for key, val in vars(args).items():
        if val is not None:
            d[key] = val
    return d


def build_parser():
    parser = argparse.ArgumentParser(
        prog="snvfilter",
        description="SNV filtering tool for tumor-normal paired samples",
    )
    parser.add_argument("--version", action="version", version=f"snvfilter {__version__}")

    # Config
    parser.add_argument("--config", "-C", help="User config YAML file")
    parser.add_argument("--rules", "-R", help="Custom rules YAML file")

    # Input
    parser.add_argument("-k", "--key_site", help="Key site file (combine.tsv)")
    parser.add_argument("-l", "--product_gene", help="Product gene list file")
    parser.add_argument("-g", "--hotspot_gene", help="Hotspot gene file")
    parser.add_argument("-t", "--gene_trans", help="Gene transcript file")
    parser.add_argument("-f", "--annovar_vcf", help="ANNOVAR VCF file")
    parser.add_argument("-c", "--gene_categorie", help="Gene categories Excel file")
    parser.add_argument("-r", "--hg19_refgene", help="hg19 refGene file")
    parser.add_argument("-r2", "--hg19_gtf", help="hg19 GTF file")
    parser.add_argument("-r3", "--modify_varannovar", help="Modified varannovar file")
    parser.add_argument("-r4", "--depth_tsv", help="bamdst depth file")
    parser.add_argument("-r5", "--depth_tsv_control", help="bamdst depth control file")
    parser.add_argument("-ga", "--gene_alias", help="Gene alias file")
    parser.add_argument("-lc", "--local_freq", help="Local frequency file (zip/txt/txt,txt)")
    parser.add_argument("-mf", "--mrd_filter", help="MRD filter file")

    # Sample
    parser.add_argument("-s", "--sample_type", choices=["blood", "tissue"], default=None, help="Sample type")
    parser.add_argument("-y", "--force_germline", action="store_true", default=None, help="Force germline mode")

    # Output
    parser.add_argument("-o", "--output_prefix", help="Output prefix (default: derived from VCF)")

    # Flags
    parser.add_argument("-d", "--debug", action="store_true", default=None, help="Debug mode")

    # Sub-commands
    sub = parser.add_subparsers(dest="command")

    # show-rules
    show_rules_p = sub.add_parser("show-rules", help="Print default filtering rules")

    # init-config
    init_cfg_p = sub.add_parser("init-config", help="Generate a default config file")
    init_cfg_p.add_argument("--output", "-o", default="config.yaml", help="Output path for generated config")

    return parser


def cmd_show_rules():
    """Print the default rules to stdout."""
    rules = Rules()
    import yaml
    print(yaml.dump({
        "common": rules.common,
        "somatic": rules.somatic,
        "germline": rules.germline,
        "mrd": rules.mrd,
    }, default_flow_style=False, allow_unicode=True))


def cmd_init_config(output):
    """Copy default config to the given path."""
    import shutil
    from .config import _DEFAULT_CONFIG_PATH
    if os.path.exists(output):
        print(f"Error: {output} already exists", file=sys.stderr)
        sys.exit(1)
    shutil.copy2(_DEFAULT_CONFIG_PATH, output)
    print(f"Default config written to {output}")


def run_filtering(args):
    """Main filtering pipeline."""
    # Load config layers
    cli_dict = {k: v for k, v in vars(args).items() if v is not None}
    config = load_config(config_file=args.config, cli_args=cli_dict)

    rules = Rules(rules_file=args.rules)

    # Extract config values (CLI overrides already merged)
    inp = config.get("input", {})
    out_cfg = config.get("output", {})
    sample = config.get("sample", {})
    flags = config.get("flags", {})

    key_site = args.key_site or inp.get("key_site")
    product_gene = args.product_gene or inp.get("product_gene")
    hotspot_gene = args.hotspot_gene or inp.get("hotspot_gene")
    gene_trans = args.gene_trans or inp.get("gene_trans")
    annovar_vcf = args.annovar_vcf or inp.get("annovar_vcf")
    gene_categorie = args.gene_categorie or inp.get("gene_categorie")
    hg19_refgene = args.hg19_refgene or inp.get("hg19_refgene")
    hg19_gtf = args.hg19_gtf or inp.get("hg19_gtf")
    modify_varannovar = args.modify_varannovar or inp.get("modify_varannovar")
    depth_tsv = args.depth_tsv or inp.get("depth_tsv")
    depth_tsv_control = args.depth_tsv_control or inp.get("depth_tsv_control")
    gene_alias = args.gene_alias or inp.get("gene_alias")
    local_freq = args.local_freq or inp.get("local_freq")
    mrd_filter_file = args.mrd_filter or inp.get("mrd_filter")

    sample_type = args.sample_type or sample.get("type", "blood")
    force_germline = args.force_germline if args.force_germline is not None else sample.get("force_germline", False)
    debug = args.debug if args.debug is not None else flags.get("debug", False)
    output_prefix = args.output_prefix or out_cfg.get("prefix")

    if sample_type not in ("blood", "tissue"):
        raise ValueError(f"invalid sample type: {sample_type}")

    # Load resources
    key_gene_tx, key_gene_changes = read_key_site(key_site)
    hotspot_gene_dict = read_hotspot_gene_dict(hotspot_gene)
    product_gene_dict = read_product_gene_dict(product_gene)
    key_gene_tx_supplement = read_gene_trans(gene_trans)
    gene_categorie_dict = load_resource_gene_categories(gene_categorie)
    refgene_exon_dict = read_refgene(hg19_refgene)
    geneid = read_gtf(hg19_gtf)
    varannovar = read_varannovar(modify_varannovar)
    snvdepth = read_depthtsv(depth_tsv)
    snvdepth_control = read_depthtsv(depth_tsv_control)
    gene_alias_dict = read_gene_trans(gene_alias)
    local_anno_tuple = get_anno_file(product_gene, local_freq)
    mrd_filter_dict = read_mrd_filter(mrd_filter_file)

    # Panel flags
    panel_flags = detect_panel_flags(product_gene)
    is_brca = panel_flags["is_brca"]
    is_lynch = panel_flags["is_lynch"]
    is_105 = panel_flags["is_105"]

    # Output paths
    if not output_prefix:
        output_prefix = os.path.basename(annovar_vcf).replace(".vcf", "")

    filter_xls = f"{output_prefix}.filter.xls"
    all_xls = f"{output_prefix}.all.xls"
    germline_filter_xls = f"{output_prefix}.filter.germline.xls"
    germline_all_xls = f"{output_prefix}.all.germline.xls"

    ctx = {
        "product_gene": product_gene,
        "product_gene_dict": product_gene_dict,
        "sample_type": sample_type,
        "hotspot_gene_dict": hotspot_gene_dict,
        "key_gene_changes": key_gene_changes,
        "key_gene_tx": key_gene_tx,
        "key_gene_tx_supplement": key_gene_tx_supplement,
        "gene_categorie_dict": gene_categorie_dict,
        "refgene_exon_dict": refgene_exon_dict,
        "geneid": geneid,
        "varannovar": varannovar,
        "snvdepth": snvdepth,
        "snvdepth_control": snvdepth_control,
        "gene_alias_dict": gene_alias_dict,
        "local_anno_tuple": local_anno_tuple,
        "mrd_filter_dict": mrd_filter_dict,
        "force_germline": force_germline,
        "is_brca": is_brca,
        "is_lynch": is_lynch,
        "is_105": is_105,
    }

    with open(annovar_vcf) as fh, \
         open(filter_xls, "w") as fh_o_filter, \
         open(all_xls, "w") as fh_o_all, \
         open(germline_filter_xls, "w") as fh_o_germline_filter, \
         open(germline_all_xls, "w") as fh_o_germline_all:

        write_header = False
        for line in fh:
            if line.startswith("#"):
                continue

            v = VariantRow(line, rules, **ctx)

            if not write_header:
                hdr = v.debug_somatic_header if debug else v.header
                for f in [fh_o_filter, fh_o_all, fh_o_germline_filter, fh_o_germline_all]:
                    f.write(hdr)
                write_header = True

            if debug:
                if v.somatic_row:
                    fh_o_filter.write(v.debug_somatic_row)
                if v.germline_row:
                    fh_o_germline_filter.write(v.debug_germline_row)
                if v.is_somatic:
                    fh_o_all.write(v.debug_somatic_row)
                if v.is_germline:
                    fh_o_germline_all.write(v.debug_germline_row)
            else:
                if v.somatic_row:
                    fh_o_filter.write(v.somatic_row)
                if v.germline_row:
                    fh_o_germline_filter.write(v.germline_row)
                if v.is_somatic:
                    fh_o_all.write(v.somatic_row_all)
                if v.is_germline:
                    fh_o_germline_all.write(v.germline_row_all)

    print(f"Filtering done. Output prefix: {output_prefix}")


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "show-rules":
        cmd_show_rules()
        return
    if args.command == "init-config":
        cmd_init_config(args.output)
        return

    # Check required inputs for run mode
    required = ["key_site", "product_gene", "hotspot_gene", "gene_trans", "annovar_vcf"]
    merged = load_config(config_file=args.config)
    inp = merged.get("input", {})
    sample = merged.get("sample", {})
    flags = merged.get("flags", {})

    missing = []
    for r in required:
        val = getattr(args, r, None) or inp.get(r)
        if not val:
            missing.append(r)
    if missing:
        parser.error(f"missing required arguments: {', '.join(missing)}")

    run_filtering(args)


if __name__ == "__main__":
    main()
