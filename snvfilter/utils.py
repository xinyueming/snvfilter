"""Utility functions for data loading and parsing."""

import os
import re
import zipfile

import pandas as pd


def read_key_site(file):
    """Parse key site TSV, return (key_gene_tx, key_gene_changes)."""
    df = pd.read_csv(file, delimiter="\t")
    key_gene_tx = {}
    key_gene_changes = {}
    for row in df.iloc:
        resistance = row.get("pathogenicity")
        gene = row.get("gene")
        tx_id = row.get("transcript_19")
        key_gene_tx.setdefault(gene, tx_id)
        if resistance not in ["Oncogenic"]:
            continue
        aa_change = row.get("aa_change")
        if not key_gene_changes.get(gene):
            key_gene_changes[gene] = [aa_change]
        else:
            key_gene_changes[gene] += [aa_change]
    return key_gene_tx, key_gene_changes


def read_product_gene_dict(file):
    """Read a plain-text gene list, one gene per line."""
    gene_dict = {}
    with open(file) as fh:
        for line in fh:
            gene_symbol = line.strip()
            if gene_symbol:
                gene_dict[gene_symbol] = 1
    return gene_dict


def read_hotspot_gene_dict(file):
    """Read hotspot gene TSV, keyed as 'gene|aa_change'."""
    hotspot_gene_dict = {}
    df = pd.read_csv(file, sep="\t")
    for row in df.iloc:
        gene = row.get("Gene.refGeneWithVer")
        aa_change = row.get("HGVSp3 <VEP>")
        if not aa_change or str(aa_change) == "nan" or str(aa_change) == "-":
            aa_change = row.get("HGVSc2 <VEP>")
        hotspot_gene_dict[f"{gene}|{aa_change}"] = 1
    return hotspot_gene_dict


def read_gene_trans(file):
    """Read gene-transcript mapping TSV."""
    key_supplement_gene_tx = {}
    if file:
        df = pd.read_csv(file, delimiter="\t", header=None)
        for row in df.iloc:
            gene = row.get(0)
            tx = row.get(1)
            if not key_supplement_gene_tx.get(gene):
                key_supplement_gene_tx.setdefault(gene, tx)
    return key_supplement_gene_tx


def read_refgene(file):
    """Read refGene file, return {transcript_id: exon_count}."""
    refgene_exon_dict = {}
    with open(file, "r") as reffile:
        for line in reffile.readlines():
            lines = re.split("[\t]+", line.strip())
            if lines[0][0] == "#":
                pass
            else:
                trans_na = lines[1].strip()
                trans_id = re.split("[.]+", trans_na)[0]
                exon_n = lines[8].strip()
                refgene_exon_dict[trans_id] = int(exon_n)
    return refgene_exon_dict


def read_varannovar(file):
    """Read modified varannovar TSV."""
    varannovar = {}
    if file:
        with open(file, "r", encoding="utf-8") as varfile:
            for line in varfile.readlines()[1:]:
                lines = re.split("[\t]+", line.strip())
                chr_na = lines[0].strip()
                pos_na = lines[1].strip()
                ref_na = lines[2].strip()
                alt_na = lines[3].strip()
                key_info = f"{chr_na}_{pos_na}_{ref_na}_{alt_na}"
                genena = lines[4].strip()
                func = lines[5].strip()
                exonicfunc = lines[6].strip()
                aaChange = lines[7].strip()
                genedetails = lines[8].strip()
                tup = (genena, func, exonicfunc, aaChange, genedetails)
                varannovar[key_info] = tup
    return varannovar


def read_gtf(f_gtf):
    """Read GTF file, extract {gene_name: gene_id}."""
    geneid = {}
    with open(f_gtf, "r") as f2:
        gi = f2.readlines()
        for line in gi:
            lines = re.split("[\t]+", line.strip())
            if lines[0][0] == "#":
                pass
            else:
                anno = lines[8].strip().split(";")
                if lines[2].strip() == "gene" or lines[2].strip() == "pseudogene":
                    gene_na = anno[0].strip().split("ID=gene-")[1]
                    gene_id = re.split("[:,]+", anno[1].strip())[1]
                    geneid[gene_na] = gene_id
    return geneid


def read_depthtsv(filterfile):
    """Read depth TSV, return {chr_pos: depth}."""
    snvdepth = {}
    if filterfile:
        cf = pd.read_table(filterfile, encoding="utf-8", encoding_errors="ignore")
        cf_list = cf.to_dict(orient="records")
        if "_depth.tsv" in filterfile:
            for li in cf_list:
                if li["#Chr"] and str(li["#Chr"]) != "nan":
                    chr_na = li["#Chr"].strip()
                    pos_na = str(int(li["Pos"])).strip()
                    raw_depth = str(int(li["Raw Depth"])).strip()
                    key_info = f"{chr_na}_{pos_na}"
                    snvdepth[key_info] = raw_depth
        else:
            for li in cf_list:
                if li["chrom"] and str(li["chrom"]) != "nan":
                    chr_na = li["chrom"].strip()
                    pos_na = str(int(li["pos"])).strip()
                    raw_depth = str(int(li["coverage"])).strip()
                    key_info = f"{chr_na}_{pos_na}"
                    snvdepth[key_info] = raw_depth
    return snvdepth


def read_mrd_filter(file):
    """Read MRD filter TSV."""
    mrd_filter_dict = {}
    if file:
        mrd_filter_dict = {
            "position": {"all": [], "general": [], "unique": []},
            "cHGVS": {"all": [], "general": [], "unique": []},
            "pHGVS": {"all": [], "general": [], "unique": []},
        }
        mrdfile = pd.read_csv(file, delimiter="\t")
        for _, row in mrdfile.iterrows():
            key_info1 = f"{row.get('Chr').strip()}:{int(row.get('Start'))}:{row.get('Ref').strip()}:{row.get('Alt').strip()}"
            key_info2 = f"{row.get('Gene').strip()}:{row.get('cHGVS').strip()}"
            key_info3 = f"{row.get('Gene').strip()}:{row.get('pHGVS').strip()}"
            mrd_filter_dict["position"]["all"].append(key_info1)
            mrd_filter_dict["cHGVS"]["all"].append(key_info2)
            mrd_filter_dict["pHGVS"]["all"].append(key_info3)
            type_val = row.get("Type").strip()
            if type_val in mrd_filter_dict["position"]:
                mrd_filter_dict["position"][type_val].append(key_info1)
                mrd_filter_dict["cHGVS"][type_val].append(key_info2)
                mrd_filter_dict["pHGVS"][type_val].append(key_info3)
    return mrd_filter_dict


def load_resource_gene_categories(gene_categories_file):
    """Load gene category Excel workbook into nested dicts."""
    genetic_gene_df = pd.read_excel(gene_categories_file, engine="openpyxl", sheet_name="遗传整理")
    drug_gene_df = pd.read_excel(gene_categories_file, engine="openpyxl", sheet_name="用药相关基因")
    immune_positive_gene_df = pd.read_excel(gene_categories_file, engine="openpyxl", sheet_name="免疫正向基因")
    immune_negative_gene_df = pd.read_excel(gene_categories_file, engine="openpyxl", sheet_name="免疫负向基因")
    immune_hyperprogressive_gene_df = pd.read_excel(gene_categories_file, engine="openpyxl", sheet_name="免疫超进展基因")
    chemical_gene_df = pd.read_excel(gene_categories_file, engine="openpyxl", sheet_name="化疗相关基因")
    mmr_gene_df = pd.read_excel(gene_categories_file, engine="openpyxl", sheet_name="MMR相关基因")

    def _to_dict(series):
        return {g: 1 for g in series.dropna().tolist()}

    gene_categorie_dict = {
        "genetic": _to_dict(genetic_gene_df["遗传基因列表"]),
        "drug": _to_dict(drug_gene_df["用药相关基因"]),
        "immune_positive": _to_dict(immune_positive_gene_df["免疫正向相关基因"]),
        "immune_negative": _to_dict(immune_negative_gene_df["免疫负向相关基因"]),
        "immune_hyperprogressive": _to_dict(immune_hyperprogressive_gene_df["免疫超进展基因"]),
        "chemical": _to_dict(chemical_gene_df["化疗基因"]),
        "mmr": _to_dict(mmr_gene_df["MMR相关基因"]),
    }
    return gene_categorie_dict


def tr_content_to_anno_dict(content):
    """Parse local frequency annotation text into (header, anno_dict)."""
    anno_dict = {}
    all_row = content.strip().split("\n")
    header = all_row[0].split("\t")
    for row in all_row[1:]:
        row = row.split("\t")
        if row[2] == ".":
            key = tuple(row[:2] + row[3:5])
        else:
            key = tuple(row[2])
        anno_dict[key] = row[5:]
    return header, anno_dict


def get_anno_file(product_gene, local_freq_input):
    """Open annotation files: anno_db, germline_anno_db, header.

    Supports 3 input modes for local_freq_input:
      1. Comma-separated two files: <fn.txt>,<germline_fn.txt>
      2. Single file: <fn.txt>
      3. ZIP archive (legacy compatibility)
    """
    anno_dict, germline_anno_dict, header = {}, {}, []
    if not local_freq_input:
        return (anno_dict, germline_anno_dict, header)

    # Legacy ZIP mode
    if zipfile.is_zipfile(local_freq_input):
        if product_gene:
            filename = os.path.basename(product_gene)
            try:
                numt = int(filename.split(".")[1])
            except Exception:
                try:
                    numt = int(filename.split(".")[2])
                except Exception:
                    numt = filename.split(".")[1]
            g180 = ["180", "lynch", "15"]
            g680 = ["680", "105", "1081", "1100", "1008", "19", "210"]
            gWES = ["wes", "WES"]
            if str(numt) in g180:
                num = "1100"
            elif str(numt) in g680:
                num = "1100"
            elif str(numt) in gWES:
                num = "WES"
            elif "kszy_BloodTumor" in filename:
                num = "624"
            elif "hg38" in filename and "kszy_84panel" in filename:
                num = "84"
            elif "hg38" in filename and "CML206panel" in filename:
                num = "cml206"
            else:
                num = "120"
            if num in ("120", "1100", "WES", "624", "84", "cml206"):
                try:
                    with zipfile.ZipFile(local_freq_input, "r") as zip_fn:
                        fn = zip_fn.read(f"{num}.mutation_frequency.txt").decode("utf-8")
                        germline_fn = zip_fn.read(f"{num}.mutation_frequency_germline.txt").decode("utf-8")
                        header, anno_dict = tr_content_to_anno_dict(fn)
                        germlineheader, germline_anno_dict = tr_content_to_anno_dict(germline_fn)
                except KeyError as e:
                    print(f"注释文件缺失: {e}")
        return (anno_dict, germline_anno_dict, header)

    # Comma-separated or single file
    if "," in local_freq_input:
        parts = [p.strip() for p in local_freq_input.split(",")]
        fn_path = parts[0]
        germline_path = parts[1] if len(parts) > 1 else ""
    else:
        fn_path = local_freq_input
        germline_path = ""

    if fn_path and os.path.exists(fn_path):
        with open(fn_path, "r", encoding="utf-8") as f:
            header, anno_dict = tr_content_to_anno_dict(f.read())
    if germline_path and os.path.exists(germline_path):
        with open(germline_path, "r", encoding="utf-8") as f:
            germline_header, germline_anno_dict = tr_content_to_anno_dict(f.read())
    return (anno_dict, germline_anno_dict, header)


def match_gene_tx(changes: list, gene_tx_dict: dict, default=False):
    """Match gene transcript from changes list."""
    for change in changes:
        if change not in [".", "UNKNOWN"]:
            arr = change.split(":")
            gene, tx = arr[:2]
            if tx and isinstance(tx, str):
                tx = tx.split(".")[0]
            map_tx = gene_tx_dict.get(gene)
            if map_tx and isinstance(map_tx, str):
                map_tx = map_tx.split(".")[0]
            if map_tx == tx:
                return change
    if default:
        return change
    return "."


def detect_panel_flags(product_gene):
    """Detect panel flags from product_gene path."""
    basename = os.path.basename(product_gene) if product_gene else ""
    return {
        "is_brca": "brca" in basename,
        "is_lynch": "lynch" in basename,
        "is_105": "zy.105" in basename or "zy.genetics.15" in basename,
        "is_hrr": "zy.hrr" in basename,
        "is_wes": "FullRNA" in basename or "zy.wes" in basename,
    }
