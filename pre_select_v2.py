# keep core site
# SNV 过滤脚本
# eg: python3 pre_select.py -k ../hotspot/combine.tsv -l ../hotspot/zy.120.gene.list -t ../hotspot/sorted.gene.tx.txt -f ../1206-D3.var.hg19_multianno.vcf -o t.tsv
# eg: python3.8 .\docker\annovar\scripts\pre_select.py -k hotspot/combine.tsv -l hotspot/zy.120.gene.list -t hotspot/sorted.gene.tx.txt -f .\docker\annovar\scripts\test\H1E442G102G0.var.hg19_multianno.vcf -o .\docker\annovar\scripts\test\H1E442G102G0.var.hg19_multianno.filter.xls
# eg: python3 pre_select_v2.py -k hotspot/combine.tsv -l hotspot/zy.120.gene.list -t hotspot/sorted.gene.tx.txt -f H1E442G102G0.var.hg19_multianno.vcf -o .\docker\annovar\scripts\test\H1E442G102G0.var.hg19_multianno.filter.xls

import re
import os
import json
import argparse
import pandas as pd
import csv
import zipfile


# 标记基因类别
# 资源：/data/hotspot/tumor-gene-20230216.xlsx
def load_resource_gene_categories(gene_categories_file):
    genetic_gene_df = pd.read_excel(gene_categories_file, engine='openpyxl', sheet_name='遗传整理')
    drug_gene_df = pd.read_excel(gene_categories_file, engine='openpyxl', sheet_name='用药相关基因')
    immune_positive_gene_df = pd.read_excel(gene_categories_file, engine='openpyxl', sheet_name='免疫正向基因')
    immune_negative_gene_df = pd.read_excel(gene_categories_file, engine='openpyxl', sheet_name='免疫负向基因')
    immune_hyperprogressive_gene_df = pd.read_excel(gene_categories_file, engine='openpyxl', sheet_name='免疫超进展基因')
    chemical_gene_df = pd.read_excel(gene_categories_file, engine='openpyxl', sheet_name='化疗相关基因')
    mmr_gene_df = pd.read_excel(gene_categories_file, engine='openpyxl', sheet_name='MMR相关基因')

    genetic_gene_list = genetic_gene_df['遗传基因列表'].to_list()
    drug_gene_list = drug_gene_df['用药相关基因'].to_list()
    immune_positive_gene_list = immune_positive_gene_df['免疫正向相关基因'].to_list()
    immune_negative_gene_list = immune_negative_gene_df['免疫负向相关基因'].to_list()
    immune_hyperprogressive_gene_list = immune_hyperprogressive_gene_df['免疫超进展基因'].to_list()
    chemical_gene_list = chemical_gene_df['化疗基因'].to_list()
    mmr_gene_list = mmr_gene_df['MMR相关基因'].to_list()

    genetic_gene_dict = dict()
    for gene in genetic_gene_list:
        genetic_gene_dict[gene] = 1

    drug_gene_dict = dict()
    for gene in drug_gene_list:
        drug_gene_dict[gene] = 1

    immune_positive_gene_dict = dict()
    for gene in immune_positive_gene_list:
        immune_positive_gene_dict[gene] = 1

    immune_negative_gene_dict = dict()
    for gene in immune_negative_gene_list:
        immune_negative_gene_dict[gene] = 1

    immune_hyperprogressive_gene_dict = dict()
    for gene in immune_hyperprogressive_gene_list:
        immune_hyperprogressive_gene_dict[gene] = 1

    chemical_gene_dict = dict()
    for gene in chemical_gene_list:
        chemical_gene_dict[gene] = 1

    mmr_gene_dict = dict()
    for gene in mmr_gene_list:
        mmr_gene_dict[gene] = 1

    gene_categorie_dict = {
        'genetic': genetic_gene_dict,
        'drug': drug_gene_dict,
        'immune_positive': immune_positive_gene_dict,
        'immune_negative': immune_negative_gene_dict,
        'immune_hyperprogressive': immune_hyperprogressive_gene_dict,
        'chemical': chemical_gene_dict,
        'mmr': mmr_gene_dict,
    }

    return gene_categorie_dict


def parse_geneinfo(geneinfo):
    part1, part2 = geneinfo.split('|')
    strand, position = part2[0], part2[1:]
    gene_transcript, func_name, func_idx = part1.split(':')
    func = f'{func_name}{func_idx}'
    gene, transcript = gene_transcript.split('_')
    return gene, transcript, func, position


def read_key_site(file):
    df = pd.read_csv(file, delimiter='\t')
    key_gene_tx = dict()
    key_gene_changes = dict()
    for row in df.iloc:
        resistance = row.get('pathogenicity')
        gene = row.get('gene')
        # entrez_id = row.get('entrez_id')
        tx_id = row.get('transcript_19')
        key_gene_tx.setdefault(gene, tx_id)  # 说明：此处必须用第一个赋值，因为有重复行，且可能出现空覆原来值，MSH2 基因。
        if resistance not in ['Oncogenic']:  # 特例： chr2	29445213	A	T 注释为ALK:NM_004304.5:exon22:c.T3512A:p.I1171N。在关键位点库中存在，且为Likely Oncogenic，但是没有报道
            # 特例2:
            continue
        aa_change = row.get('aa_change')
        if not key_gene_changes.get(gene):
            key_gene_changes[gene] = [aa_change]
        else:
            key_gene_changes[gene] += [aa_change]
    return key_gene_tx, key_gene_changes


def read_product_gene_dict(file):
    gene_dict = dict()
    with open(file) as fh:
        for line in fh:
            gene_symbol = line.strip()
            if gene_symbol:
                gene_dict[gene_symbol] = 1
    return gene_dict


def read_hotspot_gene_dict(file):
    hotspot_gene_dict = dict()
    df = pd.read_csv(file, sep='\t')
    for row in df.iloc:
        gene = row.get('Gene.refGeneWithVer')
        aa_change = row.get('HGVSp3 <VEP>')
        if not aa_change or str(aa_change) == 'nan' or str(aa_change) == '-':
            aa_change = row.get('HGVSc2 <VEP>')
        hotspot_gene_dict[f'{gene}|{aa_change}'] = 1
    return hotspot_gene_dict


def read_gene_trans(file):
    key_supplement_gene_tx = dict()
    if file:
        df = pd.read_csv(file, delimiter='\t', header=None)
        for row in df.iloc:
            gene = row.get(0)
            tx = row.get(1)
            if not key_supplement_gene_tx.get(gene):
                key_supplement_gene_tx.setdefault(gene, tx)
    return key_supplement_gene_tx


def match_gene_tx(changes: list, gene_tx_dict: dict, default=False):
    for change in changes:
        if change not in ['.', 'UNKNOWN']:
            arr = change.split(':')
            gene, tx = arr[:2]
            if tx and isinstance(tx,str):
                tx = tx.split('.')[0]
            map_tx = gene_tx_dict.get(gene)
            if map_tx and isinstance(map_tx,str):
                map_tx = map_tx.split('.')[0]
            if map_tx == tx:
                return change
    if default:
        return change
    return '.'

def read_refgene(file):
    refgene_exon_dict = dict()
    with open(file,'r') as reffile:
        for line in reffile.readlines():
            lines = re.split('[\t]+', line.strip())
            if lines[0][0] == '#':
                pass
            else:
                trans_na = lines[1].strip()
                trans_id = re.split('[.]+', trans_na)[0]
                exon_n = lines[8].strip()
                refgene_exon_dict[trans_id] = int(exon_n)
    return refgene_exon_dict

def read_varannovar(file):
    varannovar = dict()
    if file:
        with open(file,'r',encoding='utf-8') as varfile:
            for line in varfile.readlines()[1:]:
                lines = re.split('[\t]+', line.strip())
                chr_na = lines[0].strip()
                pos_na = lines[1].strip()
                ref_na = lines[2].strip()
                alt_na = lines[3].strip()
                key_info = f'{chr_na}_{pos_na}_{ref_na}_{alt_na}'
                genena = lines[4].strip()
                func = lines[5].strip()
                exonicfunc = lines[6].strip()
                aaChange = lines[7].strip()
                genedetails = lines[8].strip()
                tup = (genena, func, exonicfunc, aaChange, genedetails)
                varannovar[key_info] = tup
    return varannovar

def read_gtf(f_gtf):
    '''读取基因组注释文件信息，提取geneid'''
    geneid = {}
    with open(f_gtf,'r') as f2:
        gi = f2.readlines()
        for line in gi:
            lines = re.split('[\t]+',line.strip())
            if lines[0][0] == '#':
                pass
            else:
                anno = lines[8].strip().split(';')
                if lines[2].strip() == 'gene' or lines[2].strip() == 'pseudogene':
                    gene_na = anno[0].strip().split('ID=gene-')[1]
                    gene_id = re.split('[:,]+',anno[1].strip())[1]
                    geneid[gene_na] = gene_id
    return geneid

def read_depthtsv(filterfile):
    snvdepth = dict()
    if filterfile:
        cf = pd.read_table(filterfile, encoding="utf-8", encoding_errors="ignore")
        cf_list = cf.to_dict(orient="records")
        if '_depth.tsv' in filterfile:
            for li in cf_list:
                if li["#Chr"] and str(li["#Chr"]) != "nan" :
                    chr_na = li["#Chr"].strip()
                    pos_na = str(int(li["Pos"])).strip()
                    raw_depth = str(int(li["Raw Depth"])).strip()
                    key_info = f'{chr_na}_{pos_na}'
                    snvdepth[key_info] = raw_depth
        else:
            for li in cf_list:
                if li["chrom"] and str(li["chrom"]) != "nan" :
                    chr_na = li["chrom"].strip()
                    pos_na = str(int(li["pos"])).strip()
                    raw_depth = str(int(li["coverage"])).strip()
                    key_info = f'{chr_na}_{pos_na}'
                    snvdepth[key_info] = raw_depth
    return snvdepth

def extract_exon_transvar(exoninfo):
    num_exon = re.findall('\d+',exoninfo)
    if len(num_exon) == 2 and 'intron' in exoninfo:
        return f'intron{num_exon[0]}'
    if len(num_exon) == 1 and 'cds' in exoninfo:
        return f'exon{num_exon[0]}'
    return exoninfo

def add_transvar_result(file_in, file_out, f_transvar, tx_gene):
    # 预读取 transvar 文件到 DataFrame
    transvar_df = pd.read_csv(f_transvar, sep='\t', header=None)
    transvar_df.columns = [f'col{i}' for i in range(transvar_df.shape[1])]

    def find_transvar_row(tmp_input, tsx):
        # 先按 tmp_input 匹配，再按 tsx 匹配
        df = transvar_df[transvar_df['col0'] == tmp_input]
        if tsx != '.':
            df = df[df['col1'].str.startswith(tsx)]
        else:
            # 如果 tsx 是 '.', 只取第一个匹配的结果
            df = df.head(1)
        if not df.empty:
            return df.iloc[0].tolist()
        return None
    
    xls_file = file_in
    out_file = file_out
    print(out_file)
    fho = open(out_file, 'w')
    no_change = ['p.K373Rfs*49']  # TP53  NM_000546.5	c.1113del	p.K373Rfs*105
    with open(xls_file) as fh:
        header = fh.readline().strip('\n')
        header_arr = header.split('\t')
        idx_new_gene, idx_new_exon, idx_new_transcript, idx_new_cdna, idx_new_amid, idx_tmp_transvar_input = header_arr.index('new_gene'), header_arr.index('new_exon'), header_arr.index(
            'new_transcript'), header_arr.index('new_cdna'), header_arr.index('new_amid'), header_arr.index(
            'tmp_transvar_input')
        fho.write('\t'.join(header_arr[:idx_tmp_transvar_input] + header_arr[idx_tmp_transvar_input + 1:]) + '\n')
        for line in fh:
            tsaa,tstx,tsexon = '','.','.'
            arr = line.strip('\n').split('\t')
            value_new_gene, value_new_exon, value_new_transcript, val_new_cdna, val_new_amid, val_tmp_transvar_input = arr[idx_new_gene], arr[idx_new_exon], arr[idx_new_transcript], arr[
                idx_new_cdna], arr[idx_new_amid], arr[idx_tmp_transvar_input]
            tsxq = value_new_transcript
            tsx = tsxq.split('.')[0]
            if not tsx:
                tsxq = tx_gene.get(value_new_gene,'.')
                tsx = tsxq.split('.')[0]
            if val_tmp_transvar_input:
                try:
                    #cmd = f'grep {val_tmp_transvar_input} {f_transvar} |grep {tsx} |cut -f5'
                    if tsx:
                        # cmd = "awk '{if($1 == val_tmp_transvar_input666)print}' f_transvar666 |grep tsx666"
                        # cmd = cmd.replace('val_tmp_transvar_input666', f'"{val_tmp_transvar_input}"').replace('f_transvar666', f_transvar).replace('tsx666', tsx)
                        # tsaa_list = os.popen(cmd).read().strip().split('\t')
                        tsaa_list = find_transvar_row(val_tmp_transvar_input, tsx)
                        tsaa = tsaa_list[4] if tsaa_list and len(tsaa_list) > 4 else ''
                        tstx = tsaa_list[1].split(' ')[0] if tsaa_list and len(tsaa_list) > 1 else '.'
                        exoninfo = tsaa_list[5] if tsaa_list and len(tsaa_list) > 5 else '.'
                        tsexon = extract_exon_transvar(exoninfo)
                    if not tsaa and value_new_transcript == '.':
                        # cmd = "awk '{if($1 == val_tmp_transvar_input666)print}' f_transvar666 |head -n 1"
                        # cmd = cmd.replace('val_tmp_transvar_input666', f'"{val_tmp_transvar_input}"').replace('f_transvar666', f_transvar)
                        # tsaa_list = os.popen(cmd).read().strip().split('\t')
                        tsaa_list = find_transvar_row(val_tmp_transvar_input, '.')
                        tsaa = tsaa_list[4]
                        tstx = tsaa_list[1].split(' ')[0]
                        exoninfo = tsaa_list[5]
                        tsexon = extract_exon_transvar(exoninfo)
                # 输出报错信息
                except Exception as e:
                    print(f'Error transvar processing line: {tsaa_list} with error: {e}')
                    tsaa = ''
                if tsaa:
                    tsaas = tsaa.strip().split('/')
                    cts = tsaas[1]
                    if 'del' in cts and 'ins' not in cts:
                        cts = cts.split('del')[0] + 'del'
                    elif 'dup' in cts:
                        cts = cts.split('dup')[0] + 'dup'
                    pts = tsaas[2]
                    if 'del' in pts and 'ins' not in pts:
                        pts = pts.split('del')[0] + 'del'
                    elif 'dup' in pts:
                        pts = pts.split('dup')[0] + 'dup'
                    new_cdna, new_amid = cts, pts
                    if new_cdna != '.':
                        arr[idx_new_cdna] = new_cdna
                        if value_new_transcript == '.':
                            arr[idx_new_transcript] = tstx
                            arr[idx_new_exon] = tsexon
                    if new_amid != '.':
                        arr[idx_new_amid] = [pi.replace('*fs*1', '*') if pi.endswith('*fs*1') else pi.replace('fs*1', '*') if pi.endswith('fs*1') else pi for pi in [new_amid]][0]
            fho.write('\t'.join(arr[:idx_tmp_transvar_input] + arr[idx_tmp_transvar_input + 1:]) + '\n')
    fho.close()

def tr_content_to_anno_dict(content):
    """ 把snv本地频率注释文件转注释字典
    """
    anno_dict = {}  # 创建一个字典来存储注释文件中的四列信息
    all_row = content.strip().split('\n')
    header = all_row[0].split('\t')  # 第一行表头
    for row in all_row[1:]:
        row = row.split('\t')
        if row[2] == '.':
            key = tuple(row[:2] + row[3:5])  # 使用四列作为键
        else:
            key = tuple(row[2])
        anno_dict[key] = row[5:]
    return header, anno_dict

def get_anno_file(product_gene, local_freq_input):
    """ 打开注释文件 anno_db、anno_germline_db、注释表头。
    支持 3 种 -lc 输入：
      1. 逗号分隔的两个直接文件：<fn.txt>,<germline_fn.txt>
      2. 单个直接文件：<fn.txt>（无胚系频率，germline_anno_dict 为空）
      3. zip 压缩包（旧版兼容）：按 product_gene 推断 panel 号提取
    """
    anno_dict, germline_anno_dict, header = {}, {}, []
    if not local_freq_input:
        return (anno_dict, germline_anno_dict, header)

    # ── 情况 3：zip 兼容（原逻辑） ──
    if zipfile.is_zipfile(local_freq_input):
        if product_gene:
            filename = os.path.basename(product_gene)
            try:
                numt = int(filename.split('.')[1])
            except:
                try:
                    numt = int(filename.split('.')[2])
                except:
                    numt = filename.split('.')[1]
            g180 = ['180', 'lynch', '15']
            g680 = ['680', '105', '1081', '1100', '1008', '19', '210']
            gWES = ['wes', 'WES']
            if str(numt) in g180:
                num = '1100'
            elif str(numt) in g680:
                num = '1100'
            elif str(numt) in gWES:
                num = 'WES'
            elif 'kszy_BloodTumor' in filename:
                num = '624'
            elif 'hg38' in filename and 'kszy_84panel' in filename:
                num = '84'
            elif 'hg38' in filename and 'CML206panel' in filename:
                num = 'cml206'
            else:
                num = '120'
            if num == '120' or num == '1100' or num == 'WES' or num == '624' or num == '84' or num == 'cml206':
                try:
                    with zipfile.ZipFile(local_freq_input, 'r') as zip_fn:
                        fn = zip_fn.read(f'{num}.mutation_frequency.txt').decode('utf-8')
                        germline_fn = zip_fn.read(f'{num}.mutation_frequency_germline.txt').decode('utf-8')
                        header, anno_dict = tr_content_to_anno_dict(fn)
                        germline_header, germline_anno_dict = tr_content_to_anno_dict(germline_fn)
                except KeyError as e:
                    print(f"注释文件缺失: {e}")
        return (anno_dict, germline_anno_dict, header)

    # ── 情况 1：逗号分隔的两个直接文件 ──
    if ',' in local_freq_input:
        parts = [p.strip() for p in local_freq_input.split(',')]
        fn_path = parts[0]
        germline_path = parts[1] if len(parts) > 1 else ''
    # ── 情况 2：单个直接文件 ──
    else:
        fn_path = local_freq_input
        germline_path = ''

    if fn_path and os.path.exists(fn_path):
        with open(fn_path, 'r', encoding='utf-8') as f:
            header, anno_dict = tr_content_to_anno_dict(f.read())
    if germline_path and os.path.exists(germline_path):
        with open(germline_path, 'r', encoding='utf-8') as f:
            germline_header, germline_anno_dict = tr_content_to_anno_dict(f.read())
    return (anno_dict, germline_anno_dict, header)

def read_mrd_filter(file):
    mrd_filter_dict = dict()
    if file:
        # 初始化字典结构
        mrd_filter_dict = {
            "position": {"all": [], "general": [], "unique": []},
            "cHGVS": {"all": [], "general": [], "unique": []},
            "pHGVS": {"all": [], "general": [], "unique": []}
        }
        mrdfile = pd.read_csv(file, delimiter='\t')
        for _, row in mrdfile.iterrows():
            key_info1 = f"{row.get('Chr').strip()}:{int(row.get('Start'))}:{row.get('Ref').strip()}:{row.get('Alt').strip()}"
            key_info2 = f"{row.get('Gene').strip()}:{row.get('cHGVS').strip()}"
            key_info3 = f"{row.get('Gene').strip()}:{row.get('pHGVS').strip()}"
            mrd_filter_dict["position"]["all"].append(key_info1)
            mrd_filter_dict["cHGVS"]["all"].append(key_info2)
            mrd_filter_dict["pHGVS"]["all"].append(key_info3)
            type_val = row.get('Type').strip()
            if type_val in mrd_filter_dict["position"]:
                mrd_filter_dict["position"][type_val].append(key_info1)
                mrd_filter_dict["cHGVS"][type_val].append(key_info2)
                mrd_filter_dict["pHGVS"][type_val].append(key_info3)
    return mrd_filter_dict

class VariantRow:
    # annovar row, each line is a variant

    def __init__(self, row: str, product_gene, product_gene_dict, sample_type, hotspot_gene_dict, key_gene_changes, key_gene_tx, key_gene_tx_supplement, gene_categorie_dict, refgene_exon_dict, geneid, varannovar, snvdepth, snvdepth_control, gene_alias_dict, local_anno_tuple, mrd_filter_dict, force_germline, is_brca, is_lynch, is_105, debug=False):
        # extra param
        self.force_germline = force_germline
        self.is_brca = is_brca
        self.is_lynch = is_lynch
        self.is_105 = is_105

        self.product_gene = product_gene
        self.product_gene_dict = product_gene_dict
        self.sample_type = sample_type
        if self.is_brca and self.sample_type == 'blood':
            self.force_germline = True
        if self.is_lynch:
            self.force_germline = True
        if self.is_105:
            self.force_germline = True
        if 'zy.hrr' in self.product_gene:
            self.hrr_panel = True
        else:
            self.hrr_panel = False
        self.hotspot_gene_dict = hotspot_gene_dict
        self.gene_categorie_dict = gene_categorie_dict
        self.refgene_exon_dict = refgene_exon_dict
        self.gene_alias_dict = gene_alias_dict

        self.info_dict = dict()
        self.info_dict2 = dict()
        self.arr = row.split('\t')
        self.alt = self.arr[4]
        self.filter_snv = self.arr[6]
        self.info = self.arr[7]
        self.parser_info(self.info)
        self.parser_info2(self.info)
        self.key_infos = f'{self.arr[0]}_{self.arr[1]}_{self.arr[3]}_{self.arr[4]}'
        self.varannovar = varannovar
        if self.key_infos in self.varannovar:
            self.info_dict['Gene.refGeneWithVer'] = self.varannovar[self.key_infos][0]
            self.info_dict['Func.refGeneWithVer'] = self.varannovar[self.key_infos][1]
            self.info_dict['ExonicFunc.refGeneWithVer'] = self.varannovar[self.key_infos][2]
            self.info_dict['AAChange.refGeneWithVer'] = self.varannovar[self.key_infos][3]
            self.info_dict['gene_details'] = self.varannovar[self.key_infos][4]

        self.snvdepth = snvdepth
        self.snvdepth_control = snvdepth_control
        self.dp = int(self.info_dict['DP'])
        self.vd = int(self.info_dict['VD'])
        self.af = float(self.info_dict['AF'])
        self.somatic_SBF = '.'
        self.somatic_ODDRATIO = '.'
        try:
            self.pvalue = float(self.info_dict['SSF'])
        except:
            self.pvalue = 'NA'
            self.somatic_SBF = self.info_dict.get('SBF', '.')
            self.somatic_ODDRATIO = self.info_dict.get('ODDRATIO', '.')

        self.info_prefix = self.arr[8]
        self.info_sample_1 = self.arr[9]

        self.info_sample_dict_1 = dict()
        for k, v in zip(self.info_prefix.split(':'), self.info_sample_1.split(':')):
            self.info_sample_dict_1[k] = v
        self.ald = self.info_sample_dict_1.get('ALD', '.').strip()
        self.ald_f = self.ald.split(',')[0].strip()
        self.ald_r = self.ald.split(',')[-1].strip()
        self.rd = self.info_sample_dict_1.get('RD', '.').strip()
        self.rd_f = self.rd.split(',')[0].strip()
        self.rd_r = self.rd.split(',')[-1].strip()
        self.gt = self.info_sample_dict_1.get('GT', '.').strip().replace('/', '|')
        try:
            if int(self.ald_f) + int(self.ald_r) > 0:
                self.ald = self.ald + ';' + str(round(int(self.ald_f)/(int(self.ald_f) + int(self.ald_r)), 2))
            else:
                self.ald = self.ald + ';' + '.'
            if int(self.rd_f) + int(self.rd_r) > 0:
                self.rd = self.rd + ';' + str(round(int(self.rd_f)/(int(self.rd_f) + int(self.rd_r)), 2))
            else:
                self.rd = self.rd + ';' + '.'
        except:
            pass

        self.germline_dp = '.'
        self.germline_vd = '.'
        self.germline_af = '.'
        self.germline_ald = '.'
        self.germline_ald_f = '.'
        self.germline_ald_r = '.'
        self.germline_rd = '.'
        self.germline_SBF = '.'
        self.germline_ODDRATIO = '.'
        self.germline_gt = '.'
        if self.is_pair:
            self.info_sample_2 = self.arr[10]
            self.info_sample_dict_2 = dict()
            for k, v in zip(self.info_prefix.split(':'), self.info_sample_2.split(':')):
                self.info_sample_dict_2[k] = v
            self.germline_ald = self.info_sample_dict_2.get('ALD', '.').strip()
            self.germline_ald_f = self.germline_ald.split(',')[0].strip()
            self.germline_ald_r = self.germline_ald.split(',')[-1].strip()
            self.germline_rd = self.info_sample_dict_2.get('RD', '.').strip()
            self.germline_rd_f = self.germline_rd.split(',')[0].strip()
            self.germline_rd_r = self.germline_rd.split(',')[-1].strip()
            self.germline_gt = self.info_sample_dict_2.get('GT', '.').strip().replace('/', '|')
            try:
                if int(self.germline_ald_f) + int(self.germline_ald_r) > 0:
                    self.germline_ald = self.germline_ald + ';' + str(round(int(self.germline_ald_f)/(int(self.germline_ald_f) + int(self.germline_ald_r)), 2))
                else:
                    self.germline_ald = self.germline_ald + ';' + '.'
                if int(self.germline_rd_f) + int(self.germline_rd_r) > 0:
                    self.germline_rd = self.germline_rd + ';' + str(round(int(self.germline_rd_f)/(int(self.germline_rd_f) + int(self.germline_rd_r)), 2))
                else:
                    self.germline_rd = self.germline_rd + ';' + '.'
            except:
                pass
            self.somatic_SBF = self.info_sample_dict_1.get('SBF', '.').strip()
            self.somatic_ODDRATIO = self.info_sample_dict_1.get('ODDRATIO', '.').strip()
            self.germline_SBF = self.info_sample_dict_2.get('SBF', '.').strip()
            self.germline_ODDRATIO = self.info_sample_dict_2.get('ODDRATIO', '.').strip()
            try:
                self.germline_dp = int(self.info_sample_dict_2.get('DP'))
                self.germline_vd = int(self.info_sample_dict_2.get('VD'))
                self.germline_af = float(self.info_sample_dict_2.get('AF'))
            except:
                pass
            if self.germline_vd == 0:
                self.pvalue = 'NA'

        if self.force_germline:
            self.vd, self.germline_vd = self.germline_vd, self.vd
            self.dp, self.germline_dp = self.germline_dp, self.dp
            self.af, self.germline_af = self.germline_af, self.af

        self.parser_common_tx(self.gene, key_gene_changes, key_gene_tx, key_gene_tx_supplement)
        self.parser_gene_detail()
        self.refgene_exon_n = '.'
        self.parser_exon(self.new_tx, self.refgene_exon_dict)
        try:
            self.geneid = geneid[self.gene]
        except:
            self.geneid = '.'
        self.local_anno_dict, self.local_germline_anno_dict, self.local_header = local_anno_tuple[0], local_anno_tuple[1], local_anno_tuple[2]
        self.comment, self.local_freq = self.anno_file(self.local_anno_dict)
        self.germline_comment, self.germline_local_freq = self.anno_file(self.local_germline_anno_dict)
        self.mrd_filter_dict = mrd_filter_dict

    def parser_info(self, info):
        for i in info.split(';'):
            if '=' not in i:
                continue
            k, v = i.split('=')
            if k == 'cosmic70':
                k = 'cosmic96'
            if k == 'gnomad312_AF_eas':
                k = 'AF_eas'
            if k == 'gnomad312_AF_popmax':
                k = 'AF_popmax'
            self.info_dict.setdefault(k, v)

    def parser_info2(self, info):
        for i2 in info.split(';'):
            if '=' not in i2:
                continue
            k, v = i2.split('=')
            if k == 'cosmic70':
                k = 'cosmic96'
            if k == 'gnomad312_AF_eas':
                k = 'AF_eas'
            if k == 'gnomad312_AF_popmax':
                k = 'AF_popmax'
            if k == 'gnomad312_AF':
                k = 'AF'
            self.info_dict2[k] = v

    def parser_exon(self, new_tx, refgene_exon_dict):
        #refgene_exon_n = '.'
        n_trans = new_tx.split('.')[0]
        #n_exon = new_exon.lstrip('exon')
        try:
            self.refgene_exon_n = refgene_exon_dict[n_trans]
        except:
            pass
        #try:
        #    e_ratio = (int(n_exon)/int(refgene_exon_n))*100
        #    self.pos_to_gene = "{:.2f}%".format(e_ratio)
        #except:
        #    pass

    @property
    def is_pair(self):
        if len(self.arr) == 10:
            return False
        return True

    @property
    def germline_af_string(self):
        if self.germline_af == '.':
            return '.'
        return f'{round(self.germline_af * 100,2)}%'

    @property
    def af_string(self):
        if self.af == '.':
            return '.'
        return f'{round(self.af * 100,2)}%'

    def parser_common_tx(self, gene, key_gene_changes, key_gene_tx, key_gene_tx_supplement):
        self.common_aa_change = '.'
        # 取经典转录本
        self.common_aa_change = match_gene_tx(self.gene_detail, key_gene_tx_supplement, default=False)
        # 判断是否核心库里
        if self.common_aa_change == '.':
            if key_gene_changes.get(gene):
                # 匹配上基因
                for key_change in key_gene_changes.get(gene):
                    # 匹配上氨基酸变化
                    for ichange in self.gene_detail:
                        p_info = ichange.split(':')[-1]  # 特殊情况：chr1	115256565	T	A  注释为: NRAS:NM_002524.5:exon3:c.A146T:p.E49V，且 A146T 在关键位点库，但是表示的是氨基酸变化
                        if key_change in p_info:
                            # 匹配
                            self.common_aa_change = ichange
                            break
        if self.common_aa_change == '.':
            # 位点 chr17	37883239	C	T 注释为 ERBB2:NM_004448.3:exon25:c.C3142T:p.R1048C ，但是在combine.tsv中ERBB2对NM_004448.2，导致找错转录本，需要补充查询转录本
            self.common_aa_change = match_gene_tx(self.gene_detail, key_gene_tx, default=True)
            # 经典转录本文件 key_gene_tx_supplement 需要维护
            # 手动修改了 PTEN    NM_000314.7
            # 手动修改了 PDGFRA  NM_006206.6
            # 手动修改了 NF1     NM_000267.3
        #需要特殊处理的p.
        p_list_es = ['ASXL1:p.G646Wfs*10']
        aas = self.common_aa_change.split(':')
        aaa = aas[-1]
        p_es = f'{gene}:{aaa}'
        if 'p.' in aaa and 'fs' in aaa:
            aaa = aaa.replace('X', '*')
            if 'fs*' in aaa:
                aaas = aaa.split('fs*')
                if len(aaas) > 1:
                    if p_es in p_list_es:
                        aaas[-1] = str(int(aaas[-1]) + 2)
                    else:
                        aaas[-1] = str(int(aaas[-1]) + 1)
                    aaa = 'fs*'.join(aaas)
        aas[-1] = aaa
        self.common_aa_change = ':'.join(aas)

    def parser_gene_detail(self):
        self.new_gene, self.new_tx, self.new_exon, self.new_cdna, self.new_amid = '.', '.', '.', '.', '.'
        if self.common_aa_change != '.':
            if self.common_aa_change != 'UNKNOWN':
                change_arr = self.common_aa_change.split(':')
                try:
                    self.new_gene, self.new_tx, self.new_exon = change_arr[:3]
                except:
                    pass
                if len(change_arr) == 5:
                    self.new_cdna, self.new_amid = change_arr[3:5]
                elif len(change_arr) == 4:
                    self.new_cdna = change_arr[3]
                elif len(change_arr) == 3:
                    self.new_exon = '.'
                    self.new_cdna = change_arr[-1]
                if '>' in self.new_cdna:
                    pass
                else:
                    tem_new_cdna = self.new_cdna.lstrip('c.')
                    a_new_cdna = re.findall('[A-Z]+',tem_new_cdna)
                    n_new_cdna = re.findall('\d+',tem_new_cdna)
                    if len(a_new_cdna) == 2:
                        try:
                            if len(n_new_cdna) == 1 and len(a_new_cdna[0]) == 1 and len(a_new_cdna[1]) == 1:
                                self.new_cdna = f'c.{n_new_cdna[0]}{a_new_cdna[0]}>{a_new_cdna[1]}'
                        except:
                            pass
                if 'substitution' in self.info_dict['ExonicFunc.refGeneWithVer'] or 'stopgain' in self.info_dict['ExonicFunc.refGeneWithVer']:
                    if self.arr[3][0] != self.arr[4][0]:
                        try:
                            res = re.match(r'^(c.\d+_\d+)([a-z,A-Z]+)',self.new_cdna)
                            self.new_cdna = f'{res.groups()[0]}delins{res.groups()[1]}'
                        except:
                            pass
                self.new_amid = self.new_amid.replace('X', '*')
            else:
                self.new_gene = self.gene
        else:
            self.new_gene = self.gene
        if self.gene in self.gene_alias_dict:
            self.new_gene = self.gene_alias_dict.get(self.gene)

    @property
    def end_pos(self):
        # 处理 End_pos
        pos = int(self.arr[1])
        len_ref = len(self.arr[3])
        len_alt = len(self.arr[4])
        if len_ref > len_alt:
            return pos + (len_ref - len_alt)
        else:
            return pos

    @property
    def transvar_c(self):
        '''c点和p点通过transvar转换为HGVS格式'''
        ts = dict()
        ts['c'] = self.new_cdna
        ts['p'] = self.new_amid
        chrom = self.arr[0]
        tsx = self.new_tx.split('.')[0]
        pos = int(self.arr[1])
        end_pos = int(self.arr[1])
        gkey = ''
        tsaa = ''
        gref = self.arr[3]
        galt = self.arr[4]
        len_ref = len(self.arr[3])
        len_alt = len(self.arr[4])
        if len_ref == 1 and len_alt > 1 and '>' not in galt:
            if gref == galt[0]:
                gkey = f'{chrom}:{str(pos)}_{str(end_pos)}ins{galt[1:]}'
            else:
                gkey = f'{chrom}:{str(pos)}delins{galt}'
        elif len_ref > 1 and len_alt == 1:
            end_pos += len_ref - 1
            if gref[0] == galt:
                pos += 1
                gkey = f'{chrom}:{str(pos)}_{str(end_pos)}del'
            else:
                gkey = f'{chrom}:{str(pos)}_{str(end_pos)}delins{galt}'
        elif len_ref > 1 and len_alt > 1 and '>' not in galt:
            end_pos += len_ref - 1
            gkey = f'{chrom}:{str(pos)}_{str(end_pos)}delins{galt}'
        elif len_ref == 1 and len_alt == 1 and '>' not in self.new_cdna:
            gkey = f'{chrom}:{str(pos)}delins{galt}'
        # if gkey and tsx and self.func.strip() != 'UTR3':
        if gkey and self.func.strip() != 'UTR3' and tsx:
            os.system(f'echo {gkey} >> transvar.tem.xls')
            return gkey
        elif gkey and self.func.strip() != 'UTR3' and not tsx and not self.hrr_panel:
            os.system(f'echo {gkey} >> transvar.tem.xls')
            return gkey
        elif gkey and self.common_aa_change == 'UNKNOWN' and ('exonic' in self.func or 'splicing' in self.func):
            os.system(f'echo {gkey} >> transvar.tem.xls')
            return gkey
        return ''
        #     try:
        #         cmd = f'transvar ganno -i {gkey} --refseq |grep {tsx} |cut -f5'
        #         tsaa = os.popen(cmd).read().strip()
        #     except:
        #         tsaa = ''
        # if tsaa:
        #     tsaas = tsaa.split('/')
        #     cts = tsaas[1]
        #     if 'del' in cts and 'ins' not in cts:
        #         cts = cts.split('del')[0] + 'del'
        #     elif 'dup' in cts:
        #         cts = cts.split('dup')[0] + 'dup'
        #     pts = tsaas[2]
        #     if 'del' in pts and 'ins' not in pts:
        #         pts = pts.split('del')[0] + 'del'
        #     elif 'dup' in pts:
        #         pts = pts.split('dup')[0] + 'dup'
        #     ts['c'] = cts
        #     ts['p'] = pts
        #     if ts['p'] == '.':
        #         ts['p'] = self.new_amid
        # return ts

    @property
    def status(self):
        if self.force_germline:
            return 'germline_data'
        return self.info_dict.get('STATUS', '.')

    @property
    def snv_depth(self):
        key_infos = f'{self.arr[0]}_{str(int(self.arr[1]))}'
        if key_infos in self.snvdepth:
            try:
                if int(self.snvdepth[key_infos]) > int(self.dp):
                    return self.snvdepth[key_infos]
            except:
                return self.dp
        try:
            if int(self.snvdepth[f'{self.arr[0]}_{str(int(self.arr[1])+1)}']) > int(self.dp):
                return self.snvdepth[f'{self.arr[0]}_{str(int(self.arr[1])+1)}']
        except:
            return self.dp
        return self.dp

    @property
    def snv_depth_control(self):
        key_infos = f'{self.arr[0]}_{str(int(self.arr[1]))}'
        if key_infos in self.snvdepth_control:
            try:
                if int(self.snvdepth_control[key_infos]) > int(self.germline_dp):
                    return self.snvdepth_control[key_infos]
            except:
                return self.germline_dp
        try:
            if int(self.snvdepth_control[f'{self.arr[0]}_{str(int(self.arr[1])+1)}']) > int(self.germline_dp):
                return self.snvdepth_control[f'{self.arr[0]}_{str(int(self.arr[1])+1)}']
        except:
            return self.germline_dp
        return self.germline_dp

    def anno_file(self, anno_db):
        """ 根据 anno_db 注释 snv 位点
        """
        comment = '\t'.join(['.'] * 9) + '\n'
        local_freq = 0
        try:
            key1 = tuple(self.arr[:2] + self.arr[3:5])
            key2 = tuple(self.common_aa_change)
            if key2 in anno_db:
                comment = '\t'.join(anno_db[key2]) + '\n'
                local_freq = float(anno_db[key2][0])
            elif key1 in anno_db:
                comment = '\t'.join(anno_db[key1]) + '\n'
                local_freq = float(anno_db[key1][0])
            else:
                comment = '\t'.join(['.'] * 9) + '\n'
                local_freq = 0
        except StopIteration:
            pass
        return comment,local_freq

    @property
    def mrd_label(self):
        if not self.mrd_filter_dict:
            return '.'
        key_info1 = f"{self.arr[0]}:{self.arr[1]}:{self.arr[3]}:{self.arr[4]}"
        key_info2 = f"{self.gene}:{self.new_cdna}" if self.new_cdna != '.' else None
        key_info3 = f"{self.gene}:{self.new_amid}" if self.new_amid != '.' else None
        for label in ['unique', 'general']:
            for k in ['position', 'cHGVS', 'pHGVS']:
                if key_info1 in self.mrd_filter_dict[k][label] or \
                   key_info2 in self.mrd_filter_dict[k][label] or \
                   key_info3 in self.mrd_filter_dict[k][label]:
                    return label
        return '.'

    @property
    def repeat_region_anno(self):
        complex_anno = ''
        if self.info_dict.get('simple_repeat', '.') != '.' or self.info_dict.get('rmsk', '.') != '.':
            complex_anno = 'LCR'
        if self.info_dict.get('genomicSuperDups', '.') != '.':
            if complex_anno:
                complex_anno += ';SDR'
            else:
                complex_anno = 'SDR'
        if not complex_anno:
            complex_anno = '.'
        return complex_anno

    @property
    def is_somatic(self):
        return self.status in ['.', 'StrongSomatic', 'LikelySomatic']

    @property
    def is_germline(self):
        return self.status in ['Germline', 'AFDiff', 'germline_data', 'StrongLOH', 'LikelyLOH']


    @property
    def gene(self):
        return self.info_dict.get('Gene.refGeneWithVer', '.')

    @property
    def gene_category(self):
        gene_category_list = list()
        for category in ['genetic', 'drug', 'immune_positive', 'immune_negative', 'immune_hyperprogressive', 'chemical', 'mmr']:
            if self.gene_categorie_dict[category].get(self.gene):
                gene_category_list.append(category)
        if len(gene_category_list) > 0:
            return ';'.join(gene_category_list)
        else:
            return '.'

    @property
    def func(self):
        return self.info_dict.get('Func.refGeneWithVer', '.')

    @property
    def gene_detail(self):
        # 取 change 变化
        if self.info_dict['AAChange.refGeneWithVer'] == '.':
            changes = self.info_dict['GeneDetail.refGeneWithVer'].split('\\x3b')
            changes = ['.' if i == '.' else f'{self.gene}:{i}' for i in changes]
        else:
            changes = self.info_dict['AAChange.refGeneWithVer'].split(',')
        return changes

    @property
    def gene_details_string(self):
        # 取 change 变化
        new_changes = list()
        changes = self.gene_detail
        p_list_es = ['ASXL1:p.G646Wfs*10']
        for aa in changes:
            aas = aa.split(':')
            aaa = aas[-1]
            p_es = f'{self.gene}:{aaa}'
            if 'p.' in aaa and 'fs' in aaa:
                aaa = aaa.replace('X','*')
                if 'fs*' in aaa:
                    aaas = aaa.split('fs*')
                    if len(aaas) > 1:
                        if p_es in p_list_es:
                            aaas[-1] = str(int(aaas[-1]) + 2)
                        else:
                            aaas[-1] = str(int(aaas[-1]) + 1)
                        aaa = 'fs*'.join(aaas)
            aas[-1] = aaa
            new_changes.append(':'.join(aas))
        return ','.join(new_changes)

    @property
    def exonic_func(self):
        return self.info_dict.get('ExonicFunc.refGeneWithVer', '.')

    @property
    def aa_change(self):
        return self.info_dict.get('AAChange.refGeneWithVer', '.')

    @property
    def clinvar_sig(self):
        return self.info_dict.get('CLNSIG', '.')
    
    @property
    def onco_sig(self):
        return self.info_dict.get('ONCSIG', '.')

    @property
    def onekg_af(self):
        return self.info_dict.get('1000g2015aug_all', '.')

    @property
    def gnomad_af(self):
        return self.info_dict.get('AF_eas', '.')

    @property
    def cosmic(self):
        return self.info_dict.get('cosmic96', '.')

    @property
    def in_cosmic_report(self):
        return self.cosmic and self.cosmic != '.'

    @property
    def in_pathogenic(self):
        # 保留致病变异
        # return self.clinvar_sig in ['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic', 'drug_response', 'Oncogenic/Likely_oncogenic', 'Oncogenic', 'Likely_oncogenic']
        return self.clinvar_sig in ['drug_response'] or self.onco_sig in ['Oncogenic/Likely_oncogenic', 'Oncogenic', 'Likely_oncogenic']

    @property
    def rule_1(self):
        # 规则一：过滤产品基因列表
        wes_black_genelist = ['USP17L15', 
                              'USP17L17', 
                              'USP17L19', 
                              'MUC2', 
                              'MUC5AC', 
                              'ZAN', 
                              'GOLGA6L2', 
                              'ZNF717', 
                              'USP17L25', 
                              'USP17L20', 
                              'USP17L11']
        if 'FullRNA.gene.list' in self.product_gene or 'zy.wes.gene.list' in self.product_gene:
            if self.gene in wes_black_genelist:
                return False
            if len(self.gene.split('\\x3b')) > 1:
                return False
            return True
        if self.product_gene_dict.get(self.gene) or self.product_gene_dict.get(self.gene_alias_dict.get(self.gene)):
            return True
        else:
            if '\\x3b' in self.gene:
                if self.product_gene_dict.get(self.gene.split('\\x3b')[0]) or self.product_gene_dict.get(self.gene.split('\\x3b')[1]):
                    return True
        return False

    @property
    def rule_2(self):
        # 规则二：过滤结构变异
        return self.alt not in ['<INV>', '<DEL>', '<DUP>']

    @property
    def rule_3(self):
        # 规则三：1.gene功能属于 exonic 或 splicing
        #       2.排除同义变异
        #       3.MET基因保留13号和14号内含子突变
        #       4.intron突变的CLNSIGCONF列有Pathogenic字眼需保留
        if self.gene == 'MET':
            if self.new_exon in ['intron13', 'intron14']:
                return True
        if 'intronic' in self.func and 'athogenic' in self.info_dict.get('CLNSIGCONF', '.'):
            return True
        if ('exonic' in self.func or 'splicing' in self.func) and 'ncRNA' not in self.func:
            if self.exonic_func != 'synonymous_SNV':
                return True
        return False

    @property
    def rule_4(self):
        # 规则四：排除良性变异, 可缩小范围
        return self.clinvar_sig not in ['Benign', 'Likely_benign', 'Benign/Likely_benign']

    # @property
    # def rule_5(self):
    #     # 规则五：680和1081过滤HLA基因
    #     if 'HLA-' not in self.gene:
    #         return True
    #     return False

    @property
    def rule_6(self):
        # 规则六：保留TERT基因特定突变
        tert_vars = ['chr5_1295228', 'chr5_1295250', 'chr5_1295161']
        tert_var = f'{self.arr[0]}_{self.arr[1]}'
        if self.gene == 'TERT' and  tert_var in tert_vars:
            return True
        return False

    @property
    def rule_7(self):
        # 规则七：MET,EGFR,HER2,ERBB2,KIT,BCOR基因除外，过滤长度>50bp的ins、del变异
        except_genes = ['MET', 'EGFR', 'HER2', 'ERBB2', 'KIT', 'BCOR']
        if (len(self.arr[3]) > 50 or len(self.arr[4]) > 50) and self.gene not in except_genes:
            return False
        return True

    @property
    def somatic_common_pass_1(self):
        if 'MRD' in self.product_gene:
            return True
        if self.sample_type == 'tissue':
            return self.af >= 0.005 and self.vd >= 5
        if self.sample_type == 'blood':
            return self.af >= 0.001 and self.vd >= 3
        return False
    
    @property
    def somatic_common_pass_2(self):
        if 'MRD' in self.product_gene:
            return True
        if self.sample_type == 'tissue':
            return self.af >= 0.01 and self.vd >= 10
        if self.sample_type == 'blood':
            return self.af >= 0.005 and self.vd >= 5
        return False
    
    @property
    def somatic_common_pass_3(self):
        return self.filter_snv == 'PASS'

    @property
    def somatic_common_pass_4(self):
        #保留体系频率/胚系频率大于3的位点
        if self.is_pair and self.germline_af != '.' and float(self.germline_af) > 0 and float(self.af)/float(self.germline_af) < 3:
            return False
        return True
    
    @property
    def somatic_common_pass_5(self):
        #LCR区域的插入或缺失突变，组织突变频率<3%或cfDNA突变频率<1%，则被过滤(BRCA1，BRCA2,MLH1,MSH6,PMS2,MSH2基因除外)
        if "LCR" in self.repeat_region_anno and self.gene not in ['BRCA1', 'BRCA2', 'MLH1', 'MSH6', 'PMS2', 'MSH2']:
            if len(self.arr[3]) != len(self.arr[4]) or len(self.arr[3]) > 1:
                if self.sample_type == 'tissue' and float(self.af) < 0.03:
                    return False
                if self.sample_type == 'blood' and float(self.af) < 0.01:
                    return False
        return True


    @property
    def somatic_rule_1(self):
        # 体系规则一：人群频率<0.01
        if self.onekg_af != '.' and float(self.onekg_af) > 0.01:
            return False
        if self.gnomad_af != '.' and float(self.gnomad_af) > 0.01:
            return False
        if self.info_dict2['AF'] != '.' and float(self.info_dict2['AF']) > 0.01:
            return False
        if self.info_dict["AF_popmax"] != '.' and float(self.info_dict["AF_popmax"]) > 0.01:
            return False
        return True

    @property
    def somatic_rule_2(self):
        # 体系规则二：带胚系项目：克隆性造血突变过滤（保留）
        # 1.体系频率＞胚系频率的位点（保留）
        # 2.体/胚突变频率fisher检验显著性差异位点（保留Fisher检验p值＜0.05）
        if self.is_pair and self.germline_af != '.' and float(self.germline_af) > 0 and float(self.af)/float(self.germline_af) < 1:
            return False
        if self.is_pair and self.germline_af != '.' and float(self.germline_af) > 0 and self.pvalue > 0.05:
            return False
        return True
    
    @property
    def somatic_rule_3(self):
        # 体系规则三：正负链比例：过滤正链占总reads数的分布比例为≥0.8或者≤0.2的位点
        return 0.2 <= float(int(self.ald_f) / self.vd) <= 0.8
    
    @property
    def somatic_rule_4_1(self):
        # 体系规则四：热点基因
        # 优化：热点基因规则，简化逻辑
        is_hotspot = self.hotspot_gene_dict.get(f'{self.gene}|{self.new_amid}') or self.hotspot_gene_dict.get(f'{self.gene}|{self.new_cdna}')
        is_pathogenic = self.in_pathogenic and self.somatic_common_pass_3
        return (is_hotspot or is_pathogenic) and self.somatic_common_pass_1

    @property
    def somatic_rule_4_2(self):
        # 体系规则四：非热点基因
        return self.somatic_common_pass_2 and self.somatic_common_pass_3 and self.local_freq <= 0.15 and self.somatic_common_pass_4 and self.somatic_common_pass_5

    @property
    def germline_rule_1(self):
        if self.germline_af != '.' and float(self.germline_af) <= 0.15:
            return False
        # 体系频率过滤
        # if self.af != '.' and float(self.af) <= 0.15:
        #     return False
        return True

    @property
    def germline_rule_2(self):
        if self.onekg_af != '.' and float(self.onekg_af) > 0.05:
            return False
        if self.gnomad_af != '.' and float(self.gnomad_af) > 0.05:
            return False
        if self.info_dict2['AF'] != '.' and float(self.info_dict2['AF']) > 0.05:
            return False
        if self.info_dict["AF_popmax"] != '.' and float(self.info_dict["AF_popmax"]) > 0.05:
            return False
        return True

    @property
    def germline_rule_3(self):
        genetic_dict = self.gene_categorie_dict.get('genetic')
        if genetic_dict.get(self.gene):
            return True
        return False

    # @property
    def mrd_filter(self,germline_filter=False):
        key_site1 = f'{self.arr[0]}:{self.arr[1]}:{self.arr[3]}:{self.arr[4]}'
        # new_cdna为"."时，无法命中
        key_site2 = f'{self.gene}:{self.new_cdna}' if self.new_cdna != '.' else None
        # new_amid为"."时，无法命中
        key_site3 = f'{self.gene}:{self.new_amid}' if self.new_amid != '.' else None
        # None无法命中, 只要命中一个即可
        #return key_site1 in self.mrd_filter_dict["position"] or key_site2 in self.mrd_filter_dict["cHGVS"] or key_site3 in self.mrd_filter_dict["pHGVS"]
        # 去除key_site3
        #if (key_site1 and key_site1 in self.mrd_filter_dict["position"]["all"]) or (key_site2 and key_site2 in self.mrd_filter_dict["cHGVS"]["all"]) or (key_site3 and key_site3 in self.mrd_filter_dict["pHGVS"]["all"]):
        if (key_site1 and key_site1 in self.mrd_filter_dict["position"]["all"]) or (key_site2 and key_site2 in self.mrd_filter_dict["cHGVS"]["all"]):
            if not germline_filter and self.is_pair and self.germline_af != '.' and float(self.germline_af) > 0 and float(self.af)/float(self.germline_af) < 1:
                return False
            if germline_filter and self.germline_af != '.' and float(self.germline_af) < 0.15:
                return False
            return True
        return False

    @property
    def header(self):
        return '\t'.join([
            'Chr',
            'Pos',
            'End_pos',
            'Ref',
            'Alt',
            'Filter',
            'GeneID',
            'Gene.refGeneWithVer',
            'Func.refGeneWithVer',
            'ExonicFunc.refGeneWithVer',
            'AAChange.refGeneWithVer',
            'gene_details',
            'simple_repeat',
            'rmsk',
            'genomicSuperDups',
            'Low_complexity_and_SuperDup_regions',
            '1000g2015aug_all',
            'AF_all',
            'AF_popmax',
            'AF_eas',
            'avsnp150',
            'CLNALLELEID',
            'CLNDN',
            'CLNDISDB',
            'CLNREVSTAT',
            'CLNSIG',
            'CLNSIGCONF',
            'InterVar',
            'ACMG',
            'ONCDN',
            'ONCDISDB',
            'ONCREVSTAT',
            'ONCSIG',
            'cosmic96',
            'SIFT_pred',
            'Polyphen2_HDIV_pred',
            'LRT_pred',
            'MutationTaster_pred',
            'FATHMM_pred',
            'PROVEAN_pred',
            'M-CAP_pred',
            'PrimateAI_pred',
            'REVEL_score',
            'REVEL_rankscore',
            'TAll_depth_raw',
            'TAll_depth',
            'TAlt_depth',
            'TAlt_depth_forward',
            'TAlt_depth_reverse',
            'TAlt_depth_ALD',
            'TREF_depth_RD',
            'Somatic_SBF',
            'Somatic_ODDRATIO',
            'genotype',
            'Status',
            'germline_depth_raw',
            'germline_depth',
            'germline_alt_depth',
            'germline_alt_depth_forward',
            'germline_alt_depth_reverse',
            'germline_alt_depth_ALD',
            'germline_alt_depth_RD',
            'germline_SBF',
            'germline_ODDRATIO',
            'germline_af',
            'germline_genotype',
            'fisher_p_value',
            'gene_category',
            'mrd_label',
            'new_transcript',
            'all_exon_num',
            'new_gene',
            'new_exon',
            'new_cdna',
            'new_amid',
            'dcs_supporting_reads',
            'scs_supporting_reads',
            'TAF',
            'level',
            'tmp_transvar_input'
        ] + self.local_header[5:]) + '\n'


    @property
    def row(self):
        return '\t'.join([str(s) for s in [
            self.arr[0],
            self.arr[1],
            self.end_pos,
            self.arr[3],
            self.arr[4],
            self.arr[6],
            self.geneid,
            self.gene,
            self.func,
            self.exonic_func,
            self.common_aa_change,
            self.gene_details_string,
            self.info_dict.get('simple_repeat', '.'),
            self.info_dict.get('rmsk', '.'),
            self.info_dict.get('genomicSuperDups', '.').replace('\\x3b',';').replace('\\x3d','='),
            self.repeat_region_anno,
            self.onekg_af,
            self.info_dict2['AF'],
            self.info_dict["AF_popmax"],
            self.gnomad_af,
            self.info_dict["avsnp150"],
            self.info_dict["CLNALLELEID"],
            self.info_dict["CLNDN"],
            self.info_dict["CLNDISDB"],
            self.info_dict["CLNREVSTAT"],
            self.info_dict["CLNSIG"],
            self.info_dict.get("CLNSIGCONF", "."),
            self.info_dict.get("InterVar", "."),
            self.info_dict.get("ACMG", "."),
            self.info_dict["ONCDN"],
            self.info_dict["ONCDISDB"],
            self.info_dict["ONCREVSTAT"],
            self.info_dict["ONCSIG"],
            self.info_dict["cosmic96"],
            self.info_dict["SIFT_pred"],
            self.info_dict["Polyphen2_HDIV_pred"],
            self.info_dict["LRT_pred"],
            self.info_dict["MutationTaster_pred"],
            self.info_dict["FATHMM_pred"],
            self.info_dict["PROVEAN_pred"],
            self.info_dict["M-CAP_pred"],
            self.info_dict["PrimateAI_pred"],
            self.info_dict["REVEL_score"],
            self.info_dict["REVEL_rankscore"],
            self.snv_depth,
            self.dp,
            self.vd,
            self.ald_f,
            self.ald_r,
            self.ald,
            self.rd,
            self.somatic_SBF,
            self.somatic_ODDRATIO,
            self.gt,
            self.status,
            self.snv_depth_control,
            self.germline_dp,
            self.germline_vd,
            self.germline_ald_f,
            self.germline_ald_r,
            self.germline_ald,
            self.germline_rd,
            self.germline_SBF,
            self.germline_ODDRATIO,
            self.germline_af_string,
            self.germline_gt,
            self.pvalue,
            self.gene_category,
            self.mrd_label,
            self.new_tx,
            self.refgene_exon_n,
            self.new_gene,
            self.new_exon,
            self.new_cdna,
            [pi.replace('*fs*1', '*') if pi.endswith('*fs*1') else pi.replace('fs*1', '*') if pi.endswith('fs*1') else pi for pi in [self.new_amid]][0],
            self.info_dict.get('DCS', '.'),
            self.info_dict.get('SCS', '.'),
            self.af_string,
            '',
        ]]) + '\t' + self.transvar_c + '\n'

    @property
    def pass_common_rule(self):
        # if 'zy.680' in self.product_gene or 'zy.1081' in self.product_gene or 'zy.1008' in self.product_gene:
        #     return (self.rule_1 and self.rule_2 and self.rule_3 and self.rule_4 and self.rule_5) or (self.rule_1 and self.rule_6)
        return (self.rule_1 and self.rule_2 and self.rule_3 and self.rule_4 and self.rule_7) or (self.rule_1 and self.rule_6 and self.rule_7)

    @property
    def is_somatic_filter(self):
        if 'MRD' in self.product_gene:
            if self.is_somatic:
                return self.mrd_filter()
        else:
            if self.pass_common_rule:
                if self.is_somatic:
                    if (self.somatic_rule_1 and self.somatic_rule_2 and self.somatic_rule_3 and (self.somatic_rule_4_1 or self.somatic_rule_4_2)) or self.rule_6:
                        return True
        return False

    @property
    def is_germline_filter(self):
        if 'MRD' in self.product_gene:
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
            return self.row.strip('\n') + '\t' + self.comment
        return ''

    @property
    def germline_row(self):
        if self.is_germline_filter:
            return self.row.strip('\n') + '\t' + self.germline_comment
        return ''

    @property
    def somatic_row_all(self):
        return self.row.strip('\n') + '\t' + self.comment

    @property
    def germline_row_all(self):
        return self.row.strip('\n') + '\t' + self.germline_comment

    @property
    def common_debug_header(self):
        return self.header.replace('\n', '\t') + '\t'.join(['rule_1', 'rule_2', 'rule_3', 'rule_4'])

    @property
    def debug_somatic_header(self):
        return self.common_debug_header + '\t' + '\t'.join(['somatic_rule_1', 'somatic_rule_2_1', 'somatic_rule_2_2', 'somatic_rule_3', 'fanal_rule']) + '\n'

    @property
    def debug_germline_header(self):
        return self.common_debug_header + '\t' + '\t'.join(['germline_rule_1', 'germline_rule_2', 'final_rule']) + '\n'

    @property
    def common_debug_row(self):
        new_row = self.row.replace('\n', '')
        return f'{new_row}\t{self.rule_1}\t{self.rule_2}\t{self.rule_3}\t{self.rule_4}'

    @property
    def debug_somatic_row(self):
        return f'{self.common_debug_row}\t{self.somatic_rule_1}\t{self.somatic_rule_2_1}\t{self.somatic_rule_2_2}\t{self.somatic_rule_3}\t{self.is_somatic_filter}\n'

    @property
    def debug_germline_row(self):
        return f'{self.common_debug_row}\t{self.germline_rule_1}\t{self.germline_rule_2}\t{self.is_germline_filter}\n'


def main(args):
    key_site = args.key_site
    product_gene = args.product_gene
    hotspot_gene = args.hotspot_gene
    gene_trans = args.gene_trans
    sample_type = args.sample_type
    annovar_vcf_file = args.annovar_vcf_file
    gene_categorie_file = args.gene_categorie_file
    hg19_refGene_file = args.hg19_refGene_file
    hg19_gtf_file = args.hg19_gtf_file
    modify_varannovar = args.modify_varannovar
    depth_tsv = args.depth_tsv
    depth_tsv_control = args.depth_tsv_control
    gene_alias = args.gene_alias
    local_freq = args.local_freq
    mrd_filter_file = args.mrd_filter_file
    force_germline = args.force_germline
    debug = args.debug

    key_gene_tx, key_gene_changes = read_key_site(key_site)
    hotspot_gene_dict = read_hotspot_gene_dict(hotspot_gene)
    product_gene_dict = read_product_gene_dict(product_gene)
    key_gene_tx_supplement = read_gene_trans(gene_trans)
    gene_categorie_dict = load_resource_gene_categories(gene_categorie_file)
    refgene_exon_dict = read_refgene(hg19_refGene_file)
    geneid = read_gtf(hg19_gtf_file)
    varannovar = read_varannovar(modify_varannovar)
    snvdepth = read_depthtsv(depth_tsv)
    snvdepth_control = read_depthtsv(depth_tsv_control)
    gene_alias_dict = read_gene_trans(gene_alias)
    local_anno_tuple = get_anno_file(product_gene, local_freq)
    mrd_filter_dict = read_mrd_filter(mrd_filter_file)

    # 做一些判断
    is_brca = False
    if 'brca' in product_gene.split('/')[-1]:
        is_brca = True
    is_lynch = False
    if 'lynch' in product_gene.split('/')[-1]:
        is_lynch = True
    is_105 = False
    if 'zy.105' in product_gene.split('/')[-1] or 'zy.genetics.15' in product_gene.split('/')[-1]:
        is_105 = True

    # 输出文件
    prefix = args.annovar_vcf_file.replace('.vcf', '')
    prefix = os.path.basename(prefix)
    filter_xls = f'{prefix}.filter.tem.xls'
    all_xls = f'{prefix}.all.tem.xls'
    germline_filter_xls = f'{prefix}.filter.germline.tem.xls'
    germline_all_xls = f'{prefix}.all.germline.tem.xls'

    filter_xls_o = f'{prefix}.filter.xls'
    all_xls_o = f'{prefix}.all.xls'
    germline_filter_xls_o = f'{prefix}.filter.germline.xls'
    germline_all_xls_o = f'{prefix}.all.germline.xls'

    if sample_type not in ['blood', 'tissue']:
        raise Exception(f'invalid sample type {sample_type}')

    fh_o_filter = open(filter_xls, 'w')
    fh_o_all = open(all_xls, 'w')
    fh_o_germline_filter = open(germline_filter_xls, 'w')
    fh_o_germline_all = open(germline_all_xls, 'w')

    # 总表，用于判断筛选条件
    with open(annovar_vcf_file) as fh:
        write_header_flag = False
        for line in fh:
            if line.startswith('#'):
                continue

            v = VariantRow(line, product_gene, product_gene_dict, sample_type, hotspot_gene_dict, key_gene_changes, key_gene_tx, key_gene_tx_supplement, gene_categorie_dict, refgene_exon_dict, geneid, varannovar, snvdepth, snvdepth_control, gene_alias_dict, local_anno_tuple, mrd_filter_dict, force_germline, is_brca=is_brca, is_lynch=is_lynch, is_105=is_105, debug=debug)

            if not write_header_flag:
                if debug:
                    fh_o_filter.write(v.debug_somatic_header)
                    fh_o_all.write(v.debug_somatic_header)
                    fh_o_germline_filter.write(v.debug_germline_header)
                    fh_o_germline_all.write(v.debug_germline_header)
                else:
                    fh_o_filter.write(v.header)
                    fh_o_all.write(v.header)
                    fh_o_germline_filter.write(v.header)
                    fh_o_germline_all.write(v.header)
                write_header_flag = True

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

    fh_o_all.close()
    fh_o_filter.close()
    fh_o_germline_filter.close()
    fh_o_germline_all.close()

    os.system(f'sort -u transvar.tem.xls > transvar.xls')
    cmd = f'transvar ganno -l transvar.xls --refseq > transvar_output.xls'
    rc = os.system(cmd)
    if not rc:
        print('transvar exec success')

    f_transvar = f'transvar_output.xls'

    if 'FullRNA' in product_gene or 'zy.wes' in product_gene or 'zy.hrr' in product_gene:
        add_transvar_result(filter_xls, filter_xls_o, f_transvar, key_gene_tx_supplement)
        os.system(f'mv {all_xls} {all_xls_o}')
        add_transvar_result(germline_filter_xls, germline_filter_xls_o, f_transvar, key_gene_tx_supplement)
        os.system(f'mv {germline_all_xls} {germline_all_xls_o}')
    else:
        add_transvar_result(filter_xls, filter_xls_o, f_transvar, key_gene_tx_supplement)
        add_transvar_result(all_xls, all_xls_o, f_transvar, key_gene_tx_supplement)
        add_transvar_result(germline_filter_xls, germline_filter_xls_o, f_transvar, key_gene_tx_supplement)
        add_transvar_result(germline_all_xls, germline_all_xls_o, f_transvar, key_gene_tx_supplement)

    os.system(f'rm -f *.tem.xls')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--key_site', required=True, help='key site file')
    parser.add_argument('-l', '--product_gene', required=True, help='product gene list file')
    parser.add_argument('-g', '--hotspot_gene', required=True, help='hot gene list file')
    parser.add_argument('-t', '--gene_trans', required=True, help='gene transcript file')
    parser.add_argument('-s', '--sample_type', required=True, help='sample_type')
    parser.add_argument('-f', '--annovar_vcf_file', required=True, help='annovar vcf file')
    parser.add_argument('-o', '--filtered_annovar_vcf_file', required=True, help='filtered annovar vcf file')
    parser.add_argument('-c', '--gene_categorie_file', required=True, help='Base germline vcf file')
    parser.add_argument('-r', '--hg19_refGene_file', required=True, help='hg19_refgene.txt')
    parser.add_argument('-r2', '--hg19_gtf_file', required=True, help='hg19_genomic.gff')
    parser.add_argument('-r3', '--modify_varannovar', required=False, help='modify_varannovar.xls')
    parser.add_argument('-r4', '--depth_tsv', required=False, help='bamdst depth file')
    parser.add_argument('-r5', '--depth_tsv_control', required=False, help='bamdst depth control file')
    parser.add_argument('-ga', '--gene_alias', required=False, help='gene alias file')
    parser.add_argument('-lc', '--local_freq', required=False, help='local freq file (zip / 单个txt / fn.txt,germline_fn.txt)')
    parser.add_argument('-mf', '--mrd_filter_file', required=False, help='mrd filter file')
    parser.add_argument('-y', '--force_germline', required=False, help='Base germline vcf file')
    parser.add_argument('-d', '--debug', required=False, help='debug on')
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
