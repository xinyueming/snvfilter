# SNV Filter

SNV (单核苷酸变异) 过滤脚本，用于处理肿瘤-正常配对样本的 ANNOVAR 注释 VCF 文件。

## 功能

- 体系突变筛选
- 胚系突变筛选
- TransVar HGVS 标准化注释
- 本地频率数据库注释
- MRD 突变标记

## 依赖

- Python 3.x
- pandas
- openpyxl
- TransVar (HGVS 标准化工具)

## 使用

```bash
python3 pre_select_v2.py \
  -k hotspot/combine.tsv \
  -l hotspot/zy.120.gene.list \
  -t hotspot/sorted.gene.tx.txt \
  -s blood \
  -f sample.var.hg19_multianno.vcf \
  -o output.filter.xls \
  -c gene_categories.xlsx \
  -r hg19_refgene.txt \
  -r2 hg19_genomic.gff
```

## 参数说明

| 参数 | 说明 |
|------|------|
| -k | 关键位点文件 |
| -l | 产品基因列表 |
| -g | 热点基因文件 |
| -t | 基因转录本文件 |
| -s | 样本类型 (blood/tissue) |
| -f | ANNOVAR 注释 VCF 文件 |
| -o | 输出文件 |
| -c | 基因分类文件 |
| -r | hg19 refGene 文件 |
| -r2 | hg19 GTF 文件 |
| -lc | 本地频率文件 |
| -mf | MRD 过滤文件 |