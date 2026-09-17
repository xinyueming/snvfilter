# SNV Filter

SNV (单核苷酸变异) 过滤工具，用于处理肿瘤-正常配对样本的 ANNOVAR 注释 VCF 文件。

## 安装

```bash
pip install -e .
```

## 使用

### 基本过滤

```bash
snvfilter \
  -k hotspot/combine.tsv \
  -l hotspot/zy.120.gene.list \
  -g hotspot/hotspot.tsv \
  -t hotspot/sorted.gene.tx.txt \
  -s blood \
  -f sample.var.hg19_multianno.vcf \
  -c gene_categories.xlsx \
  -r hg19_refgene.txt \
  -r2 hg19_genomic.gff
```

### 使用配置文件

```bash
# 生成默认配置
snvfilter init-config --output my_config.yaml

# 使用配置运行
snvfilter --config my_config.yaml
```

### 查看默认规则

```bash
snvfilter show-rules
```

### 命令行参数覆盖配置

```bash
snvfilter --config my_config.yaml -s tissue -o custom_prefix
```

参数优先级：**命令行 > 用户配置 > 默认配置**

## 参数说明

| 参数 | 说明 |
|------|------|
| `-C, --config` | 配置文件 (YAML) |
| `-R, --rules` | 自定义规则文件 (YAML) |
| `-k` | 关键位点文件 |
| `-l` | 产品基因列表 |
| `-g` | 热点基因文件 |
| `-t` | 基因转录本文件 |
| `-s` | 样本类型 (blood/tissue) |
| `-f` | ANNOVAR 注释 VCF 文件 |
| `-c` | 基因分类文件 (xlsx) |
| `-r` | hg19 refGene 文件 |
| `-r2` | hg19 GTF 文件 |
| `-r3` | modify_varannovar 文件 |
| `-r4` | bamdst depth 文件 |
| `-r5` | bamdst depth control 文件 |
| `-ga` | gene alias 文件 |
| `-lc` | 本地频率文件 |
| `-mf` | MRD 过滤文件 |
| `-y` | 强制胚系模式 |
| `-d` | 调试模式 |
| `-o` | 输出前缀 |

## 项目结构

```
snvfilter/
├── pyproject.toml        # 包配置
├── config/
│   ├── default.yaml      # 默认配置
│   └── rules.yaml        # 默认过滤规则
├── snvfilter/
│   ├── __init__.py
│   ├── cli.py            # CLI 入口
│   ├── config.py         # 配置加载
│   ├── rules.py          # 规则模块
│   ├── filters.py        # 过滤逻辑 (VariantRow)
│   ├── utils.py          # 工具函数
│   └── parser.py         # VCF 解析 (预留)
└── tests/
    ├── test_config.py
    └── test_rules.py
```

## 规则分类

- **common**: 通用规则（产品基因、结构变异、基因功能、良性变异、TERT 热点、大片段 indel）
- **somatic**: 体系突变规则（人群频率、克隆性造血、正负链比例、热点基因）
- **germline**: 胚系突变规则（频率、人群频率、基因类别）
- **mrd**: MRD 规则

## 依赖

- Python >= 3.8
- pandas
- openpyxl
- pyyaml
