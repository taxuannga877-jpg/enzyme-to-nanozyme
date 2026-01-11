# 纳米酶PDB库

## 概述

本目录包含所有纳米酶EC号的AlphaFold结构文件（PDB格式）。

## 统计信息

- **EC号总数**: 13
- **PDB文件总数**: 1344
- **成功下载**: 1344
- **下载失败**: 5

## 目录结构

```
pdb_library/
├── README.md                    # 本文件
├── library_index.json           # 详细索引文件
├── 1_1_3_4/                     # EC 1.1.3.4 (Glucose Oxidase)
├── 1_10_3_2/                    # EC 1.10.3.2 (Laccase)
├── 1_11_1_6/                    # EC 1.11.1.6 (Catalase)
├── ...
```

## EC号列表

- **1.1.3.4** (Glucose Oxidase): 8/8 成功, 0 失败
- **1.10.3.2** (Laccase): 88/88 成功, 0 失败
- **1.11.1.11** (Peroxidase): 19/19 成功, 0 失败
- **1.11.1.12** (Glutathione Peroxidase): 12/12 成功, 0 失败
- **1.11.1.21** (Peroxidase): 338/340 成功, 2 失败
- **1.11.1.6** (Catalase): 147/147 成功, 0 失败
- **1.11.1.7** (Peroxidase): 157/157 成功, 0 失败
- **1.11.1.9** (Glutathione Peroxidase): 47/47 成功, 0 失败
- **1.15.1.1** (Superoxide Dismutase): 435/438 成功, 3 失败
- **1.3.3.4** (Oxidase): 15/15 成功, 0 失败
- **1.4.3.4** (Oxidase): 23/23 成功, 0 失败
- **3.1.21.1** (DNase): 18/18 成功, 0 失败
- **3.1.3.1** (Phosphatase): 37/37 成功, 0 失败

## 使用说明

1. 所有PDB文件按EC号组织在各自的子目录中
2. 文件名格式: `AF-{alphafold_id}-F1-model_v6.pdb`
3. 详细索引信息请查看 `library_index.json`

## 更新日期

2026-01-12 01:22:46
