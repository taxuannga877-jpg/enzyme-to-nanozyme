# 新功能使用指南

## 概述

Enzyme Structure Viewer 现已集成三大核心功能：

1. **PDB 搜索与下载** - 为没有PDB文件的条目提供一键下载功能
2. **EasIFA 催化位点预测** - 使用深度学习模型预测未标注的催化位点
3. **Motif 提取与可视化** - 提取催化残基的几何信息并可视化

---

## 功能详解

### 1. PDB 搜索与下载

#### 使用场景
当查询某个EC号时，如果某些条目没有本地PDB文件，系统会自动显示"下载PDB"按钮。

#### 操作步骤
1. 在左侧选择一个EC号
2. 查看PDB列表，找到显示"PDB文件不存在"的条目
3. 点击绿色的"下载PDB"按钮
4. 等待下载完成（通常几秒钟）
5. 页面自动刷新，现在可以查看该结构了

#### API 端点
```bash
POST /api/download_pdb
Content-Type: application/json

{
  "alphafold_id": "Q84J37",
  "uniprot_id": "Q84J37"
}
```

#### 响应示例
```json
{
  "status": "success",
  "pdb_path": "/home/user/.111tangboshi/cache/pdb/AF-Q84J37-F1-model_v6.pdb",
  "message": "PDB文件下载成功: AF-Q84J37-F1-model_v6.pdb"
}
```

---

### 2. EasIFA 催化位点预测

#### 使用场景
当某个蛋白质结构**没有已知的催化位点标注**时，可以使用EasIFA模型进行预测。

#### 操作步骤
1. 点击"查看结构"按钮打开结构视图
2. 如果该结构没有活性位点标注，右上角会显示黄色的"预测催化位点"按钮
3. 点击"预测催化位点"
4. 等待模型推理完成（可能需要几分钟，取决于序列长度）
5. 预测结果会自动显示在结构视图和催化位点信息表中

#### 预测结果
- **Binding sites** (结合位点) - 标记为红色
- **Catalytic sites** (催化位点) - 标记为绿色
- **Other sites** (其他位点) - 标记为黄色

#### API 端点
```bash
POST /api/predict_active_sites
Content-Type: application/json

{
  "pdb_path": "/path/to/structure.pdb",
  "ec_number": "1.10.3.2",
  "uniprot_id": "Q84J37"
}
```

#### 响应示例
```json
{
  "status": "success",
  "predicted_sites": [
    {
      "residue_index": 77,
      "residue_name": "HIS",
      "site_type": "Catalytic",
      "confidence": 0.95
    }
  ],
  "message": "成功预测到 21 个活性位点"
}
```

#### 注意事项
- 需要EasIFA模型文件（位于 `data/models/easifa/checkpoints/`）
- 如果模型文件不存在，系统会提示错误
- 预测时间与蛋白质序列长度成正比

---

### 3. Motif 提取与可视化

#### 使用场景
当结构**已有催化位点信息**（来自UniProt标注或EasIFA预测）时，可以提取催化Motif用于纳米酶设计。

#### 操作步骤
1. 确保结构已有活性位点标注
2. 在结构视图右上角，点击蓝色的"Motif提取"按钮
3. 系统自动提取催化残基的几何信息
4. 自动跳转到Motif可视化页面

#### Motif展示页面包含：
- **催化残基信息**：每个锚点原子的详细信息（残基类型、编号、坐标等）
- **几何约束**：残基之间的距离和角度关系
- **3D可视化**：（可选）催化位点的3D结构
- **下载功能**：可下载Motif的JSON文件

#### API 端点
```bash
POST /api/extract_motif
Content-Type: application/json

{
  "pdb_path": "/path/to/structure.pdb",
  "ec_number": "1.10.3.2",
  "uniprot_id": "Q84J37",
  "nanozyme_type": "LAC"
}
```

#### 响应示例
```json
{
  "status": "success",
  "motif": {
    "motif_id": "Q84J37_1.10.3.2_LAC",
    "anchor_atoms": [...],
    "geometry_constraints": [...]
  },
  "motif_file": "/home/user/.111tangboshi/motifs/Q84J37_1.10.3.2_LAC.json",
  "message": "成功提取 21 个催化残基"
}
```

#### Motif JSON 格式
```json
{
  "motif_id": "Q84J37_1.10.3.2_LAC",
  "source_uniprot_id": "Q84J37",
  "source_ec_number": "1.10.3.2",
  "nanozyme_type": "LAC",
  "anchor_atoms": [
    {
      "atom_name": "ND1",
      "residue_name": "HIS",
      "residue_number": 77,
      "coordinates": [0.657, 9.551, -6.275],
      "is_donor": true
    }
  ],
  "geometry_constraints": [
    {
      "constraint_type": "distance",
      "atom_indices": [0, 1],
      "value": 3.85,
      "unit": "angstrom"
    }
  ]
}
```

---

## 工作流程示例

### 完整流程：从EC号到Motif提取

```
1. 选择EC号 (如 1.10.3.2)
   ↓
2. 查看PDB列表
   ↓
3. 如果没有PDB → 下载PDB
   ↓
4. 查看结构
   ↓
5. 如果没有活性位点 → 预测催化位点
   ↓
6. 提取Motif
   ↓
7. 查看并下载Motif数据
```

---

## 技术架构

### 后端 (Flask)
- `app.py`: 主应用，集成所有功能模块
- `nanozyme_mining/database/uniprot_fetcher.py`: PDB下载
- `nanozyme_mining/prediction/easifa_predictor.py`: EasIFA预测
- `nanozyme_mining/extraction/extractor.py`: Motif提取

### 前端 (HTML/JavaScript)
- `index.html`: 主界面，集成新功能按钮
- `motif_view.html`: Motif可视化页面

### 数据存储
- `cache/pdb/`: PDB结构文件
- `cache/json/`: 活性位点标注数据
- `motifs/`: 提取的Motif JSON文件

---

## 测试

运行测试脚本验证所有功能：

```bash
cd /home/tangboshi/.111tangboshi/enzyme_viewer
python3 test_new_features.py
```

---

## 常见问题

### Q: EasIFA预测失败怎么办？
A: 确保模型文件已正确放置在 `data/models/easifa/checkpoints/` 目录下。

### Q: 下载PDB失败？
A: 可能是网络问题或AlphaFold数据库中不存在该结构，请稍后重试。

### Q: Motif提取显示"未找到催化残基"？
A: 检查活性位点信息是否存在，或者尝试先运行EasIFA预测。

### Q: 支持哪些纳米酶类型？
A: 目前支持 POD (过氧化物酶), CAT (过氧化氢酶), SOD (超氧化物歧化酶), LAC (漆酶), OXD (氧化酶), GSH (谷胱甘肽过氧化物酶)。

---

## 更新日志

### v2.0.0 (2024-01-09)
- ✨ 新增 PDB搜索与下载功能
- ✨ 新增 EasIFA催化位点预测功能
- ✨ 新增 Motif提取与可视化功能
- 🎨 优化UI，添加功能入口按钮
- 📝 添加完整的API文档

---

## 联系方式

如有问题或建议，请联系开发团队。

