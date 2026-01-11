# Nanozyme Assembly Guide

## 核心理念

**关键记忆：纳米酶不是蛋白质！**

纳米酶是人工设计的具有酶活性的纳米材料，包括：
- 金属氧化物 (Fe3O4, CeO2, MnO2等)
- 单原子催化剂 (Fe-N4-C, Cu-N4等)
- 金属有机框架 (MOF)
- 碳基纳米材料 (石墨烯，碳纳米管)
- 金属团簇 (Pt, Au, Ag)

## 快速开始

### 1. 创建催化Motif

```python
from nanozyme_mining.assembly import (
    NanozymeMotif,
    AnchorAtom,
    MaterialType,
    CoordinationType,
    MetalProperties,
)

# 定义锚点原子
anchor_atoms = [
    AnchorAtom(
        atom_name="FE1",
        element="Fe",
        coordinates=[0.0, 0.0, 0.0],
        is_metal_center=True,
        oxidation_state=3,
    ),
    # ... 更多原子
]

# 创建motif
motif = NanozymeMotif(
    motif_id="FePOD_001",
    nanozyme_type="POD",
    anchor_atoms=anchor_atoms,
    material_type=MaterialType.METAL_OXIDE,
)

# 保存
motif.save("./motif_library/POD/FePOD_001.json")
```

### 2. 组装纳米酶

```python
from nanozyme_mining.assembly import NanozymeAssembler

# 初始化组装器
assembler = NanozymeAssembler(
    strategy="rule",  # 规则组装
    material_type=MaterialType.METAL_OXIDE,
    output_dir="./output",
)

# 组装
nanozyme = assembler.assemble(
    motifs=motif,
    num_metal_centers=3,
)

# 保存结构
assembler.save_structure(nanozyme, "Fe3O4_nanozyme", formats=['xyz', 'json'])
```

### 3. 验证结构

```python
from nanozyme_mining.assembly import NanozymeValidator

validator = NanozymeValidator()
result = validator.validate(nanozyme)

print(f"Valid: {result.is_valid}")
print(f"Errors: {result.errors}")
print(f"Scores: {result.scores}")
```

## 组装策略

### 策略1：规则组装 (Rule-based)

**适用场景：** 简单纳米酶，明确的材料类型

```python
assembler = NanozymeAssembler(strategy="rule")

# 不同材料类型的组装
# 1. 金属氧化物
nanozyme = assembler.assemble(
    motifs=motif,
    material_type=MaterialType.METAL_OXIDE,
    size=20,  # 原子数
)

# 2. 单原子催化剂
nanozyme = assembler.assemble(
    motifs=motif,
    material_type=MaterialType.SAC,
)

# 3. 金属配合物
nanozyme = assembler.assemble(
    motifs=motif,
    material_type=MaterialType.METAL_COMPLEX,
)

# 4. MOF
nanozyme = assembler.assemble(
    motifs=motif,
    material_type=MaterialType.MOF,
)
```

**规则组装器工作原理：**
1. 提取金属中心从motif
2. 根据配位几何放置配位原子
3. 添加骨架/支撑结构
4. 连接多个motif（如果有）

### 策略2：扩散模型组装 (Diffusion-based)

**适用场景：** 复杂纳米酶，需要灵活生成

```python
# TODO: 实现中
assembler = NanozymeAssembler(strategy="diffusion")
nanozyme = assembler.assemble(
    motifs=[motif1, motif2],
    target_size=50,
)
```

**借鉴DiffLinker的思路：**
- 固定片段：催化motif
- 生成连接器：金属团簇、配体
- 保持3D等变性

### 策略3：模板组装 (Template-based)

**适用场景：** 已知结构类型

```python
# TODO: 实现中
assembler = NanozymeAssembler(strategy="template")
nanozyme = assembler.assemble(
    template="single_atom_catalyst",
    motifs=motif,
    parameters={'support': 'graphene'},
)
```

## 支持的材料类型

### 1. 金属氧化物 (Metal Oxide)

```python
from nanozyme_mining.assembly import MaterialType

material_type = MaterialType.METAL_OXIDE

# 典型例子
# - Fe3O4: 过氧化物酶样
# - CeO2: 抗氧化
# - MnO2: 氧化酶样
```

**特点：**
- 纳米颗粒结构
- 表面活性位点
- 可调节尺寸和形貌

### 2. 单原子催化剂 (SAC)

```python
material_type = MaterialType.SAC

# 典型例子
# - Fe-N4-C: POD/CAT活性
# - Co-N4: OER/ORR
# - Cu-N4: CO2还原
```

**特点：**
- 单个金属原子
- 锚定在载体上（C, N掺杂)
- 最大原子利用效率

### 3. 金属有机框架 (MOF)

```python
material_type = MaterialType.MOF

# 典型例子
# - Fe-MOF: POD活性
# - Zr-MOF: 酸催化
# - Cu-MOF: 氧化反应
```

**特点：**
- 金属节点 + 有机连接体
- 高比表面积
- 可调孔结构

### 4. 金属配合物

```python
material_type = MaterialType.METAL_COMPLEX

# 典型例子
# - Fe-porphyrin: 模拟细胞色素P450
# - Cu-phenanthroline: 核酸酶活性
```

**特点：**
- 分子级精确结构
- 配位化学控制
- 易于修饰

## Motif库管理

### 加载Motif库

```python
from nanozyme_mining.assembly import MotifLibrary

library = MotifLibrary(library_dir="./motif_library")

print(f"Total motifs: {len(library)}")
print(library)  # 显示各类型统计
```

### 查询Motif

```python
# 按纳米酶类型
pod_motifs = library.get_by_type("POD")

# 按材料类型
metal_oxide_motifs = library.get_by_material(MaterialType.METAL_OXIDE)

# 复杂搜索
motifs = library.search(
    nanozyme_type="POD",
    material_type=MaterialType.METAL_OXIDE,
    min_confidence=0.8,
)
```

### 转换现有Motif

如果你有基于蛋白质提取的motif，可以转换为纳米酶motif：

```python
from nanozyme_mining.assembly.motif_enhanced import convert_basic_to_nanozyme_motif

nanozyme_motif = convert_basic_to_nanozyme_motif(
    basic_motif_path="./motif_library/Catalase/P49317_1.11.1.6_Catalase.json",
    material_type=MaterialType.METAL_OXIDE,
)

nanozyme_motif.save("./motif_library/CAT/P49317_nanozyme.json")
```

## 结构操作

### 查看结构信息

```python
# 基本信息
print(nanozyme)
print(f"Formula: {nanozyme.compute_formula()}")
print(f"Composition: {nanozyme.compute_composition()}")

# 活性位点
active_atoms = nanozyme.get_active_site_atoms()
print(f"Active site atoms: {len(active_atoms)}")

# 金属中心
metal_atoms = nanozyme.get_metal_atoms()
print(f"Metal atoms: {len(metal_atoms)}")
for metal in metal_atoms:
    print(f"  {metal.element} at {metal.coordinates}")
```

### 结构变换

```python
import numpy as np

# 平移
nanozyme.translate(np.array([5.0, 0.0, 0.0]))

# 旋转
from scipy.spatial.transform import Rotation
R = Rotation.from_euler('z', 45, degrees=True).as_matrix()
nanozyme.rotate(R)

# 居中
centroid = nanozyme.get_centroid()
nanozyme.translate(-centroid)
```

### 导出结构

```python
# XYZ格式
nanozyme.to_xyz("nanozyme.xyz", comment="Fe3O4 nanozyme")

# JSON格式（保留所有信息）
nanozyme.to_json("nanozyme.json")

# 读取
from nanozyme_mining.assembly import NanozymeStructure
loaded = NanozymeStructure.from_json("nanozyme.json")
```

## 验证检查项

### 化学有效性

```python
from nanozyme_mining.assembly import NanozymeValidator

validator = NanozymeValidator(
    strict=False,  # 严格模式
    check_coordination=True,  # 检查配位
    check_geometry=True,  # 检查几何
    check_accessibility=True,  # 检查可接近性
)

result = validator.validate(nanozyme)

# 检查结果
if result.is_valid:
    print("✓ Structure is valid")
else:
    print("✗ Structure has issues:")
    for error in result.errors:
        print(f"  ERROR: {error}")

for warning in result.warnings:
    print(f"  WARNING: {warning}")

# 质量评分
print(f"Completeness: {result.scores['completeness']:.2f}")
print(f"Geometry quality: {result.scores['geometry_quality']:.2f}")
print(f"Coordination quality: {result.scores['coordination_quality']:.2f}")
```

### 检查项目

1. **键长检查**
   - Fe-O: 1.7-2.3 Å
   - Fe-N: 1.8-2.4 Å
   - Cu-N: 1.8-2.3 Å
   - ...

2. **原子重叠**
   - 非成键原子最小距离 > 1.5 Å

3. **金属配位**
   - 配位数是否合理
   - 配位几何是否正确

4. **活性位点可接近性**
   - 活性原子是否被埋藏

## 批量组装

### 从多个motif批量生成

```python
# 准备motif组
motif_groups = [
    [motif1],
    [motif2],
    [motif3, motif4],  # 组合motif
]

# 批量组装
nanozymes = assembler.assemble_batch(
    motif_groups=motif_groups,
)

# 保存所有结构
for i, nanozyme in enumerate(nanozymes):
    assembler.save_structure(nanozyme, f"nanozyme_{i:03d}")
```

## 与下游计算对接

### 准备DFT计算输入

```python
# 导出XYZ格式供Gaussian/VASP使用
nanozyme.to_xyz("nanozyme.xyz")

# TODO: 添加底物并优化几何
# TODO: 过渡态搜索
# TODO: 吸附能计算
```

### 催化活性预测

```python
# TODO: 实现
from nanozyme_mining.assembly import CatalyticActivityPredictor

predictor = CatalyticActivityPredictor()

# 计算吸附能
substrate = Molecule.from_smiles("OO")  # H2O2
adsorption_energy = predictor.compute_adsorption_energy(
    nanozyme, substrate
)

# 搜索过渡态
ts = predictor.search_transition_state(
    nanozyme, substrate, product
)

# 预测活性
activity = predictor.predict_activity(nanozyme)
print(f"Predicted kcat: {activity['kcat']}")
```

## 完整工作流示例

```python
from nanozyme_mining.assembly import (
    MotifLibrary,
    NanozymeAssembler,
    NanozymeValidator,
    MaterialType,
)

# 1. 加载motif库
library = MotifLibrary("./motif_library")
pod_motifs = library.get_by_type("POD")

# 2. 选择组装策略
assembler = NanozymeAssembler(
    strategy="rule",
    material_type=MaterialType.METAL_OXIDE,
    output_dir="./output/pod_nanozymes",
)

# 3. 组装纳米酶
nanozyme = assembler.assemble(
    motifs=pod_motifs[0],
    num_metal_centers=3,
    size=50,
)

# 4. 验证结构
validator = NanozymeValidator()
result = validator.validate(nanozyme)

if result.is_valid:
    print(f"✓ Valid nanozyme: {nanozyme.compute_formula()}")
    
    # 5. 保存
    assembler.save_structure(
        nanozyme,
        "FePOD_nanozyme",
        formats=['xyz', 'json']
    )
    
    # 6. 下游计算
    # TODO: DFT, activity prediction
else:
    print("✗ Invalid structure, skipping...")
```

## 常见问题

### Q: 为什么强调"纳米酶不是蛋白质"？

A: 这是最重要的概念！纳米酶是人工纳米材料，不能直接套用蛋白质的组装方法。主要区别：

| 特性 | 蛋白质酶 | 纳米酶 |
|------|---------|--------|
| 结构单元 | 氨基酸 | 金属、氧化物、碳 |
| 连接方式 | 肽键 | 配位键、共价键、离子键 |
| 结构层次 | 一级到四级结构 | 原子/纳米尺度 |
| 活性中心 | 催化残基 | 金属中心、表面位点 |
| 设计方法 | 蛋白质工程 | 材料化学 |

### Q: 如何选择组装策略？

A:
- **规则组装**: 简单纳米酶，已知材料类型，快速原型
- **扩散模型**: 复杂结构，需要灵活生成，有训练数据
- **模板组装**: 标准拓扑，系统探索

### Q: Motif从哪里来？

A: 三种来源：
1. 从天然酶PDB提取（已实现）
2. 从纳米酶文献手动构建
3. 从实验数据逆向工程

### Q: 如何提高组装质量？

A:
1. 使用高质量motif (confidence > 0.8)
2. 添加更多几何约束
3. 多次组装选择最优
4. 用DFT优化结构

## 下一步

1. [ ] 实现扩散模型组装
2. [ ] 集成DFT计算接口
3. [ ] 催化活性预测模型
4. [ ] 更多材料类型模板
5. [ ] 可视化工具

## 参考资料

- [DiffLinker paper](https://www.nature.com/articles/s42256-024-00815-9)
- [LigandDiff repository](https://github.com/...)
- [stk documentation](https://stk.readthedocs.io)
- 纳米酶综述文章（待补充）

## 联系

如有问题，请查看 `examples/assemble_nanozyme_example.py` 中的完整示例。

