# 从片段生成到纳米酶结构组装的迁移方案

## 核心理念

**关键记忆：纳米酶不是蛋白质！**

纳米酶是人工设计的具有酶活性的纳米材料（金属氧化物、碳材料、金属有机框架等），不是天然酶蛋白质。因此：
- ❌ 不能直接套用蛋白质折叠方法
- ❌ 不能假设氨基酸序列和肽键
- ✅ 需要基于材料化学和配位化学
- ✅ 催化位点是金属中心、金属簇、或表面活性位点

## 参考项目分析

### 1. DiffLinker - 分子连接器设计

**核心思想：**
- 输入：多个断开的3D片段（fragments）
- 输出：连接这些片段的完整分子
- 方法：扩散模型（DDPM）+ 等变图神经网络（EGNN）

**关键组件：**
```python
# 数据表示
- positions: 3D坐标 [N_atoms, 3]
- one_hot: 原子类型 [N_atoms, n_atom_types]
- fragment_mask: 哪些原子是固定片段 [N_atoms]
- linker_mask: 哪些原子是要生成的连接器 [N_atoms]
- anchors: 连接点原子 [N_atoms]

# 生成过程
1. 固定片段位置和类型
2. 使用扩散模型在片段之间生成原子
3. 保持3D等变性（旋转/平移不变）
```

**可迁移特性：**
- ✅ 3D几何约束保持
- ✅ 片段组装逻辑
- ✅ 锚点原子概念
- ❌ 原子类型限制（只支持有机分子C/N/O/S/F）

### 2. LigandDiff - 配体生成

**核心思想：**
- 输入：金属中心 + 部分配体
- 输出：完整的金属配合物
- 方法：扩散模型 + 配体分解

**关键组件：**
```python
# 配位化学特性
- metal center: 金属中心固定
- ligand_breakdown: 自动识别配体单元
- context: 固定的配体部分
- ligand_diff: 要生成的配体部分

# 配位约束
- 金属配位数（coordination number）
- 配位几何（octahedral, square planar等）
- 配体齿合度（denticity）
```

**可迁移特性：**
- ✅ 金属中心处理
- ✅ 配位几何约束
- ✅ 多组分组装
- ✅ 金属-配体键生成

### 3. stk - 分子构建库

**核心思想：**
- 通用的分子组装框架
- 拓扑图（topology graphs）+ 构建块（building blocks）
- 支持各种分子类型：金属配合物、笼状物、聚合物等

**关键组件：**
```python
# 构建范式
class MetalComplex:
    def __init__(
        self,
        topology_graph,  # Octahedral, SquarePlanar等
        building_blocks,  # 金属中心 + 配体
    )

# 顶点定义
- MetalVertex: 金属中心位置
- MonoDentateLigandVertex: 单齿配体
- BiDentateLigandVertex: 双齿配体
```

**可迁移特性：**
- ✅ 拓扑图框架
- ✅ 构建块组合
- ✅ 几何优化
- ⚠️ 需要预定义拓扑（不够灵活）

## 迁移方案

### 阶段1：Motif表示增强

**目标：** 将催化motif扩展为支持纳米酶材料

```python
class NanozymeMotif:
    """纳米酶催化位点motif"""
    
    # 基础信息
    motif_id: str
    nanozyme_type: str  # POD/CAT/SOD等
    
    # 原子信息（扩展支持金属）
    anchor_atoms: List[AnchorAtom]
    # 支持的原子类型：
    # - 金属: Fe, Cu, Mn, Co, Ni, Zn, Ce等
    # - 非金属: C, N, O, S, P
    # - 氧化物: Fe3O4, CeO2, MnO2等的片段表示
    
    # 几何约束
    geometry_constraints: List[GeometryConstraint]
    # - 金属-配体距离
    # - 配位角度
    # - 多金属中心间距
    
    # 材料特性
    material_type: MaterialType  # METAL_OXIDE, CARBON, MOF, METAL_CLUSTER
    coordination_geometry: CoordinationType  # OCTAHEDRAL, TETRAHEDRAL等
    
    # 化学性质
    oxidation_state: Optional[int]  # 金属氧化态
    electron_config: Optional[str]   # 电子构型
```

### 阶段2：片段组装引擎

**架构：** 借鉴DiffLinker + LigandDiff，适配纳米酶

```python
class NanozymeAssembler:
    """纳米酶结构组装器"""
    
    def __init__(
        self,
        motif_library: MotifLibrary,
        diffusion_model: Optional[DiffusionModel] = None,
        rule_engine: ChemicalRuleEngine,
    ):
        self.motif_library = motif_library
        self.diffusion_model = diffusion_model
        self.rule_engine = rule_engine
    
    def assemble(
        self,
        nanozyme_type: str,
        material_type: MaterialType,
        num_active_sites: int = 1,
        scaffold_template: Optional[str] = None,
    ) -> NanozymeStructure:
        """
        组装策略：
        1. 选择motif片段
        2. 确定拓扑结构
        3. 片段定位
        4. 连接器生成（DiffLinker方法）
        5. 结构优化
        """
```

**组装策略：**

#### 策略A：基于规则的组装（Rule-based）

适用场景：简单纳米酶，明确的材料类型

```python
class RuleBasedAssembler:
    """基于化学规则的组装"""
    
    def assemble_metal_oxide_nanozyme(
        self,
        motif: NanozymeMotif,
        size: int,  # 纳米颗粒尺寸
    ):
        """
        组装金属氧化物纳米酶
        
        规则：
        1. 生成晶体结构单元
        2. 从motif提取表面活性位点
        3. 切割纳米颗粒
        4. 将活性位点嵌入表面
        """
        
    def assemble_mof_nanozyme(
        self,
        motif: NanozymeMotif,
        linker_library: List[OrganicLinker],
    ):
        """
        组装MOF纳米酶
        
        规则：
        1. 金属节点（从motif提取）
        2. 有机连接体（从库中选择）
        3. 拓扑匹配（pcu, sql等）
        4. 结构生成
        """
```

#### 策略B：基于扩散模型的生成（Diffusion-based）

适用场景：复杂纳米酶，需要灵活生成

```python
class DiffusionBasedAssembler:
    """基于扩散模型的组装"""
    
    def __init__(self):
        # 借鉴DiffLinker架构
        self.dynamics = NanozymeEGNN(
            in_node_nf=16,  # 扩展到支持金属类型
            n_dims=3,
            hidden_nf=256,
        )
        self.edm = InpaintingEDM(
            dynamics=self.dynamics,
            diffusion_steps=1000,
        )
    
    def assemble_from_motifs(
        self,
        motifs: List[NanozymeMotif],
        target_size: int,
    ):
        """
        从多个motif生成纳米酶
        
        过程：
        1. 将motif转换为3D片段
        2. 片段放置（初始构型）
        3. 扩散生成连接原子/团簇
        4. 化学有效性检查
        """
        
        # 准备数据
        fragments = [self.motif_to_fragment(m) for m in motifs]
        
        # 生成
        data = {
            'positions': fragments_positions,
            'atom_types': fragments_atom_types,
            'fragment_mask': mask_fixed,
            'linker_mask': mask_generate,
            'metal_constraints': metal_coordination_constraints,
        }
        
        # 扩散采样
        generated = self.edm.sample(data)
        
        return self.build_nanozyme(generated)
```

#### 策略C：模板匹配组装（Template-based）

适用场景：已知结构类型的纳米酶

```python
class TemplateBasedAssembler:
    """基于模板的组装（类似stk）"""
    
    def __init__(self):
        self.templates = {
            'single_atom_catalyst': SingleAtomTemplate(),
            'metal_cluster': MetalClusterTemplate(),
            'mof_node': MOFNodeTemplate(),
        }
    
    def assemble_from_template(
        self,
        template_name: str,
        motif: NanozymeMotif,
        parameters: Dict,
    ):
        """
        使用模板组装
        
        例子：
        - 单原子纳米酶：Fe-N4-C
        - 金属团簇：Fe3O4核壳结构
        - MOF节点：Fe3-O-COO
        """
        template = self.templates[template_name]
        return template.build(motif, parameters)
```

### 阶段3：化学有效性验证

**目标：** 确保生成的纳米酶结构化学合理

```python
class NanozymeValidator:
    """纳米酶结构验证器"""
    
    def validate(self, structure: NanozymeStructure) -> ValidationResult:
        """
        验证内容：
        1. 金属配位数
        2. 键长键角
        3. 氧化还原对
        4. 溶剂可接近性
        5. 稳定性估计
        """
        
    def check_metal_coordination(self, metal_atom, neighbors):
        """检查金属配位是否合理"""
        
    def check_oxidation_states(self, structure):
        """检查氧化态是否平衡"""
        
    def check_catalytic_pocket(self, structure, substrate):
        """检查催化口袋是否可容纳底物"""
```

### 阶段4：下游计算对接

**目标：** 与底物进行过渡态搜索和吸附能计算

```python
class CatalyticActivityPredictor:
    """催化活性预测器"""
    
    def __init__(self):
        self.dft_engine = DFTCalculator()  # ASE/Gaussian等
        self.ml_model = ActivityPredictor()  # 机器学习模型
    
    def compute_adsorption_energy(
        self,
        nanozyme: NanozymeStructure,
        substrate: Molecule,
    ) -> float:
        """计算底物吸附能"""
        
    def search_transition_state(
        self,
        nanozyme: NanozymeStructure,
        substrate: Molecule,
        product: Molecule,
    ) -> TransitionState:
        """搜索过渡态"""
        
    def predict_activity(
        self,
        nanozyme: NanozymeStructure,
    ) -> Dict[str, float]:
        """预测催化活性"""
        return {
            'kcat': ...,
            'Km': ...,
            'activation_energy': ...,
        }
```

## 实现路线图

### Phase 1: 基础设施（2周）
- [x] Motif表示扩展
- [ ] 材料类型分类系统
- [ ] 几何约束引擎
- [ ] 化学规则库

### Phase 2: 规则组装（2周）
- [ ] 金属氧化物组装器
- [ ] MOF组装器
- [ ] 单原子催化剂组装器
- [ ] 验证系统

### Phase 3: 扩散模型（4周）
- [ ] 数据集准备（纳米酶训练数据）
- [ ] 模型架构适配（EGNN扩展到金属）
- [ ] 训练pipeline
- [ ] 采样和后处理

### Phase 4: 集成测试（2周）
- [ ] 端到端测试
- [ ] 与DFT计算对接
- [ ] 案例研究（POD/CAT/SOD）
- [ ] 文档和教程

## 技术挑战与解决方案

### 挑战1：原子类型扩展
**问题：** DiffLinker只支持C/N/O/S/F，需要支持金属

**解决方案：**
- 扩展one-hot编码维度
- 添加金属特征（电负性、半径、氧化态）
- 训练金属感知的EGNN

### 挑战2：配位化学约束
**问题：** 金属配位有严格的几何要求

**解决方案：**
- 在扩散过程中添加硬约束
- 使用LigandDiff的配体分解方法
- 后处理几何优化

### 挑战3：多尺度结构
**问题：** 纳米酶从原子到纳米尺度

**解决方案：**
- 分层组装：原子团簇 → 活性位点 → 完整纳米酶
- 粗粒化表示（coarse-grained）
- 多分辨率模型

### 挑战4：训练数据稀缺
**问题：** 纳米酶结构数据远少于有机分子

**解决方案：**
- 迁移学习（从有机分子预训练）
- 数据增强（旋转、平移、扰动）
- 结合物理模拟生成数据

## 示例工作流

```python
# 完整流程示例
from nanozyme_mining import (
    MotifLibrary,
    NanozymeAssembler,
    NanozymeValidator,
    CatalyticActivityPredictor,
)

# Step 1: 加载motif库
motif_lib = MotifLibrary.from_directory("./motif_library")
pod_motifs = motif_lib.get_by_type("POD")

# Step 2: 选择组装策略
assembler = NanozymeAssembler(
    strategy="diffusion",  # 或 "rule" 或 "template"
    material_type="metal_oxide",
)

# Step 3: 组装纳米酶
nanozyme = assembler.assemble(
    motifs=pod_motifs[:3],  # 使用3个POD motif
    target_size=50,  # 50个原子
    composition={"Fe": 10, "O": 20, "C": 15, "N": 5},
)

# Step 4: 验证结构
validator = NanozymeValidator()
validation_result = validator.validate(nanozyme)
if not validation_result.is_valid:
    print(f"Validation failed: {validation_result.errors}")

# Step 5: 计算催化活性
predictor = CatalyticActivityPredictor()
substrate = Molecule.from_smiles("OO")  # H2O2
adsorption_energy = predictor.compute_adsorption_energy(
    nanozyme, substrate
)
print(f"Adsorption energy: {adsorption_energy} eV")

# Step 6: 导出结构
nanozyme.to_xyz("nanozyme_pod.xyz")
nanozyme.to_json("nanozyme_pod.json")
```

## 参考文献

1. Igashov et al. "Equivariant 3D-conditional diffusion model for molecular linker design." Nature Machine Intelligence (2024)
2. LigandDiff paper (add citation)
3. stk: "An Extendable Python Framework for Automated Molecular and Supramolecular Structure Assembly" (2021)
4. 纳米酶相关文献（待补充）

## 总结

**核心要点：**
1. 纳米酶 ≠ 蛋白质：需要材料化学视角
2. 三个项目各有所长：
   - DiffLinker: 片段连接逻辑
   - LigandDiff: 金属配位处理
   - stk: 通用组装框架
3. 混合策略：规则 + 扩散 + 模板
4. 分阶段实现：先简单后复杂
5. 持续验证：化学有效性 + 催化活性

**下一步行动：**
- 完善motif表示
- 实现规则组装器（最容易）
- 准备训练数据
- 开发扩散模型（最有潜力）

