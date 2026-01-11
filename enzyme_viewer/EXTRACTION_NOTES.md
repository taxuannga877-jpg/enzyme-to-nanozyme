# 提取说明文档

## 从 ChemEnzyRetroPlanner 项目中提取的内容

### 1. 界面组件

#### 提取的文件：
- **界面布局**: 参考 `webapp/templates/results.html` 的左右分栏布局
- **样式文件**: 参考 `webapp/static/css/results.css` 的样式设计
- **导航栏**: 参考 `webapp/templates/base.html` 的导航栏设计

#### 保留的功能：
- ✅ 左右分栏布局（左侧控制面板，右侧结构展示）
- ✅ EC号选择下拉框（改为输入框）
- ✅ PDB结构列表展示
- ✅ 酶结构3D可视化区域
- ✅ 催化位点信息表格展示

#### 移除的功能：
- ❌ 化学反应路径网络图（vis-network）
- ❌ 反应节点点击交互
- ❌ 反应式展示
- ❌ 下载结果功能
- ❌ 队列管理功能

### 2. 后端功能

#### 提取的文件：
- **结构可视化**: `webapp/images.py` → `images.py`
  - `get_structure_html_and_active_data()` 函数
  - py3Dmol 3D结构可视化
  
- **工具函数**: `webapp/utils.py` → `utils.py`
  - `convert_easifa_results()` 函数
  - `color_cell()` 和 `strong_cell()` 辅助函数

- **UniProt查询**: `retro_planner/packages/easifa/easifa/interface/utils.py`
  - `UniProtParserEC` 类
  - `query_enzyme_pdb_by_ec()` 方法

#### 新增的功能：
- ✅ `/api/query_ec` - EC号查询接口
- ✅ `/api/get_structure` - 获取结构HTML接口
- ✅ 简化的Flask应用结构

### 3. 核心功能实现

#### EC号查询流程：
1. 用户输入EC号
2. 调用 `UniProtParserEC.query_enzyme_pdb_by_ec()`
3. 从UniProt API查询数据
4. 下载PDB文件到本地
5. 返回PDB列表给前端

#### 结构展示流程：
1. 用户点击"查看结构"按钮
2. 后端读取PDB文件
3. 使用 `get_structure_html_and_active_data()` 生成3Dmol.js HTML
4. 返回HTML和活性位点数据
5. 前端渲染3D结构

### 4. 依赖关系

#### 必需依赖：
- Flask - Web框架
- py3Dmol - 3D结构可视化
- pandas - 数据处理

#### 可选依赖：
- easifa - UniProt查询和催化位点预测（如果可用）

### 5. 数据存储

#### 自动创建的文件夹：
- `data/pdb_files/` - 存储下载的PDB文件
- `data/uniprot_json/` - 存储UniProt JSON数据
- `data/uniprot_csv/` - 存储UniProt CSV数据

### 6. 与原项目的区别

| 功能 | 原项目 | 提取版本 |
|------|--------|----------|
| 化学反应规划 | ✅ | ❌ |
| EC号查询 | ✅ | ✅ |
| PDB结构展示 | ✅ | ✅ |
| 催化位点预测 | ✅ | ⚠️ (可选) |
| 反应路径可视化 | ✅ | ❌ |
| 界面风格 | 完整版 | 简化版（保持一致） |

### 7. 使用建议

1. **首次使用**：
   - 安装依赖：`pip install -r requirements.txt`
   - 运行应用：`python app.py`
   - 访问：http://localhost:5000

2. **EC号查询**：
   - 输入标准EC号格式：`X.X.X.X`（例如：1.11.1.7）
   - 系统会自动从UniProt下载数据

3. **结构查看**：
   - 点击PDB列表中的"查看结构"按钮
   - 3D结构会在右侧展示
   - 可以旋转、缩放查看结构

### 8. 扩展建议

如果需要添加催化位点预测功能：
1. 安装easifa包
2. 在 `app.py` 中初始化 `EasIFAInferenceAPI`
3. 在 `get_structure` 路由中调用预测接口
4. 将预测结果传递给 `get_structure_html_and_active_data()`


