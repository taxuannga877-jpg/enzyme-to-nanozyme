# 项目结构

```
enzyme_viewer/
├── app.py                      # Flask主应用文件
├── images.py                   # 酶结构可视化模块（从原项目提取）
├── utils.py                    # 工具函数模块（从原项目提取）
├── requirements.txt            # Python依赖列表
├── README.md                   # 项目说明文档
├── EXTRACTION_NOTES.md         # 提取说明文档
├── PROJECT_STRUCTURE.md        # 本文件
├── run.sh                      # 启动脚本
├── .gitignore                  # Git忽略文件
│
├── templates/                  # HTML模板目录
│   └── index.html              # 主页面模板
│
└── static/                     # 静态文件目录
    ├── css/
    │   ├── bootstrap.min.css   # Bootstrap样式（需要从原项目复制）
    │   └── style.css           # 自定义样式
    └── js/
        └── bootstrap.min.js    # Bootstrap脚本（需要从原项目复制）
```

## 文件说明

### 核心文件

1. **app.py** - Flask应用主文件
   - 定义路由：`/` (主页), `/api/query_ec` (EC查询), `/api/get_structure` (结构获取)
   - 初始化UniProt解析器
   - 处理PDB文件查询和结构可视化

2. **images.py** - 结构可视化模块
   - `get_structure_html_and_active_data()` - 生成3Dmol.js HTML
   - 从原项目 `webapp/images.py` 提取

3. **utils.py** - 工具函数
   - `convert_easifa_results()` - 转换预测结果为表格
   - `color_cell()` 和 `strong_cell()` - 表格样式辅助函数
   - 从原项目 `webapp/utils.py` 提取

### 前端文件

1. **templates/index.html** - 主页面
   - EC号输入框
   - PDB列表展示
   - 3D结构可视化区域
   - JavaScript交互逻辑

2. **static/css/style.css** - 自定义样式
   - 参考原项目 `webapp/static/css/results.css`
   - 保持一致的界面风格

### 配置文件

1. **requirements.txt** - Python依赖
   - Flask, py3Dmol, pandas等

2. **.gitignore** - Git忽略规则
   - 忽略数据文件、缓存等

## 数据目录（自动创建）

运行应用后会自动创建以下目录：

```
data/
├── pdb_files/          # PDB文件存储
├── uniprot_json/       # UniProt JSON数据
└── uniprot_csv/        # UniProt CSV数据
```

## 从原项目提取的代码位置

| 提取的文件 | 原项目位置 |
|-----------|-----------|
| `images.py` | `webapp/images.py` |
| `utils.py` | `webapp/utils.py` |
| 界面布局 | `webapp/templates/results.html` |
| 样式设计 | `webapp/static/css/results.css` |
| UniProt查询 | `retro_planner/packages/easifa/easifa/interface/utils.py` |

## 需要手动复制的文件

如果Bootstrap文件不存在，需要从原项目复制：

```bash
# 从原项目复制Bootstrap文件
cp 参考项目/ChemEnzyRetroPlanner-main/webapp/static/css/bootstrap.min.css static/css/
cp 参考项目/ChemEnzyRetroPlanner-main/webapp/static/js/bootstrap.min.js static/js/
```

或者使用CDN版本（修改HTML模板中的引用）。


