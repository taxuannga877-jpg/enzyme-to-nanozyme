# 酶结构查看器 (Enzyme Structure Viewer)

这是一个从 ChemEnzyRetroPlanner 项目中提取的简化版本，专注于酶结构的可视化展示。

## 功能特性

- ✅ EC号查询：输入EC号查询对应的PDB结构文件
- ✅ PDB结构列表：展示所有匹配的PDB结构文件
- ✅ 3D结构可视化：使用py3Dmol展示酶的三维结构
- ✅ 催化位点信息：显示预测的催化位点信息（如果有）

## 安装

1. 安装依赖：
```bash
pip install -r requirements.txt
```

2. 如果需要使用UniProt查询功能，需要安装easifa包：
```bash
# 从ChemEnzyRetroPlanner项目中复制easifa包
# 或者使用简化版本（不依赖easifa）
```

## 运行

```bash
python app.py
```

然后在浏览器中访问：http://localhost:5000

## 项目结构

```
enzyme_viewer/
├── app.py              # Flask主应用
├── images.py           # 酶结构可视化函数
├── utils.py            # 工具函数
├── templates/          # HTML模板
│   └── index.html
├── static/            # 静态文件
│   ├── css/
│   │   └── style.css
│   └── js/
└── data/              # 数据文件夹（自动创建）
    ├── pdb_files/      # PDB文件存储
    ├── uniprot_json/   # UniProt JSON数据
    └── uniprot_csv/    # UniProt CSV数据
```

## 使用说明

1. 在输入框中输入EC号（例如：1.11.1.7）
2. 点击"查询"按钮或按回车键
3. 左侧会显示所有匹配的PDB结构文件列表
4. 点击"查看结构"按钮查看对应的3D结构和催化位点信息

## 注意事项

- 首次查询某个EC号时，系统会从UniProt下载数据，可能需要一些时间
- PDB文件会自动下载到 `data/pdb_files/` 目录
- 如果easifa包不可用，系统会使用简化模式（不显示催化位点预测）

## 从原项目提取的文件

- `images.py`: 酶结构可视化函数（来自 `webapp/images.py`）
- `utils.py`: 工具函数（来自 `webapp/utils.py`）
- 界面样式：参考 `webapp/static/css/results.css` 和 `webapp/templates/results.html`


