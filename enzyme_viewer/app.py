import os
import json
import sys
from pathlib import Path
from flask import Flask, render_template, jsonify, request
from flask_cors import CORS
import py3Dmol
import pandas as pd
from images import get_structure_html_and_active_data

# 添加nanozyme_mining到Python路径
sys.path.insert(0, str(Path(__file__).parent.parent))
from nanozyme_mining.database.uniprot_fetcher import UniProtFetcher
from nanozyme_mining.prediction.easifa_predictor import EasIFAPredictor, LABEL_TO_SITE_TYPE
from nanozyme_mining.extraction.extractor import MotifExtractor
from enzyme_viewer.motif_db import MotifDatabase, classify_motif

app = Flask(__name__)
CORS(app)

# 配置路径 - 使用本地 cache 数据
BASE_DIR = Path(__file__).parent.parent
app.config['CACHE_DIR'] = BASE_DIR / 'cache'
app.config['JSON_CACHE_DIR'] = app.config['CACHE_DIR'] / 'json'
app.config['PDB_CACHE_DIR'] = app.config['CACHE_DIR'] / 'pdb'  # 旧缓存目录（兼容）
app.config['PDB_LIBRARY_DIR'] = BASE_DIR / 'pdb_library'  # 新的PDB库目录（按EC号组织）
app.config['MOTIF_LIBRARY_DIR'] = BASE_DIR / 'motif_library'
app.config['MOTIF_OUTPUT_DIR'] = BASE_DIR / 'motifs'
app.config['MOTIF_DB_PATH'] = Path(__file__).parent / 'motif_index.db'

# 确保文件夹存在
app.config['JSON_CACHE_DIR'].mkdir(parents=True, exist_ok=True)
app.config['PDB_CACHE_DIR'].mkdir(parents=True, exist_ok=True)
app.config['PDB_LIBRARY_DIR'].mkdir(parents=True, exist_ok=True)
app.config['MOTIF_OUTPUT_DIR'].mkdir(parents=True, exist_ok=True)

print(f"✓ 使用本地缓存数据:")
print(f"  - JSON 缓存: {app.config['JSON_CACHE_DIR']}")
print(f"  - PDB 缓存: {app.config['PDB_CACHE_DIR']} (兼容)")
print(f"  - PDB 库: {app.config['PDB_LIBRARY_DIR']} (主要)")
print(f"  - Motif 库: {app.config['MOTIF_LIBRARY_DIR']}")

# 初始化功能模块
print("✓ 初始化功能模块...")
uniprot_fetcher = UniProtFetcher(cache_dir=str(app.config['CACHE_DIR']))
print("  - UniProt Fetcher 初始化完成")

# EasIFA预测器（延迟初始化，只在需要时加载）
easifa_predictor = None
def get_easifa_predictor():
    global easifa_predictor
    if easifa_predictor is None:
        try:
            print("  - 初始化 EasIFA 预测器...")
            easifa_predictor = EasIFAPredictor(device="cpu")
            print("  - EasIFA 预测器初始化完成")
        except Exception as e:
            print(f"  ⚠️  EasIFA 预测器初始化失败: {e}")
            return None
    return easifa_predictor

# Motif提取器
motif_extractor = MotifExtractor(output_dir=str(app.config['MOTIF_OUTPUT_DIR']))
print("  - Motif Extractor 初始化完成")

# Motif数据库（延迟初始化）
motif_db = None
def get_motif_db():
    """获取Motif数据库实例（延迟初始化）"""
    global motif_db
    if motif_db is None:
        db_path = app.config['MOTIF_DB_PATH']
        if db_path.exists():
            motif_db = MotifDatabase(db_path)
            print("  - Motif数据库已加载")
        else:
            print(f"  ⚠️  Motif数据库不存在: {db_path}")
            print(f"     请运行: python enzyme_viewer/motif_db.py 构建索引")
    return motif_db

# M-CSA查询器（延迟初始化）
mcsa_query = None
def get_mcsa_query_instance():
    global mcsa_query
    if mcsa_query is None:
        try:
            print("  - 初始化 M-CSA 查询器...")
            # 使用项目根目录下的M-CSA数据库文件
            mcsa_file = BASE_DIR / 'data' / 'mcsa_database' / 'mcsa_processed.json'
            from nanozyme_mining.database.mcsa_query import MCSAQuery
            mcsa_query = MCSAQuery(mcsa_file=str(mcsa_file))
            mcsa_query.load()
            print("  - M-CSA 查询器初始化完成")
        except Exception as e:
            print(f"  ⚠️  M-CSA 查询器初始化失败: {e}")
            import traceback
            traceback.print_exc()
            return None
    return mcsa_query

@app.route('/')
def index():
    """主页面"""
    return render_template('index.html')

@app.route('/motif_library')
def motif_library():
    """Motif库浏览页面"""
    return render_template('motif_library.html')

@app.route('/api/list_ec', methods=['GET'])
def list_ec():
    """列出所有可用的EC号"""
    try:
        ec_list = []
        json_dir = app.config['JSON_CACHE_DIR']
        
        # 扫描所有 _sites.json 文件
        for json_file in json_dir.glob("*_sites.json"):
            ec_number = json_file.stem.replace("_sites", "")
            ec_list.append(ec_number)
        
        # 按EC号排序
        ec_list.sort()
        
        return jsonify({
            "status": "success",
            "ec_list": ec_list,
            "total": len(ec_list)
        })
    except Exception as e:
        import traceback
        return jsonify({
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc(),
            "ec_list": []
        }), 500

@app.route('/api/list_nanozyme_types', methods=['GET'])
def list_nanozyme_types():
    """列出所有可用的纳米酶类型（从motif库目录中扫描）"""
    try:
        nanozyme_types = []
        motif_library_dir = app.config['MOTIF_LIBRARY_DIR']
        
        if not motif_library_dir.exists():
            return jsonify({
                "status": "error",
                "error": f"Motif library directory not found: {motif_library_dir}"
            }), 404
        
        # 扫描所有目录，排除EC号格式的目录（只保留nanozyme类型目录）
        for sub_dir in motif_library_dir.iterdir():
            if not sub_dir.is_dir():
                continue
            
            dir_name = sub_dir.name
            # 判断是否为EC号格式（如 1_11_1_6）或nanozyme类型目录
            # EC号格式：数字_数字_数字_数字
            # Nanozyme类型：通常是英文单词，可能包含下划线
            is_ec_format = False
            parts = dir_name.split('_')
            if len(parts) == 4:
                try:
                    # 尝试将每部分转换为数字，如果成功则是EC号格式
                    [int(p) for p in parts]
                    is_ec_format = True
                except ValueError:
                    pass
            
            # 如果不是EC号格式，则认为是nanozyme类型目录
            if not is_ec_format:
                nanozyme_types.append(dir_name)
        
        # 按名称排序
        nanozyme_types.sort()
        
        return jsonify({
            "status": "success",
            "nanozyme_types": nanozyme_types,
            "total": len(nanozyme_types)
        })
    except Exception as e:
        import traceback
        return jsonify({
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500

@app.route('/api/query_ec', methods=['POST'])
def query_ec():
    """根据EC号查询PDB结构文件 - 从本地缓存读取"""
    if not request.is_json:
        return jsonify({"error": "Missing JSON in request"}), 400
    
    data = request.get_json()
    ec_number = data.get('ec_number', '').strip()
    
    if not ec_number:
        return jsonify({"error": "EC number is required"}), 400
    
    try:
        # 从本地 JSON 缓存读取数据
        json_file = app.config['JSON_CACHE_DIR'] / f"{ec_number}_sites.json"
        
        if not json_file.exists():
            return jsonify({
                "status": "success",
                "ec_number": ec_number,
                "pdb_list": [],
                "message": f"本地缓存中没有找到 EC {ec_number} 的数据"
            })
        
        # 读取 JSON 数据
        with open(json_file, 'r', encoding='utf-8') as f:
            enzyme_data = json.load(f)
        
        if not enzyme_data:
            return jsonify({
                "status": "success",
                "ec_number": ec_number,
                "pdb_list": [],
                "message": f"EC {ec_number} 的数据为空"
            })
        
        # 构建PDB列表
        pdb_list = []
        AFDB_VERSION = "6"  # AlphaFold DB version
        
        # 准备PDB库的EC号目录路径（将EC号中的点替换为下划线）
        ec_dir_name = ec_number.replace('.', '_')
        pdb_library_ec_dir = app.config['PDB_LIBRARY_DIR'] / ec_dir_name
        
        for idx, entry in enumerate(enzyme_data):
            uniprot_id = entry.get('uniprot_id', '')
            alphafold_id = entry.get('alphafold_id', uniprot_id)
            sequence = entry.get('sequence', '')
            
            # 构建 PDB 文件路径 - 优先从PDB库查找，然后回退到旧缓存目录
            pdb_path = None
            
            # ========== 优先从PDB库查找（按EC号组织）==========
            if pdb_library_ec_dir.exists():
                # 方法1: 标准格式 AF-{id}-F1-model_v6.pdb
                pdb_filename = f"AF-{alphafold_id}-F1-model_v{AFDB_VERSION}.pdb"
                pdb_path = pdb_library_ec_dir / pdb_filename
                
                # 方法2: 如果不存在，尝试查找任何包含该 ID 的 PDB 文件
                if not pdb_path.exists():
                    matching_pdb = list(pdb_library_ec_dir.glob(f"AF-{alphafold_id}-*.pdb"))
                    if matching_pdb:
                        pdb_path = matching_pdb[0]
                
                # 方法3: 如果还是找不到，尝试使用 uniprot_id
                if not pdb_path or not pdb_path.exists():
                    if uniprot_id and uniprot_id != alphafold_id:
                        pdb_filename = f"AF-{uniprot_id}-F1-model_v{AFDB_VERSION}.pdb"
                        pdb_path = pdb_library_ec_dir / pdb_filename
                        if not pdb_path.exists():
                            matching_pdb = list(pdb_library_ec_dir.glob(f"AF-{uniprot_id}-*.pdb"))
                            if matching_pdb:
                                pdb_path = matching_pdb[0]
                            else:
                                pdb_path = None
                    else:
                        pdb_path = None
            
            # ========== 回退到旧缓存目录（兼容性）==========
            if not pdb_path or not pdb_path.exists():
                # 方法1: 标准格式 AF-{id}-F1-model_v6.pdb
                pdb_filename = f"AF-{alphafold_id}-F1-model_v{AFDB_VERSION}.pdb"
                pdb_path = app.config['PDB_CACHE_DIR'] / pdb_filename
                
                # 方法2: 如果不存在，尝试查找任何包含该 ID 的 PDB 文件
                if not pdb_path.exists():
                    matching_pdb = list(app.config['PDB_CACHE_DIR'].glob(f"AF-{alphafold_id}-*.pdb"))
                    if matching_pdb:
                        pdb_path = matching_pdb[0]
                
                # 方法3: 如果还是找不到，尝试使用 uniprot_id
                if not pdb_path or not pdb_path.exists():
                    if uniprot_id and uniprot_id != alphafold_id:
                        pdb_filename = f"AF-{uniprot_id}-F1-model_v{AFDB_VERSION}.pdb"
                        pdb_path = app.config['PDB_CACHE_DIR'] / pdb_filename
                        if not pdb_path.exists():
                            matching_pdb = list(app.config['PDB_CACHE_DIR'].glob(f"AF-{uniprot_id}-*.pdb"))
                            if matching_pdb:
                                pdb_path = matching_pdb[0]
                            else:
                                pdb_path = None
                    else:
                        pdb_path = None
            
            pdb_info = {
                "id": idx + 1,
                "alphafolddb_id": alphafold_id,
                "uniprot_id": uniprot_id,
                "ec_number": ec_number,
                "pdb_path": str(pdb_path) if pdb_path and pdb_path.exists() else "",
                "sequence_length": len(sequence) if sequence else 0,
                "has_active_sites": len(entry.get('active_sites', [])) > 0
            }
            pdb_list.append(pdb_info)
        
        return jsonify({
            "status": "success",
            "ec_number": ec_number,
            "pdb_list": pdb_list,
            "total": len(pdb_list),
            "message": f"从本地缓存找到 {len(pdb_list)} 个结构"
        })
        
    except Exception as e:
        import traceback
        return jsonify({
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc(),
            "pdb_list": []
        }), 500

@app.route('/api/get_structure', methods=['POST'])
def get_structure():
    """获取酶结构HTML和催化位点信息 - 从本地缓存读取活性位点数据"""
    if not request.is_json:
        return jsonify({"error": "Missing JSON in request"}), 400
    
    data = request.get_json()
    pdb_path = data.get('pdb_path', '')
    ec_number = data.get('ec_number', '')
    uniprot_id = data.get('uniprot_id', '')
    
    if not pdb_path or not os.path.exists(pdb_path):
        return jsonify({"error": f"PDB file not found: {pdb_path}"}), 404
    
    try:
        # 从本地 JSON 缓存读取活性位点数据
        site_labels = None
        active_sites_data = []
        
        json_file = app.config['JSON_CACHE_DIR'] / f"{ec_number}_sites.json"
        if json_file.exists():
            with open(json_file, 'r', encoding='utf-8') as f:
                enzyme_data = json.load(f)
            
            # 找到对应的酶数据
            for entry in enzyme_data:
                if entry.get('uniprot_id') == uniprot_id or entry.get('alphafold_id') == uniprot_id:
                    active_sites = entry.get('active_sites', [])
                    if active_sites:
                        # 构建 site_labels 字典：{residue_index: site_type}
                        # site_type: 0=非活性位点, 1=Binding site, 2=Active site, 3=Other site
                        site_labels = {}
                        for site in active_sites:
                            site_type_str = site.get('type', '').lower()
                            start = site.get('start', 0)
                            end = site.get('end', start)
                            
                            # 映射类型
                            if 'active site' in site_type_str:
                                site_type = 2
                            elif 'binding site' in site_type_str:
                                site_type = 1
                            else:
                                site_type = 3
                            
                            # 标记所有残基
                            for res_idx in range(start, end + 1):
                                site_labels[res_idx] = site_type
                        
                        # 构建活性位点数据表格
                        for site in active_sites:
                            site_type = site.get('type', 'Unknown')
                            start = site.get('start', 0)
                            end = site.get('end', start)
                            description = site.get('description', '')
                            
                            # 确定颜色
                            if 'active site' in site_type.lower():
                                color = "#00B050"  # 绿色
                            elif 'binding site' in site_type.lower():
                                color = "#FF0000"  # 红色
                            else:
                                color = "#FFFF00"  # 黄色
                            
                            if start == end:
                                active_sites_data.append({
                                    "Residue Index": start,
                                    "Residue Name": "",  # 可以从 PDB 文件读取
                                    "Color": color,
                                    "Active Type": site_type,
                                    "Description": description
                                })
                            else:
                                active_sites_data.append({
                                    "Residue Index": f"{start}-{end}",
                                    "Residue Name": "",
                                    "Color": color,
                                    "Active Type": site_type,
                                    "Description": description
                                })
                    break
        
        # 生成结构 HTML
        structure_html, active_data = get_structure_html_and_active_data(
            enzyme_structure_path=pdb_path,
            site_labels=site_labels,  # 使用从缓存读取的活性位点数据
            view_size=(900, 900),
            show_active=(site_labels is not None and len(site_labels) > 0)
        )
        
        # 构建活性位点数据表格 HTML
        active_data_html = ""
        if active_sites_data:
            active_data_df = pd.DataFrame(active_sites_data)
            active_data_html = active_data_df.to_html(
                index=False,
                escape=False,
                classes="table table-striped"
            )
        elif active_data and len(active_data) > 0:
            # 如果 get_structure_html_and_active_data 返回了数据，使用它
            active_data_df = pd.DataFrame(
                active_data,
                columns=["Residue Index", "Residue Name", "Color", "Active Type"]
            )
            active_data_html = active_data_df.to_html(
                index=False,
                escape=False,
                classes="table table-striped"
            )
        else:
            active_data_html = "<p>未找到活性位点数据</p>"
        
        return jsonify({
            "status": "success",
            "structure_html": structure_html,
            "active_data_html": active_data_html,
            "ec_number": ec_number,
            "uniprot_id": uniprot_id,
            "has_active_sites": len(active_sites_data) > 0
        })
        
    except Exception as e:
        import traceback
        return jsonify({
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500

@app.route('/api/download_pdb', methods=['POST'])
def download_pdb():
    """下载PDB文件 - 从AlphaFold数据库"""
    if not request.is_json:
        return jsonify({"error": "Missing JSON in request"}), 400
    
    data = request.get_json()
    alphafold_id = data.get('alphafold_id', '').strip()
    uniprot_id = data.get('uniprot_id', '').strip()
    
    if not alphafold_id and not uniprot_id:
        return jsonify({"error": "AlphaFold ID or UniProt ID is required"}), 400
    
    # 优先使用alphafold_id，否则使用uniprot_id
    download_id = alphafold_id or uniprot_id
    
    try:
        print(f"  Downloading PDB for {download_id}...")
        pdb_path = uniprot_fetcher.download_pdb(download_id)
        
        if pdb_path and pdb_path.exists():
            return jsonify({
                "status": "success",
                "pdb_path": str(pdb_path),
                "message": f"PDB文件下载成功: {pdb_path.name}"
            })
        else:
            return jsonify({
                "status": "error",
                "error": f"无法下载PDB文件，可能AlphaFold数据库中不存在该结构"
            }), 404
            
    except Exception as e:
        import traceback
        return jsonify({
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500

@app.route('/api/predict_active_sites', methods=['POST'])
def predict_active_sites():
    """使用EasIFA模型预测催化位点"""
    if not request.is_json:
        return jsonify({"error": "Missing JSON in request"}), 400
    
    data = request.get_json()
    pdb_path = data.get('pdb_path', '')
    ec_number = data.get('ec_number', '')
    uniprot_id = data.get('uniprot_id', '')
    
    if not pdb_path or not os.path.exists(pdb_path):
        return jsonify({"error": f"PDB file not found: {pdb_path}"}), 404
    
    try:
        # 获取EasIFA预测器
        predictor = get_easifa_predictor()
        if predictor is None:
            return jsonify({
                "error": "EasIFA预测器初始化失败，请检查模型文件是否存在"
            }), 500
        
        print(f"  Predicting active sites for {uniprot_id}...")
        
        # 执行预测
        result = predictor.predict_with_details(
            pdb_path=pdb_path,
            uniprot_id=uniprot_id,
            reaction_smiles="C>>C"  # 默认反应
        )
        
        if result is None or not result.sites:
            return jsonify({
                "status": "success",
                "predicted_sites": [],
                "message": "未预测到活性位点"
            })
        
        # 构建预测结果
        predicted_sites = []
        site_labels = {}  # 用于可视化
        
        for site in result.sites:
            site_data = {
                "residue_index": site.residue_index,
                "residue_name": site.residue_name,
                "site_type": site.site_type,
                "confidence": site.confidence,
                "coordinates": site.coordinates
            }
            predicted_sites.append(site_data)
            
            # 构建site_labels (用于结构可视化)
            # 1=Binding, 2=Catalytic, 3=Other
            if site.site_type == "Catalytic":
                site_labels[site.residue_index] = 2
            elif site.site_type == "Binding":
                site_labels[site.residue_index] = 1
            else:
                site_labels[site.residue_index] = 3
        
        # 保存预测结果到JSON缓存（可选）
        json_file = app.config['JSON_CACHE_DIR'] / f"{ec_number}_sites.json"
        if json_file.exists():
            with open(json_file, 'r') as f:
                enzyme_data = json.load(f)
            
            # 更新对应条目的active_sites
            for entry in enzyme_data:
                if entry.get('uniprot_id') == uniprot_id:
                    # 转换预测结果为active_sites格式
                    entry['active_sites'] = [
                        {
                            "type": "Active site" if s['site_type'] == "Catalytic" else s['site_type'],
                            "start": s['residue_index'],
                            "end": s['residue_index'],
                            "description": f"Predicted by EasIFA (confidence: {s['confidence']:.2f})"
                        }
                        for s in predicted_sites
                    ]
                    break
            
            # 保存更新后的数据
            with open(json_file, 'w') as f:
                json.dump(enzyme_data, f, indent=2)
        
        # 生成预测结果的HTML表格
        active_data_df = pd.DataFrame([
            {
                "Residue Index": s['residue_index'],
                "Residue Name": s['residue_name'],
                "Site Type": s['site_type'],
                "Confidence": f"{s['confidence']:.2f}" if s['confidence'] > 0 else "N/A"
            }
            for s in predicted_sites
        ])
        active_data_html = active_data_df.to_html(
            index=False,
            escape=False,
            classes="table table-striped"
        )
        
        # 重新生成结构HTML（带预测的活性位点）
        structure_html, _ = get_structure_html_and_active_data(
            enzyme_structure_path=pdb_path,
            site_labels=site_labels,
            view_size=(900, 900),
            show_active=True
        )
        
        return jsonify({
            "status": "success",
            "predicted_sites": predicted_sites,
            "active_data_html": active_data_html,
            "structure_html": structure_html,
            "message": f"成功预测到 {len(predicted_sites)} 个活性位点"
        })
        
    except Exception as e:
        import traceback
        return jsonify({
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500

@app.route('/motif_view')
def motif_view():
    """Motif提取展示页面"""
    return render_template('motif_view.html')

@app.route('/api/get_motif', methods=['GET'])
def get_motif():
    """获取Motif数据"""
    motif_id = request.args.get('motif_id', '')
    
    if not motif_id:
        return jsonify({"error": "Motif ID is required"}), 400
    
    try:
        motif_file = app.config['MOTIF_OUTPUT_DIR'] / f"{motif_id}.json"
        
        if not motif_file.exists():
            return jsonify({
                "status": "error",
                "error": f"Motif文件不存在: {motif_id}"
            }), 404
        
        with open(motif_file, 'r') as f:
            motif_data = json.load(f)
        
        return jsonify({
            "status": "success",
            "motif": motif_data
        })
        
    except Exception as e:
        import traceback
        return jsonify({
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500

@app.route('/api/query_mcsa', methods=['POST'])
def query_mcsa():
    """查询M-CSA数据库中的金属位点和催化残基信息"""
    if not request.is_json:
        return jsonify({"error": "Missing JSON in request"}), 400
    
    data = request.get_json()
    ec_number = data.get('ec_number', '').strip()
    
    if not ec_number:
        return jsonify({"error": "EC number is required"}), 400
    
    try:
        # 获取M-CSA查询器
        query = get_mcsa_query_instance()
        if query is None:
            return jsonify({
                "status": "error",
                "error": "M-CSA查询器初始化失败，请检查数据库文件是否存在"
            }), 500
        
        print(f"  Querying M-CSA for EC {ec_number}...")
        
        # 查询金属位点信息
        metal_sites = query.get_metal_sites(ec_number)
        
        # 构建返回数据
        result = {
            "status": "success",
            "ec_number": ec_number,
            "has_metal": metal_sites.get("has_metal", False),
            "metal_types": metal_sites.get("metal_types", []),
            "metal_coordination": metal_sites.get("metal_coordination", []),
            "catalytic_residues": metal_sites.get("catalytic_residues", []),
            "mcsa_references": metal_sites.get("mcsa_references", []),
        }
        
        # 生成HTML表格用于显示
        html_content = generate_mcsa_html(metal_sites)
        result["html_content"] = html_content
        
        message = f"找到 {len(metal_sites.get('mcsa_references', []))} 个M-CSA条目"
        if metal_sites.get("has_metal"):
            message += f"，包含 {len(metal_sites.get('metal_types', []))} 种金属类型"
        result["message"] = message
        
        return jsonify(result)
        
    except Exception as e:
        import traceback
        return jsonify({
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500

def generate_mcsa_html(metal_sites):
    """生成M-CSA数据的HTML展示，重点突出金属位点信息"""
    html_parts = []
    
    # 检查是否有任何数据
    has_metal = metal_sites.get("has_metal", False)
    catalytic_residues = metal_sites.get("catalytic_residues", [])
    mcsa_refs = metal_sites.get("mcsa_references", [])
    
    if not has_metal and len(catalytic_residues) == 0 and len(mcsa_refs) == 0:
        return """
        <div class="alert alert-info" style="margin-top: 15px;">
            <i class="fas fa-info-circle"></i> <strong>查询结果：</strong>未在M-CSA数据库中找到该EC号的相关信息
        </div>
        """
    
    # ========== 金属位点信息（最重要，放在最前面）==========
    html_parts.append("<div class='mcsa-section' style='background-color: #fff3cd; padding: 15px; border-radius: 5px; border: 2px solid #ffc107; margin-bottom: 20px;'>")
    html_parts.append("<h4 style='color: #856404; margin-bottom: 15px;'><i class='fas fa-atom'></i> 金属位点信息</h4>")
    
    if has_metal:
        metal_types = metal_sites.get("metal_types", [])
        if metal_types:
            html_parts.append(f"""
            <div style='background-color: white; padding: 10px; border-radius: 5px; margin-bottom: 15px;'>
                <strong style='color: #856404;'><i class='fas fa-tag'></i> 金属类型：</strong>
                <span style='font-size: 1.1em; font-weight: bold; color: #d9534f;'>{', '.join(metal_types)}</span>
            </div>
            """)
        
        metal_coord = metal_sites.get("metal_coordination", [])
        if metal_coord:
            html_parts.append("<h5 style='color: #856404; margin-top: 15px;'><i class='fas fa-link'></i> 金属配位残基：</h5>")
            coord_data = []
            for coord in metal_coord:
                # 获取该残基配位的金属类型
                metal_types = coord.get("metal_types", [])
                metal_display = ", ".join(metal_types) if metal_types else "未知"
                if metal_types:
                    metal_display = f"<span style='color: #d9534f; font-weight: bold;'>{metal_display}</span>"
                
                coord_data.append({
                    "残基类型": coord.get("residue_type", "N/A"),
                    "残基编号": coord.get("residue_number", "N/A"),
                    "链": coord.get("chain", "N/A"),
                    "PDB ID": coord.get("pdb_id", "N/A"),
                    "配位金属": metal_display,  # 新增：显示配位的具体金属类型
                    "角色": ", ".join(coord.get("roles", [])) if coord.get("roles") else "N/A"
                })
            coord_df = pd.DataFrame(coord_data)
            html_parts.append(coord_df.to_html(index=False, escape=False, classes="table table-sm table-striped table-hover"))
    else:
        html_parts.append("""
        <div class="alert alert-warning" style="margin-top: 10px;">
            <i class="fas fa-exclamation-triangle"></i> <strong>该EC号不含金属位点</strong>
        </div>
        """)
    
    html_parts.append("</div>")
    
    # ========== 催化残基信息 ==========
    if catalytic_residues:
        html_parts.append("<div class='mcsa-section' style='margin-top: 20px; padding: 15px; background-color: #e7f3ff; border-radius: 5px; border: 1px solid #b3d9ff;'>")
        html_parts.append("<h5 style='color: #004085;'><i class='fas fa-flask'></i> 催化残基信息</h5>")
        
        residue_data = []
        for residue in catalytic_residues:
            is_metal_ligand = residue.get("is_metal_ligand", False)
            metal_label = "<span style='color: #d9534f; font-weight: bold;'>是 (金属配体)</span>" if is_metal_ligand else "否"
            
            # 如果是金属配体，显示配位的金属类型
            metal_types = residue.get("metal_types", [])
            metal_type_display = ", ".join(metal_types) if metal_types and is_metal_ligand else ("-" if not is_metal_ligand else "未知")
            if metal_types and is_metal_ligand:
                metal_type_display = f"<span style='color: #d9534f; font-weight: bold;'>{metal_type_display}</span>"
            
            residue_data.append({
                "残基类型": residue.get("residue_type", "N/A"),
                "残基编号": residue.get("residue_number", "N/A"),
                "链": residue.get("chain", "N/A"),
                "PDB ID": residue.get("pdb_id", "N/A"),
                "角色": residue.get("roles_summary", ", ".join(residue.get("roles", []))) if residue.get("roles") else "N/A",
                "金属配体": metal_label,
                "配位金属": metal_type_display  # 新增：显示配位的具体金属类型
            })
        residue_df = pd.DataFrame(residue_data)
        html_parts.append(residue_df.to_html(index=False, escape=False, classes="table table-sm table-striped table-hover"))
        html_parts.append("</div>")
    
    # ========== M-CSA参考条目 ==========
    if mcsa_refs:
        html_parts.append("<div class='mcsa-section' style='margin-top: 20px; padding: 15px; background-color: #f8f9fa; border-radius: 5px; border: 1px solid #dee2e6;'>")
        html_parts.append("<h5><i class='fas fa-book'></i> M-CSA参考条目</h5>")
        html_parts.append(f"<p class='text-muted' style='margin-bottom: 10px;'>共找到 <strong>{len(mcsa_refs)}</strong> 个M-CSA条目</p>")
        html_parts.append("<ul style='list-style-type: none; padding-left: 0;'>")
        for ref in mcsa_refs[:10]:  # 最多显示10个
            url = ref.get("url", "")
            enzyme_name = ref.get("enzyme_name", "Unknown")
            mcsa_id = ref.get("mcsa_id", "")
            uniprot_id = ref.get("uniprot_id", "")
            description = ref.get("description", "")
            
            link = f"<a href='{url}' target='_blank' style='color: #007bff;'><strong>MCSA-{mcsa_id}</strong></a>" if url else f"<strong>MCSA-{mcsa_id}</strong>"
            uniprot_link = f"<a href='https://www.uniprot.org/uniprot/{uniprot_id}' target='_blank' style='color: #007bff;'>{uniprot_id}</a>" if uniprot_id else "N/A"
            
            html_parts.append(f"""
            <li style='padding: 8px; margin-bottom: 8px; background-color: white; border-radius: 3px; border-left: 3px solid #007bff;'>
                {link} - <strong>{enzyme_name}</strong><br>
                <small class='text-muted'>UniProt: {uniprot_link}</small>
                {f'<br><small class="text-muted">{description}</small>' if description else ''}
            </li>
            """)
        if len(mcsa_refs) > 10:
            html_parts.append(f"<li style='padding: 8px;'><em class='text-muted'>... 还有 {len(mcsa_refs) - 10} 个条目</em></li>")
        html_parts.append("</ul>")
        html_parts.append("</div>")
    
    return "\n".join(html_parts)

@app.route('/api/extract_motif', methods=['POST'])
def extract_motif():
    """提取催化Motif"""
    if not request.is_json:
        return jsonify({"error": "Missing JSON in request"}), 400
    
    data = request.get_json()
    pdb_path = data.get('pdb_path', '')
    ec_number = data.get('ec_number', '')
    uniprot_id = data.get('uniprot_id', '')
    nanozyme_type = data.get('nanozyme_type', 'POD')  # 默认POD
    
    if not pdb_path or not os.path.exists(pdb_path):
        return jsonify({"error": f"PDB file not found: {pdb_path}"}), 404
    
    try:
        # 从JSON缓存读取活性位点信息
        active_site_indices = []
        json_file = app.config['JSON_CACHE_DIR'] / f"{ec_number}_sites.json"
        
        if json_file.exists():
            with open(json_file, 'r') as f:
                enzyme_data = json.load(f)
            
            for entry in enzyme_data:
                if entry.get('uniprot_id') == uniprot_id:
                    active_sites = entry.get('active_sites', [])
                    for site in active_sites:
                        start = site.get('start', 0)
                        end = site.get('end', start)
                        active_site_indices.extend(range(start, end + 1))
                    break
        
        print(f"  Extracting motif for {uniprot_id}...")
        print(f"  Active site indices: {active_site_indices}")
        
        # 提取Motif
        motif = motif_extractor.extract_motif(
            pdb_path=pdb_path,
            uniprot_id=uniprot_id,
            ec_number=ec_number,
            nanozyme_type=nanozyme_type,
            active_site_indices=active_site_indices if active_site_indices else None
        )
        
        if motif is None:
            return jsonify({
                "status": "error",
                "error": "未能提取到催化Motif，可能没有找到催化残基"
            }), 404
        
        # 保存Motif到JSON
        motif_dict = motif.to_dict()
        motif_file = app.config['MOTIF_OUTPUT_DIR'] / f"{motif.motif_id}.json"
        with open(motif_file, 'w') as f:
            json.dump(motif_dict, f, indent=2)
        
        # 生成Motif的详细信息HTML
        motif_info = {
            "motif_id": motif.motif_id,
            "uniprot_id": uniprot_id,
            "ec_number": ec_number,
            "nanozyme_type": nanozyme_type,
            "anchor_atoms": [
                {
                    "atom_name": atom.atom_name,
                    "residue_name": atom.residue_name,
                    "residue_number": atom.residue_number,
                    "chain_id": atom.chain_id,
                    "coordinates": atom.coordinates
                }
                for atom in motif.anchor_atoms
            ],
            "geometry_constraints": [
                {
                    "type": constraint.constraint_type,
                    "atoms": constraint.atom_indices,
                    "value": f"{constraint.value:.2f}",
                    "unit": constraint.unit
                }
                for constraint in motif.geometry_constraints
            ]
        }
        
        return jsonify({
            "status": "success",
            "motif": motif_info,
            "motif_file": str(motif_file),
            "message": f"成功提取 {len(motif.anchor_atoms)} 个催化残基"
        })
        
    except Exception as e:
        import traceback
        return jsonify({
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500

@app.route('/api/list_motifs', methods=['POST'])
def list_motifs():
    """列出指定纳米酶类型的所有Motif，按类型分类（优先使用本地数据库）"""
    if not request.is_json:
        return jsonify({"error": "Missing JSON in request"}), 400
    
    data = request.get_json()
    nanozyme_type = data.get('nanozyme_type', '')
    
    if not nanozyme_type:
        return jsonify({"error": "Missing nanozyme type"}), 400
    
    try:
        motifs_by_category = {
            'metal_sites': [],
            'catalytic_sites': [],
            'binding_sites': [],
            'other': []
        }
        total_count = 0
        
        # 优先使用本地数据库
        db = get_motif_db()
        if db:
            # 从数据库查询（需要检查数据库是否有按nanozyme_type查询的方法）
            # 如果没有，则回退到文件系统扫描
            try:
                # 尝试使用数据库的get_by_nanozyme_type方法（如果存在）
                if hasattr(db, 'get_by_nanozyme_type'):
                    db_motifs = db.get_by_nanozyme_type(nanozyme_type)
                else:
                    # 如果没有该方法，获取所有motif然后过滤
                    all_motifs = db.get_all()
                    db_motifs = [m for m in all_motifs if m.get('nanozyme_type', '').upper() == nanozyme_type.upper()]
                
                for db_motif in db_motifs:
                    category = db_motif['category']
                    motif_info = {
                        'motif_id': db_motif['motif_id'],
                        'uniprot_id': db_motif['uniprot_id'],
                        'ec_number': db_motif['ec_number'],
                        'nanozyme_type': db_motif['nanozyme_type'],
                        'anchor_atoms_count': db_motif['anchor_atoms_count'],
                        'category': category,
                        'file_path': db_motif['file_path']
                    }
                    motifs_by_category[category].append(motif_info)
                    total_count += 1
                
                print(f"  ✓ 从数据库查询到 {total_count} 个motif (nanozyme_type: {nanozyme_type})")
            except Exception as db_error:
                print(f"  ⚠️  数据库查询失败，回退到文件系统扫描: {db_error}")
                db = None  # 设置为None以触发文件系统扫描
        
        # 如果数据库查询失败或不存在，使用文件系统扫描
        if not db:
            print(f"  ⚠️  使用文件系统扫描 (nanozyme_type: {nanozyme_type})")
            motif_library_dir = app.config['MOTIF_LIBRARY_DIR']
            
            # 直接查找对应nanozyme类型的目录
            nanozyme_dir = motif_library_dir / nanozyme_type
            
            if not nanozyme_dir.exists() or not nanozyme_dir.is_dir():
                return jsonify({
                    "status": "error",
                    "error": f"Nanozyme type directory not found: {nanozyme_type}"
                }), 404
            
            # 递归查找该目录下的所有JSON文件
            for motif_file in nanozyme_dir.rglob("*.json"):
                try:
                    with open(motif_file, 'r') as f:
                        motif_data = json.load(f)
                    
                    # 验证nanozyme类型是否匹配（不区分大小写）
                    file_nanozyme_type = motif_data.get('nanozyme_type', '').upper()
                    if file_nanozyme_type != nanozyme_type.upper():
                        # 如果文件中的nanozyme_type不匹配，但目录名匹配，也接受
                        # 因为目录名就是nanozyme类型
                        pass
                    
                    # 分类motif
                    category = classify_motif(motif_data)
                    
                    # 准备motif信息
                    motif_info = {
                        'motif_id': motif_data.get('motif_id', ''),
                        'uniprot_id': motif_data.get('source_uniprot_id', ''),
                        'ec_number': motif_data.get('source_ec_number', ''),
                        'nanozyme_type': motif_data.get('nanozyme_type', nanozyme_type),
                        'anchor_atoms_count': len(motif_data.get('anchor_atoms', [])),
                        'category': category,
                        'file_path': str(motif_file)
                    }
                    
                    motifs_by_category[category].append(motif_info)
                    total_count += 1
                    
                except Exception as e:
                    print(f"  ⚠️  Error loading motif {motif_file}: {e}")
                    continue
        
        return jsonify({
            "status": "success",
            "nanozyme_type": nanozyme_type,
            "motifs": motifs_by_category,
            "total_count": total_count,
            "source": "database" if db else "filesystem"
        })
        
    except Exception as e:
        import traceback
        return jsonify({
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500

@app.route('/api/get_motif_structure', methods=['POST'])
def get_motif_structure():
    """获取指定Motif的完整结构信息（优先使用本地数据库）"""
    if not request.is_json:
        return jsonify({"error": "Missing JSON in request"}), 400
    
    data = request.get_json()
    motif_id = data.get('motif_id', '')
    # ec_number 和 nanozyme_type 都是可选的，用于辅助查找
    ec_number = data.get('ec_number', '')
    nanozyme_type = data.get('nanozyme_type', '')
    
    if not motif_id:
        return jsonify({"error": "Missing motif_id"}), 400
    
    try:
        motif_file = None
        
        # 优先使用本地数据库查找文件路径
        db = get_motif_db()
        if db:
            db_motif = db.get_by_id(motif_id)
            if db_motif:
                motif_file = Path(db_motif['file_path'])
                print(f"  ✓ 从数据库找到motif文件: {motif_file}")
        
        # 如果数据库中没有，回退到文件系统扫描
        if not motif_file or not motif_file.exists():
            motif_library_dir = app.config['MOTIF_LIBRARY_DIR']
            
            # 如果提供了nanozyme_type，优先在该目录中查找
            if nanozyme_type:
                nanozyme_dir = motif_library_dir / nanozyme_type
                if nanozyme_dir.exists() and nanozyme_dir.is_dir():
                    potential_file = nanozyme_dir / f"{motif_id}.json"
                    if potential_file.exists():
                        motif_file = potential_file
                    else:
                        # 递归查找
                        for motif_file_path in nanozyme_dir.rglob(f"{motif_id}.json"):
                            motif_file = motif_file_path
                            break
            
            # 如果还没找到，遍历所有目录查找
            if not motif_file or not motif_file.exists():
                for nanozyme_dir in motif_library_dir.iterdir():
                    if not nanozyme_dir.is_dir():
                        continue
                    
                    # 跳过EC号格式的目录（如果已经指定了nanozyme_type）
                    if nanozyme_type:
                        dir_name = nanozyme_dir.name
                        parts = dir_name.split('_')
                        is_ec_format = False
                        if len(parts) == 4:
                            try:
                                [int(p) for p in parts]
                                is_ec_format = True
                            except ValueError:
                                pass
                        if is_ec_format:
                            continue
                    
                    potential_file = nanozyme_dir / f"{motif_id}.json"
                    if potential_file.exists():
                        motif_file = potential_file
                        break
                    
                    # 递归查找
                    for motif_file_path in nanozyme_dir.rglob(f"{motif_id}.json"):
                        motif_file = motif_file_path
                        break
                    
                    if motif_file and motif_file.exists():
                        break
        
        if not motif_file or not motif_file.exists():
            return jsonify({
                "status": "error",
                "error": f"Motif file not found: {motif_id}"
            }), 404
        
        # 加载motif数据
        with open(motif_file, 'r') as f:
            motif_data = json.load(f)
        
        # 格式化数据
        formatted_motif = {
            'motif_id': motif_data.get('motif_id', ''),
            'uniprot_id': motif_data.get('source_uniprot_id', ''),
            'ec_number': motif_data.get('source_ec_number', ''),
            'nanozyme_type': motif_data.get('nanozyme_type', ''),
            'anchor_atoms': motif_data.get('anchor_atoms', []),
            'geometry_constraints': [
                {
                    'constraint_type': c.get('constraint_type', ''),
                    'atom_indices': c.get('atom_indices', []),
                    'value': f"{c.get('value', 0):.2f}",
                    'unit': c.get('unit', '')
                }
                for c in motif_data.get('geometry_constraints', [])
            ],
            'chemistry_tag': motif_data.get('chemistry_tag', ''),
            'reaction_smiles': motif_data.get('reaction_smiles', ''),
            'extraction_method': motif_data.get('extraction_method', ''),
            'confidence_score': motif_data.get('confidence_score', 0.0),
            # 新增字段：残基结构和2D结构
            'residue_structures': motif_data.get('residue_structures', {}),
            'structure_2d_svg': motif_data.get('structure_2d_svg', '')
        }
        
        return jsonify({
            "status": "success",
            "motif": formatted_motif,
            "source": "database" if db and db.get_by_id(motif_id) else "filesystem"
        })
        
    except Exception as e:
        import traceback
        return jsonify({
            "status": "error",
            "error": str(e),
            "traceback": traceback.format_exc()
        }), 500

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)

