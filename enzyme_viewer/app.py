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
from nanozyme_mining.structure.pdb_parser import PDBParser as ComprehensivePDBParser
from enzyme_viewer.motif_db import MotifDatabase, classify_motif

app = Flask(__name__)
CORS(app)

# 配置路径 - 使用 pdb_library（统一的数据目录）
BASE_DIR = Path(__file__).parent.parent
app.config['PDB_LIBRARY_DIR'] = BASE_DIR / 'pdb_library'  # PDB库目录（按EC号组织，包含JSON和PDB）
app.config['MOTIF_LIBRARY_DIR'] = BASE_DIR / 'motif_library'
app.config['MOTIF_OUTPUT_DIR'] = BASE_DIR / 'motifs'
app.config['MOTIF_DB_PATH'] = Path(__file__).parent / 'motif_index.db'

# 向后兼容：保留旧路径（但不再使用）
app.config['CACHE_DIR'] = BASE_DIR / 'cache'  # 仅用于向后兼容
app.config['JSON_CACHE_DIR'] = app.config['CACHE_DIR'] / 'json'  # 仅用于向后兼容
app.config['PDB_CACHE_DIR'] = app.config['CACHE_DIR'] / 'pdb'  # 仅用于向后兼容

# 确保文件夹存在
app.config['PDB_LIBRARY_DIR'].mkdir(parents=True, exist_ok=True)
app.config['MOTIF_OUTPUT_DIR'].mkdir(parents=True, exist_ok=True)

print(f"✓ 使用本地数据:")
print(f"  - PDB 库: {app.config['PDB_LIBRARY_DIR']} (包含JSON和PDB文件)")

# 初始化功能模块
print("✓ 初始化功能模块...")
uniprot_fetcher = UniProtFetcher(
    cache_dir=str(app.config['CACHE_DIR']),
    pdb_library_dir=str(app.config['PDB_LIBRARY_DIR'])
)
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

@app.route('/test_nanozyme')
def test_nanozyme():
    """测试纳米酶类型API页面"""
    return render_template('test_nanozyme.html')

def get_json_file_path(ec_number: str) -> Path:
    """获取JSON文件路径（优先从pdb_library，向后兼容旧路径）"""
    ec_dir_name = ec_number.replace(".", "_")
    json_file = app.config['PDB_LIBRARY_DIR'] / ec_dir_name / f"{ec_number}_sites.json"
    if not json_file.exists() and app.config['JSON_CACHE_DIR'].exists():
        # 向后兼容：尝试旧路径
        json_file = get_json_file_path(ec_number)
    return json_file

@app.route('/api/list_ec', methods=['GET'])
def list_ec():
    """列出所有可用的EC号"""
    try:
        ec_list = []
        
        # 扫描 pdb_library 中的所有 _sites.json 文件
        for json_file in app.config['PDB_LIBRARY_DIR'].glob("*/*_sites.json"):
            ec_number = json_file.stem.replace("_sites", "")
            ec_list.append(ec_number)
        
        # 向后兼容：如果pdb_library中没有，尝试旧路径
        if not ec_list and app.config['JSON_CACHE_DIR'].exists():
            for json_file in app.config['JSON_CACHE_DIR'].glob("*_sites.json"):
                ec_number = json_file.stem.replace("_sites", "")
                if ec_number not in ec_list:
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
    """列出所有可用的纳米酶类型（优先使用数据库，回退到文件系统扫描）"""
    try:
        nanozyme_types_set = set()
        
        # 优先使用数据库查询（快速）
        db = get_motif_db()
        if db:
            try:
                # 使用专门的方法获取所有纳米酶类型（更快）
                if hasattr(db, 'get_all_nanozyme_types'):
                    nanozyme_types_list = db.get_all_nanozyme_types()
                    nanozyme_types_set = set(nanozyme_types_list)
                else:
                    # 回退到获取所有motif然后提取
                    all_motifs = db.get_all()
                    for motif in all_motifs:
                        nanozyme_type = motif.get('nanozyme_type', '')
                        if nanozyme_type:
                            nanozyme_types_set.add(nanozyme_type)
                print(f"  ✓ 从数据库查询到 {len(nanozyme_types_set)} 种纳米酶类型")
            except Exception as db_error:
                print(f"  ⚠️  数据库查询失败，回退到文件系统扫描: {db_error}")
                db = None
        
        # 如果数据库查询失败或不存在，使用文件系统扫描
        if not db or len(nanozyme_types_set) == 0:
            motif_library_dir = app.config['MOTIF_LIBRARY_DIR']
            
            if not motif_library_dir.exists():
                return jsonify({
                    "status": "error",
                    "error": f"Motif library directory not found: {motif_library_dir}"
                }), 404
            
            print(f"  ⚠️  使用文件系统扫描 (目录: {motif_library_dir})")
            
            # 扫描所有EC号目录下的motif JSON文件，提取nanozyme_type
            for sub_dir in motif_library_dir.iterdir():
                if not sub_dir.is_dir():
                    continue
                
                # 检查是否为EC号格式目录（如 1_11_1_6）
                dir_name = sub_dir.name
                parts = dir_name.split('_')
                is_ec_format = False
                if len(parts) == 4:
                    try:
                        [int(p) for p in parts]
                        is_ec_format = True
                    except ValueError:
                        pass
                
                # 如果是EC号格式目录，扫描其下的所有分类子目录
                if is_ec_format:
                    # 扫描所有分类子目录（metal_sites, catalytic_sites, binding_sites, other）
                    for category_dir in sub_dir.iterdir():
                        if not category_dir.is_dir():
                            continue
                        
                        # 扫描该分类目录下的所有JSON文件
                        for motif_file in category_dir.glob("*.json"):
                            try:
                                with open(motif_file, 'r', encoding='utf-8') as f:
                                    motif_data = json.load(f)
                                    nanozyme_type = motif_data.get('nanozyme_type', '')
                                    if nanozyme_type:
                                        nanozyme_types_set.add(nanozyme_type)
                            except Exception as e:
                                print(f"  ⚠️  读取motif文件失败 {motif_file}: {e}")
                                continue
                else:
                    # 如果不是EC号格式，可能是旧的nanozyme类型目录结构，直接添加
                    nanozyme_types_set.add(dir_name)
        
        # 转换为列表并排序
        nanozyme_types = sorted(list(nanozyme_types_set))
        
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
        # 从本地 JSON 缓存读取数据（优先从pdb_library）
        json_file = get_json_file_path(ec_number)
        
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
            pdb_id = entry.get('pdb_id', '')
            sequence = entry.get('sequence', '')
            
            # 构建 PDB 文件路径 - 优先从PDB库查找，然后回退到旧缓存目录
            # 优先级：实验PDB > AlphaFold PDB
            pdb_path = None
            
            # ========== 优先从PDB库查找（按EC号组织）==========
            if pdb_library_ec_dir.exists():
                # 优先级1: 查找实验PDB (格式: {PDB_ID}.pdb)
                if pdb_id:
                    pdb_id = pdb_id.upper().strip()
                    exp_pdb_filename = f"{pdb_id}.pdb"
                    exp_pdb_path = pdb_library_ec_dir / exp_pdb_filename
                    if exp_pdb_path.exists():
                        pdb_path = exp_pdb_path
                
                # 优先级2: AlphaFold PDB - 标准格式 AF-{id}-F1-model_v6.pdb
                if not pdb_path or not pdb_path.exists():
                    pdb_filename = f"AF-{alphafold_id}-F1-model_v{AFDB_VERSION}.pdb"
                    pdb_path = pdb_library_ec_dir / pdb_filename
                
                # 优先级3: AlphaFold PDB - 查找任何包含该 ID 的 PDB 文件
                if not pdb_path or not pdb_path.exists():
                    matching_pdb = list(pdb_library_ec_dir.glob(f"AF-{alphafold_id}-*.pdb"))
                    if matching_pdb:
                        pdb_path = matching_pdb[0]
                
                # 优先级4: 如果还是找不到，尝试使用 uniprot_id
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
        
        json_file = get_json_file_path(ec_number)
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
        json_file = get_json_file_path(ec_number)
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
        json_file = get_json_file_path(ec_number)
        
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
        # 辅助函数：判断是否为金属元素
        def _is_metal_element(element: str) -> bool:
            """Check if an element symbol represents a metal."""
            if not element:
                return False
            element_upper = element.upper().strip()
            metal_elements = {
                "FE", "CU", "ZN", "MN", "CA", "MG", "CO", "NI",
                "MO", "W", "V", "CR", "TI", "AL", "PB", "HG",
                "CD", "AG", "AU", "PT", "PD", "RU", "RH", "IR",
                "OS", "RE", "TC", "NB", "TA", "HF", "ZR", "Y",
                "SC", "LA", "CE", "PR", "ND", "PM", "SM", "EU",
                "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU",
                "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM",
                "BK", "CF", "ES", "FM", "MD", "NO", "LR", "K", "NA"
            }
            return element_upper in metal_elements
        
        # 辅助函数：判断配体是否为金属
        def _is_metal_ligand(ligand: dict) -> bool:
            """Check if a ligand is a metal (by element or residue_name)."""
            element = ligand.get("element", "").strip()
            res_name = ligand.get("residue_name", "").upper().strip()
            
            # 检查元素是否为金属
            if element and _is_metal_element(element):
                return True
            
            # 检查配体名称是否为金属（如FE, ZN, MG等）
            # 注意：HEM和HEC是复合配体（包含金属的完整分子），不作为纯金属处理
            metal_residue_names = {
                "FE", "FE2", "FE3", "CU", "CU1", "CU2", "ZN", "ZN2",
                "MN", "MN2", "CO", "CO2", "NI", "NI2", "CA", "MG",
                "K", "NA", "MO", "W", "V", "SE"
            }
            if res_name in metal_residue_names:
                return True
            
            return False
        
        # 辅助函数：判断是否为酸根配体
        def _is_acid_anion_ligand(residue_name: str) -> bool:
            """Check if a ligand is an acid anion (SO4, PO4, etc.)."""
            if not residue_name:
                return False
            res_name_upper = residue_name.upper().strip()
            acid_anions = {
                "SO4", "PO4", "NO3", "CO3", "CLO4", "CLO3",
                "SO3", "PO3", "NO2", "CO2", "HCO3", "HPO4",
                "H2PO4", "H3PO4", "HSO4", "H2SO4"
            }
            return res_name_upper in acid_anions
        
        # 辅助函数：获取金属类型
        def _get_metal_type(ligand: dict) -> str:
            """Get metal type from ligand (element or residue_name)."""
            element = ligand.get("element", "").strip().upper()
            res_name = ligand.get("residue_name", "").strip().upper()
            
            # 优先使用元素符号
            if element and _is_metal_element(element):
                return element
            
            # 如果residue_name是金属名称，使用它
            if res_name in {"FE", "FE2", "FE3", "CU", "CU1", "CU2", "ZN", "ZN2",
                           "MN", "MN2", "CO", "CO2", "NI", "NI2", "CA", "MG",
                           "K", "NA", "MO", "W", "V", "SE", "HEM", "HEC"}:
                # 对于复合物如HEM，提取主要金属元素
                if res_name in {"HEM", "HEC"}:
                    return "FE"  # Heme contains Fe
                return res_name
            
            return element or res_name or "UNK"
        
        motifs_by_category = {
            'metal_sites': [],
            'catalytic_sites': [],
            'binding_sites': [],
            'ligands_cofactors': [],
            'other': []
        }
        total_count = 0
        
        # 用于金属去重的字典：metal_type -> metal_entry
        metal_deduplication = {}
        
        # 用于配体去重的字典：ligand_name -> ligand_entry（只保留第一个出现的实例）
        ligand_deduplication = {}
        
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
                    motif_id = db_motif['motif_id']
                    # 检查是否为M-CSA来源的motif
                    is_mcsa = '_mcsa_' in motif_id.lower()
                    
                    motif_info = {
                        'motif_id': motif_id,
                        'uniprot_id': db_motif['uniprot_id'],
                        'ec_number': db_motif['ec_number'],
                        'nanozyme_type': db_motif['nanozyme_type'],
                        'anchor_atoms_count': db_motif['anchor_atoms_count'],
                        'category': category,
                        'file_path': db_motif['file_path'],
                        'source': 'M-CSA' if is_mcsa else 'standard'
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
            
            # 扫描所有EC号目录，查找匹配nanozyme_type的motif
            for sub_dir in motif_library_dir.iterdir():
                if not sub_dir.is_dir():
                    continue
                
                # 检查是否为EC号格式目录（如 1_11_1_6）
                dir_name = sub_dir.name
                parts = dir_name.split('_')
                is_ec_format = False
                if len(parts) == 4:
                    try:
                        [int(p) for p in parts]
                        is_ec_format = True
                    except ValueError:
                        pass
                
                # 如果是EC号格式目录，扫描其下的所有分类子目录
                if is_ec_format:
                    # 扫描所有分类子目录（metal_sites, catalytic_sites, binding_sites, other）
                    for category_dir in sub_dir.iterdir():
                        if not category_dir.is_dir():
                            continue
                        
                        # 扫描该分类目录下的所有JSON文件
                        for motif_file in category_dir.glob("*.json"):
                            try:
                                with open(motif_file, 'r', encoding='utf-8') as f:
                                    motif_data = json.load(f)
                                
                                # 验证nanozyme类型是否匹配（不区分大小写）
                                file_nanozyme_type = motif_data.get('nanozyme_type', '').upper()
                                if file_nanozyme_type != nanozyme_type.upper():
                                    continue  # 不匹配，跳过
                                
                                # 分类motif（使用目录名作为分类，如果motif数据中没有）
                                category = classify_motif(motif_data)
                                # 如果分类失败，使用目录名
                                if category not in motifs_by_category:
                                    category = category_dir.name if category_dir.name in motifs_by_category else 'other'
                                
                                # 准备motif信息
                                motif_id = motif_data.get('motif_id', '')
                                # 检查是否为M-CSA来源的motif
                                source = motif_data.get('source', '')
                                is_mcsa = '_mcsa_' in motif_id.lower() or source == 'M-CSA'
                                
                                motif_info = {
                                    'motif_id': motif_id,
                                    'uniprot_id': motif_data.get('source_uniprot_id', ''),
                                    'ec_number': motif_data.get('source_ec_number', ''),
                                    'nanozyme_type': motif_data.get('nanozyme_type', nanozyme_type),
                                    'anchor_atoms_count': len(motif_data.get('anchor_atoms', [])),
                                    'category': category,
                                    'file_path': str(motif_file),
                                    'source': 'M-CSA' if is_mcsa else 'standard'
                                }
                                
                                motifs_by_category[category].append(motif_info)
                                total_count += 1
                                
                            except Exception as e:
                                print(f"  ⚠️  Error loading motif {motif_file}: {e}")
                                continue
                else:
                    # 如果不是EC号格式，可能是旧的nanozyme类型目录结构
                    if dir_name.upper() == nanozyme_type.upper():
                        # 递归查找该目录下的所有JSON文件
                        for motif_file in sub_dir.rglob("*.json"):
                            try:
                                with open(motif_file, 'r', encoding='utf-8') as f:
                                    motif_data = json.load(f)
                                
                                # 分类motif
                                category = classify_motif(motif_data)
                                
                                # 准备motif信息
                                motif_id = motif_data.get('motif_id', '')
                                # 检查是否为M-CSA来源的motif
                                source = motif_data.get('source', '')
                                is_mcsa = '_mcsa_' in motif_id.lower() or source == 'M-CSA'
                                
                                motif_info = {
                                    'motif_id': motif_id,
                                    'uniprot_id': motif_data.get('source_uniprot_id', ''),
                                    'ec_number': motif_data.get('source_ec_number', ''),
                                    'nanozyme_type': motif_data.get('nanozyme_type', nanozyme_type),
                                    'anchor_atoms_count': len(motif_data.get('anchor_atoms', [])),
                                    'category': category,
                                    'file_path': str(motif_file),
                                    'source': 'M-CSA' if is_mcsa else 'standard'
                                }
                                
                                motifs_by_category[category].append(motif_info)
                                total_count += 1
                            
                            except Exception as e:
                                print(f"  ⚠️  Error loading motif {motif_file}: {e}")
                                continue
        
        # 提取配体/辅因子信息
        try:
            pdb_parser = ComprehensivePDBParser()
            pdb_library_dir = app.config['PDB_LIBRARY_DIR']
            
            # 扫描所有EC号目录，查找匹配nanozyme_type的PDB文件
            for sub_dir in pdb_library_dir.iterdir():
                if not sub_dir.is_dir():
                    continue
                
                # 检查是否为EC号格式目录（如 1_11_1_6）
                dir_name = sub_dir.name
                parts = dir_name.split('_')
                is_ec_format = False
                if len(parts) == 4:
                    try:
                        [int(p) for p in parts]
                        is_ec_format = True
                    except ValueError:
                        pass
                
                if not is_ec_format:
                    continue
                
                # 查找该EC号对应的JSON文件，检查nanozyme_type
                ec_number = dir_name.replace('_', '.')
                json_file = sub_dir / f"{ec_number}_sites.json"
                if not json_file.exists():
                    continue
                
                try:
                    with open(json_file, 'r', encoding='utf-8') as f:
                        enzyme_data = json.load(f)
                    
                    # 检查是否有匹配nanozyme_type的条目
                    has_matching_entry = False
                    for entry in enzyme_data:
                        # 这里需要根据EC号推断nanozyme_type，或者从其他地方获取
                        # 暂时跳过nanozyme_type检查，直接提取所有配体
                        has_matching_entry = True
                        break
                    
                    if not has_matching_entry:
                        continue
                    
                    # 扫描该目录下的所有PDB文件
                    for pdb_file in sub_dir.glob("*.pdb"):
                        try:
                            parsed_data = pdb_parser.parse_pdb_file(pdb_file)
                            ligands = parsed_data.get("ligands", [])
                            
                            if not ligands:
                                continue
                            
                            # 按配体名称、链和残基编号分组
                            ligand_groups = {}
                            for ligand in ligands:
                                res_name = ligand.get("residue_name", "UNK")
                                
                                # 过滤掉水分子（HOH）
                                if res_name.upper() == "HOH":
                                    continue
                                
                                # 过滤掉酸根配体（SO4, PO4等）
                                if _is_acid_anion_ligand(res_name):
                                    continue
                                
                                chain = ligand.get("chain", "")
                                res_num = ligand.get("residue_number")
                                
                                key = f"{res_name}_{chain}_{res_num}"
                                if key not in ligand_groups:
                                    ligand_groups[key] = {
                                        "ligand_name": res_name,
                                        "chain": chain,
                                        "residue_number": res_num,
                                        "pdb_id": pdb_file.stem.replace(".pdb", "").upper(),
                                        "atom_count": 0,
                                        "atoms": [],
                                        "is_metal": False
                                    }
                                
                                ligand_groups[key]["atom_count"] += 1
                                ligand_groups[key]["atoms"].append(ligand)
                                
                                # 检查是否为金属配体
                                if _is_metal_ligand(ligand):
                                    ligand_groups[key]["is_metal"] = True
                            
                            # 从JSON文件中获取uniprot_id和ec_number
                            uniprot_id = ""
                            if enzyme_data:
                                uniprot_id = enzyme_data[0].get("uniprot_id", "")
                            
                            # 处理配体：金属配体整合到metal_sites，其他配体添加到ligands_cofactors
                            for ligand_key, ligand_info in ligand_groups.items():
                                res_name = ligand_info['ligand_name']
                                
                                # 如果是金属配体，整合到metal_sites分类
                                if ligand_info["is_metal"]:
                                    # 获取金属类型
                                    metal_type = _get_metal_type(ligand_info["atoms"][0])
                                    
                                    # 金属去重：每种金属类型只保留一个条目（第一个出现的实例）
                                    if metal_type not in metal_deduplication:
                                        metal_entry = {
                                            'metal_type': metal_type,
                                            'metal_name': res_name,
                                            'ligand_name': res_name,  # 和配体条目保持一致
                                            'ligand_id': f"{metal_type}_metal",  # 金属的唯一标识
                                            'uniprot_id': uniprot_id,
                                            'ec_number': ec_number,
                                            'occurrence_count': 1,  # 记录出现次数
                                            'atom_count': ligand_info['atom_count'],  # 只保留第一个实例的原子数
                                            'chain': ligand_info['chain'],  # 保存第一个实例的chain
                                            'residue_number': ligand_info['residue_number'],  # 保存第一个实例的residue_number
                                            'pdb_ids': set(),
                                            'file_paths': [],
                                            'category': 'metal_sites',
                                            'source': 'standard'
                                        }
                                        metal_deduplication[metal_type] = metal_entry
                                        metal_deduplication[metal_type]['pdb_ids'].add(ligand_info['pdb_id'])
                                        if str(pdb_file) not in metal_deduplication[metal_type]['file_paths']:
                                            metal_deduplication[metal_type]['file_paths'].append(str(pdb_file))
                                    else:
                                        # 如果已经存在，只增加出现次数计数，不累计原子数
                                        metal_deduplication[metal_type]['occurrence_count'] += 1
                                        metal_deduplication[metal_type]['pdb_ids'].add(ligand_info['pdb_id'])
                                        if str(pdb_file) not in metal_deduplication[metal_type]['file_paths']:
                                            metal_deduplication[metal_type]['file_paths'].append(str(pdb_file))
                                else:
                                    # 非金属配体：只保留第一个出现的实例
                                    ligand_name = ligand_info['ligand_name'].upper()
                                    
                                    # 如果这个配体类型还没有被处理过，创建条目
                                    if ligand_name not in ligand_deduplication:
                                        ligand_id = f"{ligand_info['pdb_id']}_{ligand_info['ligand_name']}_{ligand_info['chain']}_{ligand_info['residue_number']}"
                                        
                                        ligand_entry = {
                                            'ligand_id': ligand_id,
                                            'ligand_name': ligand_info['ligand_name'],
                                            'pdb_id': ligand_info['pdb_id'],
                                            'uniprot_id': uniprot_id,
                                            'ec_number': ec_number,
                                            'atom_count': ligand_info['atom_count'],
                                            'chain': ligand_info['chain'],
                                            'residue_number': ligand_info['residue_number'],
                                            'file_path': str(pdb_file),
                                            'category': 'ligands_cofactors',
                                            'source': 'standard',
                                            'occurrence_count': 1  # 记录出现次数
                                        }
                                        
                                        ligand_deduplication[ligand_name] = ligand_entry
                                        motifs_by_category['ligands_cofactors'].append(ligand_entry)
                                        total_count += 1
                                    else:
                                        # 如果已经存在，只增加出现次数计数
                                        ligand_deduplication[ligand_name]['occurrence_count'] += 1
                        
                        except Exception as e:
                            print(f"  ⚠️  Error parsing PDB file {pdb_file}: {e}")
                            continue
                
                except Exception as e:
                    print(f"  ⚠️  Error reading JSON file {json_file}: {e}")
                    continue
        
        except Exception as e:
            print(f"  ⚠️  Error extracting ligands/cofactors: {e}")
            import traceback
            traceback.print_exc()
        
        # 将去重后的金属条目添加到metal_sites分类
        for metal_type, metal_entry in metal_deduplication.items():
            # 转换set为list以便JSON序列化
            metal_entry['pdb_ids'] = sorted(list(metal_entry['pdb_ids']))
            # 创建metal_id和motif_id（前端需要motif_id）
            metal_id = f"{metal_type}_metal_{metal_entry['occurrence_count']}occ"
            metal_entry['metal_id'] = metal_id
            metal_entry['motif_id'] = metal_id  # 前端需要motif_id字段
            metal_entry['anchor_atoms_count'] = 0  # 金属位点没有锚点原子
            # 确保ligand_name存在（用于前端显示）
            if 'ligand_name' not in metal_entry:
                metal_entry['ligand_name'] = metal_entry.get('metal_name', metal_type)
            motifs_by_category['metal_sites'].append(metal_entry)
            total_count += 1
        
        # 过滤掉未分类（other）的条目
        if 'other' in motifs_by_category:
            # 从total_count中减去other分类的数量
            other_count = len(motifs_by_category['other'])
            total_count -= other_count
            # 清空other分类
            motifs_by_category['other'] = []
        
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

def generate_pdb_from_residues(residue_structures, anchor_atoms):
    """
    从residue_structures生成PDB格式字符串
    用于3Dmol.js可视化完整的分子结构
    
    对于配体/小分子（如HEM），只保留第一个出现的实例，不合并所有出现的原子
    """
    pdb_lines = []
    atom_serial = 1
    
    # 标准氨基酸列表
    standard_aa = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL", "SEC", "PYL"
    }
    
    # 创建锚定原子的映射，用于后续高亮
    anchor_atom_map = {}
    for anchor in anchor_atoms:
        key = (anchor['residue_name'], anchor['residue_number'], anchor['atom_name'])
        anchor_atom_map[key] = anchor
    
    # 用于跟踪已处理的配体/小分子（只保留第一个出现的实例）
    processed_ligands = set()
    
    # 先处理标准氨基酸残基（保留所有实例）
    for key_str, residue in residue_structures.items():
        # 解析key (格式: "RESNAME_RESNUM")
        parts = key_str.split('_')
        if len(parts) < 2:
            continue
        res_name = parts[0].upper()
        
        # 如果是标准氨基酸，处理所有实例
        if res_name in standard_aa:
            try:
                res_num = int(parts[1])
            except ValueError:
                continue
            
            chain_id = residue.get('chain_id', 'A')
            atoms = residue.get('atoms', [])
            
            # 添加每个原子
            for atom_info in atoms:
                atom_name = atom_info.get('atom_name', '')
                element = atom_info.get('element', '')
                coords = atom_info.get('coordinates', [0, 0, 0])
                occupancy = atom_info.get('occupancy', 1.0)
                bfactor = atom_info.get('bfactor', 0.0)
                
                if len(coords) < 3:
                    continue
                
                x, y, z = coords[0], coords[1], coords[2]
                
                # 格式化PDB ATOM行
                atom_name_padded = atom_name[:4].ljust(4)
                res_name_padded = res_name[:3].ljust(3)
                element_padded = element[:2].rjust(2) if element else "  "
                
                pdb_line = (
                    f"ATOM  {atom_serial:5d} {atom_name_padded:4s} {res_name_padded:3s} {chain_id:1s}"
                    f"{res_num:4d}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}          {element_padded:2s}  \n"
                )
                pdb_lines.append(pdb_line)
                atom_serial += 1
    
    # 再处理配体/小分子（只保留第一个出现的实例）
    for key_str, residue in residue_structures.items():
        # 解析key (格式: "RESNAME_RESNUM")
        parts = key_str.split('_')
        if len(parts) < 2:
            continue
        res_name = parts[0].upper()
        
        # 如果是配体/小分子（非标准氨基酸），只保留第一个出现的实例
        if res_name not in standard_aa:
            # 检查是否已经处理过这个配体类型
            if res_name in processed_ligands:
                continue  # 跳过，只保留第一个出现的实例
            
            processed_ligands.add(res_name)
            
            try:
                res_num = int(parts[1])
            except ValueError:
                continue
            
            chain_id = residue.get('chain_id', 'A')
            atoms = residue.get('atoms', [])
            
            # 添加每个原子（只处理第一个出现的配体实例）
            for atom_info in atoms:
                atom_name = atom_info.get('atom_name', '')
                element = atom_info.get('element', '')
                coords = atom_info.get('coordinates', [0, 0, 0])
                occupancy = atom_info.get('occupancy', 1.0)
                bfactor = atom_info.get('bfactor', 0.0)
                
                if len(coords) < 3:
                    continue
                
                x, y, z = coords[0], coords[1], coords[2]
                
                # 格式化PDB HETATM行（配体使用HETATM）
                atom_name_padded = atom_name[:4].ljust(4)
                res_name_padded = res_name[:3].ljust(3)
                element_padded = element[:2].rjust(2) if element else "  "
                
                pdb_line = (
                    f"HETATM{atom_serial:5d} {atom_name_padded:4s} {res_name_padded:3s} {chain_id:1s}"
                    f"{res_num:4d}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}          {element_padded:2s}  \n"
                )
                pdb_lines.append(pdb_line)
                atom_serial += 1
    
    # 添加CONECT记录（基于距离约束，这里简化处理）
    # 注意：完整的CONECT需要知道键连接信息，这里先不添加
    
    pdb_lines.append("END\n")
    return "".join(pdb_lines)

@app.route('/api/get_ligand_structure', methods=['POST'])
def get_ligand_structure():
    """获取指定配体的3D结构数据"""
    if not request.is_json:
        return jsonify({"error": "Missing JSON in request"}), 400
    
    data = request.get_json()
    ligand_id = data.get('ligand_id', '')
    pdb_path = data.get('pdb_path', '')
    ligand_name = data.get('ligand_name', '')
    chain = data.get('chain', '')
    residue_number = data.get('residue_number', '')
    
    if not pdb_path or not os.path.exists(pdb_path):
        return jsonify({"error": f"PDB file not found: {pdb_path}"}), 404
    
    if not ligand_name:
        return jsonify({"error": "Ligand name is required"}), 400
    
    try:
        pdb_parser = ComprehensivePDBParser()
        parsed_data = pdb_parser.parse_pdb_file(Path(pdb_path))
        ligands = parsed_data.get("ligands", [])
        
        # 筛选出匹配的配体原子（排除水分子HOH）
        matching_atoms = []
        for ligand in ligands:
            res_name = ligand.get("residue_name", "").upper()
            
            # 过滤掉水分子（HOH）
            if res_name == "HOH":
                continue
            
            if (res_name == ligand_name.upper() and
                ligand.get("chain", "") == chain and
                ligand.get("residue_number") == residue_number):
                matching_atoms.append(ligand)
        
        if not matching_atoms:
            return jsonify({
                "status": "error",
                "error": f"Ligand {ligand_name} not found in PDB file"
            }), 404
        
        # 生成PDB格式字符串（仅包含该配体的原子）
        pdb_lines = []
        atom_serial = 1
        
        for atom in matching_atoms:
            atom_name = atom.get("atom_name", "")
            element = atom.get("element", "")
            coords = atom.get("coordinates", [0, 0, 0])
            
            if len(coords) < 3:
                continue
            
            x, y, z = coords[0], coords[1], coords[2]
            
            # 格式化PDB HETATM行
            atom_name_padded = atom_name[:4].ljust(4)
            res_name_padded = ligand_name[:3].ljust(3)
            element_padded = element[:2].rjust(2) if element else "  "
            
            pdb_line = (
                f"HETATM{atom_serial:5d} {atom_name_padded:4s} {res_name_padded:3s} {chain:1s}"
                f"{residue_number:4d}   {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element_padded:2s}  \n"
            )
            pdb_lines.append(pdb_line)
            atom_serial += 1
        
        pdb_lines.append("END\n")
        pdb_string = "".join(pdb_lines)
        
        return jsonify({
            "status": "success",
            "ligand_id": ligand_id,
            "ligand_name": ligand_name,
            "pdb_string": pdb_string,
            "atom_count": len(matching_atoms)
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
    # 前端可能传递的额外信息（用于金属位点）
    file_paths_from_frontend = data.get('file_paths', [])
    ligand_name_from_frontend = data.get('ligand_name', '')
    metal_type_from_frontend = data.get('metal_type', '')
    chain_from_frontend = data.get('chain') or ''
    residue_number_from_frontend = data.get('residue_number') or ''
    
    # 处理null值（前端可能传递字符串"null"或null）
    if chain_from_frontend in [None, 'null', 'None', '']:
        chain_from_frontend = ''
    if residue_number_from_frontend in [None, 'null', 'None', '']:
        residue_number_from_frontend = ''
    
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
        
        # 如果找不到文件，检查是否是金属条目（金属条目没有对应的JSON文件）
        if not motif_file or not motif_file.exists():
            # 检查是否是金属条目ID（格式：{metal_type}_metal_{count}occ）
            if '_metal_' in motif_id:
                # 这是金属条目，需要从PDB文件中提取结构
                # 从motif_id中提取金属类型
                metal_type = motif_id.split('_metal_')[0]
                
                # 尝试从数据库或前端传递的信息获取金属条目信息
                metal_entry = None
                file_paths = file_paths_from_frontend or []
                
                # 优先从数据库获取完整信息
                chain = None
                residue_number = None
                if db:
                    db_motif = db.get_by_id(motif_id)
                    if db_motif:
                        metal_entry = db_motif
                        # 从数据库获取metal_name和ligand_name（优先使用ligand_name，因为它更准确）
                        metal_name = db_motif.get('ligand_name') or db_motif.get('metal_name') or metal_type
                        if not file_paths:
                            file_paths = db_motif.get('file_paths', [])
                        # 从数据库获取chain和residue_number（用于精确匹配）
                        chain = db_motif.get('chain', '')
                        residue_number = db_motif.get('residue_number', '')
                        print(f"  [Metal Site] Found in database: metal_name={metal_name}, chain={chain}, residue_number={residue_number}, file_paths={len(file_paths)} files")
                
                # 如果数据库中没有，使用前端传递的信息
                if not metal_entry:
                    metal_name = ligand_name_from_frontend or metal_type  # 优先使用前端传递的ligand_name
                    # 从前端获取chain和residue_number（如果数据库中没有）
                    if not chain:
                        chain = chain_from_frontend
                    if not residue_number:
                        residue_number = residue_number_from_frontend
                    print(f"  [Metal Site] Using frontend data: metal_type={metal_type}, metal_name={metal_name}, chain={chain}, residue_number={residue_number}, file_paths_from_frontend={len(file_paths)} files")
                
                print(f"  [Metal Site] Final: metal_type={metal_type}, metal_name={metal_name}, chain={chain}, residue_number={residue_number}, file_paths={len(file_paths)} files")
                
                # 如果还是没有，尝试从文件系统中查找第一个PDB文件
                # 注意：PDB文件在PDB_LIBRARY_DIR中，不在MOTIF_LIBRARY_DIR中
                if not file_paths:
                    pdb_library_dir = app.config['PDB_LIBRARY_DIR']
                    print(f"  [Metal Site] Searching for PDB files in: {pdb_library_dir}")
                    
                    # 方法1：如果有nanozyme_type，尝试在对应的EC号目录中查找
                    if nanozyme_type:
                        # 先尝试在nanozyme_type目录中查找（如果存在）
                        nanozyme_dir = pdb_library_dir / nanozyme_type
                        if nanozyme_dir.exists() and nanozyme_dir.is_dir():
                            pdb_files = list(nanozyme_dir.rglob("*.pdb"))
                            if pdb_files:
                                file_paths = [str(pdb_files[0])]
                                print(f"  [Metal Site] Found PDB file in nanozyme_type directory: {file_paths[0]}")
                    
                    # 方法2：如果还没找到，扫描所有EC号目录，查找包含该金属的PDB文件
                    if not file_paths:
                        print(f"  [Metal Site] Scanning all EC directories for PDB files...")
                        for sub_dir in pdb_library_dir.iterdir():
                            if not sub_dir.is_dir():
                                continue
                            
                            # 检查是否为EC号格式目录（如 1_11_1_6）
                            dir_name = sub_dir.name
                            parts = dir_name.split('_')
                            is_ec_format = False
                            if len(parts) == 4:
                                try:
                                    [int(p) for p in parts]
                                    is_ec_format = True
                                except ValueError:
                                    pass
                            
                            if is_ec_format:
                                pdb_files = list(sub_dir.glob("*.pdb"))
                                if pdb_files:
                                    # 使用第一个找到的PDB文件
                                    file_paths = [str(pdb_files[0])]
                                    print(f"  [Metal Site] Found PDB file in EC directory {dir_name}: {file_paths[0]}")
                                    break
                    
                    # 方法3：如果还是没找到，尝试在整个PDB_LIBRARY_DIR中递归查找
                    if not file_paths:
                        print(f"  [Metal Site] Recursively searching all PDB files...")
                        all_pdb_files = list(pdb_library_dir.rglob("*.pdb"))
                        if all_pdb_files:
                            file_paths = [str(all_pdb_files[0])]
                            print(f"  [Metal Site] Found PDB file: {file_paths[0]}")
                    
                    if not file_paths:
                        print(f"  [Metal Site] WARNING: No PDB files found in {pdb_library_dir}")
                
                # 尝试从第一个PDB文件中提取配体结构
                pdb_string = ""
                print(f"  [Metal Site] file_paths={file_paths}, len={len(file_paths) if file_paths else 0}")
                if file_paths and len(file_paths) > 0:
                        first_pdb_path = Path(file_paths[0])
                        print(f"  [Metal Site] Trying to extract from: {first_pdb_path}, exists={first_pdb_path.exists()}")
                        if first_pdb_path.exists():
                            try:
                                # 使用ComprehensivePDBParser解析PDB文件
                                pdb_parser = ComprehensivePDBParser()
                                parsed_data = pdb_parser.parse_pdb_file(first_pdb_path)
                                ligands = parsed_data.get("ligands", [])
                                print(f"  [Metal Site] Found {len(ligands)} ligands in PDB file")
                                if ligands:
                                    print(f"  [Metal Site] First ligand keys: {list(ligands[0].keys())}")
                                    print(f"  [Metal Site] First ligand residue_name: {ligands[0].get('residue_name', 'N/A')}")
                                
                                # 查找匹配的金属配体（使用ligand_name或metal_name，以及chain和residue_number进行精确匹配）
                                matching_atoms = []
                                print(f"  [Metal Site] Searching for metal_name={metal_name.upper()}, metal_type={metal_type.upper()}, chain={chain}, residue_number={residue_number}")
                                
                                # 构建搜索名称列表：包括metal_name、metal_type，以及特殊处理（如HEM->FE）
                                search_names = [metal_name.upper(), metal_type.upper()]
                                # 特殊处理：如果metal_type是FE，也尝试匹配HEM和HEC
                                if metal_type.upper() == "FE":
                                    search_names.extend(["HEM", "HEC"])
                                # 特殊处理：如果metal_name是HEM或HEC，也尝试匹配FE
                                if metal_name.upper() in ["HEM", "HEC"]:
                                    search_names.append("FE")
                                
                                search_names = list(set(search_names))  # 去重
                                print(f"  [Metal Site] Search names: {search_names}")
                                
                                # 精确匹配：如果有chain和residue_number，使用它们；否则只匹配residue_name
                                for ligand in ligands:
                                    res_name = ligand.get("residue_name", "").upper()
                                    ligand_chain = ligand.get("chain", "")
                                    ligand_res_num = ligand.get("residue_number", "")
                                    
                                    # 检查residue_name是否匹配
                                    if res_name in search_names:
                                        # 如果有chain和residue_number信息（非空），进行精确匹配（就像配体目录那样）
                                        if chain and residue_number and chain.strip() and str(residue_number).strip():
                                            # 转换residue_number为字符串进行比较（兼容不同格式）
                                            try:
                                                ligand_res_num_str = str(ligand_res_num) if ligand_res_num else ""
                                                residue_number_str = str(residue_number) if residue_number else ""
                                                if (ligand_chain == chain and 
                                                    ligand_res_num_str == residue_number_str):
                                                    matching_atoms.append(ligand)
                                                    print(f"  [Metal Site] Found exact match: {res_name}, chain={ligand_chain}, residue_number={ligand_res_num}")
                                            except Exception as e:
                                                print(f"  [Metal Site] Error comparing residue_number: {e}")
                                        else:
                                            # 如果没有chain和residue_number信息，只匹配residue_name（向后兼容）
                                            matching_atoms.append(ligand)
                                            print(f"  [Metal Site] Found matching ligand (name only): {res_name}")
                                
                                # 如果还是没找到，列出所有可用的配体名称
                                if not matching_atoms:
                                    unique_res_names = sorted(set(ligand.get("residue_name", "").upper() for ligand in ligands))
                                    print(f"  [Metal Site] No match found. Available residue names in PDB: {unique_res_names}")
                                    if chain and residue_number:
                                        print(f"  [Metal Site] Tried to match chain={chain}, residue_number={residue_number}")
                                        # 列出所有匹配residue_name但chain/residue_number不匹配的配体
                                        for ligand in ligands:
                                            res_name = ligand.get("residue_name", "").upper()
                                            if res_name in search_names:
                                                print(f"  [Metal Site] Found {res_name} but chain={ligand.get('chain')}, residue_number={ligand.get('residue_number')} (not matching)")
                                
                                # 生成PDB格式字符串（只保留第一个出现的实例）
                                if matching_atoms:
                                    print(f"  [Metal Site] Found {len(matching_atoms)} matching atoms")
                                    pdb_lines = []
                                    atom_serial = 1
                                    
                                    # 找到第一个配体实例的chain和residue_number
                                    first_atom = matching_atoms[0]
                                    first_chain = first_atom.get("chain", "")
                                    first_res_num_raw = first_atom.get("residue_number", "")
                                    # 转换residue_number为整数（如果是字符串）
                                    try:
                                        first_res_num = int(first_res_num_raw) if first_res_num_raw else 1
                                    except (ValueError, TypeError):
                                        first_res_num = 1
                                    
                                    print(f"  [Metal Site] First instance: chain={first_chain}, residue_number={first_res_num}")
                                    
                                    # 只处理第一个配体实例的所有原子
                                    atoms_processed = 0
                                    atoms_skipped = 0
                                    for atom in matching_atoms:
                                        atom_chain = atom.get("chain", "")
                                        res_num_raw = atom.get("residue_number", "")
                                        
                                        # 转换residue_number为整数（如果是字符串）
                                        try:
                                            res_num = int(res_num_raw) if res_num_raw else 1
                                        except (ValueError, TypeError):
                                            res_num = 1
                                        
                                        # 只保留第一个出现的配体实例（通过chain和residue_number匹配）
                                        if atom_chain != first_chain or res_num != first_res_num:
                                            atoms_skipped += 1
                                            continue  # 跳过后续出现的实例
                                        
                                        atom_name = atom.get("atom_name", "")
                                        element = atom.get("element", "")
                                        coords = atom.get("coordinates", [0, 0, 0])
                                        
                                        if len(coords) < 3:
                                            print(f"  [Metal Site] Warning: Atom {atom_name} has invalid coordinates: {coords}")
                                            continue
                                        
                                        x, y, z = coords[0], coords[1], coords[2]
                                        
                                        # 格式化PDB HETATM行
                                        # 使用实际的residue_name（从PDB文件中读取），而不是metal_name
                                        actual_res_name = atom.get("residue_name", metal_name)
                                        atom_name_padded = atom_name[:4].ljust(4)
                                        res_name_padded = actual_res_name[:3].ljust(3)
                                        element_padded = element[:2].rjust(2) if element else "  "
                                        
                                        pdb_line = (
                                            f"HETATM{atom_serial:5d} {atom_name_padded:4s} {res_name_padded:3s} {atom_chain:1s}"
                                            f"{res_num:4d}   {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element_padded:2s}  \n"
                                        )
                                        pdb_lines.append(pdb_line)
                                        atom_serial += 1
                                        atoms_processed += 1
                                    
                                    print(f"  [Metal Site] Processed {atoms_processed} atoms, skipped {atoms_skipped} atoms")
                                    
                                    if pdb_lines:
                                        pdb_lines.append("END\n")
                                        pdb_string = "".join(pdb_lines)
                                        # 如果找到了配体，更新metal_name为实际找到的配体名称（例如HEM而不是FE）
                                        if matching_atoms:
                                            actual_ligand_name = matching_atoms[0].get("residue_name", metal_name)
                                            if actual_ligand_name.upper() != metal_name.upper():
                                                print(f"  [Metal Site] Updating metal_name from {metal_name} to {actual_ligand_name}")
                                                metal_name = actual_ligand_name
                                        print(f"  [Metal Site] Successfully extracted {len(pdb_lines)-1} atoms, PDB string length={len(pdb_string)}")
                                    else:
                                        print(f"  [Metal Site] ERROR: No atoms processed! matching_atoms={len(matching_atoms)}, atoms_processed={atoms_processed}, atoms_skipped={atoms_skipped}")
                                else:
                                    print(f"  [Metal Site] ERROR: No matching atoms found! metal_name={metal_name}, metal_type={metal_type}, chain={chain}, residue_number={residue_number}")
                            except Exception as e:
                                print(f"Warning: Failed to extract metal ligand structure from {first_pdb_path}: {e}")
                                import traceback
                                traceback.print_exc()
                
                # 返回响应
                print(f"  [Metal Site] Returning response, pdb_string length={len(pdb_string) if pdb_string else 0}")
                return jsonify({
                    "status": "success",
                    "motif": {
                        'motif_id': motif_id,
                        'metal_type': metal_type,
                        'metal_name': metal_name,
                        'ligand_name': metal_name,
                        'category': 'metal_sites',
                        'anchor_atoms': [],
                        'anchor_atoms_count': 0,
                        'geometry_constraints': [],
                        'chemistry_tag': f'Metal site: {metal_type}',
                        'extraction_method': 'metal_deduplication',
                        'is_metal_site': True,
                        'uniprot_id': metal_entry.get('uniprot_id', '') if metal_entry else '',
                        'ec_number': metal_entry.get('ec_number', '') if metal_entry else '',
                        'nanozyme_type': nanozyme_type if nanozyme_type else ''
                    },
                    "pdb_string": pdb_string,
                    "source": "synthetic"
                })
            
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
        
        # 生成PDB格式字符串用于3D可视化
        pdb_string = ""
        if formatted_motif.get('residue_structures'):
            try:
                pdb_string = generate_pdb_from_residues(
                    formatted_motif['residue_structures'],
                    formatted_motif['anchor_atoms']
                )
            except Exception as e:
                print(f"Warning: Failed to generate PDB string: {e}")
        
        return jsonify({
            "status": "success",
            "motif": formatted_motif,
            "pdb_string": pdb_string,  # 添加PDB字符串
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

