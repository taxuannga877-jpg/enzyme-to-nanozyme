"""
EasIFA Active Site Predictor - 必须集成模块
============================================

基于 ChemEnzyRetroPlanner 的 EasIFA 模型进行活性位点预测。
此模块为必须集成，用于处理未标注的酶数据。

两手抓策略：
- 标注数据：直接使用 UniProt/M-CSA 的活性位点注释
- 未标注数据：使用 EasIFA 模型预测活性位点
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field

# 项目根目录
PROJECT_ROOT = Path(__file__).parent.parent.parent

# 本地模型检查点路径
DEFAULT_ENZYME_MODEL_PATH = str(
    PROJECT_ROOT / "data" / "models" / "models" / "easifa" / "checkpoints" /
    "enzyme_site_type_predition_model" /
    "train_in_uniprot_ecreact_cluster_split_merge_dataset_all_limit_3_at_2023-12-19-16-06-42" /
    "global_step_284000"
)

DEFAULT_REACTION_MODEL_PATH = str(
    PROJECT_ROOT / "data" / "models" / "models" / "easifa" / "checkpoints" /
    "reaction_attn_net" /
    "model-ReactionMGMTurnNet_train_in_uspto_at_2023-04-05-23-46-25"
)

# ChemEnzyRetroPlanner 路径配置 (用于导入 EasIFA 模块)
CHEMENZY_PATH = os.environ.get(
    "CHEMENZY_PATH",
    str(PROJECT_ROOT / "参考项目" / "ChemEnzyRetroPlanner-main")
)

# 添加 EasIFA 到 Python 路径
EASIFA_PACKAGE_PATH = os.path.join(CHEMENZY_PATH, "retro_planner", "packages", "easifa")
if os.path.exists(EASIFA_PACKAGE_PATH):
    sys.path.insert(0, EASIFA_PACKAGE_PATH)


# 活性位点类型映射 (来自原始 ChemEnzyRetroPlanner)
LABEL_TO_SITE_TYPE = {
    0: None,           # 非活性位点
    1: "Binding",      # 结合位点
    2: "Catalytic",    # 催化位点 (活性中心)
    3: "Other",        # 其他位点
}


@dataclass
class PredictedActiveSite:
    """预测的活性位点"""
    residue_index: int          # 残基索引 (1-based)
    residue_name: str           # 残基名称 (如 HIS, SER)
    site_type: str              # 位点类型: Binding/Catalytic/Other
    confidence: float = 0.0     # 置信度
    coordinates: List[float] = field(default_factory=list)  # CA原子坐标


@dataclass
class ActiveSiteResult:
    """活性位点结果（标注或预测）"""
    uniprot_id: str
    pdb_path: str
    source: str                 # "annotated" 或 "predicted"
    sites: List[PredictedActiveSite] = field(default_factory=list)
    raw_labels: Optional[List[int]] = None  # EasIFA 原始预测标签


class EasIFAPredictor:
    """
    EasIFA 活性位点预测器 - 必须集成

    基于 ChemEnzyRetroPlanner 的 EasIFAInferenceAPI。
    用于预测未标注酶的活性位点。
    """

    def __init__(
        self,
        device: str = "cpu",
        enzyme_model_path: Optional[str] = None,
        reaction_model_path: Optional[str] = None,
        max_sequence_length: int = 600
    ):
        """
        初始化 EasIFA 预测器。

        Args:
            device: 运行设备 ("cpu" 或 "cuda")
            enzyme_model_path: 酶位点预测模型路径，None则使用默认
            reaction_model_path: 反应注意力模型路径，None则使用默认
            max_sequence_length: 最大序列长度
        """
        self.device = device
        self.enzyme_model_path = enzyme_model_path or DEFAULT_ENZYME_MODEL_PATH
        self.reaction_model_path = reaction_model_path or DEFAULT_REACTION_MODEL_PATH
        self.max_sequence_length = max_sequence_length
        self._model = None
        self._initialized = False

        # 尝试初始化 EasIFA
        self._init_easifa()

    def _init_easifa(self):
        """初始化 EasIFA 模型"""
        # 检查模型文件是否存在
        enzyme_model_file = os.path.join(self.enzyme_model_path, "model.pth")
        reaction_model_file = os.path.join(self.reaction_model_path, "model.pth")

        if not os.path.exists(enzyme_model_file):
            raise FileNotFoundError(
                f"[EasIFA] 酶位点模型不存在: {enzyme_model_file}\n"
                f"请确保模型文件已放置在正确位置"
            )
        if not os.path.exists(reaction_model_file):
            raise FileNotFoundError(
                f"[EasIFA] 反应注意力模型不存在: {reaction_model_file}\n"
                f"请确保模型文件已放置在正确位置"
            )

        print(f"[EasIFA] 酶位点模型: {self.enzyme_model_path}")
        print(f"[EasIFA] 反应模型: {self.reaction_model_path}")

        try:
            from easifa.interface.utils import EasIFAInferenceAPI

            # 使用本地模型检查点
            # 注意：反应模型路径在 EasIFA 内部硬编码，不需要传入
            self._model = EasIFAInferenceAPI(
                device=self.device,
                model_checkpoint_path=self.enzyme_model_path,
                max_enzyme_aa_length=self.max_sequence_length
            )

            self._initialized = True
            print("[EasIFA] 模型初始化成功")

        except ImportError as e:
            raise ImportError(
                f"[EasIFA] 无法导入 EasIFA 模块: {e}\n"
                f"请确保 ChemEnzyRetroPlanner 已正确安装，路径: {CHEMENZY_PATH}"
            )
        except Exception as e:
            raise RuntimeError(f"[EasIFA] 模型初始化失败: {e}")

    def predict(
        self,
        pdb_path: str,
        reaction_smiles: str = "C>>C"
    ) -> Optional[List[int]]:
        """
        预测活性位点。

        Args:
            pdb_path: PDB 文件路径
            reaction_smiles: 反应 SMILES

        Returns:
            预测标签列表，每个残基一个标签
        """
        if not self._initialized:
            raise RuntimeError("[EasIFA] 模型未初始化")

        if not os.path.exists(pdb_path):
            print(f"[EasIFA] PDB 文件不存在: {pdb_path}")
            return None

        try:
            predictions = self._model.inference(
                rxn=reaction_smiles,
                enzyme_structure_path=pdb_path
            )
            return predictions
        except Exception as e:
            print(f"[EasIFA] 预测失败: {e}")
            return None

    def predict_with_details(
        self,
        pdb_path: str,
        uniprot_id: str,
        reaction_smiles: str = "C>>C"
    ) -> Optional[ActiveSiteResult]:
        """
        预测活性位点并返回详细结果。

        Args:
            pdb_path: PDB 文件路径
            uniprot_id: UniProt ID
            reaction_smiles: 反应 SMILES

        Returns:
            ActiveSiteResult 对象
        """
        labels = self.predict(pdb_path, reaction_smiles)
        if labels is None:
            return None

        # 解析 PDB 获取残基信息
        sites = self._parse_sites_from_pdb(pdb_path, labels)

        return ActiveSiteResult(
            uniprot_id=uniprot_id,
            pdb_path=pdb_path,
            source="predicted",
            sites=sites,
            raw_labels=labels
        )

    def _parse_sites_from_pdb(
        self,
        pdb_path: str,
        labels: List[int]
    ) -> List[PredictedActiveSite]:
        """从 PDB 文件解析活性位点信息"""
        sites = []
        residue_info = {}  # {res_idx: (res_name, coords)}

        with open(pdb_path, 'r') as f:
            first_res_idx = None
            for line in f:
                if not line.startswith("ATOM"):
                    continue

                res_idx = int(line[22:26].strip())
                if first_res_idx is None:
                    first_res_idx = res_idx

                rel_idx = res_idx - first_res_idx
                atom_name = line[12:16].strip()

                if atom_name == "CA":
                    res_name = line[17:20].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    residue_info[rel_idx] = (res_name, [x, y, z])

        # 构建活性位点列表
        for idx, label in enumerate(labels):
            if label > 0:  # 非零 = 活性位点
                site_type = LABEL_TO_SITE_TYPE.get(label, "Unknown")
                res_name, coords = residue_info.get(idx, ("UNK", []))

                sites.append(PredictedActiveSite(
                    residue_index=idx + 1,
                    residue_name=res_name,
                    site_type=site_type,
                    coordinates=coords
                ))

        return sites
