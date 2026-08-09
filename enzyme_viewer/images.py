"""
酶结构可视化模块
从ChemEnzyRetroPlanner项目中提取

PR1-6 (H3 fix): explicit encoding + file-existence check + per-line try/except
around coordinate parsing. Non-ASCII chains, truncated PDB lines, and missing
files now produce structured warnings instead of crashing the request handler.
"""
import logging
import os
import re
import py3Dmol

log = logging.getLogger("e2n.images")

_PY3DMOL_REMOTE_LOADER = re.compile(
    r"var loadScriptAsync = function\(uri\)\{.*?\n\};\s*"
    r"if\(typeof \$3Dmolpromise === 'undefined'\) \{\s*"
    r"\$3Dmolpromise = null;\s*"
    r"\$3Dmolpromise = loadScriptAsync\('[^']+/3Dmol-min\.js'\);\s*"
    r"\}",
    re.DOTALL,
)

_LOCAL_3DMOL_BOOTSTRAP = """\
if (typeof $3Dmol === 'undefined') {
  throw new Error('Local 3Dmol.js asset is required before rendering a structure.');
}
if (typeof $3Dmolpromise === 'undefined') {
  $3Dmolpromise = Promise.resolve();
}
"""

LABEL2ACTIVE_TYPE = {
    0: None,
    1: 'Binding Site',
    2: 'Catalytic Site',    # Active Site in UniProt
    3: 'Other Site'
}


def _use_preloaded_3dmol(html):
    """Remove py3Dmol's remote loader and reuse the page's local 3Dmol asset."""

    rewritten, replacements = _PY3DMOL_REMOTE_LOADER.subn(
        _LOCAL_3DMOL_BOOTSTRAP,
        html,
        count=1,
    )
    if replacements != 1:
        raise RuntimeError("unexpected py3Dmol HTML: remote loader block not found")
    return rewritten


def get_structure_html_and_active_data(
    enzyme_structure_path,
    site_labels=None,
    view_size=(900, 900),
    res_colors={
        0: '#73B1FF',   # 非活性位点
        1: '#FF0000',     # Binding Site
        2: '#00B050',     # Active Site
        3: '#FFFF00',     # Other Site
    },
    show_active=True,
    debug=False,
):
    """
    生成酶结构的3D可视化HTML和活性位点数据

    Args:
        enzyme_structure_path: PDB文件路径
        site_labels: 残基标签列表，None表示不显示活性位点
        view_size: 视图大小 (width, height)
        res_colors: 残基颜色映射
        show_active: 是否显示活性位点
        debug: 调试模式

    Returns:
        tuple: (structure_html, active_data)
            structure_html: 3Dmol.js生成的HTML字符串
            active_data: 活性位点数据列表 [(res_idx, res_name, color, active_type), ...]

    Raises:
        FileNotFoundError: PDB 文件不存在
    """
    # PR1-6 (H3 fix): file existence + readable check up front.
    if not os.path.exists(enzyme_structure_path):
        raise FileNotFoundError(f"PDB file not found: {enzyme_structure_path}")

    # PR1-6 (H3 fix): explicit encoding='utf-8' + errors='replace' so non-ASCII
    # bytes in headers (occasional in older PDB entries) don't blow up read.
    with open(enzyme_structure_path, encoding='utf-8', errors='replace') as ifile:
        system = ifile.read()

    view = py3Dmol.view(width=view_size[0], height=view_size[1])
    view.addModelsAsFrames(system)

    active_data = []

    if show_active and (site_labels is not None) and not debug:
        i = 0
        for line_no, line in enumerate(system.split("\n"), start=1):
            split = line.split()
            if len(split) == 0 or split[0] != "ATOM":
                continue
            # PR1-6 (H3 fix): wrap coordinate parsing in try/except. A truncated
            # ATOM line (fewer than 54 chars) used to raise ValueError and abort
            # the whole request; now we log and skip that single residue.
            try:
                pdb_res_number = int(line[22:26].strip())
            except (ValueError, IndexError):
                log.debug("images.py: bad residue number at line %d: %r", line_no, line[:80])
                i += 1
                continue
            # Gracefully handle residues that are not annotated; default to non-active color
            label = site_labels.get(pdb_res_number, 0)
            color = res_colors[label]
            view.setStyle({'model': -1, 'serial': i+1}, {"cartoon": {'color': color}})
            atom_name = line[12:16].strip() if len(line) >= 16 else ""
            if (atom_name == 'CA') and (label != 0):
                try:
                    residue_name = line[17:20].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except (ValueError, IndexError):
                    log.debug("images.py: bad coords at line %d: %r", line_no, line[:80])
                    i += 1
                    continue
                view.addLabel(
                    f'{residue_name} {pdb_res_number}',
                    {
                        "fontSize": 12,
                        "position": {"x": x, "y": y, "z": z},
                        "fontColor": 'black',
                        "fontOpacity": 1.0,
                        "backgroundColor": color,
                        "bold": True,
                        "backgroundOpacity": 0.6
                    }
                )
                active_data.append((
                    pdb_res_number,
                    residue_name,
                    color,
                    LABEL2ACTIVE_TYPE[label]
                ))

            i += 1
    else:
        view.setStyle({'model': -1}, {"cartoon": {'color': res_colors[0]}})

    view.zoomTo()
    view.zoom(1.6, 600)
    return _use_preloaded_3dmol(view.write_html()), active_data
