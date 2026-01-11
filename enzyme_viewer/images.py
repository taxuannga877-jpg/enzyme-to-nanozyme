"""
酶结构可视化模块
从ChemEnzyRetroPlanner项目中提取
"""
import py3Dmol

LABEL2ACTIVE_TYPE = {
    0: None,
    1: 'Binding Site',
    2: 'Catalytic Site',    # Active Site in UniProt
    3: 'Other Site'
}

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
    """
    with open(enzyme_structure_path) as ifile:
        system = ''.join([x for x in ifile])
    
    view = py3Dmol.view(width=view_size[0], height=view_size[1])
    view.addModelsAsFrames(system)
    
    active_data = []
    
    if show_active and (site_labels is not None) and not debug:
        i = 0
        for line in system.split("\n"):
            split = line.split()
            if len(split) == 0 or split[0] != "ATOM":
                continue
            # 使用PDB文件中的绝对残基编号（与site_labels的key一致）
            pdb_res_number = int(line[22:26].strip())
            # Gracefully handle residues that are not annotated; default to non-active color
            label = site_labels.get(pdb_res_number, 0)
            color = res_colors[label]
            view.setStyle({'model': -1, 'serial': i+1}, {"cartoon": {'color': color}})
            atom_name = line[12:16].strip()
            if (atom_name == 'CA') and (label != 0):
                residue_name = line[17:20].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
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
    return view.write_html(), active_data

