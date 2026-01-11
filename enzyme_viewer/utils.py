"""
工具函数模块
"""
import pandas as pd

def color_cell(str_col, color_col):
    """为表格单元格添加颜色样式"""
    return f'<div style="color:{color_col}">{str_col}</div>'

def strong_cell(str_col):
    """为单元格内容添加粗体样式"""
    return f"<strong>{str_col}</strong>"

def convert_easifa_results(site_labels, pdb_fpath, view_size=(790, 600), add_style=True):
    """
    将EasIFA预测结果转换为可视化HTML和表格
    
    Args:
        site_labels: 残基标签列表
        pdb_fpath: PDB文件路径
        view_size: 视图大小
        add_style: 是否添加样式
    
    Returns:
        tuple: (structure_html, active_data_df)
    """
    from images import get_structure_html_and_active_data
    
    structure_html, active_data = get_structure_html_and_active_data(
        pdb_fpath, site_labels=site_labels, view_size=view_size
    )

    active_data_df = pd.DataFrame(
        active_data, columns=["Residue Index", "Residue Name", "Color", "Active Type"]
    )
    if not active_data_df.empty:
        if add_style:
            active_data_df["Active Type"] = active_data_df.apply(
                lambda row: color_cell(row["Active Type"], row["Color"]), axis=1
            )
        active_data_df = active_data_df[
            ["Residue Index", "Residue Name", "Active Type"]
        ]
        if add_style:
            for col in active_data_df.columns.tolist():
                active_data_df[col] = active_data_df[col].apply(lambda x: strong_cell(x))
    return structure_html, active_data_df


