#!/bin/bash
# ============================================
# EasIFA 依赖安装脚本
# ============================================
# 用途：安装 EasIFA 活性位点预测所需的所有依赖
# 使用：bash scripts/install_dependencies.sh

set -e  # 遇到错误立即退出

echo "=========================================="
echo "开始安装 EasIFA 依赖"
echo "=========================================="

# 检查是否在正确的环境中
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "警告：未检测到 conda 环境"
    echo "建议先激活环境：source envs/bin/activate"
    read -p "是否继续？(y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

echo ""
echo "步骤 1/5: 安装 DGL (图神经网络库)"
echo "----------------------------------------"
pip install dgl -f https://data.dgl.ai/wheels/cu126/repo.html
echo "✓ DGL 安装完成"

echo ""
echo "步骤 2/5: 安装 DGLLife (化学分子图处理)"
echo "----------------------------------------"
pip install dgllife
echo "✓ DGLLife 安装完成"

echo ""
echo "步骤 3/5: 安装 RDKit (化学信息学工具)"
echo "----------------------------------------"
pip install rdkit
echo "✓ RDKit 安装完成"

echo ""
echo "步骤 4/5: 安装 fair-esm (ESM-2 蛋白质语言模型)"
echo "----------------------------------------"
pip install fair-esm
echo "✓ fair-esm 安装完成"

echo ""
echo "步骤 5/5: 安装其他依赖"
echo "----------------------------------------"
pip install py3Dmol matplotlib pandas flask flask-cors
echo "✓ 其他依赖安装完成"

echo ""
echo "=========================================="
echo "所有依赖安装完成！"
echo "=========================================="
echo ""
echo "下一步："
echo "1. ESM-2 模型权重会在首次使用时自动下载（约2.5GB）"
echo "2. 运行测试：python -c 'from nanozyme_mining.prediction import EasIFAPredictor'"
echo ""
