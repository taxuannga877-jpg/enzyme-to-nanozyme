#!/bin/bash
# 启动酶结构查看器

echo "启动酶结构查看器..."
echo "确保已安装所有依赖: pip install -r requirements.txt"
echo ""

# 检查Python环境
if ! command -v python3 &> /dev/null; then
    echo "错误: 未找到Python3"
    exit 1
fi

# 运行应用
python3 app.py


