#!/usr/bin/env python3
"""构建累计外扩 2 Å/侧的第二档骨架加粗模型。"""
from build_step01 import main


if __name__ == "__main__":
    main(offset_each_side_A=2.0, step=2)
