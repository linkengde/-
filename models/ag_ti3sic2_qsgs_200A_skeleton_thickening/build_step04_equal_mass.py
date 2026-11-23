#!/usr/bin/env python3
"""构建 Ag/Ti3SiC2 质量比维持 1:1 的最大可行小幅加粗版本。"""
from build_step01 import main


if __name__ == "__main__":
    # 0.465625 Å/侧是本次搜索到仍有足够连通 Ag 位点守住 1:1 的最大档位。
    main(offset_each_side_A=0.465625, step=4, target_mass_ratio=1.0)
