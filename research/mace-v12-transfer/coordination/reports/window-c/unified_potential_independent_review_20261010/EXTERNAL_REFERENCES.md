# 参考来源和可访问性

此任务的引用用于定义候选参考，不表示新势已通过验证。

1. **Ag 熔化温度：** NIST Chemistry WebBook, Silver, CAS 7440-22-4, phase-change data, `https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440224&Mask=4`. 常用熔化温度值约 1234.93 K；本执行环境访问 NIST WebBook 被出口 tunnel 403 拒绝，因此值作为公认表值引用、链接留下供复核，未能在本会话下载原页面。
2. **液态密度数据来源线索：** T. Iida and R. I. L. Guthrie, *The Physical Properties of Liquid Metals*, Clarendon Press, Oxford (1988), chapter/data tables for liquid-metal density. 二手汇编中银在熔点附近的常见量级约 9.3 g·cm⁻³。准确的 Ag ρ(T) 方程/表号和 1250 K 插值在当前环境未核实；检索 NIST/DOI/OpenAlex 均遇 tunnel 403。**因此 9.3 不能视作已审定的 ρref(1250 K)**，获得可复查表值及不确定度前不生成最终体积/DFT 输入。
3. **候选经典 Ag EAM 源：** S. M. Foiles, M. I. Baskes, and M. S. Daw, “Embedded-atom-method functions for the fcc metals Cu, Ag, Au, Ni, Pd, Pt, and their alloys,” *Phys. Rev. B* 33, 7983 (1986), DOI `10.1103/PhysRevB.33.7983`. 只作为候选构型采样器来源；具体参数文件、下载位置、许可和液态适用性需后续逐项核实。候选 EAM 不作为 DFT 标签。
4. **项目证据：** V18 audit commit `e455b71f6ef91a0955495001765e659aaee12174`; V19 design commit `a373d0ed7300c240006c360892808df0096f8aa1`; 三模型报告/修订 commit `2ee376322f544429dbe7923413c4dfd504e11b8b`; 上一轮 DFT plan commit `3ef38c8af17ed975ea2ceb96baeadce804320eb1`; completed state commit `81e7a218312f187957538cda322bf49b7e9dceae`; state snapshot used here `02aa85df73dd070b50e688baa42937dee66ba367`.
