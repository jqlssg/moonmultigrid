# MoonMultigrid

面向 **二维矩形节点网格上扩散–反应方程** 的 MoonBit 几何多重网格库：

`-div(k grad u) + c u = f`，矩形四边均采用给定值的 Dirichlet 边界条件。

模块：**`jqlssg/moonmultigrid@0.1.0`**，采用 MIT 许可证。这是几何多重网格求解器，不是代数多重网格（AMG）、通用稀疏求解器或 CFD 框架。项目为原创实现，不是 PyAMG 的移植版。

## 安装与使用

```sh
moon add jqlssg/moonmultigrid@0.1.0
```

在包的导入配置中添加 `"jqlssg/moonmultigrid" @mg`。下面是一个只使用公开 API 的完整示例：

```moonbit
fn main raise {
  let grid = @mg.Grid::new(32, 32)
  let op = @mg.Operator::constant(grid)
  let source = @mg.Field::sample(grid, fn(_x, _y) { 1.0 })
  let boundary = @mg.Field::zeros(grid)
  let problem = @mg.Problem::new(op, source, boundary)
  let hierarchy = @mg.Hierarchy::new(op)
  let report = hierarchy.pcg(problem)
  assert_true(report.converged())
  println(report.relative_residual())
}
```

使用 `solution()` 前，应先检查 `converged()` 或 `reason()`，不能把所有返回结果都视为成功求解。达到迭代次数或估算工作量上限时，报告会明确标记未收敛；非法输入和数值计算失效会抛出 `SolverError`。

## 三类可复用工作流

| 使用者 | 输入与处理过程 | 输出与可运行示例 |
| --- | --- | --- |
| 传热与扩散仿真开发者 | 输入正导热系数场、源项及任意边界温度，先用 FMG 构造初值，再运行 MG-PCG | 温度场、真实残差和带方向的面通量；`moon run cmd/heat --target js` |
| 图像处理开发者 | 输入源图和目标图的矩形单通道补丁，采用源梯度或混合前向梯度，并固定目标补丁边界 | 未进行数值裁剪的 Poisson 融合补丁，可贴回新的图像场；`moon run cmd/edit --target js` |
| 数值计算生态开发者 | 输入解析或离散制造解方程，对比 CG 与 MG-PCG，并导出内部 COO 矩阵 | 收敛历史、离散化误差及稀疏矩阵／向量接口；`moon run cmd/compare --target js` |

这些是基于库的集成工作流，不包含图像文件编解码器或交互式可视化编辑器。各通道需要分别处理；粘贴补丁时，网格的物理像素间距必须一致。

## 数据模型与索引约定

* `Grid::new(nx, ny, lx?, ly?)` 中的 `nx`、`ny` 是 **区间数**，因此采样点数为 `(nx+1)*(ny+1)`。节点坐标为 `(i*lx/nx,j*ly/ny)`，展平索引为 `j*(nx+1)+i`，即 x 方向优先。每轴区间数最少为 2，最多为 1024。这里的 `?` 仅用于表示参数可选。
* `Field::new` 校验全部采样值并复制输入数组。`samples()`、通量数组、矩阵条目和迭代历史均返回快照。可使用 `with_value`、`map`、`add` 创建新场；`interpolate` 在闭合矩形范围内执行双线性插值。
* 导热系数定义在 **节点上**，必须为有限正数；共享面的系数取相邻节点的谐均值。反应系数是有限常数 `c >= 0`。消去 Dirichlet 边界后，内部系统为五点对称正定系统。
* `Operator::apply` 计算内部模板时会使用输入场的边界值，输出场的边界则置零。`Problem` 仅使用源场的 **内部值** 和边界场的 **边缘值**，并对非零边界进行提升处理。返回的解会恢复物理边界值。
* 与其他线性求解器对接时，将 `Problem::effective_rhs()` 的 `Field::interior_values()` 与 `Operator::interior_matrix()` 配对使用，再通过 `Problem::restore` 将齐次边界修正解恢复为物理解。
* 网格相等性包含物理长度。不能将已有层级复用于导热系数或反应系数不同的算子，否则会被拒绝。

## 多层算法

`Hierarchy::new(op, coarse_limit?=9)` 仅在 **两个轴的区间数均为偶数** 时同时减半，并通过节点注入传递导热系数。粗网格算子采用重新离散化，而非 Galerkin 乘积。任一轴为奇数或尺寸不足时停止粗化；若终端系统仍超过 **81 个未知量**，则拒绝构造，不会悄然退化为大型稠密求解。例如，32×16 网格可以使用，128×3 网格则无法粗化至允许的粗系统规模。

终端算子只进行一次稠密 Cholesky 分解，后续重复使用其因子。在齐次边界的内部空间上，全加权限制与双线性延拓满足 `R = P^T / 4`。V/W 循环支持加权 Jacobi（默认权重为 2/3），或按字典序正向与反向扫描的对称 Gauss–Seidel。`CycleConfig::new` 允许设置 1～8 次前／后平滑扫描；SGS 的一次扫描指一对正向／反向扫描。

* `h.solve(problem, options?, config?, initial?)`：反复执行 V/W 循环，以残差控制终止。
* `h.pcg(...)`：以固定 V 循环作为预条件器的共轭梯度法；`preconditioned=false` 用于普通 CG 对比。不支持用 W 循环进行 PCG。算法检查内积正性，并在每一步重新计算真实残差。
* `h.precondition(rhs, config?)`：从零初值开始执行齐次边界循环，可复用为预条件器。
* `h.full_multigrid(problem, cycles_per_level?=1, config?)`：仅提供嵌套 FMG 初值，**不保证残差已经收敛**；应将结果交给受控求解流程继续处理。
* `solve(problem, ...)`：创建新层级的便捷函数。重复求解不同右端项时，宜保留并复用已有层级。

`SolveOptions` 默认值为：`rtol=1e-8`、`atol=1e-12`、`max_iterations=100`、`work_budget=100000000`。收敛判据在内部节点上计算：`||b_eff-Au_correction||2 <= max(atol,rtol*||b_eff||2)`。历史记录包含初始残差。当 `||b_eff||=0` 时，`relative_residual()` 返回绝对残差。边界提升可能使 `b_eff` 较大，因此绝对容差应结合物理量尺度选择。独立计算恢复后物理方程的残差时，浮点消去可能造成小幅差异。

工作量单位是对行运算量的确定性估算，**不是** 实际耗时、指令数或完整的内存分配配额。迭代预算不包含层级构造、因子分解和初始残差计算；粗层工作量按二次规模估算，并计入 W 循环的重复访问。FMG 另设每层 1～8 次循环的上限。`StopReason` 区分 `Converged`、`IterationLimit` 和 `WorkLimit`。

## 诊断与数据交换

`op.flux(u)` 返回 `-k_face * forward_difference(u)`。x 方向的面数组包含 `nx*(ny+1)` 个元素，索引为 `j*nx+i`；y 方向包含 `(nx+1)*ny` 个元素，索引为 `j*(nx+1)+i`。内部散度加上 `c*u` 与 `op.apply(u)` 一致。这些数组包含边缘行／列上沿切向的面，**不是 Neumann 边界条件接口**。`op.contrast()` 返回节点系数最大值与最小值之比，极端比值可能溢出。

`InteriorMatrix` 以从 0 开始的内部行优先编号导出 COO 条目，对角为正、非对角为负。导出采用保守的 `5*unknowns` 条目预算。这是集成接口，不替代通用稀疏矩阵库。

`Field::to_json()` 返回 `{nx,ny,lx,ly,values}`，可调用 `.stringify()` 转为文本。`field_from_json` 要求键集合恰好匹配上述字段，区间数是范围内的整数，数值有限且采样数量准确。文本长度上限为 3200 万个 UTF-16 单元，嵌套深度上限为 4。重复键沿用 core JSON 解析器的行为（最后一个值生效），因此不要将快照作为可签名的规范化文档。`to_csv()` 导出 `i,j,x,y,value`，包含全部边界节点。JSON 快照只表示场数据，不接受任意求解器执行请求。

## 复现验证

发布前的最终本地验证使用独立安装的官方稳定工具链：moon 0.1.20261012（748b63c），moonc v0.10.15+943464bfa（2026-10-12）。CI 安装官方稳定版、记录准确版本，并在干净的运行环境中测试；上游稳定版会更新，未来可能需要兼容性调整。

```sh
moon info
moon fmt --check
moon check --deny-warn
moon test --target wasm-gc --deny-warn
moon test --target wasm --deny-warn
moon test --target js --deny-warn
moon test --target native --deny-warn
moon build --target js --deny-warn
python scripts/differential.py --target js
python scripts/audit.py
```

CI 对 wasm-gc、wasm、JS、native 四个目标分别执行检查、构建、测试及全部三个下游示例；JS 还运行独立差分验证。native 需要 C 编译器，Ubuntu CI 已提供；原 Windows 环境缺少 C 编译器，因此 **未宣称本地 native 验证通过**。生成接口时会与已跟踪的 `.mbti` 文件进行差异检查。

已发布 0.1.0 的 [四目标 CI](https://github.com/jqlssg/moonmultigrid/actions/runs/39102514242) 全部通过，每个目标均有 **25 项 MoonBit 测试**。Python 独立构建了 12 个变系数矩阵（每轴分别取 4／6／8／12 与 4／6／8 个区间的矩形组合），计入非零边界，并以部分选主元高斯消元求解；JS 最大节点绝对误差为 **1.78e-12**，验收阈值为 `<1e-9`。参考求解器不使用本库导出的模板或 COO 数据。解析四次多项式的网格加密测试验证了约二阶误差收敛；另有算子／预条件器的对称性、线性性、正性及传输算子伴随关系测试。

在 **32×32 个区间的解析四次多项式示例** 上，上述 Windows 稳定工具链的 JS 目标在 `rtol=1e-10` 时，MG-PCG 与普通 CG 分别使用 **9 次和 100 次迭代**，离散化 RMS 误差为 **0.0001662418888414**。这是单个可复现实例的迭代次数观察，**不是实际耗时加速比**，也不是与 PyAMG 的普遍性能对比。

## 能力边界、来源与许可证

不支持 AMG、非结构网格、三维、Neumann／周期边界、GPU／分布式求解，也不提供非线性或瞬态 PDE 引擎。在高对比系数、各向异性或网格未能解析的系数跳变条件下，点平滑器与系数注入方案可能收敛缓慢。本库不承诺对全部允许输入均具有网格无关收敛性。即使输入有限，中间算术仍可能溢出；检测到数值失效时会拒绝继续计算。容差应处于当前问题尺度的有效数值精度范围内。

更多说明见 [架构](docs/architecture.md)、[安全约束](SECURITY.md)、[参考来源与 AI 使用披露](NOTICE.md)、[更新记录](CHANGELOG.md)、[限定范围的查重证据](docs/competition/duplicate-check.md) 和 [验收约定](docs/competition/acceptance.md)。项目不依赖第三方 MoonCakes 模块，仅使用 MoonBit core 库。分发内容不包含申报材料或凭据。
