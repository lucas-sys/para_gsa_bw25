# LCA 全局敏感性分析（GSA）工具包

基于 [Brightway2.5](https://brightway.dev/)、[SALib](https://salib.readthedocs.io/) 和 Python 多进程的模块化 **参数扰动生命周期评价（LCA）** 与 **全局敏感性分析（GSA）** 工具包。

## 功能特性

- **确定性 LCA** -- 计算多个功能单元和影响类别的基准环境影响分数
- **贡献分析** -- 递归 Tier-2 贡献分析，识别关键工艺贡献
- **公式验证** -- 验证存储的交换量与参数化公式的一致性
- **OAT 敏感性** -- 一次一因子有限差分敏感性分析，含弹性系数
- **解析方差** -- 基于不确定性元数据的方差分解（支持 13 种分布类型）
- **蒙特卡洛模拟** -- 可配置工作进程数的并行随机采样
- **Sobol GSA** -- 基于 SALib 的一阶（S1）和全局（ST）Sobol 指数分析，支持并行计算
- **FAST / RBD-FAST GSA** -- 基于傅里叶幅值灵敏度测试（FAST）及随机平衡设计变体（RBD-FAST）的全局敏感性分析，支持并行计算

## 项目结构

```
.
├── utils/                            # 核心 LCA 模块
│   ├── __init__.py                   # 公共 API 统一导出
│   ├── lca_config.py                 # 配置、日志、计时器、LCASetup
│   ├── lca_project_setup.py          # Brightway 项目初始化、MultiLCA
│   ├── lca_parameters.py             # 参数查询、公式求解
│   ├── lca_matrices.py               # 稀疏矩阵更新工具
│   ├── lca_contribution.py           # 贡献分析、OAT、解析方差
│   ├── lca_monte_carlo.py            # 蒙特卡洛采样与并行评估
│   ├── lca_sobol.py                  # Sobol GSA（构建、采样、评估、分析）
│   ├── lca_fast.py                   # FAST / RBD-FAST GSA（构建、采样、评估、分析）
│   └── lca_gsa_helper.py             # 向后兼容的重导出层
├── data/                             # 输入数据库（用户自行提供）
│   └── ...                           # 将 Excel 数据库文件放在此处
├── results/                          # 输出文件（按分析类型分类）
│   ├── deterministic/                # 确定性 LCA 结果
│   ├── contribution/                 # 贡献分析结果
│   ├── sensitivity_oat/              # OAT 敏感性结果
│   ├── sensitivity_analytical/       # 解析方差结果
│   ├── monte_carlo/                  # 蒙特卡洛结果
|   ├── sobol/                        # Sobol GSA 结果
│   └── fast/                         # Fast GSA 结果
├── gsa_workflow.ipynb                # 交互式 Notebook（完整工作流）
├── init_databases.py                 # 数据库初始化脚本
├── test_full_workflow.py             # 端到端测试
├── test_run.py   c                    # 快速验证测试
├── test_sobol_only.py                # Sobol 专项测试
├── requirements.txt
├── LICENSE
└── README.md
```

## 环境要求

- Python >= 3.10
- [Brightway2.5](https://brightway.dev/) (`bw2data`, `bw2calc`, `bw2io`, `bw2analyzer`)
- ecoinvent 数据库（需要 ecoinvent 许可证）

## 安装

```bash
# 创建并激活 conda 环境
conda create -n bw25 python=3.10
conda activate bw25

# 安装依赖
pip install -r requirements.txt
```

## 配置说明

运行任何分析前，需要配置以下参数：

| 参数                            | 说明                         | 配置位置                              |
| ----------------------------- | -------------------------- | --------------------------------- |
| `PROJECT_NAME`                | Brightway 项目名称             | `init_databases.py`、测试脚本、Notebook |
| `FOREGROUND_DB`               | 前景数据库名称                    | 测试脚本、Notebook                     |
| `GROUP_NAME`                  | ActivityParameter 参数组名称    | 测试脚本、Notebook                     |
| `FU_CODE`                     | 功能单元的活动代码                  | 测试脚本、Notebook                     |
| `METHODS`                     | 影响类别方法元组列表                 | 测试脚本、Notebook                     |
| `ECOINVENT_USERNAME/PASSWORD` | ecoinvent 凭据               | `init_databases.py` 或环境变量         |
| `CUSTOM_DATABASES`            | (数据库名, xlsx路径, 匹配数据库列表) 元组 | `init_databases.py`               |

模板文件中所有参数均使用 `your_*` 占位符，运行前请替换为实际值。

## 数据库初始化

```bash
# 设置 ecoinvent 凭据（推荐使用环境变量）
export ECOINVENT_USERNAME="your_username"
export ECOINVENT_PASSWORD="your_password"

# 或直接编辑 init_databases.py 中的配置部分

# 运行初始化脚本
python init_databases.py
```

将自定义数据库 Excel 文件放入 `data/` 目录，并在 `init_databases.py` 的 `CUSTOM_DATABASES` 列表中添加相应条目。

## 快速开始

### Python 脚本

```python
from utils import (
    set_project, initialize_multilca, deterministic_lca,
    oat_sensitivity, run_parallel_monte_carlo,
    build_sobol_problem, generate_sobol_samples,
    run_parallel_sobol_from_samples, sobol_indices_from_results,
)

import bw2data as bd

# 初始化
PROJECT_NAME = 'your_project_name'
FOREGROUND_DB = 'your_foreground_db'
GROUP_NAME = 'your_parameter_group'

set_project(PROJECT_NAME)
db = bd.Database(FOREGROUND_DB)
fu = [db.get(code='your_fu_code')]

methods = [
    ('your_ecoinvent_db', 'EF v3.1', 'climate change',
     'global warming potential (GWP100)'),
]

# 确定性 LCA
df = deterministic_lca(fu, methods)

# 蒙特卡洛（1000 次迭代，2 个工作进程）
mc_results = run_parallel_monte_carlo(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    n_samples=1000,
    n_workers=2,
)

# Sobol GSA
problem, params = build_sobol_problem(GROUP_NAME)
samples = generate_sobol_samples(problem, N=512)
sobol_results = run_parallel_sobol_from_samples(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    sobol_problem=problem,
    sobol_samples=samples,
    n_workers=2,
)

# FAST / RBD-FAST GSA
from utils import lca_fast

fast_problem, fast_params = lca_fast.build_fast_problem(GROUP_NAME)

# --- Classic FAST (S1 only, fewer samples) ---
fast_samples = lca_fast.generate_fast_samples(fast_problem, M=4, method='fast')
fast_results = lca_fast.run_parallel_fast_from_samples(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    fast_problem=fast_problem,
    fast_samples=fast_samples,
    n_workers=2,
)
Si_fast = lca_fast.fast_indices_from_results(
    fast_problem=fast_problem,
    fast_results_df=fast_results,
    functional_unit=fu[0]['name'],
    method=methods[0],
    fast_samples=fast_samples,
    M=4,
    analysis_method='fast',
)
df_fast = lca_fast.fast_indices_to_dataframe(Si_fast, fast_problem, analysis_method='fast')

# --- RBD-FAST (S1 + ST) ---
rbd_samples = lca_fast.generate_fast_samples(fast_problem, M=4, method='rbd_fast', seed=42)
rbd_results = lca_fast.run_parallel_fast_from_samples(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    fast_problem=fast_problem,
    fast_samples=rbd_samples,
    n_workers=2,
)
Si_rbd = lca_fast.fast_indices_from_results(
    fast_problem=fast_problem,
    fast_results_df=rbd_results,
    functional_unit=fu[0]['name'],
    method=methods[0],
    fast_samples=rbd_samples,
    M=4,
    analysis_method='rbd_fast',
)
df_rbd = lca_fast.fast_indices_to_dataframe(Si_rbd, fast_problem, analysis_method='rbd_fast')

# One-click workflow
result = lca_fast.run_full_fast_workflow(
    project_name=PROJECT_NAME,
    group_name=GROUP_NAME,
    functional_units=fu,
    methods=methods,
    M=4,
    n_workers=2,
    method='rbd_fast',  # or 'fast'
)
```

### Jupyter Notebook

打开 `gsa_workflow.ipynb` 并按顺序执行单元格。Notebook 覆盖了从数据库设置到 Sobol GSA 的完整工作流。运行前请更新配置单元格中的项目特定值。

## 工作流章节

| 章节  | 描述                  | 输出目录                              |
| --- | ------------------- | --------------------------------- |
| 0   | Brightway 项目设置      | --                                |
| 1   | 影响类别选择              | --                                |
| 2   | 功能单元选择              | --                                |
| 3   | 确定性 LCA             | `results/deterministic/`          |
| 4   | 贡献分析                | `results/contribution/`           |
| 5   | OAT 敏感性             | `results/sensitivity_oat/`        |
| 6   | 解析方差                | `results/sensitivity_analytical/` |
| 7   | 蒙特卡洛模拟              | `results/monte_carlo/`            |
| 8   | Sobol GSA           | `results/sobol/`                  |
| 10  | FAST / RBD-FAST GSA | `results/fast/`                   |

## 支持的不确定性分布

| 代码   | 分布类型      | 参数                                       |
| ---- | --------- | ---------------------------------------- |
| 0, 1 | 无不确定性     | --                                       |
| 2    | 对数正态      | loc (ln_mean), scale (ln_std)            |
| 3    | 正态        | loc (均值), scale (标准差)                    |
| 4    | 均匀        | minimum, maximum                         |
| 5    | 三角        | minimum, maximum, loc (众数)               |
| 6    | 伯努利       | loc (概率)                                 |
| 7    | 离散均匀      | minimum, maximum                         |
| 8    | 威布尔       | shape, loc, scale                        |
| 9    | 伽马        | shape (alpha), loc, scale                |
| 10   | 贝塔        | shape (alpha), shape2 (beta), loc, scale |
| 11   | 耿贝尔 (GEV) | shape (xi), loc, scale                   |
| 12   | t 分布      | shape (自由度), loc, scale                  |

## FAST / RBD-FAST GSA

| 变体       | 函数                  | 指数      | 样本量公式       | 说明                     |
| -------- | ------------------- | ------- | ----------- | ---------------------- |
| 经典 FAST  | `method='fast'`     | S1      | `(4M²+1)·k` | 确定性，仅一阶，样本量低           |
| RBD-FAST | `method='rbd_fast'` | S1 + ST | `≥65·k`     | 随机，一阶+全阶，SALib 拉丁超立方采样 |

其中 `k` = 不确定参数个数，`M` = 干扰因子（推荐默认值 4）。

**与 Sobol 比较**

| 方法       | 一阶 S1 | 全阶 ST | 最小推荐样本      | 适用场景        |
| -------- | ----- | ----- | ----------- | ----------- |
| Sobol    | ✓     | ✓     | `N·(k+2)`   | 标准全局 GSA    |
| 经典 FAST  | ✓     | ✗     | `(4M²+1)·k` | 仅需 S1，样本预算紧 |
| RBD-FAST | ✓     | ✓     | `65·k`      | 样本量介于二者之间   |

主要 API：

| 函数                                              | 说明                                |
| ----------------------------------------------- | --------------------------------- |
| `build_fast_problem(group)`                     | 从 ActivityParameter 构建 SALib 问题字典 |
| `sanitize_fast_problem(problem)`                | 校验并修复非法参数边界                       |
| `generate_fast_samples(problem, M, method)`     | 生成 FAST / RBD-FAST 样本矩阵           |
| `expected_fast_sample_rows(problem, M, method)` | 预测样本行数                            |
| `run_parallel_fast_from_samples(...)`           | 多进程评估样本矩阵                         |
| `fast_indices_from_results(...)`                | 计算 FAST / RBD-FAST 灵敏度指数          |
| `fast_indices_to_dataframe(Si, problem)`        | 转为整洁 DataFrame                    |
| `run_full_fast_workflow(...)`                   | 一键工作流（构建→采样→评估）                   |

## 多进程

蒙特卡洛和 Sobol 评估均使用 Python 的 `spawn` 多进程上下文以确保 Windows 兼容性。每个工作进程独立加载 Brightway 项目和数据库，因此内存占用与工作进程数成正比。建议使用 2-4 个工作进程。

## 许可证

本项目基于 MIT 许可证发布。详见 [LICENSE](LICENSE)。

## 致谢

- [Brightway2.5](https://brightway.dev/) -- LCA 框架
- [SALib](https://salib.readthedocs.io/) -- 敏感性分析库
- [ecoinvent](https://ecoinvent.org/) -- 生命周期清单数据库
