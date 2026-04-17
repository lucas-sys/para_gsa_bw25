# results/fast/

FAST / RBD-FAST 全局敏感性分析输出目录。

| 文件 | 说明 |
|------|------|
| `fast_raw_results.xlsx` | 并行评估的原始 LCA 分数（每行对应一个样本） |
| `fast_indices.xlsx` | 各功能单元 × 影响类别的 FAST 灵敏度指数（S1、ST 等） |
| `fast_sensitivity_indices.png` | 一阶（S1）和全阶（ST）灵敏度指数条形图 |
| `fast_vs_sobol_comparison.png` | FAST S1 与 Sobol S1 散点对比图（可选） |

> 由 `gsa_workflow.ipynb` Section 10 自动生成。
