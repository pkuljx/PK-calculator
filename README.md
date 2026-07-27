# PK Calculator

基于经典药代动力学公式，从临床前动物数据预测人体 PK 参数（清除率 CL、稳态分布容积 Vss、半衰期 t½）的工具集。

本仓库同时提供两种使用形态：
1. **Claude Skill**（`SKILL.md` + `references/formulas.md`）— 供 Claude / Agent 直接调用；
2. **Python 参考实现**（`pk_calculator.py`）— 可作为模块 import，也可作为命令行工具运行。

两条路径共享同一套公式，结果完全一致。

---

## 一、支持的 8 种方法

| # | 方法 | 预测参数 | 最低输入要求 |
|---|------|---------|--------------|
| 1 | **IVIVE** | CL | 体外 CLint（微粒体或肝细胞）、fu,p、cLogP 或 fu,inc |
| 2 | **Allometry CL** | CL | ≥2 个物种的 BW + CL |
| 3 | **Allometry Vss** | Vss | ≥2 个物种的 BW + Vss |
| 4 | **FCIM** | CL | ≥2 个物种的 BW + CL、参比种和人的 fu |
| 5 | **Wajima CL** | CL | 大鼠 CL、犬 CL、MW、氢键受体数 Ha |
| 6 | **Wajima Vss** | Vss | 大鼠 Vss、犬 Vss |
| 7 | **Obach Vss** | Vss | 犬 Vss、犬和人的 fu |
| 8 | **Oie–Tozer Vss** | Vss | 大鼠 / 猴 / 犬 的 Vss + fu，人的 fu |

外加一个便捷函数 `half_life(CL, Vss)`：`t½ = 0.693 × Vss / CL`。

---

## 二、实现思路

### 2.1 为什么同时做成 Skill 和 Python 脚本

- **Skill 模式**：把公式、常数、方法选择规则写成 Markdown，让 Agent 依据自然语言输入自动挑方法、准备数据、给出解释。适合"帮我算一下"的对话式场景。
- **脚本模式**：把同一套公式落到确定性的 Python 代码里。凡是需要"可复现、可测试、可嵌入流水线"的场合，都走这里，避免每次由 LLM 现算带来的算术漂移。

两者互相引用：`SKILL.md` 指向 `pk_calculator.py` 作为首选执行路径；`pk_calculator.py` 的注释指向 `references/formulas.md` 作为公式出处。

### 2.2 目录结构

```
PK-calculator/
├── SKILL.md                  # Skill 入口：workflow、方法选择指南、单位约定
├── references/
│   └── formulas.md           # 全部公式、常数表、回归步骤的权威来源
├── pk_calculator.py          # Python 参考实现（本 README 的重点）
└── README.md
```

### 2.3 `pk_calculator.py` 内部结构

单文件、零第三方依赖（仅用标准库 `math` / `argparse` / `json`），共分四层：

1. **常数层** — 把 `formulas.md` 中的所有物理化学常数直接编码为 dict：
   - `SPECIES_BW`：各物种默认体重
   - `IVIVE_DEFAULTS`：每个物种的 MPPGL / HPGL / Qh / 肝重
   - `BRAIN_SCALE`：脑重换算比例
   - `OIE_TOZER_VOL`：Oie–Tozer 用到的生理容积（Vp / Ve / Vr / R_e/i）

2. **通用辅助层**
   - `brw(species, bw)` — 脑重
   - `mlp_hours(species, bw)` — 最大寿命 MLP（小时）
   - `loglog_regression(bws, ys)` — 对数-对数空间的最小二乘线性回归，返回 `a, b, log10_a, R²`。所有异速生长类方法（Allometry CL / Vss、FCIM）都调用这一个函数，保证回归口径一致。
   - `fu_inc_from_clogp(clogp)` — 由 cLogP 估算 fu,inc

3. **方法层** — 每种预测方法是一个纯函数，输入关键字参数，输出结构化字典：
   - `ivive(...)`：well-stirred 模型；同时支持微粒体和肝细胞路径；如果给了 `CL_obs` 就顺带算 IVIVC。
   - `allometry_cl(...)`：先做初拟合拿到斜率 b，按规则自动选模型：
     - `b ≤ 0.71` → simple
     - `0.71 < b ≤ 1.0` → MLP 校正
     - `b > 1.0` → BrW 校正

     再用选中的 Y 变量二次拟合，最后把人体 Y 值反算回 CL。
   - `allometry_vss(...)`：直接 simple allometry，可选 fu 校正。
   - `fcim(...)`：**不做 fu 校正**（这是 FCIM 的定义），只用截距 a 和 `Rfu = fu_ref/fu_human`，套 `33.35 × (a/Rfu)^0.77`。
   - `wajima_cl / wajima_vss`：直接把多项回归公式写死；Vss 版注意先把 L/kg → mL/kg 再取 log。
   - `obach_vss`：`Vss_dog × fu_human/fu_dog × BW`。
   - `oie_tozer_vss`：先对大鼠/猴/犬各自算 fut，取平均后代入人的生理容积公式反推 Vdss。

4. **调度 + CLI 层**
   - `METHODS` 字典把方法名映射到函数；
   - `run(method, payload)` 是编程入口；
   - `_cli()` 用 `argparse` 提供 `list-methods` 与 `run <method> -i input.json` 两个子命令，输入输出都是 JSON，方便被别的脚本 pipe。

### 2.4 一个隐含约定：单位

代码严格遵循 `SKILL.md` 里声明的单位约定，避免"到底是 mL/min/kg 还是 L/h"之类的错乱：

- CL 输入统一 **mL/min/kg**；内部换算到 L/h 时乘以 `0.06 × BW`
- Vss 输入统一 **L/kg**
- fu 无量纲，0–1
- 体重 kg；默认值来自 `SPECIES_BW`，也可显式覆盖

---

## 三、使用示例

### 3.1 作为 Python 模块

```python
from pk_calculator import allometry_cl, wajima_vss, half_life

# 4 个物种数据推人 CL，自动选模型
r = allometry_cl([
    {"species": "mouse",  "CL": 50},
    {"species": "rat",    "CL": 30},
    {"species": "dog",    "CL": 10},
    {"species": "monkey", "CL": 15},
])
print(r["selected_model"], r["human_CL_L_per_h"])

# Wajima 法预测 Vss
v = wajima_vss(Vss_rat=2.0, Vss_dog=3.0)
print(v["human_Vss_L"])

# 半衰期
print(half_life(CL_L_per_h=15.5, Vss_L=v["human_Vss_L"]))
```

### 3.2 作为命令行工具

```bash
# 列出全部方法
python pk_calculator.py list-methods

# 通过 stdin 传入 JSON
echo '{"Vss_dog":3.0,"fu_dog":0.2,"fu_human":0.05}' \
    | python pk_calculator.py run obach-vss -i -

# 或者从文件读入，结果写到文件
python pk_calculator.py run allometry-cl -i input.json -o result.json
```

一个 `input.json` 示例（Allometry CL 带 fu 校正）：

```json
{
  "species_data": [
    {"species": "mouse", "CL": 50, "fu": 0.15},
    {"species": "rat",   "CL": 30, "fu": 0.12},
    {"species": "dog",   "CL": 10, "fu": 0.20},
    {"species": "monkey","CL": 15, "fu": 0.18}
  ],
  "fu_correction": true,
  "fu_human": 0.10,
  "human_weight": 70
}
```

### 3.3 作为 Claude Skill

在 Claude Code 里把本仓库放入 `~/.claude/skills/pk-calculator/`，然后自然语言描述数据即可：

> "我有 mouse/rat/dog/monkey 的 CL 分别是 50/30/10/15 mL/min/kg，帮我预测人体 CL。"

Agent 会读取 `SKILL.md`，选择合适方法，调用 `pk_calculator.py` 拿到确定性结果，再把过程与结论以 Markdown 表格返回。

---

## 四、方法选择建议

| 你手里的数据 | 推荐方法 |
|--------------|----------|
| 只有体外 CLint | IVIVE |
| ≥3 个物种的 CL | Allometry CL（自动选 simple/MLP/BrW）+ FCIM 对照 |
| 只有大鼠 + 犬 | Wajima CL、Wajima Vss |
| 只有犬 Vss + fu | Obach Vss |
| 大鼠/猴/犬 Vss + fu 齐全 | Oie–Tozer Vss（生理机制解释力最强） |

数据充足时，建议**同时跑多种方法并列对比**——`pk_calculator.py` 的函数式设计正是为了方便这种批量比较。

---

## 五、参考

所有公式、常数与推导步骤见 [`references/formulas.md`](references/formulas.md)。修改任何公式请以该文件为准，然后同步更新 `pk_calculator.py`。
