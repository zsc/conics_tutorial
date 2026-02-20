---
name: conics-report-generator
description: 生成/维护圆锥曲线图文报告 dist/conics_report.html（MathJax + markdown2 + Python 生成器 + 纯 JS/SVG 互动图）。当你需要：从 AGENTS.md 生成 HTML、修复公式被 Markdown 吃掉/转义、把代码块与局部段落默认折叠、添加/调试带 slider 的互动图（展示距离和/差、焦点-准线离心率、旋转去交叉项、反射性质、导圆、圆锥截线统一图等）时使用。
---

# conics_report 迭代经验教训（Skill）

## 0) 关键事实（先记住）

- 报告是“生成物”：`dist/conics_report.html` 会被 `python3 scripts/generate_conics_report.py` 覆盖；不要手改 `dist/` 来做“永久修复”。
- 内容主来源：`AGENTS.md`；结构/互动图/折叠/修复逻辑在 `scripts/generate_conics_report.py`。
- 渲染链路：`AGENTS.md` → 预处理（修 LaTeX/数学定界符）→ `markdown2` → HTML 后处理（折叠/修复）→ MathJax 渲染 + 内联 JS mount 互动图。

## 1) 最小工作流（改一次就跑一次）

1. 修改 `AGENTS.md`（文字/公式）或 `scripts/generate_conics_report.py`（生成逻辑/互动图/CSS）。
2. 运行生成：`python3 scripts/generate_conics_report.py`
3. 打开/刷新：`dist/conics_report.html`（遇到 JS/CSS 变更，用硬刷新）。
4. 用 `rg` 定位问题与产物是否生效（例如查 `祖暅原理`、`widget-rotation-abc`、`details class="codeblock"`）。

## 2) markdown2 + MathJax 的高频坑（以及有效解法）

### 2.1 `_` 会被 Markdown 当强调，破坏 TeX

- 现象：`\lambda_1` 变成 HTML `<em>`，MathJax 解析失败。
- 解法：在“数学区域”里把 `_` 替换为 `&#95;`（浏览器会还原成 `_`，MathJax 可正常读）。
- 原则：只在数学定界符内替换；不要全局替换。

### 2.2 数学里的 `<` `>` 会触发 HTML 解析

- 现象：`0<e<1` 在 HTML 里被吞/变形。
- 解法：在数学区域转义为 `&lt;` `&gt;`。

### 2.3 `\(...\)` 在 Markdown 里容易被“反斜杠规则”吃掉

- 建议：在 Python 生成的 Markdown 字符串里用 `\\(` `\\)`（双反斜杠）保留到最终 HTML 给 MathJax。

### 2.4 `pmatrix`/矩阵换行 `\\` 很脆

- 现象：`\\` 在多层（Python 字符串 → Markdown → HTML）后丢失，导致矩阵不换行或出现奇怪转义。
- 解法优先级：
  1) 能不用矩阵就不用：改成“展开写法”（更稳）。
  2) 必须矩阵：对已知坏模式做“定点修复/替换”（小范围字符串修正），并生成后在 `dist/conics_report.html` 里 grep 验证。

### 2.5 `<details>` 在 markdown2 里不是可靠“块元素”

- 现象：`<p><details>...` 结构不合法，折叠失效或排版怪。
- 解法：在 HTML 后处理阶段包裹（字符串替换/正则替换），不要在 Markdown 源里硬写 `<details>`。
- 本项目采用：对 codehilite 包一层 `<details class="codeblock">`，对指定段落（如祖暅原理）按 marker 折叠。

## 3) 互动图（纯 JS + SVG）通用套路

### 3.1 先画“能看见的几何”，再谈“公式正确”

- 每个 widget 至少要画出：曲线（path）、关键点（焦点/点 P）、关键线段（PF、到准线垂足等）、读数（数值 + 误差）。
- 经验：如果画面只有文字/轴，通常是“忘画曲线”或“JS 抛异常导致 redraw 中断”。

### 3.2 坐标系与 viewBox（别让线宽/字号失控）

- 用 world 坐标，SVG 里统一 y 翻转：`y_svg = -y_world`（线段/圆/文本都要一致）。
- 线条：加 `vector-effect: non-scaling-stroke`，并把 `stroke-width` 设小一点（避免“线宽很粗”的反馈）。
- 文本：SVG 的字体大小会跟 viewBox 缩放产生错觉；务必显式设置 `font-size`（推荐用 `px`）。
  - 典型 bug：viewBox 是 `[-6,6]` 这种“小世界”，默认 `16` 会显得“巨大字母占满屏”。

### 3.3 Slider 交互与可读性

- 让 slider 控制“一个参数一件事”（例如点参数、a/b/p、旋转角 θ）。
- readout 同时展示：
  - 距离和/差守恒：`PF1+PF2` 或 `|PF1-PF2|`
  - 焦点-准线离心率：`PF / dist(P, directrix) ≈ e`
  - 误差：`computed - target`（用来证明“守恒”是真的，不是画出来像）

### 3.4 调试建议（最快定位）

- 打开浏览器控制台看首个异常（一个 widget 抛错可能阻断后续 mount）。
- 用 `rg` 确认：
  - HTML 中有 `div id="widget-..."` 占位
  - JS 中有 `mount...()` 实现
  - `init()` 里确实调用了 `mount...()`
- 做语法检查：把 HTML 最后那段 `<script>` 提取成临时文件，用 `node --check` 检查；不要把临时文件提交进 Git。

## 4) “折叠默认收起”经验

- 代码块：基于 `markdown2` 输出的 `<div class="codehilite">...` 做 HTML 后处理包 `<details class="codeblock">`。
- 局部说明（如祖暅原理）：用稳定 marker 在 HTML 中“截取并包裹”成 `<details class="fold">`。
- 原则：折叠逻辑尽量做“幂等”（重复生成不会越包越多），并且 marker 要足够独特。

## 5) 每次迭代结束前的检查清单

- 跑生成：`python3 scripts/generate_conics_report.py` 无报错。
- `dist/conics_report.html` 里：
  - 关键公式能被 MathJax 正确渲染（尤其含 `_`、`<` `>` 的位置）。
  - 新 widget 的占位 div + JS mount + init 调用都存在。
  - 折叠生效：代码块/祖暅原理默认收起。
- 视觉层面：
  - 线宽不过粗；缩放后线宽不变（`vector-effect` 生效）。
  - SVG 文本不“巨大盖住画面”（显式 `font-size`）。
  - 互动图读数显示“守恒量 + 误差”。

