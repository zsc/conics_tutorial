from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path

import markdown2
import numpy as np


@dataclass(frozen=True)
class ExampleAResult:
    disc: float
    kind: str
    center: tuple[float, float]
    theta_deg: float
    eigenvalues: tuple[float, float]
    a: float
    b: float
    c: float
    e: float
    foci: tuple[tuple[float, float], tuple[float, float]]


@dataclass(frozen=True)
class ExampleBResult:
    p: float
    point: tuple[float, float]
    focus: tuple[float, float]
    directrix_x: float
    pf: float
    dist_directrix: float
    tangent: str


@dataclass(frozen=True)
class HyperbolaResult:
    a: float
    b: float
    c: float
    e: float
    asymptote_slope: float


def classify_conic(A: float, B: float, C: float) -> tuple[float, str]:
    disc = B * B - 4 * A * C
    if disc < 0:
        kind = "ellipse-type (Δ2 < 0)"
    elif disc == 0:
        kind = "parabola-type (Δ2 = 0)"
    else:
        kind = "hyperbola-type (Δ2 > 0)"
    return disc, kind


def _plot_example_a(out_path: Path) -> ExampleAResult:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    # 9x^2 + 4xy + 6y^2 -54x -24y +36 = 0
    A, B, C, D, E, F0 = 9.0, 4.0, 6.0, -54.0, -24.0, 36.0
    disc, kind = classify_conic(A, B, C)

    Q = np.array([[A, B / 2], [B / 2, C]], dtype=float)
    p = np.array([D, E], dtype=float)
    M = np.array([[2 * A, B], [B, 2 * C]], dtype=float)
    center = np.linalg.solve(M, -p)

    eigvals, eigvecs = np.linalg.eigh(Q)
    F_at_center = (center @ Q @ center) + (p @ center) + F0
    rhs = -F_at_center
    axes_sq = rhs / eigvals  # semi-axis^2 along each eigenvector

    idx_major = int(np.argmax(axes_sq))
    idx_minor = 1 - idx_major
    a2 = float(axes_sq[idx_major])
    b2 = float(axes_sq[idx_minor])
    a = math.sqrt(a2)
    b = math.sqrt(b2)
    c = math.sqrt(max(a2 - b2, 0.0))
    e = c / a if a != 0 else float("nan")

    v_major = eigvecs[:, idx_major]
    v_minor = eigvecs[:, idx_minor]
    V = np.column_stack([v_major, v_minor])  # columns: major, minor

    t = np.linspace(0, 2 * np.pi, 900)
    U = np.vstack([a * np.cos(t), b * np.sin(t)])  # coords in principal basis
    pts = (V @ U).T + center

    f1 = center + c * v_major
    f2 = center - c * v_major

    theta = 0.5 * math.atan2(B, A - C)  # tan(2θ)=B/(A-C)
    theta_deg = math.degrees(theta)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    ax.plot(pts[:, 0], pts[:, 1], linewidth=2.2, label="Conic A (ellipse)")
    ax.scatter([center[0]], [center[1]], s=45, zorder=3, label="center")
    ax.scatter([f1[0], f2[0]], [f1[1], f2[1]], s=45, zorder=3, label="foci")

    # Principal axes
    ax.plot(
        [center[0] - a * v_major[0], center[0] + a * v_major[0]],
        [center[1] - a * v_major[1], center[1] + a * v_major[1]],
        linewidth=2,
        label="major axis",
    )
    ax.plot(
        [center[0] - b * v_minor[0], center[0] + b * v_minor[0]],
        [center[1] - b * v_minor[1], center[1] + b * v_minor[1]],
        linewidth=2,
        label="minor axis",
    )

    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Example A: General conic → ellipse (center, axes, foci)")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)

    return ExampleAResult(
        disc=float(disc),
        kind=kind,
        center=(float(center[0]), float(center[1])),
        theta_deg=float(theta_deg),
        eigenvalues=(float(eigvals[0]), float(eigvals[1])),
        a=a,
        b=b,
        c=c,
        e=e,
        foci=((float(f1[0]), float(f1[1])), (float(f2[0]), float(f2[1]))),
    )


def _plot_example_b(out_path: Path) -> ExampleBResult:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    # y^2 = 8x  => 4p = 8 => p = 2
    p = 2.0
    tt = np.linspace(-4, 4, 500)
    x_par = p * tt**2
    y_par = 2 * p * tt

    t0 = 1.0
    P = np.array([p * t0**2, 2 * p * t0], dtype=float)  # (2,4)
    F = np.array([p, 0.0], dtype=float)  # (2,0)
    directrix_x = -p

    pf = float(np.linalg.norm(P - F))
    dist_directrix = float(abs(P[0] - directrix_x))

    x_line = np.linspace(-3, 10, 200)
    y_tan = (x_line + p * (t0**2)) / t0  # t0*y = x + p*t0^2

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    ax.plot(x_par, y_par, linewidth=2.2, label=r"parabola $y^2=8x$")
    ax.scatter([F[0]], [F[1]], s=45, zorder=3, label="focus (2,0)")
    ax.axvline(directrix_x, linestyle="--", linewidth=1.8, label="directrix x=-2")
    ax.scatter([P[0]], [P[1]], s=45, zorder=3, label="P at t=1 (2,4)")
    ax.plot(x_line, y_tan, linewidth=2, label="tangent at P: y=x+2")

    # Latus rectum endpoints (x=p, y=±2p)
    ax.scatter([p, p], [2 * p, -2 * p], s=35, zorder=3, label="latus rectum ends")

    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Example B: Parabola with focus/directrix/tangent")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)

    return ExampleBResult(
        p=p,
        point=(float(P[0]), float(P[1])),
        focus=(float(F[0]), float(F[1])),
        directrix_x=float(directrix_x),
        pf=pf,
        dist_directrix=dist_directrix,
        tangent="y = x + 2",
    )


def _plot_hyperbola(out_path: Path) -> HyperbolaResult:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    a, b = 3.0, 2.0
    c = math.sqrt(a * a + b * b)
    e = c / a

    x1 = np.linspace(-10, -a - 1e-6, 600)
    y1 = b * np.sqrt(x1**2 / a**2 - 1)
    x2 = np.linspace(a + 1e-6, 10, 600)
    y2 = b * np.sqrt(x2**2 / a**2 - 1)

    x_as = np.linspace(-10, 10, 2)
    slope = b / a
    y_as = slope * x_as

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    ax.plot(x1, y1, linewidth=2.2, label=r"hyperbola $x^2/9 - y^2/4 = 1$")
    ax.plot(x1, -y1, linewidth=2.2)
    ax.plot(x2, y2, linewidth=2.2)
    ax.plot(x2, -y2, linewidth=2.2)
    ax.plot(x_as, y_as, linestyle="--", linewidth=1.8, label="asymptotes")
    ax.plot(x_as, -y_as, linestyle="--", linewidth=1.8)
    ax.scatter([c, -c], [0, 0], s=45, zorder=3, label="foci")

    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Bonus: Hyperbola with asymptotes and foci")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)

    return HyperbolaResult(a=a, b=b, c=c, e=e, asymptote_slope=slope)


_PMATRIX_FIXES: dict[str, str] = {
    # Q matrix and p vector in §1.1
    r"Q=\begin{pmatrix}A&\frac B2\ \frac B2&C\end{pmatrix}": (
        # NOTE: Markdown will collapse backslashes, so we intentionally use
        # double-escaped row separators here (\\\\ -> \\ in final HTML/TeX).
        r"Q=\begin{pmatrix}A&\frac{B}{2}\\\\ \frac{B}{2}&C\end{pmatrix}"
    ),
    r"\mathbf{p}=\begin{pmatrix}D\E\end{pmatrix}": r"\mathbf{p}=\begin{pmatrix}D\\\\E\end{pmatrix}",
    # Column vector (D,E) appears elsewhere too (e.g. §1.5).
    r"\begin{pmatrix}D\E\end{pmatrix}": r"\begin{pmatrix}D\\\\E\end{pmatrix}",
    # Center existence (§1.5)
    r"\begin{pmatrix}2A&B\B&2C\end{pmatrix}": r"\begin{pmatrix}2A&B\\\\B&2C\end{pmatrix}",
    r"\begin{pmatrix}x_0\y_0\end{pmatrix}": r"\begin{pmatrix}x_0\\\\y_0\end{pmatrix}",
    # Example A center solve (§4A)
    r"\begin{pmatrix}18&4\4&12\end{pmatrix}": r"\begin{pmatrix}18&4\\\\4&12\end{pmatrix}",
    r"\begin{pmatrix}54\24\end{pmatrix}": r"\begin{pmatrix}54\\\\24\end{pmatrix}",
}


def _fix_latex_typos_outside_code(text: str) -> str:
    for src, dst in _PMATRIX_FIXES.items():
        text = text.replace(src, dst)
    # Remove Markdown line-break backslashes that would otherwise leak into TeX.
    text = re.sub(r",\\\s*\n", ",\n", text)
    return text


def _convert_display_math_blocks(text: str) -> str:
    lines: list[str] = []
    in_code = False
    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        stripped = line.strip()
        if stripped.startswith("```"):
            in_code = not in_code
            lines.append(line)
            continue
        if not in_code and stripped == "[":
            lines.append("$$")
            continue
        if not in_code and stripped == "]":
            lines.append("$$")
            continue
        lines.append(line)
    return "\n".join(lines) + "\n"


_INLINE_MATH_HINT = re.compile(r"[\\^_=<>\u00B1*/]|[0-9]")


def _looks_like_math(content: str) -> bool:
    s = content.strip()
    if not s:
        return False
    if "http://" in s or "https://" in s:
        return False
    if "\\" in s:
        return True
    if _INLINE_MATH_HINT.search(s) is None:
        # permit single-letter variables like (Q), (x), (y), (e)
        return bool(re.fullmatch(r"[A-Za-z]", s))
    # If it contains CJK but no TeX commands, it's likely not math.
    if re.search(r"[\u4e00-\u9fff]", s):
        return False
    return True


def _convert_inline_math_parentheses(text: str) -> str:
    out_lines: list[str] = []
    in_code = False
    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        stripped = line.strip()
        if stripped.startswith("```"):
            in_code = not in_code
            out_lines.append(line)
            continue
        if in_code:
            out_lines.append(line)
            continue

        res: list[str] = []
        i = 0
        while i < len(line):
            ch = line[i]
            if ch != "(":
                res.append(ch)
                i += 1
                continue

            # Find matching ')' with nesting.
            depth = 0
            j = i
            while j < len(line):
                if line[j] == "(":
                    depth += 1
                elif line[j] == ")":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1

            if j >= len(line) or line[j] != ")":
                res.append(ch)
                i += 1
                continue

            content = line[i + 1 : j]
            if _looks_like_math(content):
                # NOTE: Markdown treats "\(" as an escape, so we output "\\(" to
                # preserve "\(" in the final HTML for MathJax.
                res.append(r"\\(" + content + r"\\)")
            else:
                res.append(line[i : j + 1])
            i = j + 1

        out_lines.append("".join(res))

    return "\n".join(out_lines) + "\n"


def _build_extra_formula_table_markdown() -> str:
    # Keep it “checklist-like”: formula → scenario → one-liner.
    # Use \\( ... \\) here (not \( ... \)) so Markdown doesn't eat the backslashes.
    return r"""
## 7) 更冷门但可用于解题的公式清单（公式 → 适用场景 → 一行例题）

| 公式 | 适用场景 | 一行例题 |
| --- | --- | --- |
| 极线/极点（齐次）：若二次曲线为 \\(X^T A X=0\\)，点 \\(P\\) 的极线：\\(P^T A X=0\\)。 | 求切线、判共线/共点（射影味但解析也能算） | 已知 \\(A\\) 与 \\(P\\)，直接矩阵乘法得到直线系数。 |
| 极线（一般式快速写法）：对 \\(Ax^2+Bxy+Cy^2+Dx+Ey+F=0\\)，点 \\(P(x_1,y_1)\\) 的极线可由“代换规则”写出：\\(x^2\\!\\to xx_1,\\ y^2\\!\\to yy_1,\\ xy\\!\\to \\frac12(xy_1+x_1y),\\ x\\!\\to \\frac12(x+x_1),\\ y\\!\\to \\frac12(y+y_1)\\)。 | 不想上齐次坐标时的手算捷径 | 把 \\(x_1,y_1\\) 代入即可得到极线。 |
| 弦中点公式（通用）：对一般二次曲线 \\(S=0\\)，若 \\(M(x_1,y_1)\\) 是某条弦的中点，则该弦满足 \\(T=S_1\\)（\\(T\\) 为“二次项双线性化 + 一次项取平均”的表达）。 | “给中点求弦方程/方向” 常见套路 | 圆 \\(x^2+y^2=r^2\\) 退化为 \\(xx_1+yy_1=x_1^2+y_1^2\\)。 |
| 直线与圆锥曲线交点判别：设直线代入得到二次方程 \\(at^2+bt+c=0\\)，判别式 \\(b^2-4ac\\) 决定交点个数；韦达 \\(t_1+t_2=-b/a,\\ t_1t_2=c/a\\)。 | “切线判定/弦长/中点” 的统一入口 | 用韦达直接得到弦中点参数 \\(\\frac{t_1+t_2}{2}\\)。 |
| 仿射不变量（常用口径）：二次项矩阵 \\(Q\\) 的惯性指数（正/负/零特征值个数）在可逆线性变换下不变；\\(\\det(Q)=0\\) 与否对应“抛物线型 vs 有中心”。 | 判断类型在更一般的线性变换下是否改变 | \\(Q\\) 非退化且正定 \\(\\Rightarrow\\) 椭圆型（经仿射可化为圆）。 |
| 过定点弦族参数化：固定点 \\(P\\)，取过 \\(P\\) 的任意方向直线 \\(L(m)\\)，与圆锥曲线交于 \\(P,Q\\)，用代入得到关于参数的二次式，另一根即为 \\(Q(m)\\)。 | “弦族/二次曲线上的点随斜率变化” | 竞赛常用：把“过定点”变成“一元二次另一根”。 |
"""


def _build_polar_section_markdown() -> str:
    # Use \\( ... \\) so Markdown doesn't eat backslashes.
    return r"""
## 8) 极线/极点（Polar line / Pole）展开说明（更偏“可计算”版本）

### 8.1 为什么它好用：把“切线/共线判定”变成线代

极线/极点本质上是二次曲线的“对偶”结构：你给一个点 \\(P\\)，它会对应一条直线 \\(\ell\\)（极线）；反过来你给一条直线 \\(\ell\\)，也会对应一个点 \\(P\\)（极点）。在解析几何里它最常用的三个等价用法是：

1. **切线是极线的特例**：若点 \\(P\\) 在曲线上，则 \\(P\\) 的极线就是该点处切线。
2. **互伴性（对称性）**：\\(P\\) 在 \\(Q\\) 的极线上 \\(\Leftrightarrow\\) \\(Q\\) 在 \\(P\\) 的极线上（非常适合判共线/共点）。
3. **“外点两切线”的弦联络**：对外点 \\(P\\)，从 \\(P\\) 引曲线两条切线，切点连线（接触弦）就是 \\(P\\) 的极线。

### 8.2 齐次坐标的一行公式（推荐记这个）

把平面点写成齐次向量 \\(X=(x,y,1)^T\\)。一般二次曲线
$$
Ax^2+Bxy+Cy^2+Dx+Ey+F=0
$$
可以写成矩阵形式
$$
X^T A X = 0,\qquad
A=\begin{pmatrix}
A & \frac{B}{2} & \frac{D}{2} \\\\
\frac{B}{2} & C & \frac{E}{2} \\\\
\frac{D}{2} & \frac{E}{2} & F
\end{pmatrix}.
$$

给定点 \\(P=(x_1,y_1,1)^T\\)，它的极线就是
$$
\ell:\ (A P)^T X = 0.
$$
写成通常的直线式就是 \\(\ell_1 x+\ell_2 y+\ell_3=0\\)（其中 \\((\ell_1,\ell_2,\ell_3)^T=A P\\)）。

> 直观理解：\\(A\\) 提供了一个“二次型内积”。极线就是把这个“内积”对 \\(P\\) 做一次线性化。

### 8.3 不想上齐次？用“代换规则”直接写极线（手算友好）

对 \\(S(x,y)=Ax^2+Bxy+Cy^2+Dx+Ey+F\\)，点 \\(P(x_1,y_1)\\) 的极线可由下面替换得到（把 \\(S\\) 里的每一项“线性化”）：

- \\(x^2\\to xx_1\\)，\\(y^2\\to yy_1\\)
- \\(xy\\to \\frac12(xy_1+x_1y)\\)
- \\(x\\to \\frac12(x+x_1)\\)，\\(y\\to \\frac12(y+y_1)\\)
- 常数 \\(F\\) 保持不变

替换后得到的方程就是极线方程 \\(\ell=0\\)。

### 8.4 两个最经典的例子（顺便记极线=切线）

**例 1：单位圆** \\(x^2+y^2=1\\)。点 \\(P(x_1,y_1)\\) 的极线为
$$
xx_1+yy_1=1.
$$
若 \\(P\\) 在圆上（即 \\(x_1^2+y_1^2=1\\)），这条直线正是圆在 \\(P\\) 处的切线。

**例 2：椭圆** \\(\frac{x^2}{a^2}+\frac{y^2}{b^2}=1\\)。点 \\(P(x_1,y_1)\\) 的极线为
$$
\frac{xx_1}{a^2}+\frac{yy_1}{b^2}=1.
$$
同样，若 \\(P\\) 在椭圆上，这就是切线方程（你在解析几何题里经常见到的那条）。

### 8.5 一个常见技巧：用极线快速处理“外点两切线”

若 \\(P\\) 是外点，从 \\(P\\) 作二次曲线两条切线，切点为 \\(T_1,T_2\\)。那么：

- **接触弦** \\(T_1T_2\\) 的方程就是 \\(P\\) 的极线；
- 反过来，如果你先求出 \\(P\\) 的极线，再和曲线联立求交点，就能直接拿到 \\(T_1,T_2\\)。
"""


def _build_generated_section(a: ExampleAResult, b: ExampleBResult, h: HyperbolaResult) -> str:
    return rf"""
---

## 6) 自动计算结果与互动图示（由脚本生成）

### 6.1 算例 A：一般二次曲线 → 标准椭圆（数值提取）

- 判别式：\\(\Delta_2 = {a.disc:.0f}\\) → {a.kind}
- 中心：\\((x_0,y_0)=({a.center[0]:.2f},{a.center[1]:.2f})\\)
- 去交叉项旋转角：\\(\theta \approx {a.theta_deg:.3f}^\\circ\\)
- 特征值：\\(\lambda_1,\lambda_2 \approx ({a.eigenvalues[0]:.3g},{a.eigenvalues[1]:.3g})\\)
- 半轴：\\(a\approx {a.a:.4f},\ b\approx {a.b:.4f}\\)，\\(c=\\sqrt{{a^2-b^2}}\\approx {a.c:.4f}\\)，离心率 \\(e=c/a\\approx {a.e:.4f}\\)
- 焦点（原坐标系）：\\(F_1\\approx({a.foci[0][0]:.3f},{a.foci[0][1]:.3f})\\)，\\(F_2\\approx({a.foci[1][0]:.3f},{a.foci[1][1]:.3f})\\)

<div id="widget-example-a" class="widget"></div>

<details class="fallback">
  <summary>静态图（备份）</summary>
  <p><img src="assets/example_a.png" alt="Example A ellipse plot" /></p>
</details>

### 6.2 算例 B：抛物线焦点-准线性质与切线验证

- \\(\,y^2=8x \\Rightarrow p={b.p:g}\\)，焦点 \\(F=({b.focus[0]:g},{b.focus[1]:g})\\)，准线 \\(x={b.directrix_x:g}\\)
- 取参数点 \\(t=1\\)：\\(P=({b.point[0]:g},{b.point[1]:g})\\)
- 数值验证：\\(PF\\approx {b.pf:.6g}\\)，点到准线距离 \\(\\approx {b.dist_directrix:.6g}\\)
- 切线：\({b.tangent}\)

<div id="widget-parabola" class="widget"></div>

<details class="fallback">
  <summary>静态图（备份）</summary>
  <p><img src="assets/example_b.png" alt="Example B parabola plot" /></p>
</details>

### 6.3 Bonus：标准双曲线 + 渐近线 + 焦点

- \\(\frac{{x^2}}{{a^2}}-\frac{{y^2}}{{b^2}}=1\\)，取 \\(a={h.a:g},\\,b={h.b:g}\\)
- \\(c=\\sqrt{{a^2+b^2}}\\approx {h.c:.4f}\\)，离心率 \\(e=c/a\\approx {h.e:.4f}\\)
- 渐近线：\\(y=\\pm \\frac{{b}}{{a}}x=\\pm {h.asymptote_slope:.4f}x\\)

<div id="widget-hyperbola" class="widget"></div>

<details class="fallback">
  <summary>静态图（备份）</summary>
  <p><img src="assets/hyperbola.png" alt="Hyperbola plot" /></p>
</details>

### 6.4 Director circle（导圆）互动图

<div id="widget-director-circle" class="widget"></div>
"""


def _preprocess_agents_markdown(md: str) -> str:
    md = _fix_latex_typos_outside_code(md)
    md = _convert_display_math_blocks(md)
    md = _convert_inline_math_parentheses(md)
    return md


def _escape_markdown_emphasis_in_math(md: str) -> str:
    """
    markdown2 will interpret "_" as emphasis markers. That corrupts TeX like
    "\\lambda_1" into HTML "<em>" tags, which then breaks MathJax parsing.

    Workaround: inside math regions, replace "_" with "&#95;" so Markdown won't
    treat it as emphasis, while the browser will still decode it back to "_"
    before MathJax runs.
    """

    def escape_inline_math(line: str) -> str:
        open_delim = r"\\("
        close_delim = r"\\)"
        out: list[str] = []
        i = 0
        while True:
            start = line.find(open_delim, i)
            if start < 0:
                out.append(line[i:])
                break
            end = line.find(close_delim, start + len(open_delim))
            if end < 0:
                out.append(line[i:])
                break
            out.append(line[i:start])
            content = line[start + len(open_delim) : end]
            content = content.replace("_", "&#95;")
            # Avoid Markdown's raw-HTML parsing inside math (e.g. 0<e<1).
            content = content.replace("<", "&lt;").replace(">", "&gt;")
            out.append(open_delim + content + close_delim)
            i = end + len(close_delim)
        return "".join(out)

    out_lines: list[str] = []
    in_code = False
    in_display = False
    for raw_line in md.splitlines():
        line = raw_line.rstrip("\n")
        stripped = line.strip()
        if stripped.startswith("```"):
            in_code = not in_code
            out_lines.append(line)
            continue
        if in_code:
            out_lines.append(line)
            continue
        if stripped == "$$":
            in_display = not in_display
            out_lines.append(line)
            continue
        if in_display:
            escaped = line.replace("_", "&#95;")
            escaped = escaped.replace("<", "&lt;").replace(">", "&gt;")
            out_lines.append(escaped)
            continue
        out_lines.append(escape_inline_math(line))

    return "\n".join(out_lines) + "\n"


def _collapse_code_blocks(body_html: str) -> str:
    """
    Wrap markdown2/pygments code blocks in <details> so they are collapsed by default.
    """

    pattern = re.compile(r'(<div class="codehilite">.*?</div>)', re.DOTALL)

    def repl(m: re.Match[str]) -> str:
        block = m.group(1)
        return (
            '<details class="codeblock">'
            '<summary>代码（点击展开）</summary>'
            f"{block}"
            "</details>"
        )

    return pattern.sub(repl, body_html)


def _interactive_js(exa: ExampleAResult) -> str:
    # Keep the JS standalone (no extra dependencies).
    cx, cy = exa.center
    f1x, f1y = exa.foci[0]
    f2x, f2y = exa.foci[1]
    return f"""
(() => {{
  'use strict';

  const NS = 'http://www.w3.org/2000/svg';

  const EXAMPLE_A = {{
    a: {exa.a:.12g},
    b: {exa.b:.12g},
    center: {{ x: {cx:.12g}, y: {cy:.12g} }},
    f1: {{ x: {f1x:.12g}, y: {f1y:.12g} }},
    f2: {{ x: {f2x:.12g}, y: {f2y:.12g} }},
  }};

  function el(tag, attrs = {{}}, parent = null) {{
    const node = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs)) {{
      if (v === null || v === undefined) continue;
      if (k === 'text') node.textContent = String(v);
      else node.setAttribute(k, String(v));
    }}
    if (parent) parent.appendChild(node);
    return node;
  }}

  function svg(tag, attrs = {{}}, parent = null) {{
    const node = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs)) {{
      if (v === null || v === undefined) continue;
      node.setAttribute(k, String(v));
    }}
    if (parent) parent.appendChild(node);
    return node;
  }}

  function fmt(x, digits = 6) {{
    if (!Number.isFinite(x)) return 'NaN';
    const abs = Math.abs(x);
    if (abs !== 0 && (abs < 1e-4 || abs >= 1e6)) return x.toExponential(4);
    return x.toFixed(digits);
  }}

  function dist(p, q) {{
    const dx = p.x - q.x;
    const dy = p.y - q.y;
    return Math.hypot(dx, dy);
  }}

  function pathFromWorld(points) {{
    if (!points.length) return '';
    let d = `M ${{points[0].x}} ${{-points[0].y}}`;
    for (let i = 1; i < points.length; i++) d += ` L ${{points[i].x}} ${{-points[i].y}}`;
    return d;
  }}

  function worldBounds(points) {{
    let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
    for (const p of points) {{
      xmin = Math.min(xmin, p.x);
      xmax = Math.max(xmax, p.x);
      ymin = Math.min(ymin, p.y);
      ymax = Math.max(ymax, p.y);
    }}
    return {{ xmin, xmax, ymin, ymax }};
  }}

  function setViewBox(svgNode, bounds, pad) {{
    const xmin = bounds.xmin - pad;
    const xmax = bounds.xmax + pad;
    const ymin = bounds.ymin - pad;
    const ymax = bounds.ymax + pad;
    const x = xmin;
    const y = -ymax; // flip-y
    const w = xmax - xmin;
    const h = ymax - ymin;
    svgNode.setAttribute('viewBox', `${{x}} ${{y}} ${{w}} ${{h}}`);
  }}

  function lineFromPoints(line, p, q) {{
    line.setAttribute('x1', p.x);
    line.setAttribute('y1', -p.y);
    line.setAttribute('x2', q.x);
    line.setAttribute('y2', -q.y);
  }}

  function circleAt(circ, p) {{
    circ.setAttribute('cx', p.x);
    circ.setAttribute('cy', -p.y);
  }}

  function buildSliderRow(parent, opts) {{
    const row = el('div', {{ class: 'control-row' }}, parent);
    el('div', {{ class: 'control-label', text: opts.label }}, row);
    const input = el('input', {{
      type: 'range',
      min: opts.min,
      max: opts.max,
      step: opts.step,
      value: opts.value,
    }}, row);
    const val = el('div', {{ class: 'control-value' }}, row);
    const updateVal = () => {{
      val.textContent = opts.format ? opts.format(Number(input.value)) : String(input.value);
    }};
    input.addEventListener('input', () => {{
      updateVal();
      opts.onInput?.(Number(input.value));
    }});
    updateVal();
    return input;
  }}

  function mountExampleA() {{
    const container = document.getElementById('widget-example-a');
    if (!container) return;
    container.innerHTML = '';
    container.classList.add('widget');

    el('div', {{ class: 'widget-title', text: 'Example A 椭圆：验证 PF₁+PF₂ = 2a（滑块改变点 P）' }}, container);
    const controls = el('div', {{ class: 'widget-controls' }}, container);
    const svgRoot = svg('svg', {{ class: 'conic-svg', role: 'img' }}, container);
    const readout = el('div', {{ class: 'widget-readout' }}, container);

    const a = EXAMPLE_A.a;
    const b = EXAMPLE_A.b;
    const center = EXAMPLE_A.center;
    const f1 = EXAMPLE_A.f1;
    const f2 = EXAMPLE_A.f2;
    const c = dist(center, f1);
    const vMajor = {{ x: (f1.x - center.x) / c, y: (f1.y - center.y) / c }};
    const vMinor = {{ x: -vMajor.y, y: vMajor.x }};

    const pathCurve = svg('path', {{ class: 'curve' }}, svgRoot);
    const seg1 = svg('line', {{ class: 'seg' }}, svgRoot);
    const seg2 = svg('line', {{ class: 'seg' }}, svgRoot);
    const circF1 = svg('circle', {{ class: 'focus', r: 0.12 }}, svgRoot);
    const circF2 = svg('circle', {{ class: 'focus', r: 0.12 }}, svgRoot);
    const circP = svg('circle', {{ class: 'point', r: 0.14 }}, svgRoot);

    const ellipsePts = [];
    const N = 480;
    for (let i = 0; i <= N; i++) {{
      const t = (i / N) * 2 * Math.PI;
      const u = a * Math.cos(t);
      const v = b * Math.sin(t);
      ellipsePts.push({{
        x: center.x + vMajor.x * u + vMinor.x * v,
        y: center.y + vMajor.y * u + vMinor.y * v,
      }});
    }}
    pathCurve.setAttribute('d', pathFromWorld(ellipsePts));
    circleAt(circF1, f1);
    circleAt(circF2, f2);

    const bounds = worldBounds([...ellipsePts, f1, f2]);
    setViewBox(svgRoot, bounds, 0.9);

    const tSlider = buildSliderRow(controls, {{
      label: 't（角参数）',
      min: 0,
      max: (2 * Math.PI).toFixed(4),
      step: 0.01,
      value: (Math.PI / 3).toFixed(4),
      format: (t) => `${{fmt(t, 3)}} rad  (≈ ${{fmt(t * 180 / Math.PI, 1)}}°)`,
      onInput: update,
    }});

    function update() {{
      const t = Number(tSlider.value);
      const u = a * Math.cos(t);
      const v = b * Math.sin(t);
      const P = {{
        x: center.x + vMajor.x * u + vMinor.x * v,
        y: center.y + vMajor.y * u + vMinor.y * v,
      }};

      circleAt(circP, P);
      lineFromPoints(seg1, P, f1);
      lineFromPoints(seg2, P, f2);

      const d1 = dist(P, f1);
      const d2 = dist(P, f2);
      const sum = d1 + d2;
      const target = 2 * a;
      const err = sum - target;
      readout.innerHTML = `
        <div><span class="k">PF₁</span>=${{fmt(d1)}}, <span class="k">PF₂</span>=${{fmt(d2)}}, <span class="k">PF₁+PF₂</span>=${{fmt(sum)}}；常数 <span class="k">2a</span>=${{fmt(target)}}；误差=${{fmt(err, 8)}}</div>
      `;
    }}

    update();
  }}

  function mountParabola() {{
    const container = document.getElementById('widget-parabola');
    if (!container) return;
    container.innerHTML = '';
    container.classList.add('widget');

    el('div', {{ class: 'widget-title', text: '抛物线：验证 PF = 到准线距离（e=1）' }}, container);
    const controls = el('div', {{ class: 'widget-controls' }}, container);
    const svgRoot = svg('svg', {{ class: 'conic-svg', role: 'img' }}, container);
    const readout = el('div', {{ class: 'widget-readout' }}, container);

    const pathCurve = svg('path', {{ class: 'curve' }}, svgRoot);
    const directrix = svg('line', {{ class: 'directrix' }}, svgRoot);
    const segFocus = svg('line', {{ class: 'seg' }}, svgRoot);
    const segDir = svg('line', {{ class: 'seg2' }}, svgRoot);
    const circF = svg('circle', {{ class: 'focus', r: 0.12 }}, svgRoot);
    const circP = svg('circle', {{ class: 'point', r: 0.14 }}, svgRoot);

    const tMax = 3.5;
    const N = 520;

    let p = 2.0;
    let t = 1.0;

    const pSlider = buildSliderRow(controls, {{
      label: 'p',
      min: 0.5,
      max: 5.0,
      step: 0.05,
      value: p,
      format: (v) => fmt(v, 2),
      onInput: (v) => {{ p = v; redraw(); }},
    }});

    const tSlider = buildSliderRow(controls, {{
      label: 't（参数点）',
      min: -tMax,
      max: tMax,
      step: 0.01,
      value: t,
      format: (v) => fmt(v, 2),
      onInput: (v) => {{ t = v; redraw(); }},
    }});

    function redraw() {{
      p = Number(pSlider.value);
      t = Number(tSlider.value);

      const focus = {{ x: p, y: 0 }};
      const directrixX = -p;
      const P = {{ x: p * t * t, y: 2 * p * t }};
      const Q = {{ x: directrixX, y: P.y }}; // foot to directrix (horizontal segment)

      // curve
      const pts = [];
      for (let i = 0; i <= N; i++) {{
        const tt = -tMax + (2 * tMax * i) / N;
        pts.push({{ x: p * tt * tt, y: 2 * p * tt }});
      }}
      pathCurve.setAttribute('d', pathFromWorld(pts));

      // features
      circleAt(circF, focus);
      circleAt(circP, P);
      lineFromPoints(segFocus, P, focus);
      lineFromPoints(segDir, P, Q);

      // directrix line segment
      const ySpan = 2 * p * tMax;
      directrix.setAttribute('x1', directrixX);
      directrix.setAttribute('y1', -(-ySpan));
      directrix.setAttribute('x2', directrixX);
      directrix.setAttribute('y2', -(ySpan));

      // viewBox
      const bounds = worldBounds([...pts, focus, P, Q, {{ x: directrixX, y: -ySpan }}, {{ x: directrixX, y: ySpan }}]);
      setViewBox(svgRoot, bounds, 0.9);

      // invariants
      const dFocus = dist(P, focus);
      const dDir = Math.abs(P.x - directrixX);
      const err = dFocus - dDir;
      readout.innerHTML = `
        <div><span class="k">P</span>=(${{fmt(P.x, 4)}}, ${{fmt(P.y, 4)}})，<span class="k">F</span>=(${{fmt(focus.x, 3)}},0)，准线 <span class="k">x</span>=${{fmt(directrixX, 3)}}</div>
        <div><span class="k">PF</span>=${{fmt(dFocus)}}, <span class="k">dist(P, directrix)</span>=${{fmt(dDir)}}；误差=${{fmt(err, 8)}}</div>
      `;
    }}

    redraw();
  }}

  function mountHyperbola() {{
    const container = document.getElementById('widget-hyperbola');
    if (!container) return;
    container.innerHTML = '';
    container.classList.add('widget');

    el('div', {{ class: 'widget-title', text: '双曲线：验证 |PF₁−PF₂| = 2a（滑块改变 a,b 与点 P）' }}, container);
    const controls = el('div', {{ class: 'widget-controls' }}, container);
    const svgRoot = svg('svg', {{ class: 'conic-svg', role: 'img' }}, container);
    const readout = el('div', {{ class: 'widget-readout' }}, container);

    const pathR = svg('path', {{ class: 'curve' }}, svgRoot);
    const pathL = svg('path', {{ class: 'curve' }}, svgRoot);
    const as1 = svg('line', {{ class: 'asym' }}, svgRoot);
    const as2 = svg('line', {{ class: 'asym' }}, svgRoot);
    const seg1 = svg('line', {{ class: 'seg' }}, svgRoot);
    const seg2 = svg('line', {{ class: 'seg' }}, svgRoot);
    const circF1 = svg('circle', {{ class: 'focus', r: 0.12 }}, svgRoot);
    const circF2 = svg('circle', {{ class: 'focus', r: 0.12 }}, svgRoot);
    const circP = svg('circle', {{ class: 'point', r: 0.14 }}, svgRoot);

    const Umax = 2.0;
    const N = 560;

    let a = 3.0;
    let b = 2.0;
    let u = 1.0;

    const aSlider = buildSliderRow(controls, {{
      label: 'a',
      min: 0.8,
      max: 6.0,
      step: 0.05,
      value: a,
      format: (v) => fmt(v, 2),
      onInput: (v) => {{ a = v; redraw(); }},
    }});
    const bSlider = buildSliderRow(controls, {{
      label: 'b',
      min: 0.5,
      max: 5.0,
      step: 0.05,
      value: b,
      format: (v) => fmt(v, 2),
      onInput: (v) => {{ b = v; redraw(); }},
    }});
    const uSlider = buildSliderRow(controls, {{
      label: 'u（参数）',
      min: -Umax,
      max: Umax,
      step: 0.01,
      value: u,
      format: (v) => fmt(v, 2),
      onInput: (v) => {{ u = v; redraw(); }},
    }});

    function redraw() {{
      a = Number(aSlider.value);
      b = Number(bSlider.value);
      u = Number(uSlider.value);

      const c = Math.hypot(a, b);
      const f1 = {{ x: c, y: 0 }};
      const f2 = {{ x: -c, y: 0 }};

      const P = {{ x: a * Math.cosh(u), y: b * Math.sinh(u) }}; // right branch

      // curves
      const ptsR = [];
      const ptsL = [];
      for (let i = 0; i <= N; i++) {{
        const uu = -Umax + (2 * Umax * i) / N;
        const x = a * Math.cosh(uu);
        const y = b * Math.sinh(uu);
        ptsR.push({{ x, y }});
        ptsL.push({{ x: -x, y }});
      }}
      pathR.setAttribute('d', pathFromWorld(ptsR));
      pathL.setAttribute('d', pathFromWorld(ptsL));

      // asymptotes
      const slope = b / a;
      const xSpan = a * Math.cosh(Umax);
      const A1 = {{ x: -xSpan, y: -slope * xSpan }};
      const A2 = {{ x: xSpan, y: slope * xSpan }};
      const B1 = {{ x: -xSpan, y: slope * xSpan }};
      const B2 = {{ x: xSpan, y: -slope * xSpan }};
      lineFromPoints(as1, A1, A2);
      lineFromPoints(as2, B1, B2);

      // points + segments
      circleAt(circF1, f1);
      circleAt(circF2, f2);
      circleAt(circP, P);
      lineFromPoints(seg1, P, f1);
      lineFromPoints(seg2, P, f2);

      const bounds = worldBounds([...ptsR, ...ptsL, f1, f2, P, A1, A2, B1, B2]);
      setViewBox(svgRoot, bounds, 1.2);

      const d1 = dist(P, f1);
      const d2 = dist(P, f2);
      const diff = Math.abs(d2 - d1);
      const target = 2 * a;
      const err = diff - target;
      readout.innerHTML = `
        <div><span class="k">PF₁</span>=${{fmt(d1)}}, <span class="k">PF₂</span>=${{fmt(d2)}}, <span class="k">|PF₁−PF₂|</span>=${{fmt(diff)}}；常数 <span class="k">2a</span>=${{fmt(target)}}；误差=${{fmt(err, 8)}}</div>
      `;
    }}

    redraw();
  }}

  function mountDirectorCircle() {{
    const container = document.getElementById('widget-director-circle');
    if (!container) return;
    container.innerHTML = '';
    container.classList.add('widget');

    el('div', {{ class: 'widget-title', text: '导圆：从圆上点引曲线切线互相垂直（滑块调参/选点）' }}, container);
    const controls = el('div', {{ class: 'widget-controls' }}, container);
    const svgRoot = svg('svg', {{ class: 'conic-svg', role: 'img' }}, container);
    const readout = el('div', {{ class: 'widget-readout' }}, container);

    const typeRow = el('div', {{ class: 'control-row' }}, controls);
    el('div', {{ class: 'control-label', text: '类型' }}, typeRow);
    const typeSel = el('select', {{ class: 'control-select' }}, typeRow);
    el('option', {{ value: 'ellipse', text: '椭圆' }}, typeSel);
    el('option', {{ value: 'hyperbola', text: '双曲线' }}, typeSel);
    const typeVal = el('div', {{ class: 'control-value', text: '' }}, typeRow);

    let a = 4.0;
    let b = 2.5;
    let phi = 0.8;

    const aSlider = buildSliderRow(controls, {{
      label: 'a',
      min: 1.2,
      max: 7.0,
      step: 0.05,
      value: a,
      format: (v) => fmt(v, 2),
      onInput: (v) => {{ a = v; redraw(); }},
    }});
    const bSlider = buildSliderRow(controls, {{
      label: 'b',
      min: 0.6,
      max: 6.0,
      step: 0.05,
      value: b,
      format: (v) => fmt(v, 2),
      onInput: (v) => {{ b = v; redraw(); }},
    }});
    const phiSlider = buildSliderRow(controls, {{
      label: '选点角 φ',
      min: 0,
      max: (2 * Math.PI).toFixed(4),
      step: 0.01,
      value: phi,
      format: (t) => `${{fmt(t, 3)}} rad  (≈ ${{fmt(t * 180 / Math.PI, 1)}}°)`,
      onInput: (v) => {{ phi = v; redraw(); }},
    }});

    const curve1 = svg('path', {{ class: 'curve' }}, svgRoot);
    const curve2 = svg('path', {{ class: 'curve' }}, svgRoot);
    const circleDir = svg('path', {{ class: 'dircircle' }}, svgRoot);
    const tang1 = svg('line', {{ class: 'tangent' }}, svgRoot);
    const tang2 = svg('line', {{ class: 'tangent' }}, svgRoot);
    const pointS = svg('circle', {{ class: 'point', r: 0.14 }}, svgRoot);
    const t1 = svg('circle', {{ class: 'focus', r: 0.11 }}, svgRoot);
    const t2 = svg('circle', {{ class: 'focus', r: 0.11 }}, svgRoot);

    function circlePath(R) {{
      // SVG arc circle centered at origin in world coords, then flip y.
      const y = 0;
      const d = [
        `M ${{R}} ${{-y}}`,
        `A ${{R}} ${{R}} 0 1 0 ${{-R}} ${{-y}}`,
        `A ${{R}} ${{R}} 0 1 0 ${{R}} ${{-y}}`,
      ].join(' ');
      return d;
    }}

    function conicEllipsePts(a, b) {{
      const pts = [];
      const N = 520;
      for (let i = 0; i <= N; i++) {{
        const t = (i / N) * 2 * Math.PI;
        pts.push({{ x: a * Math.cos(t), y: b * Math.sin(t) }});
      }}
      return pts;
    }}

    function conicHyperbolaPts(a, b) {{
      const ptsR = [];
      const ptsL = [];
      const U = 1.9;
      const N = 520;
      for (let i = 0; i <= N; i++) {{
        const u = -U + (2 * U * i) / N;
        const x = a * Math.cosh(u);
        const y = b * Math.sinh(u);
        ptsR.push({{ x, y }});
        ptsL.push({{ x: -x, y }});
      }}
      return {{ ptsR, ptsL }};
    }}

    function lineSegmentForAxByC(A, B, C, bounds) {{
      // Solve A x + B y + C = 0 within bounding box (world coords).
      const xmin = bounds.xmin, xmax = bounds.xmax, ymin = bounds.ymin, ymax = bounds.ymax;
      const pts = [];
      // Intersections with x = xmin/xmax
      if (Math.abs(B) > 1e-12) {{
        const y1 = (-C - A * xmin) / B;
        const y2 = (-C - A * xmax) / B;
        if (y1 >= ymin && y1 <= ymax) pts.push({{ x: xmin, y: y1 }});
        if (y2 >= ymin && y2 <= ymax) pts.push({{ x: xmax, y: y2 }});
      }}
      // Intersections with y = ymin/ymax
      if (Math.abs(A) > 1e-12) {{
        const x1 = (-C - B * ymin) / A;
        const x2 = (-C - B * ymax) / A;
        if (x1 >= xmin && x1 <= xmax) pts.push({{ x: x1, y: ymin }});
        if (x2 >= xmin && x2 <= xmax) pts.push({{ x: x2, y: ymax }});
      }}
      if (pts.length < 2) return null;
      // pick farthest pair
      let best = [pts[0], pts[1]];
      let bestD = -1;
      for (let i = 0; i < pts.length; i++) {{
        for (let j = i + 1; j < pts.length; j++) {{
          const d = dist(pts[i], pts[j]);
          if (d > bestD) {{
            bestD = d;
            best = [pts[i], pts[j]];
          }}
        }}
      }}
      return best;
    }}

    function intersectCircleLineUnit(n) {{
      // Unit circle x^2+y^2=1 with line n·x = 1 (n=(nx,ny)).
      const nx = n.x, ny = n.y;
      const nn = Math.hypot(nx, ny);
      if (!(nn > 1 + 1e-12)) return null;
      const x0 = nx / (nn * nn);
      const y0 = ny / (nn * nn);
      const h = Math.sqrt(1 - (x0 * x0 + y0 * y0));
      const dx = -ny / nn;
      const dy = nx / nn;
      return [
        {{ x: x0 + h * dx, y: y0 + h * dy }},
        {{ x: x0 - h * dx, y: y0 - h * dy }},
      ];
    }}

    function tangentPointsEllipseFromExternalPoint(P, a, b) {{
      // Scale to unit circle: X=x/a, Y=y/b.
      const n = {{ x: P.x / a, y: P.y / b }};
      const sols = intersectCircleLineUnit(n);
      if (!sols) return null;
      return sols.map((T) => ({{ x: a * T.x, y: b * T.y }}));
    }}

    function tangentPointsHyperbolaFromExternalPoint(P, a, b) {{
      // Scale: X=x/a, Y=y/b. Hyperbola: X^2 - Y^2 = 1.
      const X0 = P.x / a;
      const Y0 = P.y / b;

      if (Math.abs(Y0) < 1e-10) {{
        // X0*X = 1
        if (Math.abs(X0) < 1e-12) return null;
        const X = 1 / X0;
        const Yabs = Math.sqrt(Math.max(X * X - 1, 0));
        return [
          {{ x: a * X, y: b * Yabs }},
          {{ x: a * X, y: -b * Yabs }},
        ];
      }}

      const A = Y0 * Y0 - X0 * X0;
      const B = 2 * X0;
      const C = -(1 + Y0 * Y0);
      const disc = B * B - 4 * A * C;
      if (disc < 0) return null;
      const s = Math.sqrt(disc);
      if (Math.abs(A) < 1e-12) {{
        // Linear
        const X = -C / B;
        const Y = (X0 * X - 1) / Y0;
        return [{{ x: a * X, y: b * Y }}];
      }}
      const X1 = (-B + s) / (2 * A);
      const X2 = (-B - s) / (2 * A);
      const Y1 = (X0 * X1 - 1) / Y0;
      const Y2 = (X0 * X2 - 1) / Y0;
      return [
        {{ x: a * X1, y: b * Y1 }},
        {{ x: a * X2, y: b * Y2 }},
      ];
    }}

    function redraw() {{
      a = Number(aSlider.value);
      b = Number(bSlider.value);
      phi = Number(phiSlider.value);
      const type = typeSel.value;
      typeVal.textContent = '';

      // director circle radius
      let R = NaN;
      if (type === 'ellipse') R = Math.sqrt(a * a + b * b);
      else R = a > b ? Math.sqrt(a * a - b * b) : NaN;

      const S = Number.isFinite(R)
        ? {{ x: R * Math.cos(phi), y: R * Math.sin(phi) }}
        : {{ x: 0, y: 0 }};
      circleAt(pointS, S);

      // conic + director circle
      const allPts = [];
      if (type === 'ellipse') {{
        const pts = conicEllipsePts(a, b);
        curve1.setAttribute('d', pathFromWorld(pts));
        curve2.setAttribute('d', '');
        allPts.push(...pts);
      }} else {{
        const pts = conicHyperbolaPts(a, b);
        curve1.setAttribute('d', pathFromWorld(pts.ptsR));
        curve2.setAttribute('d', pathFromWorld(pts.ptsL));
        allPts.push(...pts.ptsR, ...pts.ptsL);
      }}

      if (Number.isFinite(R) && R > 0) {{
        circleDir.setAttribute('d', circlePath(R));
      }} else {{
        circleDir.setAttribute('d', '');
      }}

      // Tangents from S (only when director circle is real)
      let tangentsOk = false;
      let dot = NaN;
      let T1 = null, T2 = null;
      if (Number.isFinite(R) && R > 0) {{
        const tangPts =
          type === 'ellipse'
            ? tangentPointsEllipseFromExternalPoint(S, a, b)
            : tangentPointsHyperbolaFromExternalPoint(S, a, b);

        if (tangPts && tangPts.length >= 2) {{
          T1 = tangPts[0];
          T2 = tangPts[1];
          circleAt(t1, T1);
          circleAt(t2, T2);

          // Line coefficients in world coords:
          // ellipse: (x1/a^2)x + (y1/b^2)y = 1
          // hyperbola: (x1/a^2)x - (y1/b^2)y = 1
          const A1 = T1.x / (a * a);
          const B1 = (type === 'ellipse' ? T1.y : -T1.y) / (b * b);
          const A2 = T2.x / (a * a);
          const B2 = (type === 'ellipse' ? T2.y : -T2.y) / (b * b);

          const bounds = worldBounds([...allPts, S]);
          const segA = lineSegmentForAxByC(A1, B1, -1, {{
            xmin: bounds.xmin - 1,
            xmax: bounds.xmax + 1,
            ymin: bounds.ymin - 1,
            ymax: bounds.ymax + 1,
          }});
          const segB = lineSegmentForAxByC(A2, B2, -1, {{
            xmin: bounds.xmin - 1,
            xmax: bounds.xmax + 1,
            ymin: bounds.ymin - 1,
            ymax: bounds.ymax + 1,
          }});
          if (segA && segB) {{
            lineFromPoints(tang1, segA[0], segA[1]);
            lineFromPoints(tang2, segB[0], segB[1]);
            const d1 = {{ x: segA[1].x - segA[0].x, y: segA[1].y - segA[0].y }};
            const d2 = {{ x: segB[1].x - segB[0].x, y: segB[1].y - segB[0].y }};
            dot = (d1.x * d2.x + d1.y * d2.y) / (Math.hypot(d1.x, d1.y) * Math.hypot(d2.x, d2.y));
            tangentsOk = true;
          }}
        }}
      }}

      if (!tangentsOk) {{
        tang1.setAttribute('x1', 0); tang1.setAttribute('y1', 0); tang1.setAttribute('x2', 0); tang1.setAttribute('y2', 0);
        tang2.setAttribute('x1', 0); tang2.setAttribute('y1', 0); tang2.setAttribute('x2', 0); tang2.setAttribute('y2', 0);
        t1.setAttribute('cx', 0); t1.setAttribute('cy', 0);
        t2.setAttribute('cx', 0); t2.setAttribute('cy', 0);
      }}

      const bounds = worldBounds([...allPts, ...(Number.isFinite(R) ? [{{x:R,y:0}},{{x:-R,y:0}},{{x:0,y:R}},{{x:0,y:-R}}] : []), S]);
      setViewBox(svgRoot, bounds, 1.1);

      if (type === 'hyperbola' && !(a > b)) {{
        readout.innerHTML = `<div>双曲线导圆：需要 <span class="k">a &gt; b</span> 才是实圆；当前 a=${{fmt(a,2)}}, b=${{fmt(b,2)}}。</div>`;
      }} else {{
        const circleEq =
          type === 'ellipse'
            ? `x²+y² = a²+b² = ${{fmt(a*a + b*b, 3)}}`
            : `x²+y² = a²−b² = ${{fmt(a*a - b*b, 3)}}`;
        const perp = tangentsOk ? `两条切线方向余弦点积 ≈ ${{fmt(dot, 6)}}（≈0 表示垂直）` : '（当前选点无法作两条实切线）';
        readout.innerHTML = `
          <div>导圆方程：<span class="k">${{circleEq}}</span>，半径 R=${{Number.isFinite(R) ? fmt(R, 4) : '—'}}。</div>
          <div>${{perp}}</div>
        `;
      }}
    }}

    typeSel.addEventListener('change', redraw);
    redraw();
  }}

  function init() {{
    mountExampleA();
    mountParabola();
    mountHyperbola();
    mountDirectorCircle();
  }}

  if (document.readyState === 'loading') {{
    document.addEventListener('DOMContentLoaded', init);
  }} else {{
    init();
  }}
}})();
"""

def _render_html(md: str, title: str, interactive_js: str) -> str:
    body = markdown2.markdown(
        md,
        extras=[
            "fenced-code-blocks",
            "tables",
            "strike",
            "task_list",
            "header-ids",
        ],
    )
    body = _collapse_code_blocks(body)

    return f"""<!doctype html>
<html lang="zh">
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>{title}</title>
    <style>
      :root {{
        color-scheme: light dark;
        --bg: #ffffff;
        --fg: #111827;
        --muted: #6b7280;
        --card: #f8fafc;
        --border: #e5e7eb;
        --code-bg: #0b1020;
        --code-fg: #e5e7eb;
        --link: #2563eb;
      }}
      @media (prefers-color-scheme: dark) {{
        :root {{
          --bg: #0b1020;
          --fg: #e5e7eb;
          --muted: #9ca3af;
          --card: #111827;
          --border: #1f2937;
          --code-bg: #0b1020;
          --code-fg: #e5e7eb;
          --link: #60a5fa;
        }}
      }}
      body {{
        margin: 0;
        background: var(--bg);
        color: var(--fg);
        font: 16px/1.65 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
          Helvetica, Arial, "Noto Sans", "Apple Color Emoji", "Segoe UI Emoji";
      }}
      main {{
        max-width: 980px;
        margin: 0 auto;
        padding: 32px 18px 56px;
      }}
      h1, h2, h3 {{
        line-height: 1.25;
        letter-spacing: -0.01em;
      }}
      h1 {{ font-size: 2rem; margin: 0 0 0.6rem; }}
      h2 {{ font-size: 1.45rem; margin: 2.2rem 0 0.8rem; border-bottom: 1px solid var(--border); padding-bottom: .35rem; }}
      h3 {{ font-size: 1.15rem; margin: 1.6rem 0 0.6rem; }}
      p {{ margin: 0.75rem 0; }}
      a {{ color: var(--link); }}
      blockquote {{
        margin: 1rem 0;
        padding: 0.1rem 1rem;
        border-left: 4px solid var(--border);
        color: var(--muted);
        background: color-mix(in srgb, var(--card) 70%, transparent);
      }}
      hr {{
        border: 0;
        border-top: 1px solid var(--border);
        margin: 2rem 0;
      }}
      pre {{
        background: var(--code-bg);
        color: var(--code-fg);
        padding: 14px 14px;
        border-radius: 10px;
        overflow: auto;
      }}
      code {{
        font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
        font-size: 0.92em;
      }}
      pre code {{
        font-size: 0.9em;
      }}
      table {{
        width: 100%;
        border-collapse: collapse;
        margin: 1rem 0;
        font-size: 0.95em;
      }}
      th, td {{
        border: 1px solid var(--border);
        padding: 10px 10px;
        vertical-align: top;
      }}
      th {{
        background: color-mix(in srgb, var(--card) 80%, transparent);
        text-align: left;
      }}
      img {{
        max-width: 100%;
        height: auto;
        border-radius: 10px;
        border: 1px solid var(--border);
        background: var(--card);
      }}
      details.fallback {{
        margin: 0.8rem 0 1.2rem;
        border: 1px solid var(--border);
        border-radius: 10px;
        background: color-mix(in srgb, var(--card) 75%, transparent);
        padding: 0.3rem 0.6rem;
      }}
      details.fallback > summary {{
        cursor: pointer;
        font-weight: 600;
        color: var(--muted);
        padding: 0.35rem 0.2rem;
      }}
      details.codeblock {{
        margin: 1rem 0 1.2rem;
        border: 1px solid var(--border);
        border-radius: 10px;
        background: color-mix(in srgb, var(--card) 75%, transparent);
      }}
      details.codeblock > summary {{
        cursor: pointer;
        font-weight: 650;
        padding: 10px 12px;
        color: var(--muted);
        user-select: none;
      }}
      details.codeblock .codehilite {{
        margin: 0;
        border-top: 1px solid var(--border);
      }}
      details.codeblock pre {{
        margin: 0;
        border-radius: 0 0 10px 10px;
      }}
      .widget {{
        border: 1px solid var(--border);
        border-radius: 12px;
        background: color-mix(in srgb, var(--card) 82%, transparent);
        padding: 12px 12px 10px;
        margin: 0.8rem 0 1.2rem;
      }}
      .widget-title {{
        font-weight: 700;
        margin: 2px 0 10px;
        letter-spacing: -0.01em;
      }}
      .widget-controls {{
        display: grid;
        grid-template-columns: 1fr;
        gap: 8px;
        margin: 0 0 10px;
      }}
      .control-row {{
        display: grid;
        grid-template-columns: 120px 1fr 160px;
        gap: 10px;
        align-items: center;
      }}
      .control-label {{
        color: var(--muted);
        font-weight: 600;
      }}
      .control-value {{
        text-align: right;
        font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
        font-size: 0.92em;
        color: var(--muted);
      }}
      .control-select {{
        width: 100%;
        padding: 6px 8px;
        border-radius: 10px;
        border: 1px solid var(--border);
        background: var(--bg);
        color: var(--fg);
      }}
      .conic-svg {{
        width: 100%;
        height: 420px;
        display: block;
        border-radius: 12px;
        border: 1px solid var(--border);
        background: var(--bg);
      }}
      .widget-readout {{
        margin: 10px 2px 2px;
        color: var(--muted);
        font-size: 0.95em;
      }}
      .widget-readout .k {{
        font-weight: 700;
        color: var(--fg);
      }}
      .curve {{
        fill: none;
        stroke: color-mix(in srgb, var(--link) 85%, var(--fg));
        stroke-width: 1.6;
        vector-effect: non-scaling-stroke;
      }}
      .seg {{
        stroke: color-mix(in srgb, #f59e0b 88%, var(--fg));
        stroke-width: 1.4;
        vector-effect: non-scaling-stroke;
        opacity: 0.9;
      }}
      .seg2 {{
        stroke: color-mix(in srgb, #10b981 88%, var(--fg));
        stroke-width: 1.4;
        vector-effect: non-scaling-stroke;
        opacity: 0.9;
      }}
      .focus {{
        fill: color-mix(in srgb, #ef4444 88%, var(--fg));
        stroke: var(--bg);
        stroke-width: 0.08;
      }}
      .point {{
        fill: color-mix(in srgb, #2563eb 92%, var(--fg));
        stroke: var(--bg);
        stroke-width: 0.08;
      }}
      .directrix {{
        stroke: color-mix(in srgb, var(--muted) 85%, var(--fg));
        stroke-width: 1.2;
        vector-effect: non-scaling-stroke;
        stroke-dasharray: 6 6;
      }}
      .asym {{
        stroke: color-mix(in srgb, var(--muted) 85%, var(--fg));
        stroke-width: 1.1;
        vector-effect: non-scaling-stroke;
        stroke-dasharray: 6 6;
      }}
      .dircircle {{
        fill: none;
        stroke: color-mix(in srgb, #a855f7 82%, var(--fg));
        stroke-width: 1.4;
        vector-effect: non-scaling-stroke;
        opacity: 0.9;
      }}
      .tangent {{
        stroke: color-mix(in srgb, #22c55e 82%, var(--fg));
        stroke-width: 1.4;
        vector-effect: non-scaling-stroke;
        opacity: 0.9;
      }}
      .meta {{
        color: var(--muted);
        font-size: 0.95em;
        margin: 0 0 1.4rem;
      }}
    </style>
    <script>
      window.MathJax = {{
        tex: {{
          inlineMath: [['\\\\(', '\\\\)'], ['$', '$']],
          displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
          processEscapes: true
        }},
        options: {{
          skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code']
        }}
      }};
    </script>
    <script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
  </head>
  <body>
    <main>
      <div class="meta">Generated by <code>scripts/generate_conics_report.py</code></div>
      {body}
    </main>
    <script>
{interactive_js}
    </script>
  </body>
</html>
"""


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate HTML report + plots from AGENTS.md.")
    parser.add_argument("--agents", default="AGENTS.md", help="Path to AGENTS.md")
    parser.add_argument("--out-dir", default="dist", help="Output directory")
    parser.add_argument("--title", default="圆锥曲线守恒式与关系式（图文报告）", help="HTML title")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    agents_path = (repo_root / args.agents).resolve()
    out_dir = (repo_root / args.out_dir).resolve()
    assets_dir = out_dir / "assets"
    assets_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots + computed results
    exa = _plot_example_a(assets_dir / "example_a.png")
    exb = _plot_example_b(assets_dir / "example_b.png")
    hyp = _plot_hyperbola(assets_dir / "hyperbola.png")

    agents_md = agents_path.read_text(encoding="utf-8")
    processed_md = _preprocess_agents_markdown(agents_md)
    processed_md += _build_generated_section(exa, exb, hyp)
    processed_md += _build_extra_formula_table_markdown()
    processed_md += _build_polar_section_markdown()
    processed_md = _escape_markdown_emphasis_in_math(processed_md)

    html = _render_html(processed_md, title=args.title, interactive_js=_interactive_js(exa))
    out_dir.mkdir(parents=True, exist_ok=True)
    out_html = out_dir / "conics_report.html"
    out_html.write_text(html, encoding="utf-8")

    print(f"Wrote: {out_html}")
    print(f"Assets: {assets_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
