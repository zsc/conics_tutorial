下面把**圆锥曲线（椭圆/双曲线/抛物线）里“守恒式（更准确叫：在平移/旋转等坐标变换下不变的量、以及由定义带来的不变量）”**和**常用/冷门数学关系式**做一个“可查表式”的总结，并给出**解析几何算例**与 **Python 画图代码**（可直接运行）。

---

## 1) 最常用的“守恒式 / 不变量”（解析几何最爱用的那批）

### 1.1 二次曲线一般式与矩阵形式（通用入口）

一般二次曲线：
[
Ax^2+Bxy+Cy^2+Dx+Ey+F=0
]
矩阵写法（便于“守恒量”一眼看出）：
[
\mathbf{x}^TQ\mathbf{x}+ \mathbf{p}^T\mathbf{x}+F=0,\quad
Q=\begin{pmatrix}A&\frac B2\ \frac B2&C\end{pmatrix},\
\mathbf{p}=\begin{pmatrix}D\E\end{pmatrix}
]

### 1.2 分类判别式（旋转不变：最常用）

[
\Delta_2 = B^2-4AC
]

* (\Delta_2<0)：椭圆型（可能是椭圆、圆、点、虚椭圆）
* (\Delta_2=0)：抛物线型（可能退化成平行线/重线）
* (\Delta_2>0)：双曲线型（可能退化成两相交直线）

> 这是最常被称作“守恒”的：**无论你怎么旋转坐标轴，(B^2-4AC) 不变**（本质是 (Q) 的特征值符号结构不变）。

### 1.3 消去交叉项的旋转角（常用）

当 (B\neq 0)，旋转角 (\theta) 满足
[
\tan 2\theta=\frac{B}{A-C}
]
旋转后 (xy) 项消失。

### 1.4 二次型的特征值（更“本质”的不变量）

(Q) 的特征值 (\lambda_1,\lambda_2)（以及它们的符号）在正交变换下不变。

* (\text{tr}(Q)=A+C=\lambda_1+\lambda_2)
* (\det(Q)=AC-\frac{B^2}{4}=\lambda_1\lambda_2)

### 1.5 中心存在性（常用但容易忘）

二次曲线有“中心”（椭圆/双曲线及其退化）当且仅当线性方程组有解：
[
\begin{pmatrix}2A&B\B&2C\end{pmatrix}\begin{pmatrix}x_0\y_0\end{pmatrix}=
-\begin{pmatrix}D\E\end{pmatrix}
]
抛物线一般**没有中心**。

---

## 2) 三大圆锥曲线的“定义型守恒式”（距离和/差、离心率）

### 2.1 离心率统一定义（超级常用）

焦点—准线定义：
[
\frac{\text{点到焦点距离}}{\text{点到准线距离}}=e
]

* 椭圆：(0<e<1)
* 抛物线：(e=1)
* 双曲线：(e>1)

### 2.2 椭圆（最常用关系式）

标准式（主轴沿 x）：
[
\frac{x^2}{a^2}+\frac{y^2}{b^2}=1,\quad a\ge b>0
]
参数式（常用）：
[
x=a\cos t,\quad y=b\sin t,\quad 0\le t<2\pi
]
焦距参数：
[
c^2=a^2-b^2,\quad e=\frac ca
]
**距离和守恒（定义）**：对任一点 (P)，
[
PF_1+PF_2=2a
]
准线：
[
x=\pm\frac{a}{e}=\pm\frac{a^2}{c}
]
半通径（latus rectum）：
[
\ell=\frac{b^2}{a}
]
面积（常用）：
[
S=\pi ab
]
祖暅原理（Cavalieri/“截面积守恒”思想）的一个推导思路（平面版）：

把椭圆域
[
E:\ \frac{x^2}{a^2}+\frac{y^2}{b^2}\le 1
]
与单位圆盘
[
D:\ x^2+y^2\le 1
]
作比较。令 (y=bt)（(t\in[-1,1])），则椭圆在高度 (y) 处的水平截线长度为
[
L_E(y)=2a\sqrt{1-\frac{y^2}{b^2}}=2a\sqrt{1-t^2}
]
而单位圆盘在高度 (t) 处的截线长度为
[
L_D(t)=2\sqrt{1-t^2}
]
于是对应“每一层截线长度成固定比例 (a)”且“高度参数再按 (b) 缩放”，面积按 (ab) 成比例：
[
S(E)=ab\,S(D)=ab\cdot\pi=\pi ab
]
切线（点 ((x_0,y_0)) 在椭圆上）：
[
\frac{xx_0}{a^2}+\frac{yy_0}{b^2}=1
]

### 2.3 双曲线（最常用关系式）

标准式：
[
\frac{x^2}{a^2}-\frac{y^2}{b^2}=1
]
参数式（常用）：
[
x=a\sec t,\quad y=b\tan t,\quad |t|<\frac{\pi}{2}
]
或（双曲函数参数）：
[
x=a\cosh u,\quad y=b\sinh u,\quad u\in\mathbb{R}
]
[
c^2=a^2+b^2,\quad e=\frac ca>1
]
**距离差守恒（定义）**：
[
|PF_1-PF_2|=2a
]
渐近线：
[
y=\pm\frac ba x
]
准线：
[
x=\pm\frac{a}{e}=\pm\frac{a^2}{c}
]
半通径：
[
\ell=\frac{b^2}{a}
]
切线（点 ((x_0,y_0)) 在双曲线上）：
[
\frac{xx_0}{a^2}-\frac{yy_0}{b^2}=1
]

### 2.4 抛物线（最常用关系式）

标准式：
[
y^2=4px\quad (p>0)
]
焦点、准线：
[
F(p,0),\quad x=-p
]
焦弦（过焦点的弦）长度常用结论：

* **通径（latus rectum）长**：(4p)
  参数式（常用）：
  [
  x=pt^2,\quad y=2pt
  ]
  切线（参数 (t_0)）：
  [
  ty = x + pt^2\quad\text{在 }t=t_0\text{ 处 } \Rightarrow\ t_0y=x+pt_0^2
  ]

---

## 3) 一些“冷门但很好用/很漂亮”的关系式（挑常见冷门）

### 3.1 共轭直径/共轭弦（椭圆/双曲线）

对椭圆参数点 (P(a\cos t, b\sin t))，其切线方向与“共轭直径”有经典对应（工程/几何优化里挺常用）。

### 3.2 反射性质（焦点反射）

* 椭圆：从一个焦点出发反射到另一焦点（入射角=反射角）
* 抛物线：来自焦点的光线反射后平行于轴（以及反过来）
* 双曲线：来自一个焦点反射后好像来自另一个焦点（“虚像”性质）

### 3.3 Director circle（导圆，较冷门）

* 椭圆 (\frac{x^2}{a^2}+\frac{y^2}{b^2}=1)：导圆 (x^2+y^2=a^2+b^2)
* 双曲线 (\frac{x^2}{a^2}-\frac{y^2}{b^2}=1)：导圆 (x^2+y^2=a^2-b^2)（若 (a>b) 才为实圆）

导圆：满足“从该圆上点引曲线切线互相垂直”等等等价性质（竞赛/解析几何题偶尔出现）。

### 3.4 极线/极点（更偏射影几何，但解析几何也能算）

对二次曲线 (X^TAX=0)，点 (P) 的极线可用 (P^TA X=0) 表达（坐标齐次形式）。对求切线、判定共线很强。

---

## 4) 解析几何算例（带完整步骤思路）

### 算例 A：把一般二次曲线化到标准形并提取几何量

给：
[
9x^2+4xy+6y^2-54x-24y+36=0
]
**(1) 分类：**
(\Delta_2=B^2-4AC=4^2-4\cdot 9\cdot 6=16-216=-200<0)
所以是**椭圆型**（或退化椭圆）。

**(2) 旋转去掉 (xy)：**
[
\tan 2\theta=\frac{B}{A-C}=\frac{4}{9-6}=\frac{4}{3}
]
取 (\theta) 满足即可（Python 里可直接算数值）。

**(3) 平移到中心：**
解
[
\begin{pmatrix}18&4\4&12\end{pmatrix}\begin{pmatrix}x_0\y_0\end{pmatrix}=
\begin{pmatrix}54\24\end{pmatrix}
]
得到中心 ((x_0,y_0))（Python 一行线代就能解）。

**(4) 得到标准式参数 (a,b) 与主轴方向：**
用特征值分解得到旋转后的二次项系数，再结合常数项得到 (a^2,b^2)。
下面 Python 会把这些全部算出来并画图（含中心/主轴）。

---

### 算例 B：抛物线切线、焦点、准线与“守恒式”验证

取 (y^2=8x)，即 (4p=8\Rightarrow p=2)。

* 焦点：(F(2,0))
* 准线：(x=-2)
* 任一点 (P(x,y)) 在抛物线上满足
  [
  PF=\text{到准线距离} \quad(\text{离心率 }e=1)
  ]
  再取参数点 (t=1)：(P(pt^2,2pt)=(2,4))。
  可直接数值验证 (PF) 与 (x) 到 (-2) 的距离相等。
  切线：(t y = x + pt^2\Rightarrow 1\cdot y=x+2)。

---

## 5) Python 画图与自动化提取参数（可直接运行）

下面代码会做三件事：

1. 对算例 A 自动：分类、求中心、求旋转角、化标准式参数，并画出椭圆
2. 对算例 B 画抛物线、焦点、准线、某点切线，并验证 (PF=)到准线距离
3. 顺带画一个标准双曲线示例（带渐近线、焦点）

```python
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Helper: classify conic by invariants
# -----------------------------
def classify_conic(A,B,C,D,E,F):
    disc = B*B - 4*A*C
    if disc < 0:
        kind = "ellipse-type (disc<0)"
    elif disc == 0:
        kind = "parabola-type (disc=0)"
    else:
        kind = "hyperbola-type (disc>0)"
    return disc, kind

# -----------------------------
# Example A: 9x^2 + 4xy + 6y^2 -54x -24y +36 = 0
# Use matrix form x^T Q x + p^T x + F = 0
# -----------------------------
A,B,C,D,E,F0 = 9,4,6,-54,-24,36
disc, kind = classify_conic(A,B,C,D,E,F0)
print("Example A discriminant B^2-4AC =", disc, "=>", kind)

Q = np.array([[A, B/2],
              [B/2, C]], dtype=float)
p = np.array([D, E], dtype=float)

# Center exists if [2A B; B 2C] is invertible and solution exists
M = np.array([[2*A, B],
              [B, 2*C]], dtype=float)
center = np.linalg.solve(M, -p)
x0, y0 = center
print("Example A center (x0,y0) =", center)

# Rotate to diagonalize Q: eigen-decomposition gives principal axes
eigvals, eigvecs = np.linalg.eigh(Q)  # Q = V diag(lam) V^T
# In rotated coords u = V^T (x - center)
lam1, lam2 = eigvals
V = eigvecs

# Compute constant after completing square: (x-center)^T Q (x-center) = -F'
# Original: x^T Q x + p^T x + F0 = 0
# At center, gradient zero, so value at center gives F'
F_at_center = (center @ Q @ center) + (p @ center) + F0
# In shifted coords z = x-center, equation: z^T Q z + F_at_center = 0
# Rotate: u^T diag(lam) u + F_at_center = 0
print("Example A value at center (should be <=0 for real ellipse):", F_at_center)

# For ellipse: lam1 u1^2 + lam2 u2^2 = -F_at_center
rhs = -F_at_center
a2 = rhs/lam1
b2 = rhs/lam2
# a is sqrt(max), b is sqrt(min) in geometric sense:
a = np.sqrt(max(a2,b2))
b = np.sqrt(min(a2,b2))
print("Example A semi-axes (a,b) =", (a,b))

# Generate ellipse points in principal coordinates, then map back
t = np.linspace(0, 2*np.pi, 500)
u1 = np.sqrt(a2)*np.cos(t)
u2 = np.sqrt(b2)*np.sin(t)
U = np.vstack([u1, u2])              # 2xN
X = (V @ U).T + center               # Nx2 in original coordinates

# Plot Example A
plt.figure()
plt.plot(X[:,0], X[:,1], label="Conic A (ellipse)")
plt.scatter([x0],[y0], s=40, label="center")

# Plot principal axes directions
# axes directions are columns of V, scale by semi-axes
axis1 = center + V[:,0]*np.sqrt(a2)
axis2 = center + V[:,1]*np.sqrt(b2)
plt.plot([center[0], axis1[0]], [center[1], axis1[1]], linewidth=2, label="principal axis 1")
plt.plot([center[0], axis2[0]], [center[1], axis2[1]], linewidth=2, label="principal axis 2")

plt.gca().set_aspect("equal", adjustable="box")
plt.title("Example A: general conic -> ellipse (center + principal axes)")
plt.grid(True)
plt.legend()

# -----------------------------
# Example B: Parabola y^2 = 8x => p=2
# -----------------------------
p_par = 2
# param t
tt = np.linspace(-4,4,400)
x_par = p_par*tt**2
y_par = 2*p_par*tt

# choose point at t0=1
t0 = 1
P = np.array([p_par*t0**2, 2*p_par*t0])  # (2,4)
F_focus = np.array([p_par, 0.0])         # (2,0)
# directrix x=-p
directrix_x = -p_par

# verify PF = distance to directrix
PF = np.linalg.norm(P - F_focus)
dist_directrix = abs(P[0] - directrix_x)
print("\nExample B check e=1 property at P:", P)
print("PF =", PF, " distance to directrix =", dist_directrix)

# tangent line at t0: t0*y = x + p*t0^2
# => y = (x + p*t0^2)/t0
# here t0=1: y = x + 2
x_line = np.linspace(-3, 10, 200)
y_tan = (x_line + p_par*(t0**2))/t0

plt.figure()
plt.plot(x_par, y_par, label="parabola y^2=8x")
plt.scatter([F_focus[0]],[F_focus[1]], s=40, label="focus (2,0)")
plt.axvline(directrix_x, linestyle="--", label="directrix x=-2")
plt.scatter([P[0]],[P[1]], s=40, label="P at t=1 (2,4)")
plt.plot(x_line, y_tan, label="tangent at P: y=x+2")

plt.gca().set_aspect("equal", adjustable="box")
plt.title("Example B: parabola with focus/directrix/tangent")
plt.grid(True)
plt.legend()

# -----------------------------
# Bonus: standard hyperbola x^2/a^2 - y^2/b^2 = 1
# -----------------------------
a_h, b_h = 3, 2
x = np.linspace(-10, -a_h, 400)
y = b_h*np.sqrt(x**2/a_h**2 - 1)
x2 = np.linspace(a_h, 10, 400)
y2 = b_h*np.sqrt(x2**2/a_h**2 - 1)
# asymptotes y = ±(b/a)x
x_as = np.linspace(-10, 10, 2)
y_as = (b_h/a_h)*x_as

plt.figure()
plt.plot(x, y, label="hyperbola branch")
plt.plot(x, -y)
plt.plot(x2, y2)
plt.plot(x2, -y2)
plt.plot(x_as, y_as, linestyle="--", label="asymptote")
plt.plot(x_as, -y_as, linestyle="--")
# foci
c_h = np.sqrt(a_h**2 + b_h**2)
plt.scatter([c_h, -c_h], [0,0], s=40, label="foci")
plt.gca().set_aspect("equal", adjustable="box")
plt.title("Bonus: hyperbola with asymptotes and foci")
plt.grid(True)
plt.legend()

plt.show()
```

---

如果你希望我再补一份“更冷门但可用于解题的公式清单”（比如：极线极点的坐标算法、圆锥曲线的弦中点公式、过定点弦族的参数化、二次曲线的仿射不变量、圆锥曲线与直线交点的判别式与韦达关系、以及一批竞赛常见结论），我也可以按“**公式 → 适用场景 → 一行例题**”的格式继续整理。
需要。
