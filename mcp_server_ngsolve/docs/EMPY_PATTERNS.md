# NGSolve Patterns from EMPY_Analysis

このドキュメントは、EMPY_Analysisレポジトリから抽出したNGSolveの実装パターンを記録しています。

## 概要

EMPY_Analysis (`S:\NGSolve\EMPY\EMPY_Analysis`) は電磁場解析のための実装集で、以下の重要なパターンを含んでいます:

- **T-Omega法**: 渦電流問題のための電流密度-磁気スカラーポテンシャル定式化
- **A-Phi法**: ベクトルポテンシャル-電気スカラーポテンシャル定式化
- **Loop Field計算**: 多重連結領域のトポロジー処理
- **カスタムソルバ統合**: ICCG solver with C++バックエンド

## 1. T-Omega Method (渦電流問題)

### 定式化

**変数:**
- T: 電流密度ベクトル (current density vector)
- Ω: 磁気スカラーポテンシャル (magnetic scalar potential)

**支配方程式:**
```
curl(1/σ curl T) + sμT + sμ grad(Ω) = 0  (in conductor)
div(μ(T + grad(Ω))) = 0                   (in air)
```

### NGSolve実装パターン

```python
from ngsolve import *

# Mixed finite element space
fesT = HCurl(mesh, order=2, nograds=True, definedon="conductor",
             dirichlet="conductorBND", complex=False)
fesOmega = H1(mesh, order=2, dirichlet="upper|lower", complex=False)
fespace = fesT * fesOmega

(T, Omega), (W, psi) = fespace.TnT()

# Bilinear form
a = BilinearForm(fespace)
# Conductor region
a += (1/sigma) * curl(T) * curl(W) * dx("conductor")
a += s * mu * T * W * dx("conductor")
a += s * mu * T * grad(psi) * dx("conductor")
# Air region
a += s * mu * grad(Omega) * grad(psi) * dx("air")
```

### 重要な設定

**1. HCurl space with nograds=True:**
```python
fesT = HCurl(mesh, order=2, nograds=True, definedon="conductor")
```
- `nograds=True`: Tが curl-free成分を持たないことを保証
- 電流密度は div T = 0 を満たす必要がある

**2. Dirichlet境界条件:**
```python
fesT = HCurl(..., dirichlet="conductorBND")        # Tangential T = 0
fesOmega = H1(..., dirichlet="upper|lower")        # Ω = 定数
```

### Loop Field処理 (多重連結領域)

**Genus計算とLoop Field生成:**

```python
def LoopFields(mesh, domain, connected=1):
    """多重連結領域のloop fieldを計算"""
    # Surface mesh from boundary
    smesh = surface_mesh_from_boundary(mesh, "conductorBND")

    # Genus (ループの数) を計算
    # genus g = (2 - χ)/2 where χ is Euler characteristic
    g = surface_genus(smesh, connected)

    fes = HCurl(mesh, order=1, nograds=True, definedon=domain)
    u, v = fes.TnT()

    loops = []
    for k in range(g):
        # ランダムなエッジを選択
        gfu = GridFunction(fes)
        edge_dofs = fes.GetDofNrs(random_edge(smesh))[0]
        gfu.vec[edge_dofs] = 1
        fes.FreeDofs()[edge_dofs] = False

        # curl T = 0 を解く
        a = BilinearForm(fes)
        a += curl(u) * curl(v) * dx
        a.Assemble()

        fr = -a.mat * gfu.vec
        gfu = ICCG_Solve(fes, gfu, a, fr)

        # Grad部分を除去
        fesPhi = H1(mesh, order=1, definedon="air")
        phi, psi = fesPhi.TnT()

        a = BilinearForm(fesPhi)
        a += grad(phi) * grad(psi) * dx
        f = LinearForm(fesPhi)
        f += grad(psi) * gfu * dx
        a.Assemble()
        f.Assemble()

        gfPhi = ICCG_Solve(fesPhi, GridFunction(fesPhi), a, f.vec)
        gfw = gfu - grad(gfPhi)

        # 既存ループとの直交化 (Gram-Schmidt)
        gft = gfw
        for kd in range(len(loops)):
            prod = Integrate(gfw * loops[kd] * dx, mesh)
            gft.vec -= prod * loops[kd].vec

        # 正規化
        norm = sqrt(Integrate(gft * gft * dx, mesh))
        gft.vec /= norm

        loops.append(gft)

    return loops
```

**Loop電流との連成:**

```python
def GetLoopCouplings(loopFields, fesTOmega, boundary):
    """Loop fieldとの連成項を計算"""
    g = len(loopFields)
    fs = []      # 右辺ベクトル
    fafs = []    # 連成行列要素

    (T, omega), (W, psi) = fesTOmega.TnT()

    for k in range(g):
        loopField = loopFields[k]
        gfTOmega = GridFunction(fesTOmega)
        SetBoundaryValue(loopField, 1, gfTOmega, boundary)
        gfT, gfOmega = gfTOmega.components

        # 右辺ベクトル
        f = LinearForm(fesTOmega)
        f += (1/sigma) * curl(gfT) * curl(W) * dx("conductor")
        f += s * mu * gfT * (W + grad(psi)) * dx("conductor")
        f += s * mu * loopField * grad(psi) * dx("air")
        f.Assemble()
        fs.append(f)

        # 連成行列 (k×k symmetric matrix)
        tmp = []
        for k2 in range(k+1):
            loopField2 = loopFields[k2]
            gfT2, _ = gfTs[k2].components

            faf = Integrate((1/sigma) * curl(gfT) * curl(gfT2) * dx("conductor"), mesh)
            faf += Integrate(s * mu * gfT * gfT2 * dx("conductor"), mesh)
            faf += Integrate(s * mu * loopField * loopField2 * dx("air"), mesh)

            tmp.append(faf)
            if k2 < k: fafs[k2].append(faf)
        fafs.append(tmp)

    return fs, fafs
```

**連成系の求解:**

```python
# システム行列 + loop連成
# [A   C^T] [x]   [f]
# [C   B  ] [I] = [V]
# where A: システム行列, C: loop連成ベクトル, B: loop self-term, I: loop電流, V: 電圧

x, current = SolveCoupled(fesTOmega, systemMatrix, loopCouplings, loopSelfTerms,
                          sourceTerms, voltages)

# x: DOF vector
# current: loop currents [I1, I2, ..., Ig]
```

## 2. A-Phi Method

### H-field Tangential Regularization

境界でのH-fieldの接線成分を正則化する手法:

```python
def Ht_regularization(H, mesh, boundary, feOrder):
    """Tangential H-field regularization at boundary"""
    normal = specialcf.normal(mesh.dim)

    # H1 space on boundary
    fesu = H1(mesh, order=feOrder, definedon=mesh.Boundaries(boundary), complex=False)
    u, v = fesu.TnT()

    # Minimize ||n × ∇u - H_t||^2
    a = BilinearForm(fesu)
    a += Cross(normal, grad(u).Trace()) * Cross(normal, grad(v).Trace()) * ds

    f = LinearForm(fesu)
    f += -H * Cross(normal, grad(v).Trace()) * ds

    a.Assemble()
    f.Assemble()

    gfu = GridFunction(fesu)
    gfu = ICCG_Solve(fesu, gfu, a, f.vec)

    # Regularized field
    hreg = H + Cross(normal, grad(gfu).Trace())

    return hreg
```

**使用場面:**
- 導体境界でのH-fieldの接線連続性を改善
- 不連続なHフィールドの正則化
- ポスト処理での精度向上

## 3. Custom Matrix Solver Integration

### ICCG Solver with scipy + C++ backend

```python
import scipy.sparse as sp
import EMPY_Solver  # C++ pybind11 module

def iccg_solve(fes, gf, A, Bvec, tol=1e-10, max_iter=1000, accel_factor=1.1):
    """ICCG solver with C++ acceleration"""
    # NGSolve matrix → scipy CSR
    asci = sp.csr_matrix(A.mat.CSR())
    Acut = asci[:, fes.FreeDofs()][fes.FreeDofs(), :]
    fcut = np.array(Bvec)[fes.FreeDofs()]
    ucut = np.zeros_like(fcut)

    # Extract sparse matrix data
    rows, cols = Acut.nonzero()
    vals = np.ravel(Acut[rows, cols])
    dim = fcut.size

    # C++ solver
    solver = EMPY_Solver.EMPY_Solver()
    solver.SetMatrix(dim, len(rows), rows, cols, vals)
    solver.SetScaling(True)
    solver.SetEps(tol)
    solver.SetShiftParameter(accel_factor)
    solver.SetDivCriterion(10.0, 10)

    ucut = solver.Solve(fcut, ucut)

    # Get convergence info
    log = solver.GetResidualLog()
    shift = solver.GetShiftParameter()

    # Update GridFunction
    np.array(gf.vec.FV(), copy=False)[fes.FreeDofs()] += ucut

    # Verify solution
    result = Acut.dot(ucut) - fcut
    norm = np.linalg.norm(result) / np.linalg.norm(fcut)
    print(f"Residual norm: {norm}")

    return gf
```

### Matrix Coupling for Loop Currents

```python
def AddCoupling(matrix, cvecs, amat):
    """Add coupling terms for loop currents

    [A   C^T]
    [C   B  ]

    where:
    - A: original matrix (dim × dim)
    - C: coupling vectors (dim × nadd)
    - B: self-term matrix (nadd × nadd)
    """
    dim = matrix.shape[0]
    nadd = len(cvecs)

    # Allocate extended sparse matrix
    rows, cols = matrix.nonzero()
    vals = np.ravel(matrix[rows, cols])
    size = vals.size

    # Estimate new size
    non_zeros = max(np.count_nonzero(cv) for cv in cvecs)
    sizep = size + non_zeros * 2 * nadd + nadd * nadd

    new_rows = np.zeros(sizep, dtype=int)
    new_cols = np.zeros(sizep, dtype=int)
    new_vals = np.zeros(sizep)

    # Copy original matrix
    k = size
    new_rows[:k] = rows
    new_cols[:k] = cols
    new_vals[:k] = vals

    # Add coupling terms C and C^T
    for n in range(nadd):
        fcut = cvecs[n]
        r = dim + n
        for col in range(dim):
            v = fcut[col]
            if v != 0:
                # C^T[n, col]
                new_rows[k] = col
                new_cols[k] = r
                new_vals[k] = v
                k += 1
                # C[col, n]
                new_rows[k] = r
                new_cols[k] = col
                new_vals[k] = v
                k += 1

    # Add self-term matrix B
    for n in range(nadd):
        for n2 in range(nadd):
            new_rows[k] = dim + n
            new_cols[k] = dim + n2
            new_vals[k] = amat[n][n2]
            k += 1

    # Create extended sparse matrix
    new_a = sp.csr_matrix((new_vals, (new_rows, new_cols)),
                          shape=(dim + nadd, dim + nadd))

    return new_a
```

## 4. 実装のベストプラクティス

### 1. Finite Element Space選択

| 問題 | 空間 | 設定 |
|------|------|------|
| 電流密度 T | HCurl | nograds=True, dirichlet=境界 |
| 磁気ポテンシャル Ω | H1 | dirichlet=上下境界 |
| ベクトルポテンシャル A | HCurl | dirichlet=境界 |
| 電気ポテンシャル Φ | H1 | dirichlet=境界 |

### 2. 境界条件設定

**Dirichlet:**
```python
# Tangential component = 0
fesT = HCurl(mesh, dirichlet="conductorBND")

# Normal component = 0
fesOmega = H1(mesh, dirichlet="upper|lower")
```

**Boundary Value設定:**
```python
def SetBoundaryValue(source_gf, factor, target_gf, boundaryId):
    """Copy boundary values from source to target"""
    mesh = source_gf.space.mesh
    for t in mesh.Boundaries(boundaryId).Elements():
        for e in t.edges:
            k_src = source_gf.space.GetDofNrs(e)
            k_tgt = target_gf.space.GetDofNrs(e)
            target_gf.vec[k_tgt[0]] = source_gf.vec[k_src[0]] * factor
```

### 3. Multi-connected Domain処理

**手順:**
1. Genus計算: `g = surface_genus(surface_mesh, connected)`
2. Loop field生成: `g`個の独立なloop fields
3. Gram-Schmidt直交化
4. 正規化: `∫ loop_i · loop_j dx = δ_ij`
5. 連成行列構築
6. 拡張システム求解

### 4. ソルバパラメータ調整

```python
# ICCG solver
solver.SetEps(1e-10)                      # 収束判定
solver.SetShiftParameter(1.1)             # 加速係数 (自動調整)
solver.SetDivCriterion(10.0, 10)         # 発散判定 (倍率, 回数)
solver.SetScaling(True)                   # 対角スケーリング

# 推奨設定:
# - tol: 1e-10 ~ 1e-12 (高精度)
# - accel_factor: 1.0 ~ 1.3 (自動調整推奨)
# - max_iter: 1000 ~ 2000
```

### 5. メモリ管理

```python
# 大きな行列は使用後すぐに解放
asci = sp.csr_matrix(A.mat.CSR())
Acut = asci[:, fes.FreeDofs()][fes.FreeDofs(), :]
asci = None  # 明示的に解放

# NumPy配列も同様
rows = None
cols = None
vals = None
```

## 5. パフォーマンス最適化

### Sparse Matrix操作

```python
# ✓ 良い: CSR formatで直接操作
rows, cols = matrix.nonzero()
vals = matrix[rows, cols]

# ✗ 悪い: dense変換
dense_matrix = matrix.toarray()  # メモリ爆発
```

### TaskManager活用

```python
with TaskManager():
    a.Assemble()  # 並列アセンブリ
    f.Assemble()
```

### FreeDofs indexing

```python
# ✓ 効率的: Boolean indexing
Acut = A[:, fes.FreeDofs()][fes.FreeDofs(), :]
fcut = f[fes.FreeDofs()]

# ✗ 非効率: Loop
for i in range(dim):
    if fes.FreeDofs()[i]:
        # ...
```

## 6. MCP Server統合の可能性

以下の機能をMCPツールとして実装できる:

1. **`ngsolve_eddy_current_t_omega`** - T-Omega法の自動セットアップ
2. **`ngsolve_compute_loop_fields`** - Loop field自動計算
3. **`ngsolve_h_field_regularization`** - 境界でのH-field正則化
4. **`ngsolve_coupled_system_solve`** - Loop電流連成系ソルバ
5. **`ngsolve_topology_analysis`** - メッシュトポロジー解析 (genus計算)

## 7. 2スカラー法 (Reduced-Total Scalar Potential Method)

### 概要

2スカラー法は、永久磁石や電流源を含む静磁場問題を効率的に解くための手法です。領域を2つに分け、それぞれ異なるスカラーポテンシャルを使用します。

### 定式化

**領域分割:**
- **Ωₛ (source region)**: 電流源または永久磁石を含む領域
- **Ωₙ (non-source region)**: 電流源のない領域（空気など）

**ポテンシャル:**
- **Ωᵣ (reduced scalar potential)** in Ωₛ: H = H₀ - ∇Ωᵣ
- **Ωₜ (total scalar potential)** in Ωₙ: H = -∇Ωₜ

ここで、H₀は既知の源磁場（コイルや永久磁石による）

### 支配方程式

**Source region (Ωₛ):**
```
∇·μ(H₀ - ∇Ωᵣ) = 0
```

**Non-source region (Ωₙ):**
```
∇·μ(-∇Ωₜ) = 0
```

**Interface boundary (Γ):**
```
μ₁(H₀ - ∇Ωᵣ)·n = μ₂(-∇Ωₜ)·n  (normal component continuity)
Ωᵣ = Ωₜ  (tangential component continuity)
```

### NGSolve実装パターン

```python
from ngsolve import *

# Finite element spaces
# Ωᵣ in source region
fesR = H1(mesh, order=2, definedon="source")
# Ωₜ in non-source region
fesT = H1(mesh, order=2, definedon="air")
# Mixed space
fespace = fesR * fesT

(Omega_r, Omega_t), (psi_r, psi_t) = fespace.TnT()

# Known source field H₀
H0 = CoefficientFunction((Hx0, Hy0, Hz0))  # From coils or PM

# Bilinear form
a = BilinearForm(fespace)

# Source region: ∇·μ(H₀ - ∇Ωᵣ) = 0
a += mu * grad(Omega_r) * grad(psi_r) * dx("source")

# Non-source region: ∇·μ(-∇Ωₜ) = 0
a += mu * grad(Omega_t) * grad(psi_t) * dx("air")

# Interface coupling (automatic with conforming FE)

# Linear form (source term from H₀)
f = LinearForm(fespace)
f += mu * H0 * grad(psi_r) * dx("source")

a.Assemble()
f.Assemble()

# Solve
gf = GridFunction(fespace)
gf.vec.data = a.mat.Inverse(fespace.FreeDofs()) * f.vec

# Extract solutions
gf_Omega_r, gf_Omega_t = gf.components

# Reconstruct H field
H_source = H0 - grad(gf_Omega_r)  # in source region
H_air = -grad(gf_Omega_t)         # in air region

# B field
B_source = mu * H_source
B_air = mu * H_air
```

### 境界条件

**Dirichlet境界 (外部境界):**
```python
fesT = H1(mesh, order=2, definedon="air",
          dirichlet="outer_boundary")
gf_Omega_t.Set(0, BND)  # Ωₜ = 0 on outer boundary
```

**Interface境界 (自動処理):**
- Conforming FE spaceを使用すると、界面での連続性は自動的に満たされる
- Ωᵣ = Ωₜ at interface (H1空間の連続性)
- μ₁∇Ωᵣ·n = μ₂∇Ωₜ·n (weak formulation)

### H₀の計算（Biot-Savart則）

**コイルによる磁場:**
```python
def BiotSavart(coil_current, coil_path, eval_points):
    """Compute H₀ from current-carrying coil"""
    mu0 = 4e-7 * np.pi

    H0 = np.zeros((len(eval_points), 3))

    for segment in coil_path:
        dl = segment.tangent * segment.length
        r_seg = segment.center

        for i, p in enumerate(eval_points):
            r = p - r_seg
            r_norm = np.linalg.norm(r)

            # Biot-Savart: H = (I/4π) ∫ (dl × r) / |r|³
            H0[i] += coil_current * np.cross(dl, r) / (4*np.pi*r_norm**3)

    return H0
```

**永久磁石による磁場:**
```python
def PermanentMagnetField(magnetization, magnet_volume, eval_points):
    """Compute H₀ from permanent magnet"""
    # Surface charge model or volume integral
    # H₀ = -∇Φₘ, where Φₘ = (1/4π) ∫ (M·∇)(1/|r-r'|) dV'

    # Simplified for uniform magnetization
    H0 = magnetization / (4*np.pi*mu0) * solid_angle_correction

    return H0
```

### Interface処理の注意点

**方法1: Conforming mesh (推奨)**
```python
# 界面でmeshが一致している場合、自動的に連続
fesR = H1(mesh, order=2, definedon="source")
fesT = H1(mesh, order=2, definedon="air")
# Interface DOFs are shared
```

**方法2: Non-conforming mesh**
```python
# Mortarメソッドが必要
# NGSolveではperiodic boundary条件を応用
```

### 利点と欠点

**利点:**
- 電流源領域を局所化できる（source regionのみ）
- 空気領域は単純なLaplace方程式
- メモリ効率が良い（領域分割）
- 永久磁石問題に適している

**欠点:**
- H₀の計算が必要（前処理）
- Interface処理が必要
- 非線形材料の場合は反復が複雑

### 実装例（完全版）

```python
from ngsolve import *
import numpy as np

def TwoScalarPotentialSolve(mesh, mu_source, mu_air, H0_source, **kwargs):
    """
    2スカラー法で静磁場問題を解く

    Parameters:
    -----------
    mesh : NGSolve mesh
        source領域とair領域を含むメッシュ
    mu_source : float or CF
        Source領域の透磁率
    mu_air : float
        空気領域の透磁率
    H0_source : CoefficientFunction
        Source領域の既知磁場（Biot-Savartなどで計算）
    """
    order = kwargs.get("order", 2)
    tol = kwargs.get("tolerance", 1e-10)

    # FE spaces
    fesR = H1(mesh, order=order, definedon="source")
    fesT = H1(mesh, order=order, definedon="air",
              dirichlet="outer_boundary")
    fespace = fesR * fesT

    (Omega_r, Omega_t), (psi_r, psi_t) = fespace.TnT()

    # System matrix
    a = BilinearForm(fespace)
    a += mu_source * grad(Omega_r) * grad(psi_r) * dx("source")
    a += mu_air * grad(Omega_t) * grad(psi_t) * dx("air")

    # Source term
    f = LinearForm(fespace)
    f += mu_source * H0_source * grad(psi_r) * dx("source")

    # Assemble and solve
    with TaskManager():
        a.Assemble()
        f.Assemble()

    gf = GridFunction(fespace)

    # Iterative solver
    inv = CGSolver(a.mat, fespace.FreeDofs(), maxsteps=1000, tol=tol)
    gf.vec.data = inv * f.vec

    gf_Omega_r, gf_Omega_t = gf.components

    # Reconstruct fields
    H_total = IfPos(mesh.MaterialCF("source", default=0) - 0.5,
                    H0_source - grad(gf_Omega_r),
                    -grad(gf_Omega_t))

    B_total = IfPos(mesh.MaterialCF("source", default=0) - 0.5,
                    mu_source * (H0_source - grad(gf_Omega_r)),
                    mu_air * (-grad(gf_Omega_t)))

    return {
        "Omega_r": gf_Omega_r,
        "Omega_t": gf_Omega_t,
        "H_field": H_total,
        "B_field": B_total,
        "solution": gf
    }
```

### H₀からΩへの変換（境界処理）

境界でのH₀をΩに変換する手法（`HtoOmega.py`パターン）:

```python
def HtoOmega(mesh, boundary, feOrder, H):
    """Convert H field to Omega on boundary"""
    fesOmega = H1(mesh, order=feOrder,
                  definedon=mesh.Boundaries(boundary))
    omega, psi = fesOmega.TnT()

    # Minimize ||∇Ω - H||² on boundary
    a = BilinearForm(fesOmega)
    a += grad(omega).Trace() * grad(psi).Trace() * ds

    f = LinearForm(fesOmega)
    f += (grad(psi).Trace() * H) * ds

    a.Assemble()
    f.Assemble()

    gfOmega = GridFunction(fesOmega)
    inv = CGSolver(a.mat, fesOmega.FreeDofs())
    gfOmega.vec.data = inv * f.vec

    return gfOmega
```

### 応用例

**1. 永久磁石モータ:**
```python
# PM領域: H₀ from magnetization
H0_PM = M / mu0  # M: magnetization vector

# Iron領域: High permeability
mu_iron = 1000 * mu0

# Solve 2-scalar
result = TwoScalarPotentialSolve(mesh, mu_iron, mu0, H0_PM)
```

**2. 電磁石:**
```python
# Coil領域: H₀ from Biot-Savart
H0_coil = BiotSavart(current=100, coil_geometry)

# Solve 2-scalar
result = TwoScalarPotentialSolve(mesh, mu0, mu0, H0_coil)
```

### ベストプラクティス

1. **H₀の精度**: Biot-Savart計算は十分な精度で（積分点数）
2. **Mesh refinement**: Interface近傍は細かく
3. **Material CF**: `IfPos()`で領域ごとに材料を切り替え
4. **Boundary conditions**: 外部境界は十分遠方に配置
5. **Non-linear materials**: Newton反復で透磁率を更新

### コイルによるΩジャンプ（Θジャンプ）

コイル電流が流れる領域を横切ると、磁気スカラーポテンシャルΩは不連続（ジャンプ）になります。これは**Θジャンプ**と呼ばれ、電流源による本質的な特性です。

#### 理論的背景

**Ampère's Law:**
```
∇×H = J
```

磁気スカラーポテンシャルを用いると H = -∇Ω ですが、コイル電流 J が存在する場合：
```
∇×(-∇Ω) = -∇×∇Ω = J  (矛盾！)
```

この矛盾は、Ωが多価関数であることで解決されます。コイル領域の切断面を横切るたびに：
```
[Ω]⁺₋ = Θ
```

ここで、Θは**コイルの巻数と電流に比例する量**：
```
Θ = N × I  (円筒コイルの場合)
```

より一般的には：
```
Θ = ∫ₛ (J·n) dS  (切断面を通る全電流)
```

#### 物理的意味

- **Θ = 起磁力 (Magnetomotive Force, MMF)**
- コイルを一周すると Ω は Θ だけ増加（または減少）
- 多価性により単純な ∇Ω では表現できない回転成分を表現

#### 2スカラー法での実装

**基本アイデア:**
- Coil内部では Ωᵣ（reduced potential）を使用し、ジャンプを許容
- 切断面（cut surface）を定義し、その面でΩが不連続であることを許可

**Method 1: 不連続FE空間（推奨）**

```python
from ngsolve import *

# Coil領域での不連続FE空間
fesR_disc = H1(mesh, order=2, definedon="coil",
               discontinuous_at="cut_surface")

# Air領域での連続FE空間
fesT = H1(mesh, order=2, definedon="air", dirichlet="outer")

fespace = fesR_disc * fesT

(Omega_r, Omega_t), (psi_r, psi_t) = fespace.TnT()

# ジャンプ項を追加
a = BilinearForm(fespace)
a += mu * grad(Omega_r) * grad(psi_r) * dx("coil")
a += mu * grad(Omega_t) * grad(psi_t) * dx("air")

# 切断面でのジャンプ条件: [Ω] = Θ
# Lagrange multiplier法で実装
```

**Method 2: Theta fieldを用いた処理（実用的）**

Θを既知の場数として与え、Ω = Ω_continuous + Θ と分解：

```python
from ngsolve import *
import numpy as np

def ComputeThetaField(mesh, coil_domain, cut_surface, NI, **kwargs):
    """
    コイル領域でのΘ field計算

    Parameters:
    -----------
    mesh : NGSolve mesh
    coil_domain : str
        コイル領域名
    cut_surface : str
        切断面境界名
    NI : float
        巻数×電流 (起磁力)
    """
    order = kwargs.get("order", 2)

    # H1空間（coil領域のみ）
    fesTheta = H1(mesh, order=order, definedon=coil_domain)
    theta, psi = fesTheta.TnT()

    # Laplace方程式（調和関数）
    a = BilinearForm(fesTheta)
    a += grad(theta) * grad(psi) * dx(coil_domain)

    # 切断面でのジャンプを境界条件として
    # 片側を0、反対側をNIに設定
    gfTheta = GridFunction(fesTheta)

    # Cut surfaceの片側でΘ=0、反対側でΘ=NI
    # Dirichlet境界条件
    gfTheta.Set(0, BND, mesh.Boundaries("cut_minus"))
    gfTheta.Set(NI, BND, mesh.Boundaries("cut_plus"))

    # 内部を解く
    f = LinearForm(fesTheta)
    f += -grad(gfTheta) * grad(psi) * dx(coil_domain)

    with TaskManager():
        a.Assemble()
        f.Assemble()

    # Solve with fixed boundary
    inv = CGSolver(a.mat, fesTheta.FreeDofs())
    gfTheta.vec.data += inv * f.vec

    return gfTheta

# 2スカラー法でΘを組み込む
def TwoScalarWithCoilJump(mesh, mu_coil, mu_air, H0, NI, **kwargs):
    """
    コイルジャンプを考慮した2スカラー法
    """
    # Θ field計算
    gfTheta = ComputeThetaField(mesh, "coil", "cut_surface", NI)

    # FE space setup
    fesR = H1(mesh, order=2, definedon="coil")
    fesT = H1(mesh, order=2, definedon="air", dirichlet="outer")
    fespace = fesR * fesT

    (Omega_r, Omega_t), (psi_r, psi_t) = fespace.TnT()

    # Bilinear form
    a = BilinearForm(fespace)
    a += mu_coil * grad(Omega_r) * grad(psi_r) * dx("coil")
    a += mu_air * grad(Omega_t) * grad(psi_t) * dx("air")

    # Linear form（Θの影響を含む）
    f = LinearForm(fespace)
    # H₀による項
    f += mu_coil * H0 * grad(psi_r) * dx("coil")
    # Θによる修正項（切断面での処理）
    f += -mu_coil * grad(gfTheta) * grad(psi_r) * dx("coil")

    with TaskManager():
        a.Assemble()
        f.Assemble()

    gf = GridFunction(fespace)
    inv = CGSolver(a.mat, fespace.FreeDofs())
    gf.vec.data = inv * f.vec

    gf_Omega_r, gf_Omega_t = gf.components

    # 真のΩᵣは Ω_continuous + Θ
    gf_Omega_r_full = gf_Omega_r + gfTheta

    # H field
    H_coil = H0 - grad(gf_Omega_r_full)
    H_air = -grad(gf_Omega_t)

    return {
        "Omega_r": gf_Omega_r_full,
        "Omega_t": gf_Omega_t,
        "Theta": gfTheta,
        "H_coil": H_coil,
        "H_air": H_air,
        "B_coil": mu_coil * H_coil,
        "B_air": mu_air * H_air
    }
```

#### 切断面（Cut Surface）の選び方

コイルの電流ループを「切断」する面を定義する必要があります：

**円筒コイル（軸対称）:**
```python
# 例：z軸まわりのコイル
# 切断面をθ=0平面に配置
cut_surface = "coil_cut"  # y≥0, x=0 の面
```

**任意形状コイル:**
```python
# コイル電流ループを完全に横切る面
# トポロジー的に「輪」を切る面
cut_surface = mesh.Boundaries("coil_cut_surface")
```

#### 実装上の注意

1. **切断面の一意性**: 切断面の選び方は任意だが、一度定めたら一貫させる
2. **Θの符号**: 電流の向きに応じて正負が決まる（右手の法則）
3. **メッシュとの整合**: 切断面はメッシュの面と一致させる必要がある
4. **H₀との関係**: H₀（Biot-Savart）には既にコイル電流の効果が含まれているため、Θと二重計上しないよう注意

#### 検証方法

**Ampèreの法則チェック:**
```python
# コイルを囲む閉路でH·dlを積分
# = N×I になるべき
contour = "loop_around_coil"
Hdl = Integrate(H * tangent * ds(contour), mesh)
print(f"Line integral H·dl = {Hdl}, Expected NI = {NI}")
assert abs(Hdl - NI) < 1e-6
```

#### 応用例

**電磁石設計:**
```python
# 鉄心入りコイル
NI = 1000  # 100turns × 10A
gfTheta = ComputeThetaField(mesh, "coil+iron", "cut", NI)

# 2スカラー法で磁場計算
result = TwoScalarWithCoilJump(
    mesh, mu_iron, mu0, H0=0, NI=NI
)

# 鉄心内の磁束密度
B_iron = result["B_coil"]
```

**モーター解析:**
```python
# 複数コイル（多極）
for pole_id in range(n_poles):
    theta_i = ComputeThetaField(
        mesh, f"coil_{pole_id}", f"cut_{pole_id}",
        NI * cos(pole_id * 2*pi/n_poles)
    )
    # 各極の寄与を重ね合わせ
```

## 参考文献

- EMPY_Analysis Repository: `S:\NGSolve\EMPY\EMPY_Analysis`
- T-Omega Method: `EddyCurrent/include/T_Omega_Method.py`
- Matrix Solver: `include/MatrixSolver.py`
- A-Phi Method: `include/A_Phi_Method.py`
- HtoOmega: `include/HtoOmega.py`
- Magnetostatic: `Magnetostatic/Omega-2Potential.ipynb`

---

**注意**: このドキュメントは EMPY_COIL は除外しています（ユーザー要求により）。
