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

## 参考文献

- EMPY_Analysis Repository: `S:\NGSolve\EMPY\EMPY_Analysis`
- T-Omega Method: `EddyCurrent/include/T_Omega_Method.py`
- Matrix Solver: `include/MatrixSolver.py`
- A-Phi Method: `include/A_Phi_Method.py`

---

**注意**: このドキュメントは EMPY_COIL は除外しています（ユーザー要求により）。
