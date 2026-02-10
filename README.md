# NGSolve (ksugahar fork)

Multi-purpose finite element library

Find the Open Source Community on https://ngsolve.org

Support & Services: https://cerbsim.com

## Fork Features

This fork includes additional features not yet in the upstream NGSolve:

### 1. SetGeomInfo API (PR #232)

Python API for high-order curving of externally imported meshes.
Enables accurate curving of meshes imported from external mesh generators (e.g., Coreform Cubit).

```python
from netgen.occ import *
from ngsolve import *

# Set UV parameters for curved surface elements
for el in mesh.Elements2D():
    for i, v in enumerate(el.vertices):
        el.SetGeomInfo(i, PointGeomInfo(u=..., v=...))

# Apply high-order curving
mesh.Curve(order=3)
```

### 2. Intel MKL Support

Build with Intel MKL for improved performance:
```bash
cmake .. -DUSE_MKL=ON
```

### 3. SparseSolv Preconditioners

Integrated preconditioners from [JP-MARs/SparseSolv](https://github.com/JP-MARs/SparseSolv):

- **ICPreconditioner**: Shifted Incomplete Cholesky preconditioner
- **SGSPreconditioner**: Symmetric Gauss-Seidel preconditioner

```python
from ngsolve import *
from ngsolve.krylovspace import CGSolver

# Create bilinear form and assemble
a = BilinearForm(fes)
a += grad(u)*grad(v)*dx
a.Assemble()

# Create IC preconditioner
pre = ICPreconditioner(a.mat, shift=1.05)
pre.Update()

# Use with CGSolver
inv = CGSolver(a.mat, pre, printrates=True, tol=1e-10)
gfu.vec.data = inv * f.vec
```

## Binary Releases

Pre-built binaries are available on GitHub Releases:
https://github.com/ksugahar/ngsolve/releases

| Version | Features |
|---------|----------|
| v6.2.2601-setgeominfo-mkl | SetGeomInfo API + Intel MKL + SparseSolv (Recommended) |
| v6.2.2601-setgeominfo | SetGeomInfo API + SparseSolv (OpenBLAS) |

## Building from Source

### Requirements

- CMake 3.16+
- C++17 compiler (MSVC, GCC, Clang)
- Python 3.x
- (Optional) Intel MKL

### Build Commands

```bash
git clone --recursive https://github.com/ksugahar/ngsolve.git
cd ngsolve
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
```

With Intel MKL:
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_MKL=ON
```

## Related Repositories

- [ksugahar/netgen](https://github.com/ksugahar/netgen) - Netgen fork with SetGeomInfo API
- [JP-MARs/SparseSolv](https://github.com/JP-MARs/SparseSolv) - Original SparseSolv library

## Contributors

- NGSolve Team (upstream)
- Kengo Sugahara (fork maintainer)
- SparseSolv Contributors (Takahiro Sato, Shingo Hiruma, Tomonori Tsuburaya)

## License

See LICENSE file for details.
