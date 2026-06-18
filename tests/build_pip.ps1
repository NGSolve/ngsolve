$ErrorActionPreference = "Stop"

# Build and upload the Windows wheels for all supported CPython versions.
#
# We build directly against the python.org installations (C:\Python3XX) rather
# than via cibuildwheel: NGSolve/Netgen's GUI needs Tcl/Tk + the tcl/tk stub
# libraries and tkinter, none of which are available in the cibuildwheel environments.
function Invoke-Step([scriptblock]$Block) {
    & $Block
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed with exit code ${LASTEXITCODE}: $Block"
    }
}

if (test-path wheelhouse) {
    cmd.exe /c rd /s /q wheelhouse
}

$env:NETGEN_ARCH = "avx2"

$env:CCACHE_BASEDIR = (Get-Location).Path
if (-not $env:CCACHE_DIR) { $env:CCACHE_DIR = "C:\ccache" }

$pythons = @(
    "C:\Python314",
    "C:\Python313",
    "C:\Python312",
    "C:\Python311",
    "C:\Python310"
)

# Netgen version this NGSolve build links against (matching wheel from PyPI)
$netgen_version = (& C:\Python312\python.exe external_dependencies\netgen\tests\utils.py --get-version --dir=.\external_dependencies\netgen)
Write-Host "Netgen version: $netgen_version"

foreach ($pydir in $pythons) {
    $python = Join-Path $pydir "python.exe"
    Invoke-Step { & $python --version }
    Invoke-Step { & $python -m pip install -U pip requests packaging }
    # skip this version if it is already published on PyPI (non-zero = skip)
    & $python tests\utils.py --check-pip --package ngsolve
    if ($LASTEXITCODE -ne 0) { continue }

    # Works around broken packages or dependencies
    & $python -m pip uninstall -y bleach requests-toolbelt

    Invoke-Step { & $python -m pip install -U build delvewheel "scikit-build-core>=0.10" pybind11-stubgen ninja cmake }
    Invoke-Step { & $python -m pip install -U netgen-occt==7.8.1 netgen-occt-devel==7.8.1 "netgen-mesher==$netgen_version" "ngsolve-openblas==$env:OPENBLAS_VERSION_PIP" }

    if (test-path dist) { cmd.exe /c rd /s /q dist }

    $prefix = $pydir -replace '\\', '/'
    $env:CMAKE_ARGS = "-DCMAKE_PREFIX_PATH=$prefix -DPython3_ROOT_DIR=$prefix"

    Invoke-Step { & $python -m build --wheel --no-isolation --outdir dist . }

    if (Get-Command ccache -ErrorAction SilentlyContinue) { & ccache -s }

    New-Item -ItemType Directory -Force -Path wheelhouse | Out-Null
    Get-ChildItem dist\*.whl | ForEach-Object {
        Invoke-Step { & $python -m delvewheel repair --ignore-existing --namespace-pkg netgen;netgen.libs;netgen.include --exclude "nglib.dll" --exclude "ngcore.dll" --exclude "TK*.dll" --exclude "libopenblas*.dll" --exclude "openblas*.dll" -w wheelhouse $_.FullName }
    }

    Invoke-Step { & $python -m pip install -U twine requests-toolbelt urllib3 }
    Invoke-Step { & $python -m twine upload --skip-existing wheelhouse\*.whl }
}
