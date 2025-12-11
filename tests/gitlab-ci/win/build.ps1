$ErrorActionPreference = "Stop"

Write-Host "PATH in build"
Write-Host $env:PATH

python -m pip install --upgrade pybind11-stubgen
if ($LASTEXITCODE -ne 0) {
    Write-Error "pip install failed with exit code $LASTEXITCODE"
    exit $LASTEXITCODE
}

git submodule update --init --recursive
if ($LASTEXITCODE -ne 0) {
    Write-Error "git submodule update failed with exit code $LASTEXITCODE"
    exit $LASTEXITCODE
}

if (Test-Path $env:NETGEN_BUILD_DIR) {
    Remove-Item $env:NETGEN_BUILD_DIR -Recurse -Force
}

New-Item -ItemType Directory -Force -Path $env:NETGEN_BUILD_DIR | Out-Null
if (-not (Test-Path $env:NETGEN_BUILD_DIR)) {
    Write-Error "Failed to create build directory '$env:NETGEN_BUILD_DIR'"
    exit 1
}

Set-Location $env:NETGEN_BUILD_DIR

$cmakeArgs = @()
$cmakeArgs += $env:SRC_DIR
if ($env:CMAKE_GENERATOR) {
    $cmakeArgs += "-G$($env:CMAKE_GENERATOR)"
}
$cmakeArgs += "-DCMAKE_INSTALL_PREFIX=$($env:CMAKE_INSTALL_PREFIX)"
if ($env:CMAKE_CONFIG) {
    $cmakeArgs += $env:CMAKE_CONFIG -split '\s+'
}

$cmakeArgs += "-DUSE_NATIVE_ARCH=$($env:NG_USE_NATIVE_ARCH)"
$cmakeArgs += "-DUSE_CGNS=ON"
$cmakeArgs += "-DUSE_OCC=ON"
$cmakeArgs += "-DUSE_CCACHE=ON"
$cmakeArgs += "-DUSE_MKL=ON"
$cmakeArgs += "-DMKL_STATIC=ON"
$cmakeArgs += "-DMKL_ROOT=C:/Intel/compilers_and_libraries/windows/mkl"
$cmakeArgs += "-DUSE_UMFPACK=ON"
$cmakeArgs += "-DENABLE_UNIT_TESTS=ON"
$cmakeArgs += "-DCPACK_PACKAGE_NAME=NGSolve"
$cmakeArgs += "-DCMAKE_BUILD_TYPE=Release"

& cmake @cmakeArgs
if ($LASTEXITCODE -ne 0) {
    Write-Error "cmake configure failed with exit code $LASTEXITCODE"
    exit $LASTEXITCODE
}

& cmake --build . --target install --config Release
if ($LASTEXITCODE -ne 0) {
    Write-Error "cmake build/install failed with exit code $LASTEXITCODE"
    exit $LASTEXITCODE
}
