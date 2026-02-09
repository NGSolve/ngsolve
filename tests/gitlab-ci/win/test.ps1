$ErrorActionPreference = "Stop"

$testDir = Join-Path $env:NETGEN_BUILD_DIR "ngsolve"
Set-Location $testDir

python -m pip install scipy
ctest -C Release --output-on-failure
