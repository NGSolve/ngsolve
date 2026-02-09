$ErrorActionPreference = "Stop"

Set-Location $env:CI_PROJECT_DIR
if (Test-Path $env:CI_DIR) {
    Remove-Item $env:CI_DIR -Recurse -Force
}
