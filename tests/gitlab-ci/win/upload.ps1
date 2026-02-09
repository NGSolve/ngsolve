$ErrorActionPreference = "Stop"

Set-Location $env:NETGEN_BUILD_DIR
if ($LASTEXITCODE -ne 0) {
    Write-Error "Failed to cd to NETGEN_BUILD_DIR '$env:NETGEN_BUILD_DIR' (exit code $LASTEXITCODE)"
    exit $LASTEXITCODE
}

& cmake --build . --target package --config Release
if ($LASTEXITCODE -ne 0) {
    Write-Error "cmake --build ... --target package failed (exit code $LASTEXITCODE)"
    exit $LASTEXITCODE
}

$env:HOME = "/cygdrive/c/tools"

$pipelineId = $env:CI_PIPELINE_ID
$rsyncPath  = "mkdir -p deploy/builds/$pipelineId/windows && rsync"
$remoteDest = "gitlab-runner@vector.asc.tuwien.ac.at:deploy/builds/$pipelineId/windows/"

& rsync `
    -ztrl `
    --del `
    -e "ssh -o StrictHostKeyChecking=no" `
    "--rsync-path=$rsyncPath" `
    *.msi `
    $remoteDest

if ($LASTEXITCODE -ne 0) {
    Write-Error "rsync upload failed with exit code $LASTEXITCODE"
    exit $LASTEXITCODE
}

Write-Host "Upload finished successfully."
