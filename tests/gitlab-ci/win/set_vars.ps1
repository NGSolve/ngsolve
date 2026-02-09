$ErrorActionPreference = "Stop"

$vcvarsCmd = '"C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvars64.bat"'

cmd /c "$vcvarsCmd && set" |
    ForEach-Object {
        if ($_ -match '^(.*?)=(.*)$') {
            $name  = $matches[1]
            $value = $matches[2]

            # Make vars visible both to this PS session and child processes
            Set-Item -Path "Env:$name" -Value $value -ErrorAction SilentlyContinue
            [System.Environment]::SetEnvironmentVariable($name, $value)
        }
    }

$env:CMAKE_GENERATOR    = "Ninja"
$env:NETGEN_BUILD_DIR   = Join-Path $env:CI_PROJECT_DIR "build"
$env:CMAKE_INSTALL_PREFIX = Join-Path $env:CI_PROJECT_DIR "install"
$env:SRC_DIR            = $env:CI_PROJECT_DIR
$env:NETGENDIR          = Join-Path $env:CMAKE_INSTALL_PREFIX "bin"

$env:PATH = @(
    $env:NETGENDIR
    "C:\python312"
    "C:\python312\bin"
    "C:\python312\Scripts"
    "C:\tools"
    $env:PATH
) -join ";"

$env:PYTHONPATH = Join-Path $env:CMAKE_INSTALL_PREFIX "lib\site-packages"

$env:HOME = $env:USERPROFILE

& ccache -M 20G
& ccache -s

Write-Host "PATH in set_vars"
Write-Host $env:PATH
