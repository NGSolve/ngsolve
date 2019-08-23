cd %CI_PROJECT_DIR%
rd /s /q %CI_DIR%
ccache . -c
ccache . -s
