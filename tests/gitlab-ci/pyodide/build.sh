docker build -t ngsolve_pyodide -f tests/gitlab-ci/pyodide/Dockerfile .
mkdir -p dist
docker run -u `id -u` -v `pwd`/dist:/dist ngsolve_pyodide bash -c 'cd /root/output && cp ngsolve*.tar.bz2 /dist/'
scp dist/*.tar.bz2 mhochsteger@aurora:public_html/ngsolve/
