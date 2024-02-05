docker build -t ngsolve_pyodide -f tests/gitlab-ci/pyodide/Dockerfile .
mkdir -p dist
USERID=`id -u`
GROUPID=`id -g`
# docker run -u `id -u` -v `pwd`/dist:/dist ngsolve_pyodide bash -c 'cd /root/output && cp ngsolve*.tar.bz2 /dist/'
docker run -v `pwd`/dist:/dist ngsolve_pyodide bash -c "cd /root/output && cp ngsolve*.tar.bz2 /dist/ && chown -R ${USERID}:${GROUPID} /dist/"
# scp dist/*.tar.bz2 mhochsteger@aurora:public_html/ngsolve/
