docker build -t ngsolve_pyodide -f tests/gitlab-ci/pyodide/Dockerfile .
rm -rf dist
mkdir -p dist
USERID=`id -u`
GROUPID=`id -g`
# docker run -u `id -u` -v `pwd`/dist:/dist ngsolve_pyodide bash -c 'cd /root/output && cp ngsolve*.tar.bz2 /dist/'
docker run -v `pwd`/dist:/dist ngsolve_pyodide bash -c "cd /root/output && cp ngsolve*.tar.bz2 /dist/ && chown -R ${USERID}:${GROUPID} /dist/"

# push only if $CI_COMMIT_REF_NAME is "master"
if [ "$CI_COMMIT_REF_NAME" != "master" ]; then
  exit 0
fi
scp -P 692 dist/*.tar.bz2 deploy@ngsolve.org:/opt/data/files/
