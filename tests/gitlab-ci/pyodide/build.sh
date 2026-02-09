set -e
docker build -t ngsolve_pyodide -f tests/gitlab-ci/pyodide/Dockerfile .
rm -rf dist
mkdir -p dist
USERID=`id -u`
GROUPID=`id -g`
# docker run -u `id -u` -v `pwd`/dist:/dist ngsolve_pyodide bash -c 'cd /root/output && cp ngsolve*.tar.bz2 /dist/'
docker run -v `pwd`/dist:/dist ngsolve_pyodide bash -c "cd /root/output && cp ngsolve*.tar.bz2 /dist/ && chown -R ${USERID}:${GROUPID} /dist/"

# append branch name to file name if it's other than master
cd dist
if [ "$CI_COMMIT_REF_NAME" != "master" ]; then
  rename "s/ngsolve_pyodide_/ngsolve_pyodide_${CI_COMMIT_REF_NAME}_/" *.tar.bz2
fi

# do this only in CI pipeline (i.e. if $CI_COMMIT_REF_NAME is set)
if [ -n "${CI_COMMIT_REF_NAME:-}" ]; then
  scp -P 692 *.tar.bz2 deploy@ngsolve.org:/opt/data/files/
  tar xvf *.tar.bz2
  ssh -p 692 deploy@ngsolve.org "mkdir -p /opt/data/files/pyodide-0.29.1/${CI_COMMIT_REF_NAME}"
  scp -P 692 pyodide/libngcore.zip pyodide/libnglib.zip pyodide/pyngcore.zip pyodide/netgen.zip pyodide/ngsolve.zip deploy@ngsolve.org:/opt/data/files/pyodide-0.29.1/${CI_COMMIT_REF_NAME}/
fi
