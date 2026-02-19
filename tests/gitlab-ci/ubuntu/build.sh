sed -e "s/{UBUNTU_VERSION}/$UBUNTU_VERSION/" tests/gitlab-ci/ubuntu/docker_template > docker_file
docker build -t ngsolve_${CI_PIPELINE_ID}:${IMAGE_NAME} -f docker_file .
trap 'docker rmi -f ngsolve_${CI_PIPELINE_ID}:${IMAGE_NAME} || true' EXIT

mkdir -p logs
docker run \
      --rm \
      -e MKLROOT=/opt/intel/mkl \
      -e CI=$CI \
      -e CI_BUILD_REF=$CI_BUILD_REF \
      -e CI_PIPELINE_ID=$CI_PIPELINE_ID \
      -e IMAGE_NAME=$IMAGE_NAME \
      -e CMAKE_ARGS="-DCPACK_DEBIAN_PACKAGE_NAME=ngsolve${PACKAGE_NAME_SUFFIX}" \
      -e CCACHE_DIR=/ccache \
      -e SSH_PRIVATE_KEY="$SSH_PRIVATE_KEY" \
      -e MKL_NUM_THREADS=1 \
      -e MKL_DOMAIN_THREADS="MKL_DOMAIN_PARDISO=6" \
      -e LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64 \
      -e RUN_SLOW_TESTS="$RUN_SLOW_TESTS" \
      -e SHOW_LOGS="$SHOW_LOGS" \
      -e NETGENDIR=/opt/netgen/bin \
      -e NGS_NUM_THREADS=6 \
      -e PYTHONPATH=/opt/netgen/lib/python3/dist-packages \
      -e NG_BACKTRACE=1 \
      -e RUN_TESTS_AFTER_BUILD="${RUN_TESTS_AFTER_BUILD:-1}" \
      -v /opt/intel:/opt/intel \
      -v /mnt/ccache:/ccache \
      -v "$PWD/logs:/logs" \
      ngsolve_${CI_PIPELINE_ID}:${IMAGE_NAME} \
      bash -c '/root/src/ngsolve/tests/gitlab-ci/ubuntu/build_in_docker.sh && if [ "$RUN_TESTS_AFTER_BUILD" = "1" ]; then cd /root/build/ngsolve && make test_ngsolve ARGS="--output-on-failure"; fi'
