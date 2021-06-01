sed -e "s/{UBUNTU_VERSION}/$UBUNTU_VERSION/" tests/gitlab-ci/ubuntu/docker_template >> docker_file
docker build -t ngsolve_${CI_PIPELINE_ID}:${IMAGE_NAME} -f docker_file .

rm -f ngsolve_${CI_PIPELINE_ID}_${IMAGE_NAME}.id

mkdir logs

docker run \
      --cidfile ngsolve_${CI_PIPELINE_ID}_${IMAGE_NAME}.id \
      -e MKLROOT=/opt/intel/mkl \
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
      -v /opt/intel:/opt/intel \
      -v `pwd`/logs:/logs \
      -v /mnt/ccache:/ccache ngsolve_${CI_PIPELINE_ID}:${IMAGE_NAME} \
      bash /root/src/ngsolve/tests/gitlab-ci/ubuntu/build_in_docker.sh

docker commit `cat ngsolve_${CI_PIPELINE_ID}_${IMAGE_NAME}.id` ngsolve_${CI_PIPELINE_ID}_installed:${IMAGE_NAME}
