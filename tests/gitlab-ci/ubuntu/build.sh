sed -e "s/{UBUNTU_VERSION}/$UBUNTU_VERSION/" tests/gitlab-ci/ubuntu/docker_template >> docker_file
docker build -t ngsolve_${CI_PIPELINE_ID}:$IMAGE_VERSION -f docker_file .

rm -f ngsolve_${CI_PIPELINE_ID}_$IMAGE_VERSION.id

docker run \
      --cidfile ngsolve_${CI_PIPELINE_ID}_$IMAGE_VERSION.id \
      -e MKLROOT=/opt/intel/mkl \
      -e UBUNTU_VERSION_NAME=$UBUNTU_VERSION_NAME \
      -e IMAGE_VERSION=$IMAGE_VERSION \
      -e CI_BUILD_REF=$CI_BUILD_REF \
      -e CI_PIPELINE_ID=$CI_PIPELINE_ID \
      -e CMAKE_ARGS="-DCPACK_DEBIAN_PACKAGE_NAME=ngsolve${PACKAGE_NAME_SUFFIX}" \
      -e CCACHE_DIR=/ccache \
      -e SSH_PRIVATE_KEY="$SSH_PRIVATE_KEY" \
      -v /opt/intel:/opt/intel \
      -v /mnt/ccache:/ccache ngsolve_${CI_PIPELINE_ID}:$IMAGE_VERSION \
      bash /root/src/ngsolve/tests/gitlab-ci/ubuntu/build_in_docker.sh

docker commit `cat ngsolve_${CI_PIPELINE_ID}_$IMAGE_VERSION.id` ngsolve_${CI_PIPELINE_ID}_installed:$IMAGE_VERSION
