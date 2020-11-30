docker run \
      -e NETGENDIR=/opt/netgen/bin \
      -e NGS_NUM_THREADS=6 \
      -e PYTHONPATH=/opt/netgen/lib/python3/dist-packages \
      -e MKLROOT=/opt/intel/mkl \
      -e MKL_NUM_THREADS=1 \
      -e MKL_DOMAIN_THREADS="MKL_DOMAIN_PARDISO=6" \
      -e LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64 \
      -e RUN_SLOW_TESTS="$RUN_SLOW_TESTS" \
      -e NG_BACKTRACE=1 \
      -v /opt/intel:/opt/intel \
      ngsolve_${CI_PIPELINE_ID}_installed:$IMAGE_NAME \
      bash -c 'cd /root/build/ngsolve && make test_ngsolve ARGS="--output-on-failure"'
