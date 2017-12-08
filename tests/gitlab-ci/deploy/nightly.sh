if [ -n "${RUN_DEPLOY}" ];
then
  /home/gitlab-runner/deploy/deploy_nightly.sh $CI_PIPELINE_ID;
fi
