rsync -ztrl --del -e ssh \
      --rsync-path="mkdir -p deploy/builds/$CI_PIPELINE_ID/macos && rsync" \
      *.dmg \
      gitlab-runner@vector.asc.tuwien.ac.at:deploy/builds/$CI_PIPELINE_ID/macos/
