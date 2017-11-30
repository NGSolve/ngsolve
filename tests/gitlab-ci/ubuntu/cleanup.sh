## delete all containers belonging to images of the current pipeline
docker ps -a | awk '{ print $1,$2 }' | grep ngsolve_${CI_PIPELINE_ID} | awk '{print $1 }' | xargs -I {} docker rm {} || true

## list of images to delete
declare -a arr=("16.04" "16.10" "17.04" "17.10" "debug" "avx" "avx512")

for version in "${arr[@]}"
do
    docker rmi -f ngsolve_${CI_PIPELINE_ID}:$version || true
    docker rmi -f ngsolve_${CI_PIPELINE_ID}_installed:$version || true
done
