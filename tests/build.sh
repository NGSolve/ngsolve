cd 
cd src/ngsolve
git submodule update --init --recursive
cd
mkdir -p build/ngsolve
cd build/ngsolve
cmake ../../src/ngsolve -DUSE_CCACHE=ON -DUSE_MKL=ON -DUSE_UMFPACK=ON -DINSTALL_PROFILES=ON -DMKL_STATIC=ON -DCMAKE_BUILD_TYPE=Release
make -j12
make install
make package

# Run ssh-agent (inside the build environment)
eval $(ssh-agent -s)
# Add the SSH key stored in SSH_PRIVATE_KEY variable to the agent store
ssh-add <(echo "$SSH_PRIVATE_KEY")
mkdir -p ~/.ssh
[[ -f /.dockerenv ]] && echo -e "Host *\n\tStrictHostKeyChecking no\n\n" > ~/.ssh/config

ssh tester@vector.asc.tuwien.ac.at "mkdir -p /home/tester/deploy/ubuntu/$CI_BUILD_REF/${UBUNTU_VERSION_NAME}_amd64"
scp *.deb tester@vector.asc.tuwien.ac.at:/home/tester/deploy/ubuntu/$CI_BUILD_REF/${UBUNTU_VERSION_NAME}_amd64/
        
