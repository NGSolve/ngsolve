FROM ubuntu:23.10
ENV DEBIAN_FRONTEND=noninteractive
ENV PYODIDE_VERSION=0.26.2
ENV EMSDK_VERSION=3.1.58
ENV PYTHON_VERSION=3.12.1
RUN apt-get update && apt-get -y upgrade
RUN apt-get update && apt-get -y install \
        nodejs \
        ccache \
        python3-pip \
        python3 \
        git \
        make \
        pkg-config \
        g++ \
        lbzip2 \
        xz-utils \
        autoconf \
        libtool \
        unzip \
        xxd \
        wget

WORKDIR /root
RUN git clone https://github.com/emscripten-core/emsdk.git
RUN cd emsdk && ./emsdk install ${EMSDK_VERSION} && ./emsdk activate ${EMSDK_VERSION}
RUN wget https://github.com/pyodide/pyodide/releases/download/${PYODIDE_VERSION}/xbuildenv-${PYODIDE_VERSION}.tar.bz2
RUN tar xvf xbuildenv-${PYODIDE_VERSION}.tar.bz2
RUN ls
RUN mkdir /root/output
WORKDIR /root/output
RUN wget https://github.com/pyodide/pyodide/releases/download/${PYODIDE_VERSION}/pyodide-core-${PYODIDE_VERSION}.tar.bz2
RUN tar xvf pyodide-core-${PYODIDE_VERSION}.tar.bz2
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN rm /usr/lib/python3.11/EXTERNALLY-MANAGED
RUN pip install cmake
ADD ./tests/gitlab-ci/pyodide/build_opencascade.sh /root/
RUN bash /root/build_opencascade.sh
ADD . /root/ngsolve
RUN mkdir /root/ngsolve/build
WORKDIR /root/ngsolve/build
RUN --mount=type=cache,target=/ccache/ bash /root/ngsolve/tests/gitlab-ci/pyodide/build_in_docker.sh
WORKDIR /root/output/pyodide
RUN rm /opt/netgen/python/*/*.a
RUN apt install -y zip
WORKDIR /opt/netgen/python
RUN zip -r pyngcore.zip pyngcore
RUN zip -r netgen.zip netgen
RUN zip -r ngsolve.zip ngsolve
RUN cp /opt/netgen/python/*.zip /root/output/pyodide
WORKDIR /root/output/pyodide
RUN cp /root/ngsolve/tests/gitlab-ci/pyodide/generate_repodata.js .
RUN node generate_repodata.js
RUN rm *.whl
RUN python3 /root/ngsolve/tests/gitlab-ci/pyodide/merge.py
RUN cd /root/output/ && tar -cvjSf ngsolve_pyodide_${PYODIDE_VERSION}.tar.bz2 pyodide
