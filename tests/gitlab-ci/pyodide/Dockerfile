# syntax=docker/dockerfile:1.7-labs
FROM ubuntu:24.10 AS builder
ENV DEBIAN_FRONTEND=noninteractive
ENV PYODIDE_VERSION=0.27.5
ENV EMSDK_VERSION=3.1.58
ENV PYTHON_VERSION=3.12.7
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
        zip \
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
RUN wget https://github.com/pyodide/pyodide/releases/download/${PYODIDE_VERSION}/pyodide-${PYODIDE_VERSION}.tar.bz2
RUN tar xvf pyodide-${PYODIDE_VERSION}.tar.bz2
RUN cd pyodide && unzip openblas*.zip
RUN cd pyodide && rm -r *.whl* *.tar fonts
RUN ls pyodide/*.zip | grep -v python_stdlib.zip | xargs rm -f --
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN pip install --break-system-packages "cmake<4"

FROM builder AS build_occ
ADD ./tests/gitlab-ci/pyodide/build_opencascade.sh /root/
RUN bash /root/build_opencascade.sh

FROM builder as build_ngsolve
COPY --from=build_occ /opt/opencascade /opt/opencascade
RUN mkdir /root/ngsolve_build
WORKDIR /root/ngsolve_build
COPY . /root/ngsolve
COPY tests/gitlab-ci/pyodide/build_in_docker.sh /root/ngsolve/tests/gitlab-ci/pyodide/
RUN --mount=type=cache,target=/ccache/ bash /root/ngsolve/tests/gitlab-ci/pyodide/build_in_docker.sh
RUN rm /opt/netgen/python/*/*.a
WORKDIR /opt/netgen/python
RUN zip -r pyngcore.zip pyngcore
RUN zip -r netgen.zip netgen
RUN zip -r ngsolve.zip ngsolve

FROM builder
COPY --from=build_ngsolve /opt/netgen/python/*.zip /root/output/pyodide
WORKDIR /root/output/pyodide
COPY tests/gitlab-ci/pyodide/generate_repodata.js .
RUN node generate_repodata.js
RUN rm *.whl libopenblas.so
COPY tests/gitlab-ci/pyodide/merge.py .
RUN python3 merge.py
RUN cd /root/output/ && tar -cvjSf ngsolve_pyodide_${PYODIDE_VERSION}.tar.bz2 pyodide
