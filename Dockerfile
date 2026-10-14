# The plot-only target deliberately contains no qbasis, compiler, or solver.
FROM python:3.12-slim-bookworm@sha256:34386ef0cb081344d7ec1c103ba398e6e9f64e9ab3a1509accc92a4e24a07258 AS plot-runtime
RUN apt-get update && apt-get install -y --no-install-recommends \
 ca-certificates nodejs npm fonts-dejavu-core \
 && rm -r /var/lib/apt/lists/*
ENV OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
COPY requirements-repro.txt /tmp/requirements-repro.txt
RUN python -m pip install --no-cache-dir -r /tmp/requirements-repro.txt
COPY reproducibility/slides/package.json reproducibility/slides/package-lock.json /opt/mmhk/slides/
RUN cd /opt/mmhk/slides && npm ci --ignore-scripts --no-audit --no-fund
ENV MMHK_SLIDES_ROOT=/opt/mmhk/slides
WORKDIR /workspace
ENTRYPOINT ["python", "reproduce.py"]

FROM plot-runtime AS full
RUN apt-get update && apt-get install -y --no-install-recommends \
 g++ gfortran cmake make libboost-dev libopenblas-dev liblapacke-dev git \
 && rm -r /var/lib/apt/lists/*
WORKDIR /opt/mmhk
COPY CMakeLists.txt /opt/mmhk/
COPY cmake /opt/mmhk/cmake
COPY reproducibility/compat /opt/mmhk/reproducibility/compat
COPY create_code/chain_Hubbard_tbc_response.cc create_code/chain_Hubbard_tbc_spectral.cc /opt/mmhk/create_code/
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/mmhk/solver \
 && cmake --build build -j2 && cmake --install build \
 && dpkg-query -W > /opt/mmhk/system-packages.tsv
WORKDIR /workspace
ENV MMHK_SOLVER_DIR=/opt/mmhk/solver/bin
