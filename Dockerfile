# syntax=docker/dockerfile:1
#
# Image karpos-surfex : SURFEX OFFLINE + assimilation (SODA) + orchestration Python.
# Multi-stage : builder (compile SURFEX + wheel) → runtime slim (exécutables +
# libsurfex.a + package, sans toolchain de build).
#
#   docker build -t karpos-surfex .
#   docker run --rm karpos-surfex where
#   docker run --rm -v $PWD/run:/work karpos-surfex run --workdir /work --steps pgd

# --- Stage 1 : builder -------------------------------------------------------
FROM ubuntu:22.04 AS builder
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      gfortran make perl libnetcdff-dev \
      python3 python3-pip python3-venv ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

# Build SURFEX (libsurfex.a + PGD/PREP/OFFLINE/SODA) via le wrapper portable.
RUN ./build.sh --compiler gfortran --jobs "$(nproc)"

# Wheel du package Python, lié à libsurfex.a.
RUN python3 -m pip install --no-cache-dir build \
 && cd python \
 && python3 -m build --wheel -Csetup-args=-Dsurfex_lib=/src/exe/libsurfex.a

# --- Stage 2 : runtime -------------------------------------------------------
FROM ubuntu:22.04 AS runtime
ENV DEBIAN_FRONTEND=noninteractive
# Libs runtime : NetCDF-Fortran (+C), libgfortran, OpenMP.
RUN apt-get update && apt-get install -y --no-install-recommends \
      libnetcdff7 libgfortran5 libgomp1 \
      python3 python3-venv \
 && rm -rf /var/lib/apt/lists/*

# Exécutables + bibliothèque + données ECOCLIMAP.
COPY --from=builder /src/exe /opt/surfex/exe
COPY --from=builder /src/MY_RUN/ECOCLIMAP/*.bin /opt/surfex/ecoclimap/
COPY --from=builder /src/python/dist/*.whl /tmp/

# venv isolé + install du package (+ netCDF4 pour get_field/set_obs).
RUN python3 -m venv /opt/venv \
 && /opt/venv/bin/pip install --no-cache-dir /tmp/*.whl netCDF4 \
 && rm -f /tmp/*.whl

ENV PATH=/opt/venv/bin:$PATH
ENV SURFEX_EXE_DIR=/opt/surfex/exe
ENV ECOCLIMAP_DIR=/opt/surfex/ecoclimap

WORKDIR /work
ENTRYPOINT ["python3", "-m", "karpos_surfex"]
CMD ["where"]
