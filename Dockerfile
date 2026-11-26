# Build the solver without shipping the compiler, Autotools or source tree.
FROM ubuntu:24.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive
ARG BUILD_PROFILE=optimized

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        autoconf automake gfortran make \
        libnetcdff-dev libnetcdf-dev liblapack-dev libblas-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/IMEX_SfloW2D
COPY . .

RUN set -eu; \
    autoreconf -fi; \
    case "$BUILD_PROFILE" in \
        optimized) ./configure --with-netcdf=/usr FC=gfortran ;; \
        strict-debug) ./configure --with-netcdf=/usr FC=gfortran --enable-strict-debug ;; \
        *) echo "Unsupported BUILD_PROFILE: $BUILD_PROFILE" >&2; exit 2 ;; \
    esac; \
    make -j2; \
    test -x src/IMEX_SfloW2D

# Ubuntu 24.04 runtime shared libraries, with no build tools or test suites.
FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libnetcdff7 libnetcdf19t64 liblapack3 libblas3 libgomp1 libgfortran5 \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

LABEL org.opencontainers.image.title="IMEX-SfloW2D" \
      org.opencontainers.image.description="Depth-averaged finite-volume solver for geophysical flows" \
      org.opencontainers.image.source="https://github.com/demichie/IMEX_SfloW2D_dev" \
      org.opencontainers.image.licenses="GPL-3.0-only"

COPY --from=builder /opt/IMEX_SfloW2D/src/IMEX_SfloW2D /usr/local/bin/IMEX_SfloW2D

# Fail the image build if any runtime shared library is missing.
RUN set -eu; \
    ldd /usr/local/bin/IMEX_SfloW2D > /tmp/imex-linked-libraries.txt; \
    cat /tmp/imex-linked-libraries.txt; \
    ! grep -q "not found" /tmp/imex-linked-libraries.txt; \
    rm /tmp/imex-linked-libraries.txt

ENV OMP_NUM_THREADS=1 \
    OMP_DYNAMIC=FALSE

WORKDIR /work
ENTRYPOINT ["/usr/local/bin/IMEX_SfloW2D"]
