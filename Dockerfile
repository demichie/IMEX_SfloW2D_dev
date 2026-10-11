# IMEX_SfloW2D_dev - reproducible Linux build/test environment.
# Builds the checked-out repository supplied as the Docker build context.
# For long-term reproducibility, pin the Ubuntu image to a tested digest.
FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
ARG BUILD_PROFILE=optimized

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       autoconf automake gfortran make \
       liblapack-dev libblas-dev libnetcdff-dev libnetcdf-dev \
       python3 python3-numpy unzip ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/IMEX_SfloW2D
COPY . .

# Match .github/workflows/solver-validation.yml: generate, configure, build.
RUN autoreconf -fi \
    && case "$BUILD_PROFILE" in \
         optimized) ./configure --with-netcdf=/usr FC=gfortran ;; \
         strict-debug) ./configure --with-netcdf=/usr FC=gfortran --enable-strict-debug ;; \
         *) echo "Unsupported BUILD_PROFILE: $BUILD_PROFILE" >&2; exit 2 ;; \
       esac \
    && make -j2 \
    && test -x src/IMEX_SfloW2D

ENV OMP_NUM_THREADS=1 \
    OMP_DYNAMIC=FALSE \
    IMEX_AUDIT_FORCE_THREADS=1 \
    PYTHONDONTWRITEBYTECODE=1

# The same four test stages executed by the existing solver CI job.
# Override this CMD to use the compiled executable in a mounted case directory.
CMD ["/bin/sh", "-ec", "make check; TESTS/TEST_BINARY_RESTART/run_test.sh src/IMEX_SfloW2D; TESTS/TEST_STOCHASTIC_RESTART/run_test.sh src/IMEX_SfloW2D; TESTS/TEST_PCCU_GATE_H/run_test.sh src/IMEX_SfloW2D"]
