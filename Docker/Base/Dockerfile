FROM ubuntu:16.04
RUN apt update && apt full-upgrade -y
RUN apt install -y build-essential \
                automake \
                autoconf \
                libtool \
                gfortran \
                cmake \
                wget \
                git \
                openmpi-bin \
                libopenmpi-dev \
                libglib2.0-dev \
                libnetpbm10-dev \
                libproj-dev \
                libgsl0-dev \
                libnetcdf-dev \
                libfftw3-dev \
                libgtkglext1-dev \
                libstartup-notification0-dev \
                ffmpeg
RUN git clone https://eessmann:2A10e20s@github.com/eessmann/GISS.git
RUN cd GISS/Third_Party/ode && bash bootstrap && \
    bash configure --enable-shared=yes --enable-double-precision \
                  --enable-builtin-threading-impl --enable-libccd \
                && make -j8 && make install
RUN cd GISS/Third_Party/hypre && mkdir hypre-build && cd hypre-build && \
    cmake -G "Unix Makefiles" ../src -DHYPRE_SHARED=ON -DHYPRE_BIGINT=ON \
          -HYPRE_USING_OPENMP=ON && make -j8 && make install
RUN cd GISS/Third_Party/gts && bash configure && make -j8 && make install
RUN cd GISS \
    && wget http://gerris.dalembert.upmc.fr/gerris/gerris-snapshot.tar.gz \
    && tar -xvf gerris-snapshot.tar.gz && cd gerris-snapshot \
    && bash configure && make -j8 && make install
