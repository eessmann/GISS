FROM ubuntu:16.04
RUN apt update && apt full-upgrade -y
RUN apt install -y build-essential \
                automake \
                autoconf \
                libtool \
                gfortran \
                wget \
                git \
                darcs \
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
                libgtk2.0-dev \
                libstartup-notification0-dev \
                libftgl-dev \
                libosmesa6-dev \
                ffmpeg

RUN git clone https://eessmann:2A10e20s@github.com/eessmann/GISS.git \
    && cd /GISS && darcs get http://gerris.dalembert.upmc.fr/darcs/gts-stable \
    && darcs get http://gerris.dalembert.upmc.fr/darcs/gerris-stable \
    && darcs get http://gerris.dalembert.upmc.fr/darcs/gfsview-stable \
    && ln -s /usr/bin/libtoolize /usr/bin/libtool

RUN cd /GISS/Third_Party/ode && bash bootstrap && \
    bash configure --enable-shared=yes --enable-double-precision \
                  --enable-builtin-threading-impl --enable-libccd \
                && make -j8 && make install

RUN export HYPRE=/usr/local/hypre &&  \
    export LD_LIBRARY_PATH=${HYPRE}/lib:$LD_LIBRARY_PATH && \
    cd /GISS/Third_Party/hypre/src && \
    ./configure --prefix=$HYPRE --enable-shared && \
    make -j8 && make install && \
    /sbin/ldconfig

RUN cd /GISS/gts-stable && bash autogen.sh && ./configure && make -j8 && make install
RUN cd /GISS/gerris-stable && export HYPRE=/usr/local/hypre && bash autogen.sh \
    && CPPFLAGS=-I$HYPRE/include LDFLAGS="-L$HYPRE/lib" ./configure \
    && make -j8 && make install
RUN cd /GISS/gfsview-stable && bash autogen.sh && ./configure && make -j8 && \
    make install && rm -rfv /GISS

ENTRYPOINT ["gerris3D"]
