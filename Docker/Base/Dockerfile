FROM fedora:24
RUN dnf -y update && \
    dnf -y groupinstall "Development Tools" && \
    dnf -y groupinstall "C Development Tools and Libraries" && \
    dnf -y install gcc-gfortran \
                    environment-modules \
                    openmpi \
                    openmpi-devel \
                    glib2 \
                    glib2-devel \
                    netpbm \
                    netpbm-devel \
                    proj \
                    proj-devel \
                    gsl \
                    gsl-devel \
                    netcdf \
                    netcdf-openmpi \
                    netcdf-devel \
                    netcdf-fortran \
                    netcdf-openmpi-devel \
                    netcdf-fortran-openmpi \
                    netcdf-fortran-devel \
                    netcdf-fortran-openmpi-devel \
                    fftw \
                    fftw-devel \
                    gtkglext-devel \
                    startup-notification \
                    startup-notification-devel && \
               
rpm -Uvh http://download1.rpmfusion.org/free/fedora/rpmfusion-free-release-stable.noarch.rpm http://download1.rpmfusion.org/nonfree/fedora/rpmfusion-nonfree-release-stable.noarch.rpm && \

dnf -y update && \
dnf -y install ffmpeg \
               ffmpeg-devel

RUN git clone https://eessmann:2A10e20s@github.com/eessmann/GISS.git

RUN bash /GISS/Docker/GISS_PATH.sh
