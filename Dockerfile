#
# Dockerfile for CI and Flexpart container images
# Options:
#	- JASPER = 0/1 (include jasper libs?)
#	- COMMIT = ... (use git log --pretty=format:'%h %D (%s %at)' -n 1)
#			   with this option Flexpart is build inside the container
# Build Examples:
#	- podman build -t flexpart:fp11 -f Dockerfile --build_arg JASPER=0
#
ARG JASPER
FROM almalinux:8-minimal 
#
# Build development image (with/without jasper)
# jasper was used in FP 10.4 (eccodes, emoslib)
#
RUN microdnf install -y epel-release
RUN if [ -z $JASPER ]; then \
	microdnf install -y --enablerepo=powertools make netcdf-fortran-devel.x86_64 netcdf.x86_64 cmake tar gcc-c++ perl; \
	else \
	microdnf install -y --enablerepo=powertools make netcdf-fortran-devel.x86_64 netcdf.x86_64 cmake tar gcc-c++ perl jasper-libs.x86_64; \
	fi;
#
# Download ECCODES Version 
# note 2.30.0 has an issue!!!
#
RUN curl https://confluence.ecmwf.int/download/attachments/45757960/eccodes-2.31.0-Source.tar.gz | tar xz
RUN mkdir build && \
	cd build && \
	cmake -DENABLE_ECCODES_OMP_THREADS=ON ../eccodes-*/ && \
	make -j8 && \
	make install
#
# set environment variables
#
ENV FC=gfortran
ENV LIBRARY_PATH=/usr/lib64:/usr/local/lib64
ENV CPATH=/usr/include:/usr/local/include:/usr/lib64/gfortran/modules
