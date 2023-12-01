#
# Dockerfile for CI and Flexpart container images
# Build Examples:
#	- podman build -t harbor.wolke.img.univie.ac.at/flexpart/almalinux8 -f Dockerfile
#   - podman build -t harbor.wolke.img.univie.ac.at/flexpart/rockylinux9 -f Dockerfile
#
# FROM rockylinux:8-minimal
FROM rockylinux:9-minimal
#
# Build development image (with/without jasper)
# jasper was used in FP 10.4 (eccodes, emoslib)
#
# 8: --enablerepo=powertools
# since 9, there is 2.31 eccodes in epel
RUN microdnf install -y epel-release && \
	microdnf install -y --enablerepo=crb make netcdf-fortran-devel.x86_64 netcdf.x86_64 eccodes eccodes-devel cmake tar gcc-c++ perl git && \
	microdnf clean all -y && \
	rm -rf /var/cache/yum

#
# Download ECCODES Version 
# note 2.30.0 has an issue!!!
#
# RUN curl https://confluence.ecmwf.int/download/attachments/45757960/eccodes-2.31.0-Source.tar.gz | tar xz
# RUN mkdir build && \
# 	cd build && \
# 	cmake -DENABLE_ECCODES_OMP_THREADS=ON ../eccodes-*/ && \
# 	make -j8 && \
# 	make install
#
# set environment variables
#
ENV FC=gfortran
ENV LIBRARY_PATH=/usr/lib64:/usr/local/lib64
ENV CPATH=/usr/include:/usr/local/include:/usr/lib64/gfortran/modules
