FROM rockylinux:9-minimal
RUN microdnf install -y epel-release && \
	microdnf install -y --enablerepo=crb make netcdf-fortran-devel.x86_64 netcdf.x86_64 eccodes eccodes-devel cmake tar gcc-c++ perl git wget && \
	microdnf clean all -y && \
	rm -rf /var/cache/yum
#
# Download ECCODES Version 
# note 2.30.0 has an issue!!!
# ENV ECCODES=2.31.0
# RUN curl https://confluence.ecmwf.int/download/attachments/45757960/eccodes-${ECCODES}-Source.tar.gz | tar xz
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