# can build that image locally before if needed
# podman build -t harbor.wolke.img.univie.ac.at/flexpart/rockylinux9:latest -f ../tests/Dockerfile
FROM harbor.wolke.img.univie.ac.at/flexpart/rockylinux9:latest
# cmd Arguments for COMMIT ID/tag
ARG COMMIT=0
ENV COMMIT=$COMMIT
COPY ./src /src
# since the new data set is larger, we download the data for the container
# EC 2009 01 01 00 - 06
RUN wget  -r -np -nd -A 'EC*' -P /inputs "https://webdata.wolke.img.univie.ac.at/flexpart/"
WORKDIR /src
#
# compile using a standardized 
# need to fix march to x86_64_v3 , remove mtune
# here generic means core_avx2
RUN make -f makefile_gfortran clean \
	&& make -f makefile_gfortran eta=no arch=generic \
	&& make -f makefile_gfortran arch=generic \
	&& mkdir -p /output \
	&& echo -e "/options/\n/output/\n/inputs/\n/inputs/AVAILABLE" > /pathnames

COPY ./AVAILABLE /inputs/AVAILABLE
# Entrypoint script (user tips)
COPY ./containers/entrypoint.sh /entrypoint.sh
COPY ./tests/default_options /options
# In the Dockerfile, the ENTRYPOINT command defines the executable, 
# while CMD sets the default parameter. 
# so here our entrypoint run script and default FLEXPART_ETA or FLEXPART
ENTRYPOINT ["/entrypoint.sh"]
# Could be FLEXPART_ETA pathnames_custom
CMD ["FLEXPART_ETA"]
