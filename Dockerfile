FROM alonsojasl/ubuntu-cemrg:base
# FROM ubuntu-cemrg:lite

LABEL maintainer="Jose Alonso Solis-Lemus <jose.solislemus@kcl.ac.uk>"
LABEL Description="MeshTools3D executables"

# CGAL4.8.2 in: /installs/cgal4.8.2-build/release
# Eigen headers in: /installs/Eigen-install
# Boost installed in: /installs/boost_1_65_0/build/
# CMake installed in: /installs/cmake-master/bin/cmake

COPY . /dependencies/meshtools3D_release

RUN mkdir -p /installs/meshtools3D_build

RUN cd /installs/meshtools3D_build && \
    /installs/cmake-master/bin/cmake \
          -DCGAL_DIR=/installs/cgal4.8.2-build/release \
          -DBOOST_ROOT=/installs/boost_1_65_0/build/ \
          -DBoost_NO_BOOST_CMAKE=ON \
          /dependencies/meshtools3D_release
RUN cd /installs/meshtools3D_build && make

RUN rm -rf /dependencies

COPY ./dockerM3D.sh /installs/meshtools3D_build/dockerM3D.sh

RUN chmod +x /installs/meshtools3D_build/dockerM3D.sh

# ENTRYPOINT ["/installs/meshtools3D_build/meshtools3d"]
ENTRYPOINT ["/installs/meshtools3D_build/dockerM3D.sh"]

WORKDIR /data
