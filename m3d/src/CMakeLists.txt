set(MESHTOOLS3D_INCLUDES
  "../include/CGALDataType.hpp"
  "../include/Chrono.hpp"
  "../include/GetPot.hpp"
  "../include/INRreader.hpp"
  "../include/LaplaceSolver.hpp"
  "../include/Mesh.hpp"
  "../include/mgmres.hpp"
  "../include/MyImageWrapper.hpp"
  "../include/My_Labeled_image_mesh_domain_3.hpp"
  "../include/ThicknessEvaluation.hpp"
  "../include/VtkWriter.hpp"
)

set(MESHTOOLS3D_SRCS
  Chrono.cpp
  GetPot.cpp
  INRreader.cpp
  LaplaceSolver.cpp
  Mesh.cpp
  mgmres.cpp
  ThicknessEvaluation.cpp
  VtkWriter.cpp
)

add_library(m3d ${MESHTOOLS3D_SRCS} ${MESHTOOLS3D_INCLUDES})
