# meshtools3d

This project contains tools for creating 3D atrial mesh using CGAL 4.6 library. The actual code take as input
a segmentation in .inr file is required. Possible outputs are:

* INRIA mesh files (.mesh), evaluated directly by CGAL
* Carp mesh file (.elem and .pts), together with a set of point lists belonging to extracted regions(.vtx files); the output can be rescaled according to carp units
* vtk files (but without labels)

Some functionality are implemented but not used; in particular:

* CARP output in binary format (can be used by setting out_carp_binary = 1 inside output section but not tested till now)
* triangle re-orientation with outward normals (when boudary is extracted)

### Surface region points
This functionality determines endo, epi,  mitral valves nodes, etc  on the boundary, when these region are constrained by another one.
Endocardium and Epicardium surfaces are given as output and determined as the regions with the largest number of points (epicardium)
and the second region with the largest number of points (endocardium) within the extracted ones.

## Installation

It is assumed that CGAL-4.6 is installed; if not, download it and follows instructions at http://doc.cgal.org/latest/Manual/installation.html.
Compile requires cmake. Once cmake and CGAL are installed, type:

```sh
cmake -DCGAL_DIR=<path_to_CGAL>
```
in the same folder this file is located, specifing the right CGAL path. Then type make.

## Run

Copy the `data-template` file to data, and edit the contents as appropriate:

```sh
cp data-template data
vim data
```

The segmentation section tells where input is located, the meshing section sets
the mesh generation parameters (rescaleFactor conversely is only for output and
affects only output of vtk and CARP) and the output section sets the output
directory, name and formats (medit, carp, vtk).

## Parallel run with TBB library

If Intel TBB library is installed, the code is linked with TBB. TBB is a shared memory library, similar to OpenMP but more effcient. In this code the default number of core is setted to 1:
this to avoid to take all the cores found by the library if no number of core is specified. The implementation of this code will check that the environment variable 'TBB_NUM_THREADS' is
defined; for using for example 4 cores, on a shell window type:

```sh
  export TBB_NUM_THREADS = 4
```

Remark: TBB_NUM_THREADS is NOT a environment variable of TBB; it will affects only meshtools3d code

## Notes
* To disable TBB when compile, comment (with a # at the beginning of each line) the following lines in CMakeLists.txt:
```sh
if( TBB_FOUND )
  include(${TBB_USE_FILE})
  list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
endif()
```
* To create a debug version of the code, first generate the Makefile with:
```sh
cmake -DCMAKE_BUILD_TYPE=Debug -DCGAL_DONT_OVERRIDE_CMAKE_FLAGS=FALSE -DCGAL_DISABLE_ROUNDING_MATH_CHECK=ON
```
rounding mat checks have to be disabled if one want to use valgrind (look in https://github.com/openscad/openscad/issues/1340 for a similar problem)
and overriding of cmake flags has to be allowed

## TO DO

* check and in case correct binary output for carp meshes
* implement output for .mesh file also in Mesh class (could be useful?)
* implement output also for triangles (boundary elements)

## What's new

12 may 2015

* added a class "Chrono" for time profiling
* binary output for vtk (with correct endianess)

21 May 2015

* Laplace Solver for evaluating armonic extension
* Mesh class can evaluate Tetra and Tria centroids
* Laplace solver can evaluate Tetra Gradients for post-processing
* vtk visualize also output regions
* added a routine to write Laplace solution in a VTK file

22 May 2015

* added some examples
* added a routine to the thickness in a VTK file

11 Aug 2015

* Added parallel run facilities

27 Oct 2016

* deleted some variables when no longer needed to free memory before the evaluation of the potential
* implemented a class for read and handle .inr files
* implemented a search tree to localise boundary triangles within a specified bounding box
* implemented some CGAL extensions to use wrap the segmentation in a function that tells if a point is within the segmented domain or not
* created a new mode for mesh switching off segmentation labels during meshing and then labeling triangles only
* created a new example (sphereMultilabel) for meshing and then re-labeling

Dec 2019

* Created Dockerfile to create a docker image for `meshtools3d` to be used via docker (no compilation required)
