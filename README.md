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

Specify options inside a data file as data. Segmentation section tells
where input is located; meshing section the mesh generation parameters (rescaleFactor conversely is only for output and affects only output of vtk and CARP); output section tells where to write the output, the output name and the formats to write as output (medit, carp, vtk)

## Parallel run with TBB library

If Intel TBB library is installed, the code is linked with TBB. TBB is a shared memory library, similar to OpenMP but more effcient. In this code the default number of core is setted to 1:
this to avoid to take all the cores found by the library if no number of core is specified. The implementation of this code will check that the environment variable 'TBB_NUM_THREADS' is
defined; for using for example 4 cores, on a shell window type:

```sh
  export TBB_NUM_THREADS = 4
```

Remark: TBB_NUM_THREADS is NOT a environment variable of TBB; it will affects only meshtools3d code

## TO DO

* check and in case icorrect binary output for carp meshes
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
<<<<<<< HEAD
* added a routine to write Laplace solution in a VTK file

22 May 2015

* Added some examples
=======
* added a routine to write Laplace solution in a VTK file

11 Aug 2015
* Added parallel run facilities



