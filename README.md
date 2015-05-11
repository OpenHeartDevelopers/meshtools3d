# meshtools3d

This project contains tools for creating 3D atrial mesh using CGAL 4.6 library. The actual code take as input
a segmentation in .inr file is required. Possible outputs are:

* INRIA mesh files (.mesh), evaluated directly by CGAL
* Carp mesh file (.elem and .pts), together with a set of point lists belonging to extracted regions(.vtx files); the output can be rescaled according to carp units
* vtk files (but without labels)

Some functionality are implemented but not used; in particular:

* evaluation and output of surface region points
* output in binary formats (CARP, vtk)
* triangle re-orientation with outward normals (when boudary is extracted)

### Surface region points
This functionality, actually disabled, determine endo, epii,  mitral valves nodes, etc  on the boundary, when these region are constrained by another one.

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




## TO DO

* check and in case icorrect binary output for carp meshes
* implement output for .mesh file also in Mesh class (could be useful?)
* implement output also for triangles (boundary elements)

