//  -*- c++ -*-
//  MeshingPipeline.hpp
//
//  Templated post-meshing helpers shared by both CGAL meshing branches in
//  main.cpp.  Included only by main.cpp; template instantiation happens there
//  with the concrete C3t3 / C3t3_manualseg types from CGALDataType.hpp.

#ifndef MESHINGPIPELINE_HPP
#define MESHINGPIPELINE_HPP

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include "Mesh.hpp"

// Write c3t3 to <out_dir>/<out_name>.mesh in Medit format.
template<typename C3t3>
void writeMeditFile(const C3t3& c3t3,
                    const std::string& out_dir,
                    const std::string& out_name)
{
    std::string path = out_dir + "/" + out_name + ".mesh";
    std::ofstream f(path.c_str());
    c3t3.output_to_medit(f);
    f.close();
}

// Check that the triangulation has enough vertices and cells.
// If degenerate: writes the Medit file (unless already written via out_medit),
// prints a diagnostic, and exits.  Otherwise prints vertex/tet counts.
template<typename C3t3>
void validateTriangulation(const C3t3& c3t3,
                           const std::string& out_dir,
                           const std::string& out_name,
                           bool out_medit)
{
    const auto& triang = c3t3.triangulation();
    const size_t nbPt  = static_cast<size_t>(triang.number_of_vertices());
    const size_t nbTet = static_cast<size_t>(c3t3.number_of_cells_in_complex());

    if (nbPt < 3 || nbTet < 1)
    {
        std::cerr << "Problem with Triangulations, only "
                  << nbPt << " Vertices and " << nbTet << " Thetraedra" << std::endl;
        if (!out_medit)
            writeMeditFile(c3t3, out_dir, out_name);
        exit(1);
    }

    std::cout << "Vertices: " << nbPt << std::endl;
    std::cout << "Tetra:    " << nbTet << std::endl;
}

// Populate carpMesh vertices and tetrahedra from c3t3.
// getLabel is called with the cell iterator and must return an int region label.
template<typename C3t3, typename RegionLabelFn>
void populateCarpMeshFromC3t3(const C3t3& c3t3, Mesh& carpMesh, RegionLabelFn getLabel)
{
    const auto& triang = c3t3.triangulation();
    const size_t nbPt  = static_cast<size_t>(triang.number_of_vertices());
    const size_t nbTet = static_cast<size_t>(c3t3.number_of_cells_in_complex());

    carpMesh.initializePtVector(nbPt);
    size_t icount = 0;
    std::map<typename C3t3::Triangulation::Vertex_handle, size_t> vertices;

    for (auto it = triang.finite_vertices_begin();
              it != triang.finite_vertices_end(); ++it)
    {
        vertices[it] = icount;
        carpMesh.Pt(icount).x() = it->point().x();
        carpMesh.Pt(icount).y() = it->point().y();
        carpMesh.Pt(icount).z() = it->point().z();

        for (unsigned char ic = 0; ic < 3; ++ic)
        {
            carpMesh.Info().baricenter[ic] += carpMesh.Pt(icount).coord[ic];
            carpMesh.Info().checksum[ic]   += static_cast<float>(carpMesh.Pt(icount).coord[ic]);
            if (carpMesh.Info().bbox[ic][0] > carpMesh.Pt(icount).coord[ic])
                carpMesh.Info().bbox[ic][0] = carpMesh.Pt(icount).coord[ic];
            if (carpMesh.Info().bbox[ic][1] < carpMesh.Pt(icount).coord[ic])
                carpMesh.Info().bbox[ic][1] = carpMesh.Pt(icount).coord[ic];
        }
        ++icount;
    }

    for (unsigned char ic = 0; ic < 3; ++ic)
        carpMesh.Info().baricenter[ic] /= static_cast<double>(nbPt);

    carpMesh.initializeTetraVector(nbTet);
    icount = 0;
    for (auto itc = c3t3.cells_in_complex_begin();
              itc != c3t3.cells_in_complex_end(); ++itc)
    {
        for (unsigned char iv = 0; iv < 4; ++iv)
            carpMesh.Tet(icount).vertex[iv] = vertices.at(itc->vertex(iv));
        carpMesh.Tet(icount).regionLabel = getLabel(itc);
        ++icount;
    }
}

#endif // MESHINGPIPELINE_HPP
