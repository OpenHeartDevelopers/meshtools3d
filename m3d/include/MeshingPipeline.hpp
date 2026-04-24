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
#include "Chrono.hpp"
#include "CGALDataType.hpp"

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

// Build criteria, run CGAL::make_mesh_3, write Medit, validate, and populate the
// CARP mesh. Templated over a Traits bundle (see CGALDataType.hpp) so the labeled
// and manual-segmentation branches of main.cpp share one body.
template<typename Traits, typename LabelFn>
void runCGALMeshing(typename Traits::Domain& domain,
                    const MeshingParams& p,
                    Mesh& carpMesh,
                    const std::string& out_dir,
                    const std::string& out_name,
                    bool out_medit,
                    LabelFn getLabel)
{
    typename Traits::FacetCriteria facet_criteria(p.facet_angle, p.facet_size, p.facet_distance);
    typename Traits::CellCriteria  cell_criteria(p.cell_rad_edge_ratio, p.cell_size);
    typename Traits::MeshCriteria  criteria(facet_criteria, cell_criteria);

    std::cout << "MESHING...";
    Chrono chrono;
    chrono.start();
    typename Traits::C3t3 c3t3 = CGAL::make_mesh_3<typename Traits::C3t3>(
        domain, criteria,
        CGAL::parameters::lloyd<bool>(), CGAL::parameters::odt<bool>(),
        CGAL::parameters::perturb<bool>(), CGAL::parameters::exude<bool>());
    chrono.stop();
    std::cout << " done in " << chrono << std::endl;

    if (out_medit)
        writeMeditFile(c3t3, out_dir, out_name);
    validateTriangulation(c3t3, out_dir, out_name, out_medit);
    populateCarpMeshFromC3t3(c3t3, carpMesh, getLabel);
}

#endif // MESHINGPIPELINE_HPP
