//  -*- c++ -*-
//  LaplacePipeline.hpp
//
//  Structs and pipeline function for running the Laplace solver and optional
//  thickness evaluation on a pre-built Mesh.  No CGAL dependency.

#ifndef _LAPLACEPIPELINE_HPP
#define _LAPLACEPIPELINE_HPP

#include <set>
#include <string>
#include <vector>
#include "Mesh.hpp"

// Solver tolerances and iteration limits.
struct LaplaceParams {
    double abs_toll    = 1e-6;
    double rel_toll    = 1e-6;
    int    itr_max     = 500;
    int    dimKrilovSp = 150;
    int    verbose     = 0;
};

// Boundary condition specification.
// Files (VTX format) and/or in-memory node sets may be mixed freely.
// swap_regions controls the thickness propagation direction in
// ThicknessEvaluation (mirrors the legacy swapregions data-file option).
struct LaplaceBCConfig {
    std::vector<std::string>          zero_bc_files;
    std::vector<std::string>          one_bc_files;
    std::vector<std::set<long int> >  zero_bc_sets;
    std::vector<std::set<long int> >  one_bc_sets;
    bool swap_regions = false;
};

// Output paths and format flags.
struct OutputConfig {
    std::string out_dir;
    std::string out_name;
    bool eval_thickness      = false;
    int  thickness_algorithm = 1;
    bool out_vtk             = false;
    bool out_vtk_binary      = false;
    bool out_potential       = false;
};

// Read a VTX file (count header + one node index per line) into a node set.
std::set<long int> readVtxFile(const std::string& path);

// Run the Laplace solve (and optionally thickness evaluation) on a mesh.
// The mesh must already have boundary labels evaluated if thickness is
// requested with algorithm 1 (Bishop), since that algorithm uses endoTria
// and epiTria.  unsetBoundaryLabels() may be called before this function
// as it preserves _Endo and _Epi.
void runLaplacePipeline(Mesh& mesh,
                        const LaplaceBCConfig& bc,
                        const LaplaceParams& params,
                        const OutputConfig& out);

#endif
