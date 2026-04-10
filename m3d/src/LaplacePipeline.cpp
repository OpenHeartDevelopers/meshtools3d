#include "../include/LaplacePipeline.hpp"
#include "../include/ThicknessEvaluation.hpp"
#include "../include/VtkWriter.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

std::set<long int> readVtxFile(const std::string& path)
{
    std::ifstream f(path.c_str());
    if (!f)
        throw std::runtime_error("Cannot open VTX file: " + path);

    std::set<long int> nodes;
    size_t count;
    f >> count;

    long int idx;
    while (f >> idx)
        nodes.insert(idx);

    return nodes;
}

void runLaplacePipeline(Mesh& mesh,
                        const LaplaceBCConfig& bc,
                        const LaplaceParams& params,
                        const OutputConfig& out)
{
    ThicknessEvaluation solver(&mesh);
    solver.set_abs_toll(params.abs_toll);
    solver.set_rel_toll(params.rel_toll);
    solver.set_max_it(params.itr_max);
    solver.set_Krilov_dim(params.dimKrilovSp);
    solver.set_verbosity(params.verbose);
    solver.set_algorithm(static_cast<unsigned char>(out.thickness_algorithm));
    solver.set_swapregions(bc.swap_regions);

    for (size_t i = 0; i < bc.zero_bc_files.size(); ++i)
        solver.setBCValue(readVtxFile(bc.zero_bc_files[i]), 0.0);
    for (size_t i = 0; i < bc.zero_bc_sets.size(); ++i)
        solver.setBCValue(bc.zero_bc_sets[i], 0.0);
    for (size_t i = 0; i < bc.one_bc_files.size(); ++i)
        solver.setBCValue(readVtxFile(bc.one_bc_files[i]), 1.0);
    for (size_t i = 0; i < bc.one_bc_sets.size(); ++i)
        solver.setBCValue(bc.one_bc_sets[i], 1.0);

    solver.solve();

    const std::string outBase = out.out_dir + "/" + out.out_name;

    if (out.eval_thickness) {
        if (solver.algorithm() != static_cast<unsigned char>(2))
            mesh.initializeConnectivities();
        solver.evalThickness();
        solver.writeElementGradient(outBase);
        mesh.writeTetraCentroids(outBase);
        mesh.writeTris(outBase);
    }

    const bool write_vtk_vars = out.out_vtk &&
                                 (out.out_potential || out.eval_thickness);

    VtkWriter writerVTK(&mesh, out.out_vtk_binary);
    if (write_vtk_vars) {
        writerVTK.setOutputDir(out.out_dir);
        writerVTK.setPrefixName(out.out_name);
        writerVTK.openFileForOutput();
    }

    if (out.out_potential) {
        if (write_vtk_vars)
            writerVTK.writeVariable(solver.sol(), "potential_func", VtkWriter::Scalar);
        else
            solver.writeSolution(outBase);
    }

    if (out.eval_thickness) {
        if (write_vtk_vars)
            writerVTK.writeVariable(solver.thickness(), "Thickness", VtkWriter::Scalar);
        else
            solver.writeThickness(outBase);
    }

    if (write_vtk_vars)
        writerVTK.CloseFile();
}
